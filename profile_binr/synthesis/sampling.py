"""
Module to simulate genes from learnt criteria
"""
import functools
import multiprocessing

from typing import Tuple, Callable, List

import scipy.stats as ss
import numpy as np
import pandas as pd


def __sim_zero_inf(
    _lambda: float, size: int, ignore_deprecation: bool = False
) -> np.ndarray:
    """DEPRECATED . Simulate zero-inflated genes
    This function samples from an exponential distribution with parameter
    lambda
    """
    if not ignore_deprecation:
        __err_message = ["Error: ", "ZeroInf genes cannot be directly simulated"]
        raise DeprecationWarning("".join(__err_message))
    return np.random.exponential(scale=1 / _lambda, size=size)


def _sim_unimodal(mean: float, std_dev: float, size: int) -> np.ndarray:
    """Simulate the expression of unimodal genes using
    a normal (gaussian) distribution."""
    return np.random.normal(loc=mean, scale=std_dev, size=size)


def _sim_bimodal(
    mean1: float, mean2: float, std_dev: float, weights: Tuple[float], size: int
):
    """Simulate bimodal genes. The modelling is achieved using a gaussian mixture.
    The variance is assumed to be tied i.e. both components of the mixture have the same
    variance.

    Parametres
    ----------

    Weights should be a tuple (or eventually a list) containing
    """
    # Parameters of the mixture components
    norm_params = np.array([[mean1, std_dev], [mean2, std_dev]])

    # A stream of indices from which to choose the component
    mixture_idx = np.random.choice(len(weights), size=size, replace=True, p=weights)
    # y is the mixture sample
    return np.fromiter(
        (ss.norm.rvs(*(norm_params[i])) for i in mixture_idx), dtype=np.float64
    )


def _dropout_mask(dropout_rate: float, size: int) -> np.ndarray:
    """ Dropout mask to obtain the same dropout_rate as originally estimated"""
    return np.random.choice(
        (0, 1), size=size, replace=True, p=(dropout_rate, 1.0 - dropout_rate)
    )


def simulate_gene(
    criterion: pd.Series, n_samples: int, enforce_dropout_rate: bool = True
) -> pd.Series:
    """Simulate the expression of a gene, using the information provided by
    the criteria dataframe of a profile_binr.ProfileBin class.

    Parametres
    ----------

    criterion : an entry (row) of the criteria dataframe of a ProfileBin class,
    trained on the dataset which you want to simulate.

    n_samples : number of samples to generate

    enforce_dropout_rate : should random entries of the gene be set to zero
    in order to preserve the dropout rate estimated whilst computing the criteria
    for the original expression dataset ?
    """
    _data: np.array

    if criterion["Category"] == "Discarded":
        _data = np.full(n_samples, np.nan)
    elif criterion["Category"] == "Unimodal":
        _data = _sim_unimodal(
            mean=criterion["mean"],
            std_dev=np.sqrt(criterion["variance"]),
            size=n_samples,
        )
        # TODO : discuss the validity of this approach
        if enforce_dropout_rate:
            _data *= _dropout_mask(
                dropout_rate=criterion["DropOutRate"], size=n_samples
            )
    elif criterion["Category"] == "Bimodal":
        _data = _sim_bimodal(
            mean1=criterion["gaussian_mean1"],
            mean2=criterion["gaussian_mean2"],
            std_dev=np.sqrt(criterion["gaussian_variance"]),
            weights=(criterion["gaussian_prob1"], criterion["gaussian_prob2"]),
            size=n_samples,
        )
        if enforce_dropout_rate:
            _data *= _dropout_mask(
                dropout_rate=criterion["DropOutRate"], size=n_samples
            )
    elif criterion["Category"] == "ZeroInf":
        _data = __sim_zero_inf(_lambda=criterion["lambda"], size=n_samples)
    else:
        raise ValueError(f"Unknown category `{criterion['Category']}`, aborting")

    return pd.Series(data=_data, name=criterion.name)


def _simulate_sequential(
    df: pd.DataFrame, n_samples: int, enforce_dropout_rate: bool = True
) -> pd.DataFrame:
    """ Simulate samples from a criteria dataframe, sequentially """
    return (
        df.apply(
            lambda y: simulate_gene(
                y, n_samples=n_samples, enforce_dropout_rate=enforce_dropout_rate
            ),
            axis=1,
        ).T.dropna(how="all", axis=1),
    )[0]


def simulate_from_criteria(
    criteria: pd.DataFrame,
    n_samples: int,
    enforce_dropout_rate: bool = True,
    n_workers: int = multiprocessing.cpu_count(),
) -> pd.DataFrame:
    """ Create a new expression dataframe from a criteria dataframe using multithreading  """
    _partial_simulation_function: Callable = functools.partial(
        _simulate_sequential,
        n_samples=n_samples,
        enforce_dropout_rate=enforce_dropout_rate,
    )

    _df_splitted_ls: List[pd.DataFrame] = np.array_split(criteria, n_workers)
    with multiprocessing.Pool(n_workers) as pool:
        ret_list = pool.map(_partial_simulation_function, _df_splitted_ls)

    return pd.concat(ret_list, axis=1)
