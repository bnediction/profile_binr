"""
Module to sample genes from learnt criteria
"""
import warnings
import functools
import multiprocessing

from typing import Tuple, Callable, List, Optional, Union

import scipy.stats as ss
import numpy as np
import pandas as pd

_RandType = Union[np.random.Generator, int]
_GLOBAL_RNG = np.random.default_rng()


def set_module_rng_seed(seed: int) -> np.random.Generator:
    """Fix the module-wide random number generator
    as well as the R seed.

    parameters
    ----------
            seed : any seed accepted by numpy.random.default_rng
                   see `help(numpy.random.default_rng)`
                   for more details.

    returns
    -------
            The package-wide global random number generator
            (result of calling `numpy.random.default_rng(seed)`)

    Caveats:
        TL;DR
            If the seed is not an integer, it will not be set on
            the embedded R instance.

        Seed has a type hint 'int' because this is the recommended
        type. This function attempts to set the seed both on the
        python (numpy) and R sides. To our knowledge, there is no
        explicit conversion rule provided by rpy2 for other types
        of seeds accepted by numpy.random.default_rng().
        If a numpy object is provided to seed the global generator
        (such as SeedSequence, BitGenerator, Generator), nothing
        will happen on the R side and a warning message will be
        raised printed.
    """
    # declare global to override the generator
    # that has been imported on different sections
    global _GLOBAL_RNG
    _GLOBAL_RNG = np.random.default_rng(seed)
    if isinstance(seed, int):
        # import is performed within the function,
        # to prevent an import error due to circular import
        from ..wrappers.probinr import ProfileBin

        _probin = ProfileBin(pd.DataFrame())
        _probin.r(f"set.seed({seed})")
    else:
        _no_seed_warning_message_ls = [
            f"seed of type {type(seed)} can not be passed to the embedded R instance",
            "the seed was not set on the R side.",
        ]
        warnings.warn(" ".join(_no_seed_warning_message_ls))

    return _GLOBAL_RNG


def __sim_zero_inf(
    _lambda: float,
    size: int,
    ignore_deprecation: bool = False,
    rng: Optional[_RandType] = None,
) -> np.ndarray:
    """DEPRECATED . Simulate zero-inflated genes
    This function samples from an exponential distribution with parameter
    lambda
    """
    rng = rng or _GLOBAL_RNG
    rng = np.random.default_rng(rng) if isinstance(rng, int) else rng
    if not ignore_deprecation:
        __err_message = ["Error: ", "ZeroInf genes cannot be directly sampled"]
        raise DeprecationWarning("".join(__err_message))
    return rng.exponential(scale=1 / _lambda, size=size)


def _sim_unimodal(
    mean: float, std_dev: float, size: int, rng: Optional[_RandType] = None
) -> np.ndarray:
    """Simulate the expression of unimodal genes using
    a normal (gaussian) distribution."""
    rng = rng or _GLOBAL_RNG
    rng = np.random.default_rng(rng) if isinstance(rng, int) else rng
    return rng.normal(loc=mean, scale=std_dev, size=size)


def _sim_bimodal(
    mean1: float,
    mean2: float,
    std_dev: float,
    weights: Tuple[float],
    size: int,
    rng: Optional[_RandType] = None,
):
    """Simulate bimodal genes. The modelling is achieved using a gaussian mixture.
    The variance is assumed to be tied i.e. both components of the mixture have the same
    variance.

    Parametres
    ----------

    Weights should be a tuple (or eventually a list) containing
    """
    rng = rng or _GLOBAL_RNG
    rng = np.random.default_rng(rng) if isinstance(rng, int) else rng
    # Parameters of the mixture components
    norm_params = np.array([[mean1, std_dev], [mean2, std_dev]])
    # A stream of indices from which to choose the component
    mixture_idx = rng.choice(len(weights), size=size, replace=True, p=weights)
    # y is the mixture sample
    return np.fromiter(
        (ss.norm.rvs(*(norm_params[i]), random_state=rng) for i in mixture_idx),
        dtype=np.float64,
    )


def _dropout_mask(
    dropout_rate: float, size: int, rng: Optional[_RandType] = None
) -> np.ndarray:
    """ Dropout mask to obtain the same dropout_rate as originally estimated"""
    rng = rng or _GLOBAL_RNG
    rng = np.random.default_rng(rng) if isinstance(rng, int) else rng
    return rng.choice(
        (0, 1), size=size, replace=True, p=(dropout_rate, 1.0 - dropout_rate)
    )


def sample_gene(
    criterion: pd.Series,
    n_samples: int,
    enforce_dropout_rate: bool = True,
    rng: Optional[_RandType] = None,
) -> pd.Series:
    """Simulate the expression of a gene, using the information provided by
    the criteria dataframe of a profile_binr.ProfileBin class.

    Parametres
    ----------

    criterion : an entry (row) of the criteria dataframe of a ProfileBin class,
    trained on the dataset which you want to sample.

    n_samples : number of samples to generate

    enforce_dropout_rate : should random entries of the gene be set to zero
    in order to preserve the dropout rate estimated whilst computing the criteria
    for the original expression dataset ?
    """
    rng = rng or _GLOBAL_RNG
    rng = np.random.default_rng(rng) if isinstance(rng, int) else rng
    _data: np.array

    if criterion["Category"] == "Discarded":
        _data = np.full(n_samples, np.nan)
    elif criterion["Category"] == "Unimodal":
        _data = _sim_unimodal(
            mean=criterion["mean"],
            std_dev=np.sqrt(criterion["variance"]),
            size=n_samples,
            rng=rng,
        )
        # TODO : discuss the validity of this approach
        if enforce_dropout_rate:
            _data *= _dropout_mask(
                dropout_rate=criterion["DropOutRate"], size=n_samples, rng=rng
            )
    elif criterion["Category"] == "Bimodal":
        _data = _sim_bimodal(
            mean1=criterion["gaussian_mean1"],
            mean2=criterion["gaussian_mean2"],
            std_dev=np.sqrt(criterion["gaussian_variance"]),
            weights=(criterion["gaussian_prob1"], criterion["gaussian_prob2"]),
            size=n_samples,
            rng=rng,
        )
        if enforce_dropout_rate:
            _data *= _dropout_mask(
                dropout_rate=criterion["DropOutRate"], size=n_samples, rng=rng
            )
    elif criterion["Category"] == "ZeroInf":
        _data = __sim_zero_inf(_lambda=criterion["lambda"], size=n_samples, rng=rng)
    else:
        raise ValueError(f"Unknown category `{criterion['Category']}`, aborting")

    return pd.Series(data=_data, name=criterion.name)


def _sample_sequential(
    df: pd.DataFrame,
    n_samples: int,
    enforce_dropout_rate: bool = True,
    rng: Optional[_RandType] = None,
) -> pd.DataFrame:
    """ Simulate samples from a criteria dataframe, sequentially """
    rng = rng or _GLOBAL_RNG
    rng = np.random.default_rng(rng) if isinstance(rng, int) else rng
    return (
        df.apply(
            lambda y: sample_gene(
                y,
                n_samples=n_samples,
                enforce_dropout_rate=enforce_dropout_rate,
                rng=rng,
            ),
            axis=1,
        ).T.dropna(how="all", axis=1),
    )[0]


def sample_from_criteria(
    criteria: pd.DataFrame,
    n_samples: int,
    enforce_dropout_rate: bool = True,
    n_workers: int = multiprocessing.cpu_count(),
) -> pd.DataFrame:
    """
    Create a new expression dataframe from a criteria dataframe using multithreading

    WARNING : this function has no seeding mechanism and will use
    the module's numpy.random.Generator which is not thread-safe.
    Expect undefined behaviour
    """
    _partial_simulation_function: Callable = functools.partial(
        _sample_sequential,
        n_samples=n_samples,
        enforce_dropout_rate=enforce_dropout_rate,
    )

    _df_splitted_ls: List[pd.DataFrame] = np.array_split(criteria, n_workers)
    with multiprocessing.Pool(n_workers) as pool:
        ret_list = pool.map(_partial_simulation_function, _df_splitted_ls)

    return pd.concat(ret_list, axis=1)
