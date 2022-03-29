"""
    Module to simulate genes from learnt criteria
    and/or sample their distribution.
"""
import warnings
import functools
import multiprocessing as mp

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
        from .wrappers.probinr import ProfileBin

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
    criteria: pd.DataFrame,
    n_samples: int,
    enforce_dropout_rate: bool = True,
    rng: Optional[_RandType] = None,
) -> pd.DataFrame:
    """ Simulate samples from a criteria dataframe, sequentially """
    rng = rng or _GLOBAL_RNG
    rng = np.random.default_rng(rng) if isinstance(rng, int) else rng
    return (
        criteria.apply(
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
    n_workers: int = mp.cpu_count(),
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
    with mp.Pool(n_workers) as pool:
        ret_list = pool.map(_partial_simulation_function, _df_splitted_ls)

    return pd.concat(ret_list, axis=1)


def random_nan_binariser(binary_df: pd.DataFrame, rng: Optional[_RandType] = None):
    """Assuming df is a binary matrix produced by calling
    profile_binr.ProfileBin.binarise() this function should
    randomly resolve all NaN entries to either 1 or 0

    This function operates on a deep copy of the df argument
    """
    rng = rng or _GLOBAL_RNG
    rng = np.random.default_rng(rng) if isinstance(rng, int) else rng
    binary_df = binary_df.copy(deep=True)
    for gene in binary_df.columns:
        mask = binary_df[gene].isna()
        binary_df.loc[mask, gene] = rng.choice((1.0, 0.0), mask.sum())

    return binary_df


def simulate_unimodal_distribution(
    value: float,
    mean: float,
    std_dev: float,
    size: int,
    rng: Optional[_RandType] = None,
) -> pd.Series:
    """
    Basically a wrapper for scipy.stats.halfnorm

    params:
    ------
            * value   : 0 or 1.
                * 0   : The left half of the normal distribution
                * 1   : The right half of the normal distribution
            * mean    : Mean of the normal distribution whose half should be sampled
            * std_dev : The standard deviation of the normal distribution ... ^
            * size    : The number of random variates to be sampled from the
                        half normal distribution.
    """
    rng = rng or _GLOBAL_RNG
    rng = np.random.default_rng(rng) if isinstance(rng, int) else rng
    # random variate right hand side
    _rv_rhs = ss.halfnorm.rvs(loc=mean, scale=std_dev, size=size, random_state=rng)
    # random variate left hand side
    _rv_lhs = mean - (_rv_rhs - mean)
    return _rv_rhs if np.isclose(value, 1.0) else _rv_lhs


def simulate_bimodal_gene(
    binary_gene: pd.Series,
    criterion: pd.Series,
    _verbose: bool = False,
    rng: Optional[_RandType] = None,
) -> pd.Series:
    """
    R
    """
    rng = rng or _GLOBAL_RNG
    rng = np.random.default_rng(rng) if isinstance(rng, int) else rng
    assert binary_gene.name == criterion.name, "Criterion and gene mismatch"
    # binary_gene = fully_bin.loc[:, _unimod_gene]
    # allocate array for simulated data
    simulated_normalised_expression = pd.Series(
        np.nan, index=binary_gene.index, name=criterion.name, dtype=float
    )
    # Create masks
    # apply is marginally faster, and more idiomatic so
    # there is no need for :
    # (1 - np.array([True, False, False])).astype(bool)
    one_mask = binary_gene > 0
    zero_mask = one_mask.apply(lambda x: not x)

    assert sum(one_mask) + sum(zero_mask) == len(
        binary_gene
    ), "Floating point comparison error."
    if _verbose:
        print(
            f"Lengths\n\tOnes({sum(one_mask)}) + Zeros({sum(zero_mask)}) = ArrayLength({len(binary_gene)}) *"
        )
        print(
            f"Proportions\n\tOnes({sum(one_mask)/len(binary_gene)}) + Zeros({sum(zero_mask)/len(binary_gene)})"
        )

    simulated_from_ones = ss.norm.rvs(
        loc=criterion["gaussian_mean2"],
        scale=np.sqrt(criterion["gaussian_variance"]),
        size=sum(one_mask),  # change for one_mask.sum() ?
        random_state=rng,
    )
    simulated_from_zeros = ss.norm.rvs(
        loc=criterion["gaussian_mean1"],
        scale=np.sqrt(criterion["gaussian_variance"]),
        size=sum(zero_mask),  # change for one_mask.sum() ?
        random_state=rng,
    )
    # First approach to simulating the DropOutRate,
    # put all negative simulated values to zero
    simulated_negative = simulated_from_zeros < 0.0
    simulated_from_zeros[simulated_negative] = 0.0
    natural_dor = sum(simulated_negative) / len(binary_gene)
    # This does not seem to be sufficient for most cases!

    simulated_normalised_expression[one_mask] = simulated_from_ones
    simulated_normalised_expression[zero_mask] = simulated_from_zeros

    # Correct by randomly putting values smaller than the bin threshold to zero :
    if _verbose:
        print(f"natural DropOutRate = {natural_dor}")
    if natural_dor < criterion["DropOutRate"]:
        if _verbose:
            print("Correcting DropOutRate...", end="\t")
        # check how many values do we need to put to zero
        _correction_dor = criterion["DropOutRate"] - natural_dor

        # calculate a correction mask
        # (vector of 1 and 0 at random indices, according to the correction DOR)
        simulated_normalised_expression *= _dropout_mask(
            dropout_rate=_correction_dor,
            size=len(simulated_normalised_expression),
            rng=rng,
        )
        corrected_dor = np.isclose(simulated_normalised_expression, 0).mean()
        if _verbose:
            print("Done")
            print(f"Corrected DropOutRate : {corrected_dor}")

    return simulated_normalised_expression


def simulate_unimodal_gene(
    binary_gene: pd.Series,
    criterion: pd.Series,
    _verbose: bool = False,
    rng: Optional[_RandType] = None,
) -> pd.Series:
    """
    R
    """
    rng = rng or _GLOBAL_RNG
    rng = np.random.default_rng(rng) if isinstance(rng, int) else rng
    assert binary_gene.name == criterion.name, "Criterion and gene mismatch"
    # binary_gene = fully_bin.loc[:, _unimod_gene]
    # allocate array for simulated data
    simulated_normalised_expression = pd.Series(
        0.0, index=binary_gene.index, name=criterion.name, dtype=float
    )
    # Create masks
    one_mask = binary_gene > 0
    # change negation to 1 - one_mask ? (check if performance gain is worth it)
    zero_mask = one_mask.apply(lambda x: not x)

    assert sum(one_mask) + sum(zero_mask) == len(
        binary_gene
    ), "Floating point comparison error."
    if _verbose:
        print(
            f"Lengths\n\tOnes({sum(one_mask)}) + Zeros({sum(zero_mask)}) = ArrayLength({len(binary_gene)}) *"
        )
        print(
            f"Proportions\n\tOnes({sum(one_mask)/len(binary_gene)}) + Zeros({sum(zero_mask)/len(binary_gene)})"
        )

    simulated_from_ones = simulate_unimodal_distribution(
        1.0,
        mean=criterion["mean"],
        std_dev=np.sqrt(criterion["variance"]),
        size=sum(one_mask),  # change for one_mask.sum() ?
        rng=rng,
    )

    simulated_from_zeros = simulate_unimodal_distribution(
        0.0,
        mean=criterion["mean"],
        std_dev=np.sqrt(criterion["variance"]),
        size=sum(zero_mask),
        rng=rng,
    )

    # First approach to simulating the DropOutRate,
    # put all negative simulated values to zero
    simulated_negative = simulated_from_zeros < 0.0
    simulated_from_zeros[simulated_negative] = 0.0
    natural_dor = sum(simulated_negative) / len(binary_gene)
    # This does not seem to be sufficient for most cases!

    simulated_normalised_expression[one_mask] = simulated_from_ones
    simulated_normalised_expression[zero_mask] = simulated_from_zeros

    # Correct by randomly putting values smaller than the bin threshold to zero :
    if _verbose:
        print(f"natural DropOutRate = {natural_dor}")
    if natural_dor < criterion["DropOutRate"]:
        if _verbose:
            print("Correcting DropOutRate...", end="\t")
        # check how many values do we need to put to zero
        _correction_dor = criterion["DropOutRate"] - natural_dor

        # calculate a correction mask (vector of 1 and 0 at random indices,
        # according to the correction DOR)
        simulated_normalised_expression *= _dropout_mask(
            dropout_rate=_correction_dor,
            size=len(simulated_normalised_expression),
            rng=rng,
        )
        corrected_dor = np.isclose(simulated_normalised_expression, 0).mean()
        if _verbose:
            print("Done")
            print(f"Corrected DropOutRate : {corrected_dor}")

    return simulated_normalised_expression


def simulate_gene(
    binary_gene: pd.Series,
    criterion: pd.Series,
    rng: Optional[_RandType] = None,
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
    rng = rng or _GLOBAL_RNG
    rng = np.random.default_rng(rng) if isinstance(rng, int) else rng

    if criterion["Category"] == "Discarded":
        return pd.Series(
            np.nan, index=binary_gene.index, name=criterion.name, dtype=float
        )
    elif criterion["Category"] == "Unimodal":
        return simulate_unimodal_gene(binary_gene, criterion, rng=rng)
    elif criterion["Category"] == "Bimodal":
        return simulate_bimodal_gene(binary_gene, criterion, rng=rng)
    else:
        raise ValueError(f"Unknown category `{criterion['Category']}`, aborting")


def _simulate_subset(
    binary_df: pd.DataFrame,
    simulation_criteria: pd.DataFrame,
    rng: Optional[_RandType] = None,
) -> pd.Series:
    """ helper function, wrapper for apply """
    rng = rng or _GLOBAL_RNG
    rng = np.random.default_rng(rng) if isinstance(rng, int) else rng
    return binary_df.apply(
        lambda x: simulate_gene(x, simulation_criteria.loc[x.name, :], rng=rng)
    )


def biased_simulation_from_binary_state(
    binary_df: pd.DataFrame,
    simulation_criteria: pd.DataFrame,
    n_threads: Optional[int] = None,
    n_samples: Optional[int] = None,
    seed: Optional[int] = None,
) -> pd.DataFrame:
    """n_threads defaults to multiprocessing.cpu_count()"""

    n_threads = min(abs(n_threads), mp.cpu_count()) if n_threads else mp.cpu_count()
    # verify binarised genes are contained in the simulation criteria index
    if not all(x in simulation_criteria.index for x in binary_df.columns):
        raise ValueError(
            "'binary_df' contains at least one gene for which there is no simulation criterion."
        )
    if n_samples is not None:
        binary_df = pd.concat(n_samples * [binary_df], ignore_index=False)
    # match the order
    simulation_criteria = simulation_criteria.loc[binary_df.columns, :]

    # Generate independent random number generators, with good quality seeds
    # according to :
    # https://numpy.org/doc/stable/reference/random/parallel.html
    _seed_sequence_generator = np.random.SeedSequence(seed)
    child_seeds = _seed_sequence_generator.spawn(n_threads)
    streams = [np.random.default_rng(s) for s in child_seeds]

    _criteria_ls = np.array_split(simulation_criteria, n_threads)
    _binary_ls = np.array_split(binary_df, n_threads, axis=1)
    with mp.Pool(n_threads) as pool:
        ret_list = pool.starmap(
            _simulate_subset, zip(_binary_ls, _criteria_ls, streams)
        )

    return pd.concat(ret_list, axis=1)