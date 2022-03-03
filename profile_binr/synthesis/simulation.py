"""
    Docstring
"""

import copy
import multiprocessing as mp
from typing import Optional, Union

# TODO : change import scheme
from .sampling import _dropout_mask, pd, np, ss

# TODO : add a module-wide rng ?
_DEFAULT_RNG_SEED: int = 1928374650
_MODULE_SEED_SEQUENCE_GENERATOR = np.random.SeedSequence(_DEFAULT_RNG_SEED)
# TODO : figure a way to remove the magic number ?
_MODULE_CHILD_SEEDS = _MODULE_SEED_SEQUENCE_GENERATOR.spawn(5)
# we need 5, one for each function besides biased_simulation_from_binary_state()
_MODULE_RNGS = [np.random.default_rng(s) for s in _MODULE_CHILD_SEEDS]
for m in _MODULE_RNGS:
    print(m.__getstate__()["state"]["state"])
## TODO : define this function
# def mean_binariser(binary_df, normalized_counts_df):
#    pass

# TODO : give this a better names ?
def random_nan_binariser(df):
    """Assuming df is a binary matrix produced by calling
    profile_binr.ProfileBin.binarise() this function should
    randomly resolve all NaN entries to either 1 or 0

    This function operates on a deep copy of the df argument
    """

    df = copy.deepcopy(df)
    for gene in df.columns:
        mask = df[gene].isna()
        df.loc[mask, gene] = np.random.choice((1.0, 0, 0), mask.sum())

    return df


def simulate_unimodal_distribution(
    value: float,
    mean: float,
    std_dev: float,
    size: int,
    rng: Union[np.random._generator.Generator, int] = _MODULE_RNGS[0],
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
    rng: Union[np.random._generator.Generator, int] = _MODULE_RNGS[1],
) -> pd.Series:
    """
    R
    """
    rng = np.random.default_rng(rng) if isinstance(rng, int) else rng
    assert binary_gene.name == criterion.name, "Criterion and gene mismatch"
    # binary_gene = fully_bin.loc[:, _unimod_gene]
    # allocate array for simulated data
    simulated_normalised_expression = pd.Series(
        np.nan, index=binary_gene.index, name=criterion.name, dtype=float
    )
    # Create masks
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

        # calculate a correction mask (vector of 1 and 0 at random indices, according to the correction DOR)
        simulated_normalised_expression *= _dropout_mask(
            dropout_rate=_correction_dor,
            size=len(simulated_normalised_expression),
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
    rng: Union[np.random._generator.Generator, int] = _MODULE_RNGS[2],
) -> pd.Series:
    """
    R
    """
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

        # calculate a correction mask (vector of 1 and 0 at random indices, according to the correction DOR)
        simulated_normalised_expression *= _dropout_mask(
            dropout_rate=_correction_dor, size=len(simulated_normalised_expression)
        )
        corrected_dor = np.isclose(simulated_normalised_expression, 0).mean()
        if _verbose:
            print("Done")
            print(f"Corrected DropOutRate : {corrected_dor}")

    return simulated_normalised_expression


def simulate_gene(
    binary_gene: pd.Series,
    criterion: pd.Series,
    rng: Union[np.random._generator.Generator, int] = _MODULE_RNGS[3],
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
    rng: Union[np.random._generator.Generator, int] = _MODULE_RNGS[4],
) -> pd.Series:
    """ helper function, wrapper for apply """
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

    n_threads = min(n_threads, mp.cpu_count()) if n_threads else mp.cpu_count()
    # verify binarised genes are contained in the simulation criteria index
    if not all(x in simulation_criteria.index for x in binary_df.columns):
        raise ValueError(
            "'binary_df' contains at least one gene for which there is no simulation criterion."
        )
    if n_samples is not None:
        binary_df = pd.concat(n_samples * [binary_df], ignore_index=False)
    # match the order
    simulation_criteria = simulation_criteria.loc[binary_df.columns, :]

    # Generate independent rngs, with good quality seeds
    # according to :
    # https://numpy.org/doc/stable/reference/random/parallel.html
    # seed = seed or _DEFAULT_RNG_SEED
    _seed_sequence_generator = np.random.SeedSequence(seed)
    child_seeds = _seed_sequence_generator.spawn(n_threads)
    streams = [np.random.default_rng(s) for s in child_seeds]

    _criteria_ls = np.array_split(simulation_criteria, n_threads)
    _binary_ls = np.array_split(binary_df, n_threads, axis=1)
    with mp.Pool(n_threads) as pool:
        # args = ((_bin, _crit) for _bin, _crit in zip(_binary_ls, _criteria_ls))
        ret_list = pool.starmap(
            _simulate_subset, zip(_binary_ls, _criteria_ls, streams)
        )

    # return _binary_ls, _criteria_ls
    return pd.concat(ret_list, axis=1)
