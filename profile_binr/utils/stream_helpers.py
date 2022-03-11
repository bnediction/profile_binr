"""
    Temporarily storing helper functions for branching event
    reconstruction.
"""

from typing import Iterator, Dict, Tuple, Union, List, Optional
from functools import partial
import numpy as np
import pandas as pd
from scipy.spatial.distance import jaccard as jaccard_distance

from ..simulation import biased_simulation_from_binary_state

RandomWalkGenerator = Union[Iterator[Dict[str, int]], List[Dict[str, int]]]
IntListOrTuple = Union[Tuple[int, int], List[int]]


def jaccard_similarity(u, v, w=None) -> float:
    """Compute the Jaccard similarity index
    defined as the cardinal of the intersection
    divided by the cardinal of the union.

    returns : 1.0 - scipy.spatial.distance.jaccard(u, v, w)
    """
    return 1.0 - jaccard_distance(u, v, w)


def random_walk_to_data_frame(
    random_walk_generator: RandomWalkGenerator,
) -> pd.DataFrame:
    """Build and return a DataFrame from a random walk.

    A random walk generator is obtained as follows (minimal example) :
    >>> bn: colomoto.minibn.BooleanNetwork
    >>> dyn = colomoto.minibn.FullyAsynchronousDynamics(bn)
    >>> initial_state: Dict[str,int] = {'gene1': 0, 'gene2': 1}
    >>> n_steps = 100
    >>> _random_walk = dyn.random_walk(initial_state, steps=n_steps)

    The `_random_walk` objet should be passed to this function.

    Returns:
        A DataFrame containing the genes (variables) as columns and one row
        for each observed state in the random walk. The index is simply
        the ordered set of integers in the range [0, n_steps).

    """
    _rnd_wlks = list(random_walk_generator)
    _rnd_walk_df = pd.concat(
        [pd.DataFrame(walk, index=[i]) for i, walk in enumerate(_rnd_wlks)], axis="rows"
    )
    _rnd_walk_df.index.name = "step"
    return _rnd_walk_df


def merge_random_walks(
    walk_a: pd.DataFrame,
    walk_b: pd.DataFrame,
    label_a: str,
    label_b: str,
    branching_point_query: str,
    attractor_a_query: str,
    attractor_b_query: str,
    distinguish_attractors: bool = True,
    confound_observations: bool = True,
) -> pd.DataFrame:
    """

    Merge random walks in order to simulate a bifurcation process.

    Parameters
    ----------
        walk_a : A trajectory ending in a fixed-point attractor.
        walk_b : Idem, but on a different attractor.

            note : These two should have at least one shared point


    TODO : allow params [
        branching_point_query,
        attractor_[1,2]_query
    ] to be either strings or dictionnaries
    """
    # deep-copy our frames to avoid overwritting them
    _walk_a = walk_a.copy(deep=True)
    _walk_b = walk_b.copy(deep=True)

    def trajectory_tagger(idx, branching_point, attractor, attractor_tag):
        """Tag a bifurcation process. Use a single branching point as
        reference
            Tags include:
                "common" <=> before the branching point
                "split"  <=> the branching point
                "branch" <=> states after the branching point and before the attractor
                "attractor" <=> the attractor
        """
        if idx < branching_point:
            return f"common_{idx}"
        elif idx == branching_point:
            return f"split_{idx}"
        elif idx == attractor:
            return (
                f"attractor_{idx}_{attractor_tag}"
                if distinguish_attractors
                else f"attractor_{idx}"
            )
        elif idx > branching_point:
            return f"branch_{idx}"
        else:
            raise ValueError(
                f"Undecidable index {idx}, bp: {branching_point}, attr: {attractor}"
            )

    # Tag the first walk's index, according to the branching point and attractor queries
    _walk_a = _walk_a.set_index(
        _walk_a.index.map(
            lambda x: trajectory_tagger(
                idx=x,
                branching_point=_walk_a.query(branching_point_query).index[
                    0
                ],  # safe, only one branching point
                attractor=_walk_a.query(attractor_a_query).index[
                    0
                ],  # safe, only one attractor
                attractor_tag=label_a,
            )
        )
    )

    # Tag the second walk's index, according to the branching point and attractor queries
    _walk_b = _walk_b.set_index(
        _walk_b.index.map(
            lambda x: trajectory_tagger(
                idx=x,
                branching_point=_walk_b.query(branching_point_query).index[
                    0
                ],  # idem as for _walk_b
                attractor=_walk_b.query(attractor_b_query).index[0],
                attractor_tag=label_b,
            )
        )
    )

    f_walk_index_formatter = lambda frame, label: frame.set_index(
        frame.index.astype(str).map(lambda x: f"{label}_{x}")
    )

    _walk_a = f_walk_index_formatter(_walk_a, label_a)
    _walk_b = f_walk_index_formatter(_walk_b, label_b)

    def index_masker(unique_idx, distinguish_attractors: bool = True):
        """Confound uniquely labelled samples within a trajectory."""
        if "common" in unique_idx:
            return "common"
        elif "split" in unique_idx:
            return "split"
        elif "branch" in unique_idx:
            return "branch"
        elif "attractor" in unique_idx:
            return unique_idx if distinguish_attractors else "attractor"
        else:
            raise ValueError(f"Unknown tag `{unique_idx}` found on index")

    _result = pd.concat([_walk_a, _walk_b], axis="rows")

    _partial_index_masker = partial(
        index_masker, distinguish_attractors=distinguish_attractors
    )
    if confound_observations:
        _result = _result.set_index(_result.index.map(_partial_index_masker))

    return _result


def simulate_from_boolean_trajectory(
    boolean_trajectory_df: pd.DataFrame,
    criteria_df: pd.DataFrame,
    n_samples_per_state: IntListOrTuple = 300,
    rng_seed: Optional[int] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Generate `n_samples_per_state`, for each one of the
    states found in `boolean_trajectory_df`.
    The biased simulation from the binary state is performed
    according to `criteria_df`.

    Parameter `n_samples_per_state` is 100%

    If specified, parameter `rng_seed` allows 100% reproductible results,
    which means that given the same set of parameters with the given seed
    will always produce the same pseudo random values.

    The aim of this function is generating synthetic data to be
    used as an input to STREAM, in order to evaluate the performance
    of the PROFILE-based simulation method we have developped.

    Returns
    -------
        A tuple : (simulated_expression_dataframe, metadata)

    """
    # for all runs to obtain the same results the seeds of each run should be fixed
    _rng = np.random.default_rng(rng_seed)
    _simulation_seeds = _rng.integers(
        123, rng_seed, size=len(boolean_trajectory_df.index)
    )

    # TODO : discuss range of sample sizes
    if isinstance(n_samples_per_state, int):
        sample_sizes = [n_samples_per_state] * len(boolean_trajectory_df.index)
    elif isinstance(n_samples_per_state, (list, tuple)):
        sample_sizes = _rng.integers(
            n_samples_per_state[0],
            n_samples_per_state[1],
            size=len(boolean_trajectory_df.index),
        )
    else:
        raise TypeError(
            f"Invalid type `{type(n_samples_per_state)}` for parameter n_samples_per_state"
        )

    # generate synthetic samples
    synthetic_samples = []
    for size, rnd_walk_step, _rng_seed in zip(
        sample_sizes, boolean_trajectory_df.iterrows(), _simulation_seeds
    ):
        _idx, binary_state = rnd_walk_step

        synthetic_samples.append(
            biased_simulation_from_binary_state(
                binary_state.to_frame().T,
                criteria_df,
                n_samples=size,
                seed=_rng_seed,
            )
            .reset_index()
            .rename(columns={"index": "kind"})
        )

    # merge all experiments into a single frame
    synthetic_single_cell_experiment = pd.concat(
        synthetic_samples, axis="rows", ignore_index=True
    )

    # create an informative, artificial, and unique index
    synthetic_single_cell_experiment = synthetic_single_cell_experiment.set_index(
        synthetic_single_cell_experiment.kind
        + "_"
        + synthetic_single_cell_experiment.reset_index().index.map(
            lambda x: f"obs{str(x)}"
        )
    )

    ## TODO : add a parameter to use this section of the code or remove it
    ## create an informative, artificial, and unique index
    # synthetic_single_cell_experiment = synthetic_single_cell_experiment
    #   .set_index(synthetic_single_cell_experiment.kind)

    # Create a colour map for different cell types
    _RGB_values = list("0123456789ABCDEF")
    color_map = {
        i: "#" + "".join([_rng.choice(_RGB_values) for j in range(6)])
        for i in boolean_trajectory_df.index.unique().to_list()
    }

    # Create a metadata frame
    cell_colours = (
        synthetic_single_cell_experiment.kind.apply(lambda x: color_map[x])
        .to_frame()
        .rename(columns={"kind": "label_color"})
    )
    metadata = pd.concat(
        [synthetic_single_cell_experiment.kind, cell_colours], axis="columns"
    )
    metadata = metadata.rename(columns={"kind": "label"})
    # drop the number of activated genes from the synthetic expression frame
    synthetic_single_cell_experiment = synthetic_single_cell_experiment[
        synthetic_single_cell_experiment.columns[1:]
    ]

    return synthetic_single_cell_experiment, metadata
