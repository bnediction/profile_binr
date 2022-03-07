"""
    Docstring
"""

import multiprocessing as mp
from typing import Optional, Union

# TODO : change import scheme
from .sampling import _dropout_mask, pd, np, ss, _GLOBAL_RNG

_RandType = Union[np.random.Generator, int]

# for m in _MODULE_RNGS:
#    print(m.__getstate__()["state"]["state"])
## TODO : define this function
# def mean_binariser(binary_df, normalized_counts_df):
#    pass