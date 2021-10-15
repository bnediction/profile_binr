"""
    Tools for assessing the results of ProfileBin's simulation data.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def summarise_by_criteria(data: pd.DataFrame, criteria: pd.DataFrame) -> pd.DataFrame:
    """ Summarise by criteria to get the needed mean vs std plots """
    return pd.DataFrame(
        {"mean": data.mean(), "std": data.std(), "category": criteria.Category}
    )


def expression_profile_scatterplot(
    name: str, frame: pd.DataFrame, group_col: str = "category"
):
    """plot the expression profile for a frame
    the frame should be the result of calling

    >>> frame = summarise_by_criteria(data, criteria)
    """
    sns.lmplot(
        x="mean",
        y="std",
        hue=group_col,
        data=frame.sort_values(group_col),
        fit_reg=False,
        scatter_kws={"alpha": 0.4},
    )
    plt.title(name)
    plt.ion()
    plt.show()
