"""
    ProfileBin: The binarization method of RNA expression
    as defined in the paper:

    "Personalization of Logical Models With Multi-Omics Data Allows
    Clinical Stratification of Patients"

    by authors :

    Beal, Jonas
    Montagud, Arnau
    Traynard, Pauline
    Barillot, Emmanuel
    Calzone, Laurence

    at the Computational Systems Biology of Cancer group at Institut Curie
    website : https://sysbio.curie.fr/
      email :  contact-sysbio@curie.fr

    The repository containing the original implementation in
    Rmarkdown notebooks : https://github.com/sysbio-curie/PROFILE

    Python wrappers and parallel implementation of the
    R function compute_criteria() written by Gustavo Magaña López. 
    GitHub profile : https://github.com/gmagannaDevelop
             email : gustavo.magana-lopez@u-psud.fr
"""

__all__ = ["ProfileBin"]

from typing import NoReturn, Any, Callable, List, Dict, Tuple, Optional, Union
from pathlib import Path
import logging

import random
import string

import multiprocessing

# rpy2 conversion
import rpy2.robjects as r_objs
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import globalenv as GLOBALENV
from rpy2.rinterface_lib.embedded import RRuntimeError
from rpy2.rinterface_lib.callbacks import logger as rpy2_logger

# data management
import numpy as np
import pandas as pd

# R source code locations :
__PROBINR_DIR__ = Path(__file__).parent.absolute().parent.absolute().joinpath("_R")
__PROBINR_SRC__ = __PROBINR_DIR__.joinpath("PROFILE_source.R").absolute()
__PROBINR_BOOTSTRAP__ = __PROBINR_DIR__.joinpath("install_deps.R").absolute()

rpy2_logger.setLevel(logging.ERROR)  # will display errors, but not warnings


class ProfileBin(object):
    """
    Class used to perform the binarization of RNA expression data,
    according to the methodology presented by Beal Jonas et al.
    doi = {10.3389/fphys.2018.01965}.

    Objects of this class should be instantiated by passing a
    pandas.DataFrame which contains samples as rows and genes
    as columns.

    Here is an example of the intended usage of this class.

    import the data :
    >>> exp_data = pd.read_csv("example_gene_expression.csv")

    NOTE :
        The dataframe should contain the GENES as COLUMNS
        and SAMPLES as ROWS. The individual (sample) identifier
        should be the index.

    create a ProfileBin instance
    >>> probin = ProfileBin(exp_data)

    compute "criteria" used to determine the proper binarization rule
    for each gene in the dataset
    >>> probin.fit()

    both methods "binarize" and "normalize" are now available:
    >>> probinr.binarize() # by default will binarized the expression
                           # data used to compute the criteria.

    >>> probinr.binarize(new_data) # binarize new observations
                                   # of the same genes (or a subset).
    """

    @staticmethod
    def _r_bool(value: bool) -> str:
        """Return a string representation of a Python boolean
        which is valid in R syntax. To be used when generating
        string literals which can be directly passed to self.r
        (an R instance invoked by rpy2 at initialization)."""
        if not isinstance(value, bool):
            raise TypeError(f"R boolean representation of {type(value)} undefined")
        return "T" if value else "F"

    @staticmethod
    def _r_none() -> str:
        """Return "NULL", R's equivalent of Python's None
        This is necessary as rpy2 provides no conversion rule for Python's None."""
        return "NULL"

    @staticmethod
    def _random_string(n: int = 10) -> str:
        """ Return a random string og length `n` containing ascii lowercase letters and digits. """
        return "".join(
            random.choice(string.ascii_lowercase + string.digits) for i in range(n)
        )

    def __init__(self, data: pd.DataFrame):
        # self.__addr will be used to keep track of R objects related to the instance :
        self.__addr: str = str(hex(id(self)))
        self.r = r_objs.r
        self.r_globalenv = GLOBALENV
        self._criteria: pd.DataFrame
        self._zero_inf_criteria: pd.DataFrame
        self._zero_inf_idx: pd.core.indexes.base.Index
        self._zero_inf_df: pd.DataFrame
        # try loading all packages and functions, installing them upon failure
        try:
            with open(__PROBINR_SRC__, "r") as f:
                self.r("".join(f.readlines()))
        except RRuntimeError:
            print("\nERROR : one or more R dependencies are not installed")
            print("Trying to automatically satisfy missing dependencies\n")
            try:
                # install dependencies :
                with open(__PROBINR_BOOTSTRAP__, "r") as f:
                    self.r("".join(f.readlines()))
                print("\n Missing dependencies successfully installed \n")
                # re-import the R source as functions were not saved because
                # of the previous RRuntimeError
                with open(__PROBINR_SRC__, "r") as f:
                    self.r("".join(f.readlines()))
            except RRuntimeError as _rer:
                print("Bootstrapping the installation of R dependencies failed:")
                raise _rer from None

        # sanitise inputs :
        if not isinstance(data, pd.DataFrame):
            raise TypeError(
                f"Parameter 'data' must be of type 'pandas.DataFrame' not {type(data)}"
            )
        else:
            self._data: pd.DataFrame = data

    def __repr__(self):
        _lines = [
            f"PROFILE instance at {self.__addr}",
            f"with R instance at {self.r.__repr__().split(' ')[-1][:-1]}",
            f"is_trained : {self._is_trained}",
        ]
        return "\n".join(_lines)

    def r_ls(self):
        """ Return a list containing all the names in the main R environment. """
        return list(self.r("ls()"))

    @property
    def _is_trained(self) -> bool:
        """ Determine if the criteria has been calculated by inspecting the R environment. """
        return f"criteria_{self.__addr}" in self.r_ls()

    @property
    def _data_in_r(self) -> bool:
        """ Determine if the data to fit the criteria is present in the R environment. """
        return f"META_RNA_{self.__addr}" in self.r_ls()

    def r_instantiate_data(
        self, data: Union[pd.DataFrame, Any], identifier: str
    ) -> NoReturn:
        """
        Instantiate a DataFrame within the R embedded process,
        bind it to the provided identifier on the main environment.

        Parameters:
                  data: Any data, although the original converter is intended for pandas DataFrames
            identifier: A string, the desired variable name for the object in the R env.
        """
        try:
            with localconverter(r_objs.default_converter + pandas2ri.converter):
                self.r_globalenv[identifier] = r_objs.conversion.py2rpy(data)
        except RRuntimeError as _rer:
            print(
                "An error ocurred while instantiating the data within the R session",
                "This is likely due to the type of the objects contained in the DataFrame",
                "(rpy2 may not have implemented the needed conversion rule).",
                "",
                sep="\n",
            )
            raise RRuntimeError(str(_rer)) from None
            # change this to display R's original error ?

    # TODO: verify the function's docstring
    def fit(
        self,
        n_threads: Optional[int] = multiprocessing.cpu_count(),
        dor_threshold: Optional[float] = 0.95,
        mask_zero_entries: Optional[bool] = False,
    ) -> NoReturn:
        """
        Compute the criteria needed to decide which binarization rule
        will be applied to each gene. This is performed by calling
        the corresponding R function `compute_criteria()` via rpy2.

        Arguments :
            n_threads : The number of parallel processes (threads) to be used.
             log_time : Should the call to the R function be timed ?
              logfile : The file in which logs should be saved.

        Returns :
            None

        "Side effects" : a pandas.DataFrame containing the criteria and the label
        for each gene will be stored within the class.

        It is accessible via :
        >>> probin_instance.criteria

        criteria's columns are :
            Dip, BI, Kurtosis, DropOutRate, MeanNZ, DenPeak, Amplitude, Category

        IMPORTANT : `compute_criteria()` uses descriptors and memory
        mappings, which are temporarily stored in files containing the
        current datetime and a small hash to ensure the unicity of names.

        Do not manipulate or erase files which look like this:
            Wed Apr 14 13:47:21 2021 MFNVP2486W.bin
            Wed Apr 14 13:47:21 2021 MFNVP2486W.desc

        They will be automatically removed one the criteria is calculated
        for all genes in the dataset.

        """

        # the data must be instantiated before performing the R call :
        if not self._data_in_r:
            self.r_instantiate_data(self.data, f"META_RNA_{self.__addr}")

        # call compute_criteria only once
        if not self._is_trained:
            params = [
                f"exp_dataset = META_RNA_{self.__addr}",
                f"n_threads = {n_threads}",
                f"dor_threshold = {dor_threshold}",
                f"mask_zero_entries = {self._r_bool(mask_zero_entries)}",
            ]
            try:
                with localconverter(r_objs.default_converter + pandas2ri.converter):
                    self._criteria = r_objs.conversion.rpy2py(
                        self.r(
                            f"criteria_{self.__addr} <- compute_criteria({', '.join(params)})"
                        )
                    )
            except RRuntimeError as _rer:
                _err_ls = (
                    "",
                    str(_rer),
                    "There was an error while calling compute_criteria() in R",
                    "Some likely causes are:",
                    "\t * insufficient RAM",
                    "\t * the descriptor files might have been deleted or corrupted",
                    "\t * the data.frame contains non-numerical entries",
                )
                raise RRuntimeError("\n".join(_err_ls)) from None

    def simulation_fit(
        self,
        n_threads: Optional[int] = multiprocessing.cpu_count(),
        dor_threshold: Optional[float] = 0.75,
        mask_zero_entries: Optional[bool] = True,
    ) -> NoReturn:
        """Re compute criteria for genes classified as zero-inflated,
        in order to better estimate simulation parameters."""
        if not self._is_trained:
            raise ValueError(
                "\n".join(
                    [
                        "Cannot compute simulation fit because self.criteria does not exist.",
                        "Call self.fit() before calling this method.",
                    ]
                )
            )

        self._zero_inf_idx = self._criteria[self._criteria.Category == "ZeroInf"].index
        self._zero_inf_df = self._data.loc[:, self._zero_inf_idx]
        self.r_instantiate_data(self._zero_inf_df, f"zero_inf_RNA_{self.__addr}")

        params = [
            f"exp_dataset = zero_inf_RNA_{self.__addr}",
            f"n_threads = {n_threads}",
            f"dor_threshold = {dor_threshold}",
            f"mask_zero_entries = {self._r_bool(mask_zero_entries)}",
        ]
        try:
            with localconverter(r_objs.default_converter + pandas2ri.converter):
                self._zero_inf_criteria = r_objs.conversion.rpy2py(
                    self.r(
                        f"zero_inf_criteria_{self.__addr} <- compute_criteria({', '.join(params)})"
                    )
                )

            self._zero_inf_criteria.loc[
                self._zero_inf_criteria.Category == "ZeroInf", "Category"
            ] = "Unimodal"
        except RRuntimeError as _rer:
            _err_ls = (
                "",
                str(_rer),
                "There was an error while calling compute_criteria() in R",
                "Some likely causes are:",
                "\t * insufficient RAM",
                "\t * the descriptor files might have been deleted or corrupted",
                "\t * the data.frame contains non-numerical entries",
            )
            raise RRuntimeError("\n".join(_err_ls)) from None

    def _binarize_or_normalize(
        self,
        action: str,
        data: Optional[pd.DataFrame] = None,
        gene: Optional[str] = None,
    ) -> pd.DataFrame:
        """
        Not intended to be called directly.
        Handle calls to self.binarize and self.normalize by dynamically constructing the R call.
        """

        assert action in [
            "binarize",
            "normalize",
        ], f"No method defined for action {action}"

        if not self._is_trained:
            raise AttributeError(
                f"Cannot {action} without the criteria DataFrame. Call self.fit() first."
            )

        _df_name: str = ""
        _rm_df: bool = False
        if data is not None:
            _df_name = f"bin_{self.__addr}_{self._random_string()}"
            self.r_instantiate_data(data, _df_name)
            _rm_df = True
        else:
            _df_name = f"META_RNA_{self.__addr}"
            data = self.data

        params = [
            f"exp_dataset = {_df_name}",
            f"ref_dataset = META_RNA_{self.__addr}",
            f"ref_criteria = criteria_{self.__addr}",
        ]

        if gene is not None:
            params.append(f"gene = {gene}")

        try:
            with localconverter(r_objs.default_converter + pandas2ri.converter):
                _df: pd.DataFrame = r_objs.conversion.rpy2py(
                    self.r(
                        f"""
                            {action}_exp({', '.join(params)}) %>%
                                apply(2, function(x) tidyr::replace_na(x, replace = -1)) %>%
                                    as.data.frame()
                    """
                    )
                )
                return _df.replace(-1.0, np.nan)
        except RRuntimeError as _rer:
            _err_help_ls = [
                "",
                f"{_rer}\n",
                "Please verify that your genes (columns) are a subset",
                "of the genes used to compute the criteria.\n",
                "Difference between the dataframes' column sets",
                "(genes contained in data for which there is no established criteria):",
                f"{set(data.columns).difference(set(self.data.columns))}",
            ]
            raise RRuntimeError("\n".join(_err_help_ls)) from None
        finally:
            if _rm_df:
                _ = self.r(f"rm({_df_name})")

    def binarize(
        self, data: Optional[pd.DataFrame] = None, gene: Optional[str] = None
    ) -> pd.DataFrame:
        """
        Binarize expression data.

        Parameters:
            data: a (optional) pandas.DataFrame containing genes as columns and rows as measurements
                  Defaults to self.data (the data used to compute the criteria).

            gene: an (optional) string determining the gene name to binarize. It must be contained
                  the dataframe to binarize (and the previously computed criteria).
        """
        return self._binarize_or_normalize("binarize", data, gene)

    def binarise(self, *args, **kwargs) -> pd.DataFrame:
        """ alias for self.binarize. See help(ProfileBin.binarize) """
        return self.binarize(*args, **kwargs)

    def normalise(self, *args, **kwargs) -> pd.DataFrame:
        """ alias for self.normalize. See help(ProfileBin.normalize) """
        return self.normalize(*args, **kwargs)

    def normalize(
        self, data: Optional[pd.DataFrame] = None, gene: Optional[str] = None
    ) -> pd.DataFrame:
        """
        Normalize expression data.

        Parameters:
            data: a (optional) pandas.DataFrame containing genes as columns and rows as measurements
                  Defaults to self.data (the data used to compute the criteria).

            gene: an (optional) string determining the gene name to normalize. It must be contained
                  the dataframe to binarize (and the previously computed criteria).
        """
        return self._binarize_or_normalize("normalize", data, gene)

    def plot_zeroinf_diagnostics(
        self,
        df_plot_kwargs: Dict[str, Any] = {
            "kind": "scatter",
            "x": "DropOutRate",
            "y": "zero_inf_thresh",
        },
    ):
        """Plot the ZeroInflated threshold as a function of the dropout rate
        in order to visualise the DropOutRate's value to use as new discarding
        threshold."""
        _fig = self.criteria.plot(**df_plot_kwargs)
        return _fig

    @property
    def data(self) -> pd.DataFrame:
        """The expression data used to compute the criteria.
        WARNING : this returns a reference to the DataFrame and not
        a copy. If you modify the expression data before calling
        self.fit() the binarization might be biased."""
        return self._data

    @data.deleter
    def data(self):
        raise AttributeError(
            "ProfileBin wrapper cannot operate without attribute 'data'. Aborting deletion."
        )

    @property
    def criteria(self) -> pd.DataFrame:
        """ Computed criteria to choose the binarization algorithm """
        if hasattr(self, "_criteria"):
            return self._criteria
        else:
            raise AttributeError(
                "'criteria' has not been calculated. Call self.fit() to define it"
            )

    @property
    def criteria_zero_inf(self) -> pd.DataFrame:
        """ Computed criteria to simulate synthetic data """
        if hasattr(self, "_zero_inf_criteria"):
            return self._zero_inf_criteria
        else:
            raise AttributeError(
                "'criteria_zero_inf' has not been calculated. Call self.simulation_fit() to define it"
            )

    @criteria.deleter
    def criteria(self):
        raise AttributeError(
            "Cannot delete 'criteria' as it is necessary to perform the binarization, aborting."
        )

    def clear_r_envir(self):
        """Remove an all R objects that have been created by the ProfileBin
        instance."""
        _instance_objs = [f"'{obj}'" for obj in self.r_ls() if self.__addr in obj]
        _ = self.r(f"rm(list = c({', '.join(_instance_objs)}))")
