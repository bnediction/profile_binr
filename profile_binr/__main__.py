"""
Main module, to define a simple command line interface.
"""

# TODO : export parser logic to CLI module ?

import argparse
import warnings
import multiprocessing as mp
from pathlib import Path
from _io import TextIOWrapper

import toml
import pandas as pd
from plotnine import ggplot, geom_point, aes, stat_smooth, facet_wrap

from .wrappers.probinr import ProfileBin
from .simulation import random_nan_binariser
from .utils.diagnostics import (
    expression_profile_scatterplot,
    compare_profiles,
    summarise_by_criteria,
)

__ACTIONS = [
    "binarize",
    "py_binarize",
    "simulate",
]
__HELP_MESSAGES = {
    "main": """Process a """,
    "in_file": """A csv file containing a column for each gene and a line for each observation.
            Expression data must be normalized before using this tool.""",
    "publish_dir": """Where to store the results of applying 'action'(s) to 'in_file'.""",
    "config_file": """A TOML configuration file, containing all the parameters here described.
        If specified, the contents of THE CONFIG FILE WILL OVERRIDE ALL THE OTHER PARAMETERS
        PASSED VIA THE COMMAND LINE.""",
    "n_threads": "The number of threads for parallel execution (bounded by available cores)",
    "mask_zeros": "Should zero entries be masked when computing the criteria ?",
    "action": "The action to be performed.",
}

__DEFAULT_VALUES = {
    "action": "binarize",
    "n_threads": mp.cpu_count(),
    "mask_zeros": True,
}


def main():
    """ Command line interface to rna2bool"""
    # ARG PARSING
    try:
        # IMPORTANT :
        # The PROG argument must match the name in the pyproject.toml file.
        # IMPORTANT :
        # help messages are defined within the top-level dict __HELP_MESSAGE
        parser = argparse.ArgumentParser(
            prog="rna2bool", description=__HELP_MESSAGES["main"]
        )
        parser.add_argument(
            "-a",
            "--action",
            choices=__ACTIONS,
            nargs="+",
            dest="actions",
            help=__HELP_MESSAGES["action"],
            default=__DEFAULT_VALUES["action"],
        )
        parser.add_argument(
            "-f",
            "--expression_file",
            # type=argparse.FileType("r"),
            type=lambda p: Path(p).resolve(),
            help=__HELP_MESSAGES["in_file"],
            dest="in_file",
        )
        parser.add_argument(
            "-p",
            "--publish_dir",
            type=lambda p: Path(p).resolve(),
            help=__HELP_MESSAGES["publish_dir"],
            dest="publish_dir",
            default=Path(".").resolve(),
        )
        parser.add_argument(
            "-c",
            "--config_file",
            type=argparse.FileType("r"),
            help=__HELP_MESSAGES["config_file"],
        )
        parser.add_argument(
            "-n",
            "--n_threads",
            type=int,
            default=__DEFAULT_VALUES["n_threads"],
            dest="n_threads",
            help=__HELP_MESSAGES["n_threads"],
        )
        parser.add_argument(
            "-m",
            "--mask_zeros",
            type=bool,
            default=__DEFAULT_VALUES["mask_zeros"],
            help=__HELP_MESSAGES["mask_zeros"],
        )
        # parser.add_argument(
        #    "--pandas_kw",
        #    dest="the_dict",
        #    action=StoreDictKeyPair,
        #    nargs="+",
        #    metavar="KEY=VAL",
        # )

        args = parser.parse_args()
        if not any((bool(args.in_file), args.config_file)):
            parser.print_help()
            raise SystemExit() from None

        if args.config_file:
            _config = toml.load(args.config_file)
            args.config_file.close()

            # check config file integrity
            # change .issubset() to intersection()
            # TODO : check this condition, as it might be insufficient in the future
            # i.e., what if they define it as none or something ?
            # Should all these integrity checks be performed ?
            if not set(_config.keys()).intersection(set(args.__dict__.keys())):
                # extra_args = set(_config.keys()).difference(set(args.__dict__.keys()))
                parser.print_help()
                print(
                    f"config file {args.config_file.name} has: {_config.keys()}",
                    f"params expected: {args.__dict__.keys()}",
                    "Aborting",
                    sep="\n",
                )
                raise SystemExit() from None
            # override all other params
            args.__dict__.update(_config)
            # If more filters are needed consider iterating over args.__dict__
            if isinstance(args.publish_dir, str):
                args.publish_dir = Path(args.publish_dir).resolve()
            if isinstance(args.in_file, str):
                args.in_file = Path(args.in_file).resolve()
            elif isinstance(args.in_file, TextIOWrapper):
                args.in_file = Path(args.in_file.name).resolve()
            else:
                raise TypeError(
                    f"Unexpected type for args.in_file {type(args.in_file)}"
                )

        args.n_threads = min(abs(args.n_threads), mp.cpu_count())
        try:
            args.publish_dir.mkdir(parents=True, exist_ok=True)
        except FileExistsError as _f_e_err:
            parser.print_help()
            print("\nERROR : provided `publish_dir` exists and is not a directory")
            raise SystemExit from None
    except (argparse.ArgumentError, ValueError) as _arg_err:
        parser.print_help()
        print(_arg_err)
        raise SystemExit() from None
    # END ARG PARSING

    # CLI LOGIC
    print("Flags:")
    for k, v in sorted(vars(args).items()):
        print(f"\t{k}: {v}")
    results = {}
    print("import csv...", end="\t")
    in_counts = pd.read_csv(args.in_file, index_col=0)
    print("Done.")
    print("create instance and call fit...", end="\t")
    probin_instance = ProfileBin(in_counts)
    probin_instance.fit(n_threads=args.n_threads, mask_zero_entries=args.mask_zeros)
    _criteria_file = (
        (args.publish_dir / f"{args.in_file.name.replace('.csv', '')}_criteria.csv")
        .resolve()
        .as_posix()
    )
    probin_instance.criteria.to_csv(_criteria_file)
    print("Done.")
    if "simulate" in args.actions:
        print("call simulation fit...", end="\t")
        probin_instance.simulation_fit(
            n_threads=args.n_threads, mask_zero_entries=args.mask_zeros
        )
        _simulation_criteria_file = _criteria_file.replace(
            "_criteria.csv", "_simulation_criteria.csv"
        )
        probin_instance.simulation_criteria.to_csv(_simulation_criteria_file)
        print("Done.")

    # use a results dict !
    # iterate over the different instructions given
    # check dependencies ?
    action_calls = {  # Use lambdas for lazy evaluation
        "binarize": lambda: probin_instance.binarize(),
        "py_binarize": lambda: probin_instance.py_binarize(),
        "simulate": lambda x: probin_instance.simulate(random_nan_binariser(x)),
    }
    # Sort the actions according to the specified order
    print(f"args.actions pre sort : {args.actions}")
    args.actions.sort(key=__ACTIONS.index)
    print(f"args.actions post sort : {args.actions}")
    results = {}
    for action in args.actions:
        print(f"action: {action}...", end="\t")
        if action != "simulate":
            results[action] = action_calls[action]()
        else:
            _bin = results.get("binarize")
            _bin = results.get("py_binarize") if _bin is None else _bin
            _bin = probin_instance.py_binarize() if _bin is None else _bin
            results[action] = action_calls[action](_bin)
            _summary = compare_profiles(
                probin_instance.criteria,
                in_counts,
                results[action],
                "Original",
                "Simulated",
            )
            # TODO : move this into diagnostics ?
            plot_1 = (
                ggplot(_summary, aes("Mean", "Std"))
                + geom_point(aes(color="Category", fill="Category"), alpha=0.55)
                + facet_wrap("~Data")
            )
            plot_2 = (
                ggplot(_summary, aes("Mean", "Std"))
                + geom_point(aes(color="Data", fill="Data"), alpha=0.55)
                + facet_wrap("~Category")
            )
            plot_1_out_file = (
                (
                    args.publish_dir
                    / f"{args.in_file.name.replace('.csv', '')}_{action}d_by_Data.png"
                )
                .resolve()
                .as_posix()
            )
            plot_2_out_file = plot_1_out_file.replace(
                "_by_Data.png", "_by_Category.png"
            )
            # print(plot_1_out_file)
            # print(plot_2_out_file)
            warnings.filterwarnings("ignore")
            plot_1.save(plot_1_out_file, height=25, width=40, units="cm", verbose=False)
            plot_2.save(plot_2_out_file, height=25, width=40, units="cm", verbose=False)
            warnings.filterwarnings("default")
        # print(args.in_file.name.replace(".csv", ""))
        results[action].to_csv(
            args.publish_dir / f"{args.in_file.name.replace('.csv', '')}_{action}d.csv"
        )
        print("Done.")

    # END CLI LOGIC

    # return vars()
    return 0


if __name__ == "__main__":
    # global main_env
    # main_env = main()
    # [globals().update({k: main_env[k]}) for k in main_env.keys()]
    main()
