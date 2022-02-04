"""
Main module, to define a simple command line interface.
"""
import toml
import argparse
import multiprocessing as mp
from .wrappers.probinr import ProfileBin

__HELP_MESSAGES = {
    "main": """Process a """,
    "in_file": """A csv file containing a column for each gene and a line for each observation.
            Expression data must be normalized before using this tool.""",
    "config_file": """A TOML configuration file, containing all the parameters here described.
        If specified, the contents of THE CONFIG FILE WILL OVERRIDE ALL THE OTHER PARAMETERS
        PASSED VIA THE COMMAND LINE.""",
    "n_cores": "The number of threads for parallel execution (bounded by available cores)",
    "mask_zeros": "Should zero entries be masked when computing the criteria ?",
}

# TODO : remove this class
class DictNamespace:
    """Just like argparse.Namespace, but it can be used as a dict"""

    def __init__(self, namespace):
        # only copy instance attributes from parents
        # and make a deepcopy to avoid unwanted side-effects
        # EDIT : TypeError: cannot pickle '_io.TextIOWrapper' object
        for k, v in namespace.__dict__.items():
            setattr(self, k, v)

    def __repr__(self):
        __args = ", ".join([f"{k}={v}" for k, v in self.items()])
        return f"DictNamespace({__args})"

    def __getitem__(self, name):
        return getattr(self, name)

    def __setitem__(self, key, new_value):
        setattr(self, key, new_value)

    def __iter__(self):
        yield from self.__dict__

    def keys(self):
        return self.__dict__.keys()

    def values(self):
        return self.__dict__.values()

    def items(self):
        return self.__dict__.items()


class StoreDictKeyPair(argparse.Action):
    """Created this argparse action to save kwargs
    to a dict, to be passed to pandas.read_csv()
    this functionality will be developed in the future.
    """

    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        self._nargs = nargs
        super(StoreDictKeyPair, self).__init__(
            option_strings, dest, nargs=nargs, **kwargs
        )

    def __call__(self, parser, namespace, values, option_string=None):
        my_dict = {}
        for kv in values:
            k, v = kv.split("=")
            my_dict[k] = v
        setattr(namespace, self.dest, my_dict)


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
            "-f",
            "--expression_file",
            type=argparse.FileType("r"),
            help=__HELP_MESSAGES["in_file"],
            dest="in_file",
        )
        parser.add_argument(
            "--config_file",
            type=argparse.FileType("r"),
            help=__HELP_MESSAGES["config_file"],
        )
        parser.add_argument(
            "--n_cores",
            type=int,
            default=mp.cpu_count(),
            dest="n_cores",
            help=__HELP_MESSAGES["n_cores"],
        )
        parser.add_argument(
            "--mask_zeros", type=bool, default=True, help=__HELP_MESSAGES["mask_zeros"]
        )
        # parser.add_argument(
        #    "--pandas_kw",
        #    dest="the_dict",
        #    action=StoreDictKeyPair,
        #    nargs="+",
        #    metavar="KEY=VAL",
        # )

        args = DictNamespace(parser.parse_args())
        if not any((args.in_file, args.config_file)):
            parser.print_help()
            raise SystemExit() from None

        if args.config_file:
            _config = toml.load(args.config_file)

            if not set(args.keys()) == set(_config.keys()):
                raise KeyError(
                    " ".join(
                        [
                            f"Config file {args.config_file.name} does not have all necessary parameters:",
                            f"{list(args.keys())}",
                        ]
                    )
                ) from None

            for param_name in _config:
                if param_name != "config_file":
                    _new_param_value = _config.get(param_name)
                    if _new_param_value is None:
                        pass
                    args[param_name] = _new_param_value

        args.n_cores = min(abs(args.n_cores), mp.cpu_count())
    except (argparse.ArgumentError, ValueError) as _arg_err:
        raise _arg_err from ValueError("no me la creo")
    # END ARG PARSING

    # CLI LOGIC
    print("Flags:")
    for k, v in sorted(vars(args).items()):
        print(f"\t{k}: {v}")

    # END CLI LOGIC
    return args, parser


if __name__ == "__main__":
    global args
    global parser
    args, parser = main()
