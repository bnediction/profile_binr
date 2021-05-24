"""
    development utility script intended to be run
    from the project's root in order to update the 
    string literal containing the R code 
    for the PROFILE binarisation pipeline.

"""

from pathlib import Path


if __name__ == "__main__":
    
    try: 
        r_source = Path("R/PROFILE_source.R")
        py_destination = Path("profile_binr/core/probinr.py")

        with open(r_source.absolute().as_posix(), "r") as f:
            print("Parsing R source...", end="\t")
            txt = "".join(f.readlines())
        print("Done")

        with open(py_destination.absolute().as_posix(), "w") as f:
            print("Writting to Python... ", end="\t")
            f.write(f"PROBINR_TXT_LITERAL = \"\"\"{txt}\"\"\"\n")
        print("Done")

    except FileNotFoundError as _fnfe:
        print("Please verify that the `r_source` and `py_destination` variables")
        print("are set to existing paths in your filesystem")
        print("Raised Exception :")
        print(_fnfe)
        exit()

