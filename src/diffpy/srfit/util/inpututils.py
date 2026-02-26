#!/usr/bin/env python
##############################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2010 The Trustees of Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################
"""Input utilities."""

__all__ = ["inputToString"]

import os.path
from pathlib import Path


def inputToString(input):
    """Convert input from various modes to a string.

    This is useful when you want a method to accept a string, open file object
    or file name.

    Attributes
    ----------
    input
        An open file-like object, name of a file
        or a string containing the input.


    Returns the input in a string
    Raises IOError if the input is supected to be a file name, but the file
    cannot be found.
    """
    # Get the input into a string
    inptstr = ""
    if hasattr(input, "read"):
        inptstr = input.read()
    # TODO remove handling of string input accept only file or filename
    # FIXME check for typos in the file name
    elif os.path.exists(input) or (len(input) < 80 and input.count("\n") == 0):
        with open(input, "r") as infile:
            inptstr = infile.read()
    else:
        inptstr = input

    return inptstr


def get_dict_from_results_file(
    results_filepath: Path | str,
) -> dict[str, float]:
    """Get a dictionary of parameter names and values from a results
    file.

    The file should have lines in the format:
    "parameter_name value +/- uncertainty". Lines that do not match this
    format will be ignored.

    Parameters
    ----------
    results_filepath : pathlib.Path or str
        The path to the results file.

    Returns
    -------
    parsed_results_dict : dict
        The dictionary where keys are parameter names and values are the
        corresponding parameter values as floats.
    """
    with open(results_filepath, "r") as f:
        results_string = f.read()
    parsed_results_dict = {}
    for raw_line in results_string.splitlines():
        line = raw_line.strip()
        # skip blank lines and lines that are just dashes
        if not line or set(line) == {"-"}:
            continue
        line_items = line.split()
        if len(line_items) < 2:
            continue
        if len(line_items) >= 4 and line_items[2] == "+/-":
            try:
                parsed_results_dict[line_items[0]] = float(line_items[1])
            except ValueError:
                pass
    return parsed_results_dict


# End of file
