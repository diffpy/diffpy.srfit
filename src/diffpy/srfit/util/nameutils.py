#!/usr/bin/env python
##############################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2008 The Trustees of Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################
"""Name utilities."""

__all__ = ["is_identifier", "validate_name", "isIdentifier", "validateName"]

import re

from diffpy.utils._deprecator import build_deprecation_message, deprecated

reident = re.compile(r"^[a-zA-Z_]\w*$")

nameutils_base = "diffpy.srfit.util.nameutils"
removal_version = "4.0.0"
isIdentifier_dep_msg = build_deprecation_message(
    nameutils_base, "isIdentifier", "is_identifier", removal_version
)
validateName_dep_msg = build_deprecation_message(
    nameutils_base, "validateName", "validate_name", removal_version
)


def is_identifier(name):
    """Check whether a string is a valid name.

    A valid name starts with a letter or underscore and contains only
    letters, digits and underscores.

    Parameters
    ----------
    name : str
        The string to check.

    Returns
    -------
    bool
        The flag that is True if ``name`` is a valid name.
    """
    return reident.match(name) is not None


def validate_name(name):
    """Raise an error if a string is not a valid name.

    Parameters
    ----------
    name : str
        The string to check.

    Raises
    ------
    ValueError
        If ``name`` is not a valid name, as defined by ``is_identifier``.
    """
    if not is_identifier(name):
        raise ValueError(
            f"Name '{name}' is not a valid identifier. "
            "Please use only letters, digits and underscores, "
            "and do not start the name with a digit."
        )
    return


@deprecated(isIdentifier_dep_msg)
def isIdentifier(s):
    """This function has been deprecated and will be removed in version
    4.0.0.

    Please use diffpy.srfit.util.nameutils.is_identifier instead.
    """
    return is_identifier(s)


@deprecated(validateName_dep_msg)
def validateName(name):
    """This function has been deprecated and will be removed in version
    4.0.0.

    Please use diffpy.srfit.util.nameutils.validate_name instead.
    """
    return validate_name(name)


# End of file
