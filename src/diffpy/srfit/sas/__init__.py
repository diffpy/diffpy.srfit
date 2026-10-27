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
"""SAS calculation tools."""

__all__ = [
    "SASGenerator",
    "SASParser",
    "SASProfile",
    "PrCalculator",
    "CFCalculator",
]

import warnings

from diffpy.utils._deprecator import build_deprecation_message

from .prcalculator import CFCalculator, PrCalculator
from .sasgenerator import SASGenerator
from .sasparser import SASParser
from .sasprofile import SASProfile

removal_version = "4.0.0"

warnings.warn(
    build_deprecation_message(
        "diffpy.srfit",
        "sas",
        "cmisas",
        removal_version,
        new_base="diffpy",
    ),
    DeprecationWarning,
    stacklevel=2,
)

# End of file
