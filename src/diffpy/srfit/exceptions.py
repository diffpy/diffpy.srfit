#!/usr/bin/env python
##############################################################################
#
# diffpy.srfit      Complex Modeling Initiative
#                   (c) 2016 Brookhaven Science Associates,
#                   Brookhaven National Laboratory.
#                   All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
##############################################################################
"""
Exceptions used for SrFit - specific errors.
"""


class SrFitError(Exception):
    """Generic error in SrFit expressions or recipe."""

    pass


class ParseError(Exception):
    """Exception used by ProfileParsers."""

    pass
