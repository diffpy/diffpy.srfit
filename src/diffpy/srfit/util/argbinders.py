#!/usr/bin/env python
##############################################################################
#
# diffpy.srfit      Complex Modeling Initiative
#                   (c) 2014 Brookhaven Science Associates,
#                   Brookhaven National Laboratory.
#                   All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
##############################################################################
"""Functions for binding arguments of callable objects."""


class bind2nd(object):
    """Freeze second argument of a callable object to a given constant."""

    def __init__(self, func, arg1):
        """Freeze the second argument of function func to arg1."""
        self.func = func
        self.arg1 = arg1
        return

    def __call__(self, *args, **kwargs):
        boundargs = (args[0], self.arg1) + args[1:]
        return self.func(*boundargs, **kwargs)
