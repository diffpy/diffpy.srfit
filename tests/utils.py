#!/usr/bin/env python
"""Utilities for testing."""

import diffpy.srfit.equation.literals as literals

def _makeArgs(num):
    args = []
    for i in xrange(num):
        j=i+1
        args.append(literals.Argument(name="v%i"%j, value=j))
    return args
