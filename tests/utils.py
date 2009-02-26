#!/usr/bin/env python
"""Utilities for testing."""

import diffpy.srfit.equation.literals as literals

def _makeArgs(num):
    args = []
    for i in xrange(num):
        args.append(literals.Argument())
        args[-1].name = "v%i"%(i+1)
    return args
