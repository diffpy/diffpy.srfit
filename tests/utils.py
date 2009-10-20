#!/usr/bin/env python
"""Utilities for testing."""

import diffpy.srfit.equation.literals as literals

def _makeArgs(num):
    args = []
    for i in xrange(num):
        j=i+1
        args.append(literals.Argument(value=j))
    return args

def _makeNodes(num):
    nodes = []
    ns = {}
    for i in xrange(num):
        j=i+1
        name = "v%i"%j
        ns[name] = literals.Argument(value=j)
        nodes.append(literals.Node(name, ns))
    return nodes

