#!/usr/bin/env python
########################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2008 Trustees of the Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
########################################################################
"""Utilities for creating Equation objects."""

#FIXME - Need to rename constants in the evaluated equation
#TODO - Need to handle assignments and a shared namespace between

from .Equation import Equation
import diffpy.srfit.equation.literals as literals
import numpy

class LiteralBuilder(object):
    """Class for building a literal from an equation."""

    def __init__(self):
        self.literal = None
        return

    def __eval(self, other, OperatorClass):
        op = OperatorClass()
        op.addLiteral(self.literal)
        op.addLiteral(other.literal)
        lb = LiteralBuilder()
        lb.literal = op
        return lb

    def __add__(self, other):
        return self.__eval(other, literals.AdditionOperator)

    def __sub__(self, other):
        return self.__eval(other, literals.SubtractionOperator)

    def __mul__(self, other):
        return self.__eval(other, literals.MultiplicationOperator)

    def __div__(self, other):
        return self.__eval(other, literals.DivisionOperator)

    def __pow__(self, other):
        return self.__eval(other, literals.ExponentiationOperator)

    def __mod__(self, other):
        return self.__eval(other, literals.RemainderOperator)

class ArgumentBuilder(LiteralBuilder):

    def __init__(self, value=None, name=None):
        self.literal = literals.Argument(value=value, name=name)
        return

class OperatorBuilder(object):
    """Acts like a numpy operator, but helps build a Literal tree."""

    def __init__(self, name):
        op = getattr(numpy, name)
        self.literal = literals.UfuncOperator(op)
        return

    def __call__(self, *args):
        for arg in args:
            self.literal.addLiteral(arg.literal)
        return self

def makeNamespace(consts, args, ops):
    """Make a namespace for building an equation from a string equation."""

    ns = {}

    for const in consts.values():
        # Create a class with the name of the const within the namespace
        if const in ns: continue
        c = ArgumentBuilder(value=float(const), name="const_"+const)
        ns[const] = c

    for arg in args.values():
        # Create a class with the name of the arg within the namespace
        if arg in ns: continue
        c = ArgumentBuilder(name=arg)
        ns[arg] = c

    for op in ops.values():
        if op in ns: continue
        c = OperatorBuilder(name=op)
        ns[op] = c

    return ns


def parseEquation(eqstr):
    """Scan an equation string for tokens.

    Raises: ValueError if the equation cannot be tokenized.
    """
    import tokenize
    import token
    import cStringIO
    import sys

    interface = cStringIO.StringIO(eqstr).readline
    # output is a 5-tuple
    # token[0] = token type
    # token[1] = token string
    # token[2] = (srow, scol) - row and col where the token begins
    # token[3] = (erow, ecol) - row and col where the token ends
    # token[4] = line where the token was found
    tokens = tokenize.generate_tokens(interface)
    tokens = list(tokens)
    #for val in tokens:
    #    print val[0], val[1]

    symbols = ("+", "-", "*", "/", "**", "%")

    # Scan for argumens and operators. This will be wrong on the first pass
    # since variables like "A" and "x" will scan as operators.

    args = {} # keys are index from the tokens list
    consts = {}
    ops = {}
    for i, tok in enumerate(tokens):
        if tok[0] in (token.NAME, token.OP):
            ops[i] = tok[1]
        elif tok[0] == token.NUMBER:
            consts[i] = tok[1]

    ignore = ("(", ",", ")")

    # Now scan the operators for variables. This will find variables that are
    # not in the numpy or builtin namespace (for now).

    poplist = []
    for i, op in ops.items():
        if op not in dir(numpy) and \
           op not in dir(sys.modules['__builtin__']) and \
           op not in symbols and \
           op not in ignore:
               args[i] = op
               poplist.append(i)
        if op in ignore or op in symbols:
           poplist.append(i)
    map(ops.pop, poplist)

    return consts, args, ops


def makeEquation(eqstr):

    consts, args, ops = parseEquation(eqstr)

    #print "args", args.values()
    #print "consts", consts.values()
    #print "ops", ops.values()

    ns = makeNamespace(consts, args, ops)
    #for key, val in ns.items():
    #    print key, val

    exec("_rootwrap = " + eqstr, ns)
    rootwrap = ns["_rootwrap"]
    root = rootwrap.literal
    #print "root", root
    eq = Equation(root)
    return eq

# version
__id__ = "$Id$"

#
# End of file
