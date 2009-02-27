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

#TODO - Need to handle assignments and a shared namespaces

from .Equation import Equation
import diffpy.srfit.equation.literals as literals
import numpy

class LiteralBuilder(object):
    """Class for building a literal from an equation."""

    def __init__(self):
        self.literal = None
        return

    # The infix operators
    def __eval(self, other, OperatorClass):
        #print OperatorClass.__name__, self.literal, other.literal
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

    def __neg__(self):
        op = literals.NegationOperator()
        op.addLiteral(self.literal)
        lb = LiteralBuilder()
        lb.literal = op
        return lb

class ArgumentBuilder(LiteralBuilder):

    def __init__(self, value=None, name=None):
        self.literal = literals.Argument(value=value, name=name)
        return

class OperatorBuilder(LiteralBuilder):
    """Acts like a numpy operator, but helps build a Literal tree."""

    def __init__(self, name):
        self.name = name
        self.literal = None
        self.op = getattr(numpy, name)
        return

    def __call__(self, *args):
        newobj = OperatorBuilder(self.name)
        newobj.literal = literals.UfuncOperator(self.op)
        #print self.literal,
        for arg in args:
            #print arg.literal,
            newobj.literal.addLiteral(arg.literal)
        #print
        return newobj

def makeNamespace(constmap, args, ops):
    """Make a namespace for building an equation from a string equation."""

    ns = {}

    for arg in args:
        # Create a class with the name of the arg within the namespace
        if arg in ns: continue
        value = constmap.get(arg)
        c = ArgumentBuilder(name=arg, value=value)
        ns[arg] = c

    for op in ops:
        if op in ns: continue
        c = OperatorBuilder(op)
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
    constants = ("e", "pi")

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
        if op in constants:
            consts[i] = op
            poplist.append(i)
        if op in ignore or op in symbols:
            poplist.append(i)


    map(ops.pop, poplist)

    # Now rename the constants
    constmap = {}
    usedconsts = {}
    tokranges = [ (tokens[i][2][1], tokens[i][3][1]) for i in consts ]
    tokranges.sort(reverse=True)
    idx = 0
    neweqstr = eqstr
    for lb, ub in tokranges:
        key = eqstr[lb:ub]
        if key in usedconsts:
            cname = usedconsts[key]
        else:
            cname = "_const%i"%idx
            idx += 1
            usedconsts[key] = cname
            # Evaluate the constant in the numpy namespace
            constmap[cname] = eval(eqstr[lb:ub], dict(numpy.__dict__))
            args[cname] = cname
        neweqstr = neweqstr[:lb] + cname + neweqstr[ub:]

    return neweqstr, constmap, args.values(), ops.values()


def makeEquation(eqstr):

    neweqstr, constmap, args, ops = parseEquation(eqstr)

    #print neweqstr
    #print "args", args
    #print "ops", ops
    #print "constmap", constmap

    ns = makeNamespace(constmap, args, ops)
    #for key, val in ns.items():
    #    print key, val

    exec("_rootwrap = " + neweqstr, ns)
    rootwrap = ns["_rootwrap"]
    root = rootwrap.literal
    #print "root", root
    eq = Equation(root)
    return eq

# version
__id__ = "$Id$"

#
# End of file
