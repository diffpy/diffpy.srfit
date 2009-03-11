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
"""EquationBuilder classes and utilities for creating Equations objects.
"""

# FIXME - the builder cannot handle numpy arrays on the left of a binary
# operation because this will automatically loop the operator of the right-side
# over the arguments of the array.

from .Equation import Equation
import diffpy.srfit.equation.literals as literals

import numpy

import sys

def makeEquation(eqstr, consts = {}):
    """Make an equation from an equation string.

    Arguments
    eqstr   --  An equation in string form using standard python syntax. The
                equation string can use any function that exists as a numpy
                ufunc or that has been created with the wrapFunction method. In
                the latter case, the name of a custom function used in the
                equation string must match the name used in the wrapFunction
                method.
    consts  --  A dictionary of named constants used in the equation.

    Returns an Equation instance representing the equation string.
    """
    ns = _makeNamespace(eqstr, consts)
    builder = eval(eqstr, ns)
    return builder.getEquation()

def _makeNamespace(eqstr, consts):
    """Build an evaluation namespace from an equation string.
    
    Arguments
    eqstr   --  An equation in string form using standard python syntax. The
                equation string can use any function that exists as a numpy
                ufunc or that has been created with the wrapFunction method. In
                the latter case, the name of a custom function used in the
                equation string must match the name used in the wrapFunction
                method.
    consts  --  A dictionary of named constants used in the equation.

    Returns a dictionary of the name, ArgumentBuilder pairs.
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

    # Scan for argumens and operators. This will be wrong on the first pass
    # since variables like "a" and "x" will appear as operators to the
    # tokenizer.
    eqargs = {}
    eqconsts = {}
    eqops = {}
    symbols = ("+", "-", "*", "/", "**", "%")
    ignore = ("(", ",", ")")

    consts = dict(consts)
    consts.update( { "pi" : numpy.pi, "e" : numpy.e } )

    for i, tok in enumerate(tokens):
        if tok[0] in (token.NAME, token.OP):
            eqops[i] = tok[1]
        #elif tok[0] == token.NUMBER:
        #    eqconsts[i] = tok[1]

    # Scan the tokens in ops for arguments and constants.
    poplist = []
    for i, tok in eqops.items():
        # Move genuine varibles to the eqargs dictionary
        if (
            # Check local namespace
            tok not in dir(sys.modules[__name__]) and
            # Check symbols
            tok not in symbols and
            # Check ignored characters
            tok not in ignore and
            # Check constants
            tok not in consts
            ):
            eqargs[i] = tok
            poplist.append(i)
        # Move constants to the eqconsts dictionary
        elif tok in consts:
            eqconsts[i] = tok
            poplist.append(i)
        # Discard it if it is is in the ignore or symbol list
        elif tok in ignore or tok in symbols:
            poplist.append(i)

    # Discard the tokens that were moved or ignored
    map(eqops.pop, poplist)

    # Now start making the namespace
    ns = {}
    # Get the operators
    for opname in eqops.values():
        #print opname
        opbuilder = getattr(sys.modules[__name__], opname)
        #print opbuilder
        ns[opname] = opbuilder
    # Make the arguments
    for argname in eqargs.values():
        #print argname
        argbuilder = ArgumentBuilder(name = argname)
        #print argbuilder
        ns[argname] = argbuilder
    # Get the constants
    for constname in eqconsts.values():
        #print constname
        val = eval(constname, consts)
        constbuilder = ArgumentBuilder(name = constname, value = val, 
            const = True)
        ns[constname] = constbuilder

    return ns


class EquationBuilder(object):
    """Class for building Equation objects.

    Equation builder objects can be composed like a normal function where the
    arguments can be other EquationBuilder instances or constants.
    """

    def __init__(self):
        """Initialize."""
        self.literal = None
        return

    def getEquation(self):
        """Get the equation built by this object."""
        eq = Equation(root = self.literal)
        return eq

    def __evalBinary(self, other, OperatorClass, onright = False):
        """Evaluate a binary function.

        Other can be an EquationBuilder or a constant.

        onright --  Indicates that the operator was passed on the right side
                    (defualt false).
        """
        # Create the Operator
        op = OperatorClass()

        # Reverse takes care of non-commutative operators, and assures that the
        # ordering is perserved.
        if not onright:
            # Add the literals to the operator
            op.addLiteral(self.literal)

        # Deal with constants
        if not isinstance(other, EquationBuilder):
            literal = literals.Argument(value=other, const=True)
        else:
            literal = other.literal
        #print "value:", OperatorClass, self.literal, literal
        op.addLiteral(literal)

        if onright:
            # Add the literals to the operator
            op.addLiteral(self.literal)

        # Create a new OperatorBuilder for the Operator
        opbuilder = OperatorBuilder(op.name)
        opbuilder.literal = op
        return opbuilder

    def __evalUnary(self, OperatorClass):
        """Evaluate a unary function."""
        op = OperatorClass()
        op.addLiteral(self.literal)
        opbuilder = OperatorBuilder(op.name)
        opbuilder.literal = op
        return opbuilder

    def __add__(self, other):
        return self.__evalBinary(other, literals.AdditionOperator)

    def __radd__(self, other):
        return self.__evalBinary(other, literals.AdditionOperator, True)

    def __sub__(self, other):
        return self.__evalBinary(other, literals.SubtractionOperator)

    def __rsub__(self, other):
        return self.__evalBinary(other, literals.SubtractionOperator, True)

    def __mul__(self, other):
        return self.__evalBinary(other, literals.MultiplicationOperator)

    def __rmul__(self, other):
        return self.__evalBinary(other, literals.MultiplicationOperator, True)

    def __div__(self, other):
        return self.__evalBinary(other, literals.DivisionOperator)

    def __rdiv__(self, other):
        return self.__evalBinary(other, literals.DivisionOperator, True)

    def __pow__(self, other):
        return self.__evalBinary(other, literals.ExponentiationOperator)

    def __rpow__(self, other):
        return self.__evalBinary(other, literals.ExponentiationOperator, True)

    def __mod__(self, other):
        return self.__evalBinary(other, literals.RemainderOperator)

    def __rmod__(self, other):
        return self.__evalBinary(other, literals.RemainderOperator, True)

    def __neg__(self):
        return self.__evalUnary(literals.NegationOperator)

## These are used by the class.

class ArgumentBuilder(EquationBuilder):

    def __init__(self, value=None, name=None, const=False):
        EquationBuilder.__init__(self)
        self.literal = literals.Argument(value=value, name=name, const=const)
        return

# end class ArgumentBuilder

class OperatorBuilder(EquationBuilder):
    """Acts like an operator, but helps build a Literal tree."""

    def __init__(self, name):
        EquationBuilder.__init__(self)
        self.name = name
        return

    def __call__(self, *args):
        newobj = OperatorBuilder(self.name)
        if self.literal is not None:
            op = literals.Operator()
            op.name = self.literal.name
            op.symbol = self.literal.name
            op.nin = self.literal.nin
            op.nout = self.literal.nout
            op.operation = self.literal.operation
            newobj.literal = op
        else:
            ufunc = getattr(numpy, self.name)
            newobj.literal = literals.UfuncOperator(ufunc)
        #print "value:", self.name,
        for arg in args:
            # Wrap the argument if it is not already
            if not isinstance(arg, EquationBuilder):
                arg = ArgumentBuilder(value=arg, const=True)
            newobj.literal.addLiteral(arg.literal)
            #print arg.literal, arg.literal.value
        #print
        return newobj

# end class OperatorBuilder

def wrapFunction(name, func, nin = 2, nout = 1):
    """Wrap a function in an OperatorBuilder instance.

    This will register the OperatorBuilder instance as an attribute of this
    module so it can be recognized in an equation string when parsed with the
    makeEquation method.
    
    name    --  The name of the funciton
    func    --  A callable python object
    nin     --  The number of input arguments
    nout    --  The number of return values

    Returns the OperatorBuilder instance that wraps the function.
    """
    op = literals.Operator()
    op.name = name
    op.symbol = name
    op.nin = nin
    op.nout = nout
    op.operation = func
    
    # Create the OperatorBuilder
    opbuilder = OperatorBuilder(name)
    opbuilder.literal = op

    # Register it with the module
    setattr(sys.modules[__name__], name, opbuilder)

    # Return it
    return opbuilder

# Export all numpy operators as OperatorBuilder instances.
for name in dir(numpy):
    op = getattr(numpy, name)
    if isinstance(op, numpy.ufunc):
        setattr(sys.modules[__name__], name, OperatorBuilder(name))

# version
__id__ = "$Id$"

#
# End of file
