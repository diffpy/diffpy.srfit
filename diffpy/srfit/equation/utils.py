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


def makeEquation(eqstr, consts = {}):
    """Make an equation from an equation string.

    Arguments
    eqstr   --  An equation in string form e.g. "e**x" or "A*sin(a*x + b)". The
                equation string can use any function that exists as a numpy
                ufunc.
    consts  --  A dictionary of constants indexed by name. Note that 'e' and
                'pi' are already considered, so don't need to be included.
    
    """
    builder = EquationBuilder(consts = consts)
    eq = builder.makeEquation(eqstr)
    return eq


class EquationBuilder(object):
    """This is a class for building an Equation instance from an equation string."""

    # Infix operators recognized by the builder
    symbols = ("+", "-", "*", "/", "**", "%")
    # Grouping operators recognized by the builder
    ignore = ("(", ",", ")")
    # Default constants recognized by the builder
    constants = {"e" : numpy.e, "pi" : numpy.pi}

    def __init__(self, consts = {}):
        """Initialize attribues.

        Arguments:
        eqstr   --  A string containing the equation
        name    --  An optional name for this equation.
        consts  --  A dictionary of user-defined constants. Keys are names and
                    values are the value of the constant (float, array, etc.)
        """
        self.consts = dict(consts)

        self._neweqstr = ""
        self._args = []
        self._ops = []
        return

    def parseEquationString(self, eqstr):
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


        # Scan for argumens and operators. This will be wrong on the first pass
        # since variables like "a" and "x" will appear as operators to the
        # tokenizer.
        args = {}
        consts = {}
        ops = {}
        
        for i, tok in enumerate(tokens):
            if tok[0] in (token.NAME, token.OP):
                ops[i] = tok[1]
            elif tok[0] == token.NUMBER:
                consts[i] = tok[1]

        # Scan the tokens in ops for variables and constants and move them to
        # the proper dictionary.
        poplist = []
        for i, tok in ops.items():
            # Move genuine varibles to the args dictionary
            if (
                # Check numpy
                tok not in dir(numpy) and
                # Check builtins
                tok not in dir(sys.modules['__builtin__']) and
                # Check symbols
                tok not in EquationBuilder.symbols and
                # Check ignored characters
                tok not in EquationBuilder.ignore
                ):
                args[i] = tok
                poplist.append(i)
            # Move constants to the const dictionary
            elif tok in EquationBuilder.constants or tok in self.constants:
                consts[i] = tok
                poplist.append(i)
            # Discard it if it is is in the ignore or symbol list
            elif tok in EquationBuilder.ignore or tok in EquationBuilder.symbols:
                poplist.append(i)

        # Discard the tokens that were moved or ignored
        map(ops.pop, poplist)


        # Rename numerical constants
        # this maps constants in the equation to their value
        self.consts.update( EquationBuilder.constants )
        # This gives us a namespace in which to evaluate constants, numerical
        # and otherwise.
        namespace = dict(numpy.__dict__)
        namespace.update( self.consts )
        # This maps initial constant names to their "safe" name
        constnames = dict( [(key, key) for key in self.consts] )
        constnames.update( 
                dict( [ (key, key) for key in EquationBuilder.constants] ) )
        # The ranges for each token in the original equation
        tokranges = [ (tokens[i][2][1], tokens[i][3][1]) for i in consts ]
        tokranges.sort(reverse=True)
        idx = 0
        neweqstr = eqstr
        # For each numerical constant, we need to change its name in the
        # equation, and record its new name and value in self.consts.
        for lb, ub in tokranges:
            # Get the name of the constant
            key = eqstr[lb:ub]
            # Get its new name if it has one
            if key in constnames:
                cname = constnames[key]
            # Otherwise, give it a new name and record its value
            else:
                cname = "_const%i"%idx
                idx += 1
                constnames[key] = cname
                # Evaluate the constant in the prepared namespace
                self.consts[cname] = eval(eqstr[lb:ub], namespace)

            # Make sure that the constant will be treated like an argument
            args[cname] = cname
            # Change the name of the constant in the equation string. This will
            # have no effect for named constants.
            neweqstr = neweqstr[:lb] + cname + neweqstr[ub:]

        # Now record all of the necesary info
        self._neweqstr = neweqstr
        self._args = args.values()
        self._ops = ops.values()

        return

    def makeEquation(self, eqstr):
        """Make an Equation from the parsed equation string."""

        self.parseEquationString(eqstr)
        ns = self.__makeBuilderNamespace()

        exec("_rootwrap = " + self._neweqstr, ns)
        rootwrap = ns["_rootwrap"]
        root = rootwrap.literal
        eq = Equation(root)
        return eq

    def __makeBuilderNamespace(self):
        """Make a namespace for holding LiteralBuilders.
        
        This takes the operators and arguments (literals) from the equation string and
        makes a instance of a builder class for each one (_ArgumentBuilder or
        _OperatorBuilder). These instances are stored in the namespace under the
        names of the literals for which they are built. When the equation is
        evaluated in this namespace, the builder instances will do the work of
        creating the Literal tree around which the Equation is built.
        """

        ns = {}

        # Create an _ArgumentBuilder for every argument and constant
        for arg in self._args:
            # Create a class with the name of the arg within the namespace
            if arg in ns: continue
            value = self.consts.get(arg)
            const = arg in self.consts
            c = _ArgumentBuilder(name=arg, value=value, const=const)
            ns[arg] = c

        # Create an _OperatorBuilder for each operator in the equation
        for op in self._ops:
            if op in ns: continue
            c = _OperatorBuilder(op)
            ns[op] = c

        return ns

# end class EquationBuilder

## These are used by the EquationBuilder class.

class _LiteralBuilder(object):
    """Class for building a literal from an equation."""

    def __init__(self):
        self.literal = None
        return

    # The infix operators

    def __evalBinary(self, other, OperatorClass):
        op = OperatorClass()
        op.addLiteral(self.literal)
        op.addLiteral(other.literal)
        lb = _LiteralBuilder()
        lb.literal = op
        return lb

    def __evalUnary(self, OperatorClass):
        op = OperatorClass()
        op.addLiteral(self.literal)
        lb = _LiteralBuilder()
        lb.literal = op
        return lb

    def __add__(self, other):
        return self.__evalBinary(other, literals.AdditionOperator)

    def __sub__(self, other):
        return self.__evalBinary(other, literals.SubtractionOperator)

    def __mul__(self, other):
        return self.__evalBinary(other, literals.MultiplicationOperator)

    def __div__(self, other):
        return self.__evalBinary(other, literals.DivisionOperator)

    def __pow__(self, other):
        return self.__evalBinary(other, literals.ExponentiationOperator)

    def __mod__(self, other):
        return self.__evalBinary(other, literals.RemainderOperator)

    def __neg__(self):
        return self.__evalUnary(literals.NegationOperator)

# end class _LiteralBuilder

class _ArgumentBuilder(_LiteralBuilder):

    def __init__(self, value=None, name=None, const=False):
        self.literal = literals.Argument(value=value, name=name, const=const)
        return

# end class _ArgumentBuilder

class _OperatorBuilder(_LiteralBuilder):
    """Acts like a numpy operator, but helps build a Literal tree."""

    def __init__(self, name):
        self.name = name
        self.literal = None
        self.op = getattr(numpy, name)
        return

    def __call__(self, *args):
        newobj = _OperatorBuilder(self.name)
        newobj.literal = literals.UfuncOperator(self.op)
        #print self.literal,
        for arg in args:
            #print arg.literal,
            newobj.literal.addLiteral(arg.literal)
        #print
        return newobj

# end class _OperatorBuilder

# version
__id__ = "$Id$"

#
# End of file
