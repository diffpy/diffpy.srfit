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
"""Classes and utilities for creating evaluation networks.

The EquationFactory class is used to create an evaluation network isolated
namespace. A string equation can be given to an EquationFactory to be turned
into an evaluation network. User-defined literals can be registered with the
factory so that they are used in the equation.

An example using the EquationFactory:
The makeEquation method turns the string-representation of an equation into
an evaluation network.
> factory = EquationFactory()
> eq = factory.makeEquation("A*sin(a*x)")
will create an Equation with Arguments A, a and x that evaluates as
"A*sin(a*x)". The Arguments A, a and x are created by the factory.

Custom Arguments and constants can be included in the equation:
> factory.registerConstant("offset", 3)
> A = Argument(name = "A", value = 1.0)
> factory.registerArgument("A", A)
> eq = factory.makeEquation("A*sin(a*x) + offset")
This includes a constant offset in the equation and makes sure that the
user-defined Argument is in the equation. This is can be used to assure that
the same instance of an Argument appears in multiple equation. Other literals
can be registered in a similar fashion.

The EquationBuilder class does the hard work of making an evaluation network
from a string in EquationFactory.makeEquation. EquationBuilder can be used
directly to create evaluation networks. EquationBuilder instances overload
normal arithmetic operations so that they build an evaluation network instead.
EquationBuilder is specified in the ArgumentBuilder and OperatorBuilder
classes.  All of the numpy ufunc operators are overloaded within this module as
OperatorBuilder instances. You can create builders from literals by using the
"wrap" methods within this module or by using the builder classes directly.

With a collection of EquationBuilder objects, one can create the evaluation
network using normal python syntax:
> A = ArgumentBuilder(name = "A")
> a = ArgumentBuilder(name = "a")
> x = ArgumentBuilder(name = "x")
> # sin is defined in this module as an OperatorBuilder
> beq = A*sin(a*x)
> eq = beq.getEquation()

The equation builder can also handle scalar constants. Staring with the above
setup:
> beq2 = A*sin(a*x) + 3
> eq2 = beq2.getEquation()
Here, the object returned by 'A*sin(a*x)' knows how to add itself to the scalar
'3' to return another EquationBuilder object.  Non scalars, constant or
otherwise, must be wrapped as ArgumentBuilders in order to be used in this way.

EquationBuilder can make use of user-defined functions.  Any callable python
object can be wrapped as an OperatorBuilder with the wrapFunction method. For
example.
> _f = lambda a, b : (a-b)/(a+b)
> f = wrapFunction("f", _f)
> # Using EquationBuilder
> a = ArgumentBuilder(name = "a")
> b = ArgumentBuilder(name = "b")
> c = ArgumentBuilder(name = "c")
> beq = c*f(a,b)
> eq = beq.makeEquation()

"""

# NOTE - the builder cannot handle numpy arrays on the left of a binary
# operation because the array will automatically loop the operator of the
# right-side over its arguments. This results in an array of EquationBuilder
# instances, not an EquationBuilder that contains an array.

_builders = {}

import sys
import numpy

import diffpy.srfit.equation.literals as literals

class EquationFactory(object):
    """A Factory for Equation classes.

    builders    --  A dictionary of EquationBuilders registered with the
                    factory, indexed by name.
    newargs     --  A set of new arguments created by makeEquation. This is
                    redefined whenever makeEquation is called.
    
    """

    symbols = ("+", "-", "*", "/", "**", "%")
    ignore = ("(", ",", ")")

    def __init__(self, builders = {}):
        """Initialize.

        This registers "pi" and "e" as constants within the factory.

        builders    --  A dictionary of builders to start out with (default
                        {}).

        """
        self.builders = dict(builders)
        self.registerConstant("pi", numpy.pi)
        self.registerConstant("e", numpy.e)
        self.newargs = set()
        return

    def makeEquation(self, eqstr, buildargs = True, argclass =
            literals.Argument, argkw = {}):
        """Make an equation from an equation string.

        Arguments
        eqstr       --  An equation in string form using standard python
                        syntax.  The equation string can use any function
                        registered literal or function, including numpy ufuncs
                        that are automatically registered.
        buildargs   --  A flag indicating whether missing arguments can be
                        created by the Factory (default True). If False, then
                        the a ValueError will be raised if there are undefined
                        arguments in the eqstr. Built arguments will be of type
                        argclass.
        argclass    --  Class to use when creating new Arguments (default
                        diffpy.srfit.equation.literals.Argument). The class
                        constructor must accept the 'name' key word.
        argkw       --  Key word dictionary to pass to the argclass constructor
                        (default {}).

        Returns an Equation instance representing the equation string.

        """
        ns, args = self._makeNamespace(eqstr, buildargs, argclass, argkw)
        self.newargs = args
        beq = eval(eqstr, ns)
        return beq.getEquation()

    def registerConstant(self, name, value):
        """Register a named constant with the factory."""
        arg = literals.Argument(name=name, value=value)
        arg.const =  True
        self.registerArgument(name, arg)
        return

    def registerArgument(self, name, arg):
        """Register a named Argument with the factory."""
        argbuilder = wrapArgument(name, arg)
        self.registerBuilder(name, argbuilder)
        return

    def registerFunction(self, name, func, nin = 2, nout = 1):
        """Register a named function with the factory.

        This will register the OperatorBuilder instance as an attribute of this
        module so it can be recognized in an equation string when parsed with
        the makeEquation method.
        
        name    --  The name of the funciton
        func    --  A callable python object
        nin     --  The number of input arguments (default 2)
        nout    --  The number of return values (default 1)

        """
        opbuilder = wrapFunction(name, func, nin, nout)
        self.registerBuilder(name, opbuilder)
        return

    def registerBuilder(self, name, builder):
        """Register builder in this module so it can be used in makeEquation."""
        # Check if we already have a builder with this name. If so, then we
        # need to change the target of that object and add the new builder.
        if name in self.builders:
            self.builders[name].literal.setTarget(builder.literal)
        self.builders[name] = builder
        return

    def deRegisterBuilder(self, name):
        """De-register a builder by name."""
        if name in self.builders:
            del self.builders[name]
        return

    def _makeNamespace(self, eqstr, buildargs, argclass, argkw):
        """Create an evaluation namespace from an equation string.
        
        Arguments
        eqstr       --  An equation in string as specified in the makeEquation
                        method.
        buildargs   --  A flag indicating whether missing arguments can be
                        created by the Factory. If False, then the a ValueError
                        will be raised if there are undefined arguments in the
                        eqstr.
        argclass    --  Class to use when creating new Arguments. The class
                        constructor must accept the 'name' key word.
        argkw       --  Key word dictionary to pass to the argclass
                        constructor.

        Returns a dictionary of the name, EquationBuilder pairs.

        """
        import sys

        eqops, eqargs = self._getOpsAndArgs(eqstr)

        # Raise an error if there are arguments that need to be created, but
        # this is disallowed.
        if not buildargs and eqargs:
            msg = "The equation contains undefined arguments %s"%eqargs
            raise ValueError(msg)

        # Now start making the namespace
        ns = dict(self.builders)

        # Get the operators, partitions and generators
        for opname in eqops:
            if opname not in self.builders:
                opbuilder = getBuilder(opname)
                ns[opname] = opbuilder

        # Make the arguments
        newargs = set()
        for argname in eqargs:
            if argname not in self.builders:
                arg = argclass(name = argname, **argkw)
                argbuilder = ArgumentBuilder(name = argname, arg = arg)
                ns[argname] = argbuilder
                newargs.add(arg)

        return ns, newargs

    def _getOpsAndArgs(self, eqstr):
        """Get the Operator and Argument names from an equation."""
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
        # We have overloaded list and other python types. Thus, we go to the
        # source.
        tokens = list(tokens)

        # Scan for argumens and operators. This will be wrong on the first pass
        # since variables like "a" and "x" will appear as operators to the
        # tokenizer.
        eqargs = set()
        eqops = set()

        for i, tok in enumerate(tokens):
            if tok[0] in (token.NAME, token.OP):
                eqops.add(tok[1])

        # Scan the tokens for names that are not defined in the module or in
        # self.builders. These will be treated as Arguments that need to be
        # generated.
        for tok in set(eqops):
            # Move genuine varibles to the eqargs dictionary
            if (
                # Check module namespace
                tok not in _builders and
                # Check custom builders
                tok not in self.builders and
                # Check symbols
                tok not in EquationFactory.symbols and
                # Check ignored characters
                tok not in EquationFactory.ignore
                ):
                eqargs.add(tok)
                eqops.remove(tok)
            # Discard it if it is is in the ignore or symbol list
            elif tok in EquationFactory.ignore\
                    or tok in EquationFactory.symbols:
                eqops.remove(tok)

        return eqops, eqargs

# End class EquationFactory

class EquationBuilder(object):
    """Class for building Equation objects.

    Equation builder objects can be composed like a normal function where the
    arguments can be other EquationBuilder instances or constants.

    Attributes
    literal     --  The root of the Equation being built by this instance.

    """

    def __init__(self):
        """Initialize."""
        self.literal = None
        return

    def __call__(self, *args):
        """Raises exception for easier debugging."""
        m = "%s (%s) cannot accept arguments"%\
            (self.literal.name, self.__class__.__name__)
        raise TypeError(m)

    def getEquation(self):
        """Get the equation built by this object."""
        return self.literal

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
            literal = literals.Argument(value=other)
            literal.const = True
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
    """EquationBuilder wrapper around an Argument literal.

    Equation builder objects can be composed like a normal function where the
    arguments can be other EquationBuilder instances or constants.

    Attributes
    literal     --  The Argument wrapped by this instance.

    """

    def __init__(self, value = None, name = None, const = False, arg = None):
        """Create an ArgumentBuilder instance, containing a new Argument.

        Arguments
        value   --  The value of the wrapped Argument (float, default None)
        name    --  The name of the wrapped Argument (string, default None)
        const   --  Flag indicating whether the Argument is constant (bool,
                    default False)
        arg     --  A pre-defined Argument to use. If this is None (default),
                    then a new Argument will be created from value, name and
                    const.

        """
        EquationBuilder.__init__(self)
        if arg is None:
            self.literal = literals.Argument(value=value, name=name)
            self.literal.const = const
        else:
            self.literal = arg
        return

# end class ArgumentBuilder


class OperatorBuilder(EquationBuilder):
    """EquationBuilder wrapper around an Operator literal.

    Equation builder objects can be composed like a normal function where the
    arguments can be other EquationBuilder instances or constants.

    Attributes
    literal     --  The Operator wrapped by this instance.
    name        --  The name of the operator to be wrapped

    """

    def __init__(self, name, op = None):
        """Wrap an Operator or a function by name.

        Arguments
        name    --  The name of the wrapped Operator
        op      --  If specified, this sets the literal attribute as this
                    operator (default None). Otherwise, the name is assumed to
                    be that of a numpy ufunc, which is used to specify the
                    Operator.

        """
        EquationBuilder.__init__(self)
        self.name = name
        self.literal = op
        return

    def __call__(self, *args, **kw):
        """Call the operator builder.

        This creates a new builder that encapsulates the operation.
        
        args    --  Arguments of the operation.

        """
        newobj = OperatorBuilder(self.name)

        # If all we have is a name, then we assume that it is the name of a
        # numpy operator, and we use the corresponding Operator.
        if self.literal is None:
            ufunc = getattr(numpy, self.name)
            newobj.literal = literals.UFuncOperator(ufunc)
            self.literal = newobj.literal
        # If the Operator is already specified, then copy its attributes to a
        # new Operator inside of the new OperatorBuilder.
        else:
            nin = self.literal.nin
            if nin < 0:
                nin = len(args)
            op = literals.Operator()
            op.name = self.literal.name
            op.symbol = self.literal.name
            op.nin = nin
            op.nout = self.literal.nout
            op.operation = self.literal.operation
            newobj.literal = op

        # Do a quick check. This is the only thing we can require, since args
        # can hold tags for the equation.
        nin = newobj.literal.nin
        if len(args) != nin:
            m = "%i arguments required, %i recieved"%(nin, len(args))
            raise ValueError(m)

        # Wrap scalar arguments and process tags
        for arg in args:
            # Wrap the argument if it is not already
            if not isinstance(arg, EquationBuilder):
                arg = ArgumentBuilder(value=arg, const=True)
            newobj.literal.addLiteral(arg.literal)
        return newobj

# end class OperatorBuilder

# Utility functions

def wrapArgument(name, arg):
    """Wrap an Argument as a builder."""
    argbuilder = ArgumentBuilder(arg = arg)
    return argbuilder

def wrapFunction(name, func, nin = 2, nout = 1):
    """Wrap a function in an OperatorBuilder instance.

    name    --  The name of the funciton
    func    --  A callable python object
    nin     --  The number of input arguments (default 2)
    nout    --  The number of return values (default 1)

    Returns the OperatorBuilder instance that wraps the function.

    """
    op = literals.Operator()
    op.name = name
    op.symbol = name
    op.nin = nin
    op.nout = nout
    op.operation = func
    
    # Create the OperatorBuilder
    opbuilder = OperatorBuilder(name, op)

    # Return it
    return opbuilder

def __wrapNumpyOperators():
    """Export all numpy operators as OperatorBuilder instances in the module
    namespace.

    """
    for name in dir(numpy):
        op = getattr(numpy, name)
        if isinstance(op, numpy.ufunc):
            _builders[name] = OperatorBuilder(name)

__wrapNumpyOperators()

# Register other functions as well

def __wrapSrFitOperators():
    """Export all non-base operators from the
    diffpy.srfit.equation.literals.operators module as OperatorBuilder
    instances in the module namespace.

    """
    import inspect
    opmod = literals.operators
    for opname in dir(opmod):
        opclass = getattr(opmod, opname)
        if inspect.isclass(opclass) \
            and issubclass(opclass, opmod.Operator) \
            and opclass is not opmod.Operator \
            and opclass is not opmod.UFuncOperator:

            op = opclass()
            _builders[op.name] = OperatorBuilder(op.name, op)

    return
__wrapSrFitOperators()

def getBuilder(name):
    """Get an operator from the global builders dictionary."""
    return _builders[name]

# version
__id__ = "$Id$"

#
# End of file
