#!/usr/bin/env python
##############################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2008 The Trustees of Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################

"""Classes and utilities for creating equations.

The EquationFactory class is used to create an equation (an Equation instance)
from a string equation.  User-defined Literals can be registered with the
factory so that they are used in the equation. Registered Literals are
referenced by name, and when a new Literal of the same name is registered, the
factory will swap out the old Literal for the new one in all equations built by
the factory.

An example using the EquationFactory:
The makeEquation method turns the string-representation of an equation into
a callable object.
> factory = EquationFactory()
> eq = factory.makeEquation("A*sin(a*x)")
will create an equation that evaluates as "A*sin(a*x)". The equation takes no
arguments.

Custom Arguments and constants can be included in the equation:
> factory.registerConstant("offset", 3)
> A = Argument(name = "A", value = 1.0)
> factory.registerArgument("A", A)
> eq = factory.makeEquation("A*sin(a*x) + offset")
This includes a constant offset in the equation and makes sure that the
user-defined Argument is in the equation. This can be used to assure that
the same instance of an Argument appears in multiple equations. Other literals
can be registered in a similar fashion.

The BaseBuilder class does the hard work of making an equation from a string in
EquationFactory.makeEquation. BaseBuilder can be used directly to create
equations. BaseBuilder is specified in the ArgumentBuilder and OperatorBuilder
classes.  You can create builders from Literals or equations by using the
"wrap" methods within this module or by using the builder classes directly.

With a collection of BaseBuilder objects, one can simply write the equation
using normal python syntax:
> A = ArgumentBuilder(name = "A")
> a = ArgumentBuilder(name = "a")
> x = ArgumentBuilder(name = "x")
> # sin is defined in this module as an OperatorBuilder
> sin = getBuilder("sin")
> beq = A*sin(a*x)
> eq = beq.getEquation()

The equation builder can also handle scalar constants. Staring with the above
setup:
> beq2 = A*sin(a*x) + 3
> eq2 = beq2.getEquation()
Here, we didn't have to wrap '3' in an ArgumentBuilder. Non scalars, constant
or otherwise, must be wrapped as ArgumentBuilders in order to be used in this
way.

BaseBuilder can make use of user-defined functions.  Any callable python
object can be wrapped as an OperatorBuilder with the wrapFunction method. For
example.
> _f = lambda a, b : (a-b)/(a+b)
> f = wrapFunction("f", _f)
> # Using BaseBuilder
> a = ArgumentBuilder(name = "a")
> b = ArgumentBuilder(name = "b")
> c = ArgumentBuilder(name = "c")
> beq = c*f(a,b)
> eq = beq.makeEquation()
"""

__all__ = ["EquationFactory", "BaseBuilder", "ArgumentBuilder",
           "OperatorBuilder", "wrapArgument", "wrapOperator", "wrapFunction",
           "getBuilder"]

# NOTE - the builder cannot handle numpy arrays on the left of a binary
# operation because the array will automatically loop the operator of the
# right-side over its arguments. This results in an array of BaseBuilder
# instances, not an BaseBuilder that contains an array.

_builders = {}


import inspect
import numbers
import numpy

import diffpy.srfit.equation.literals as literals
from diffpy.srfit.equation.literals.literal import Literal
from diffpy.srfit.equation.equationmod import Equation


class EquationFactory(object):
    """A Factory for equations.

    builders    --  A dictionary of BaseBuilders registered with the
                    factory, indexed by name.
    newargs     --  A set of new arguments created by makeEquation. This is
                    redefined whenever makeEquation is called.
    equations   --  Set of equations that have been built by the EquationFactory.
    """

    symbols = ("+", "-", "*", "/", "**", "%", "|")
    ignore = ("(", ",", ")")

    def __init__(self):
        """Initialize.

        This registers "pi" and "e" as constants within the factory.
        """
        self.builders = dict(_builders)
        self.newargs = set()
        self.equations = set()
        self.registerConstant("pi", numpy.pi)
        self.registerConstant("e", numpy.e)
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

        Returns a callable Literal representing the equation string.
        """
        self._prepareBuilders(eqstr, buildargs, argclass, argkw)
        beq = eval(eqstr, {}, self.builders)
        # handle scalar numbers or numpy arrays
        if isinstance(beq, (numbers.Number, numpy.ndarray)):
            lit = literals.Argument(value=beq, const=True)
            eq = Equation(name='', root=lit)
        else:
            eq = beq.getEquation()
            self.equations.add(eq)
        return eq

    def registerConstant(self, name, value):
        """Register a named constant with the factory.

        Returns the registered builder.
        """
        arg = literals.Argument(name=name, value=value, const=True)
        return self.registerArgument(name, arg)

    def registerArgument(self, name, arg):
        """Register a named Argument with the factory.

        Returns the registered builder.
        """
        argbuilder = wrapArgument(name, arg)
        return self.registerBuilder(name, argbuilder)

    def registerOperator(self, name, op):
        """Register an Operator literal with the factory.

        Operators can be used with or without arguments (or parentheses) in an
        equation string.  If used with arguments, then the Operator will use
        the passed arguments as arguments for the operation. If used without
        arguments, it is assumed that the operator is already populated with
        arguments, and those will be used.

        Returns the registered builder.
        """
        opbuilder = wrapOperator(name, op)
        return self.registerBuilder(name, opbuilder)

    def registerFunction(self, name, func, argnames):
        """Register a named function with the factory.

        This will register a builder for the function.

        name    --  The name of the function
        func    --  A callable python object
        argnames--  The argument names for func. If these names do not
                    correspond to builders, then new constants with value 0
                    will be created for each name.

        Returns the registered builder.
        """
        for n in argnames:
            if n not in self.builders:
                self.registerConstant(n, 0)
        opbuilder = wrapFunction(name, func, len(argnames))
        for n in argnames:
            b = self.builders[n]
            l = b.literal
            opbuilder.literal.addLiteral(l)

        return self.registerBuilder(name, opbuilder)

    def registerBuilder(self, name, builder):
        """Register builder in this module so it can be used in makeEquation.

        If an extant builder with the given name is already registered, this
        will replace all instances of the old builder's literal in the
        factory's equation set with the new builder's literal. Note that this
        may lead to errors if one of the replacements causes a self-reference.

        Raises ValueError if the new builder's literal causes a self-reference
        in an existing equation.
        """
        if not isinstance(name, basestring):
            raise TypeError("Name must be a string")
        if not isinstance(builder, BaseBuilder):
            raise TypeError("builder must be a BaseBuilder instance")
        # Swap out the old builder's literal, if necessary
        newlit = builder.literal
        swapbyname = isinstance(builder, ArgumentBuilder)
        bloldlits = set()
        if name in self.builders:
            bloldlits.add(self.builders[name].literal)
        for eq in self.equations:
            eqoldlits = bloldlits
            if swapbyname and name in eq.argdict:
                eqoldlits = bloldlits.union((eq.argdict[name],))
            for oldlit in eqoldlits:
                if oldlit is newlit:
                    continue
                eq.swap(oldlit, newlit)
        # Now store the new builder
        self.builders[name] = builder
        return builder

    def deRegisterBuilder(self, name):
        """De-register a builder by name.

        This does not change the equations that use the Literal wrapped by the
        builder.
        """
        if name in self.builders:
            del self.builders[name]
        return


    def wipeout(self, eq):
        """Invalidate the specified equation and remove it from the factory.

        This will remove the equation from the purview of the factory and
        also change its formula to return NaN.  This ensures that eq does
        not observe any object in the factory and thus prevents its indirect
        pickling with the factory because of observer callback function.

        No return value.
        """
        if eq is None:
            assert eq not in self.equations
            return
        self.equations.discard(eq)
        # invalidate this equation to clean up any observer relations of
        # objects in the factory towards its literals tree.
        nan = literals.Argument('nan', value=numpy.nan, const=True)
        eq.setRoot(nan)
        return


    def _prepareBuilders(self, eqstr, buildargs, argclass, argkw):
        """Prepare builders so that equation string can be evaluated.

        This method checks the equation string for errors and missing
        arguments, and creates new arguments if allowed. In the process it
        rebuilds the newargs attribute.

        Arguments
        eqstr       --  An equation in string as specified in the makeEquation
                        method.
        buildargs   --  A flag indicating whether missing arguments can be
                        created by the factory. If False, then the a ValueError
                        will be raised if there are undefined arguments in the
                        eqstr.
        argclass    --  Class to use when creating new Arguments. The class
                        constructor must accept the 'name' key word.
        argkw       --  Key word dictionary to pass to the argclass
                        constructor.

        Raises ValueError if new arguments must be created, but this is
        disallowed due to the buildargs flag.
        Raises SyntaxError if the equation string uses invalid syntax.

        Returns a dictionary of the name, BaseBuilder pairs.
        """

        eqargs = self._getUndefinedArgs(eqstr)

        # Raise an error if there are arguments that need to be created, but
        # this is disallowed.
        if not buildargs and eqargs:
            eqargsstr = ", ".join(eqargs)
            msg = "The equation contains undefined arguments: %s"%eqargsstr
            raise ValueError(msg)

        # Make the arguments
        newargs = set()
        for argname in eqargs:
            arg = argclass(name = argname, **argkw)
            argbuilder = ArgumentBuilder(name = argname, arg = arg)
            newargs.add(arg)
            self.registerBuilder(argname, argbuilder)

        self.newargs = newargs

        return

    def _getUndefinedArgs(self, eqstr):
        """Get the undefined arguments from eqstr.

        This tokenizes eqstr and extracts undefined arguments. An undefined
        argument is defined as any token that is not a special character that
        does not correspond to a builder.

        Raises SyntaxError if the equation string uses invalid syntax.
        """
        import tokenize
        import token
        import cStringIO

        interface = cStringIO.StringIO(eqstr).readline
        # output is an iterator. Each entry (token) is a 5-tuple
        # token[0] = token type
        # token[1] = token string
        # token[2] = (srow, scol) - row and col where the token begins
        # token[3] = (erow, ecol) - row and col where the token ends
        # token[4] = line where the token was found
        tokens = tokenize.generate_tokens(interface)

        # Scan for tokens. Throw a SyntaxError if the tokenizer chokes.
        args = set()

        try:
            for tok in tokens:
                if tok[0] in (token.NAME, token.OP):
                    args.add(tok[1])
        except tokenize.TokenError:
            m = "invalid syntax: '%s'"%eqstr
            raise SyntaxError(m)

        # Scan the tokens for names that do not correspond to registered
        # builders. These will be treated as arguments that need to be
        # generated.
        for tok in set(args):
            # Move genuine varibles to the eqargs dictionary
            if (
                # Check registered builders
                tok in self.builders or
                # Check symbols
                tok in EquationFactory.symbols or
                # Check ignored characters
                tok in EquationFactory.ignore
                ):
                args.remove(tok)

        return args

# End class EquationFactory

class BaseBuilder(object):
    """Class for building equations.

    Equation builder objects can be composed like a normal function where the
    arguments can be other BaseBuilder instances or constants.

    Attributes
    literal     --  The equation Literal being built by this instance.
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
        """Get the equation built by this object.

        The equation will given the name "_eq_<root>" where "<root>" is the
        name of the root node.
        """
        # We need to make a name for this, so we name it after its root
        name = "_eq_%s"%self.literal.name
        eq = Equation(name, self.literal)
        return eq

    def __evalBinary(self, other, OperatorClass, onleft = True):
        """Evaluate a binary function.

        Other can be an BaseBuilder or a constant.

        onleft  --  Indicates that the operator was passed on the left side
                    (defualt True).
        """
        # Create the Operator
        op = OperatorClass()

        # onleft takes care of non-commutative operators, and assures that the
        # ordering is perserved.
        if onleft:
            # Add the literals to the operator
            op.addLiteral(self.literal)

        # Deal with constants
        if isinstance(other, BaseBuilder):
            literal = other.literal
        elif isinstance(other, Literal):
            literal = other
        else:
            literal = literals.Argument(value=other, const=True)
        op.addLiteral(literal)

        if not onleft:
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
        return self.__evalBinary(other, literals.AdditionOperator, False)

    def __sub__(self, other):
        return self.__evalBinary(other, literals.SubtractionOperator)

    def __rsub__(self, other):
        return self.__evalBinary(other, literals.SubtractionOperator, False)

    def __mul__(self, other):
        return self.__evalBinary(other, literals.MultiplicationOperator)

    def __rmul__(self, other):
        return self.__evalBinary(other, literals.MultiplicationOperator, False)

    def __div__(self, other):
        return self.__evalBinary(other, literals.DivisionOperator)

    def __rdiv__(self, other):
        return self.__evalBinary(other, literals.DivisionOperator, False)

    def __pow__(self, other):
        return self.__evalBinary(other, literals.ExponentiationOperator)

    def __rpow__(self, other):
        return self.__evalBinary(other, literals.ExponentiationOperator, False)

    def __mod__(self, other):
        return self.__evalBinary(other, literals.RemainderOperator)

    def __rmod__(self, other):
        return self.__evalBinary(other, literals.RemainderOperator, False)

    def __neg__(self):
        return self.__evalUnary(literals.NegationOperator)

## These are used by the class.

class ArgumentBuilder(BaseBuilder):
    """BaseBuilder wrapper around an Argument literal.

    Equation builder objects can be composed like a normal function where the
    arguments can be other BaseBuilder instances or constants.

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
        BaseBuilder.__init__(self)
        if arg is None:
            self.literal = literals.Argument(value=value, name=name,
                    const=const)
        else:
            self.literal = arg
        return

# end class ArgumentBuilder

class OperatorBuilder(BaseBuilder):
    """BaseBuilder wrapper around an Operator literal.

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
        BaseBuilder.__init__(self)
        self.name = name
        self.literal = op
        return

    def __call__(self, *args):
        """Call the operator builder.

        This creates a new builder that encapsulates the operation.

        args    --  Arguments of the operation.

        Raises ValueError if self.literal.nin >= 0 and len(args) != op.nin
        """
        newobj = OperatorBuilder(self.name)

        # If all we have is a name, then we assume that it is the name of a
        # numpy operator, and we use the corresponding Operator.
        if self.literal is None:
            ufunc = getattr(numpy, self.name)
            self.literal = literals.UFuncOperator(ufunc)
        # Here the Operator is already specified.  We can copy its attributes
        # to a new Operator inside of the new OperatorBuilder.
        op = literals.makeOperator(name=self.literal.name,
                                   symbol=self.literal.symbol,
                                   nin=self.literal.nin,
                                   nout=self.literal.nout,
                                   operation=self.literal.operation)
        newobj.literal = op

        # Now that we have a literal, let's check our inputs
        literal = newobj.literal
        if literal.nin >= 0 and len(args) != literal.nin:
            raise ValueError("%s takes %i arguments (%i given)"%\
                    (self.literal, self.literal.nin, len(args)))

        # Wrap scalar arguments
        for i, arg in enumerate(args):
            # Wrap the argument if it is not already
            if not isinstance(arg, BaseBuilder):
                name = self.name + "_%i"%i
                arg = ArgumentBuilder(value = arg, name = name, const = True)
            newobj.literal.addLiteral(arg.literal)

        return newobj

# end class OperatorBuilder

# Utility functions

def wrapArgument(name, arg):
    """Wrap an Argument as a builder."""
    argbuilder = ArgumentBuilder(arg = arg)
    return argbuilder

def wrapOperator(name, op):
    """Wrap an Operator as a builder."""
    opbuilder = OperatorBuilder(name, op)
    return opbuilder

def wrapFunction(name, func, nin=2, nout=1):
    """Wrap a function in an OperatorBuilder instance.

    name    --  The name of the function
    func    --  A callable python object
    nin     --  The number of input arguments (default 2)
    nout    --  The number of return values (default 1)

    Returns the OperatorBuilder instance that wraps the function.
    """
    op = literals.makeOperator(name=name, symbol=name,
                               nin=nin, nout=nout,
                               operation=func)

    # Create the OperatorBuilder
    opbuilder = OperatorBuilder(name, op)

    return opbuilder

def getBuilder(name):
    """Get an operator from the global builders dictionary."""
    return _builders[name]

def __wrapNumpyOperators():
    """Export all numpy operators as OperatorBuilder instances in the module
    namespace.
    """
    for name in dir(numpy):
        op = getattr(numpy, name)
        if isinstance(op, numpy.ufunc):
            _builders[name] = OperatorBuilder(name)
    return
__wrapNumpyOperators()

# Register other functions as well
def __wrapSrFitOperators():
    """Export all non-base operators from the
    diffpy.srfit.equation.literals.operators module as OperatorBuilder
    instances in the module namespace.
    """
    opmod = literals.operators
    excluded_types = set((opmod.CustomOperator, opmod.UFuncOperator))
    # check if opmod member should be wrapped as OperatorBuilder
    is_exported_type = lambda cls : (
        inspect.isclass(cls) and issubclass(cls, opmod.Operator) and
        not inspect.isabstract(cls) and
        not cls in excluded_types)
    # create OperatorBuilder objects
    for nm, opclass in inspect.getmembers(opmod, is_exported_type):
        op = opclass()
        assert op.name, "Unnamed Operator should never appear here."
        _builders[op.name] = OperatorBuilder(op.name, op)
    return
__wrapSrFitOperators()

# End of file
