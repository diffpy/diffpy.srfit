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
The make_equation method turns the string-representation of an equation into
a callable object.
> factory = EquationFactory()
> eq = factory.make_equation("A*sin(a*x)")
will create an equation that evaluates as "A*sin(a*x)". The equation takes no
arguments.

Custom Arguments and constants can be included in the equation:
> factory.register_constant("offset", 3)
> A = Argument(name = "A", value = 1.0)
> factory.register_argument("A", A)
> eq = factory.make_equation("A*sin(a*x) + offset")
This includes a constant offset in the equation and makes sure that the
user-defined Argument is in the equation. This can be used to assure that
the same instance of an Argument appears in multiple equations. Other literals
can be registered in a similar fashion.

The BaseBuilder class does the hard work of making an equation from a string in
EquationFactory.make_equation. BaseBuilder can be used directly to create
equations. BaseBuilder is specified in the ArgumentBuilder and OperatorBuilder
classes.  You can create builders from Literals or equations by using the
"wrap" functions within this module or by using the builder classes directly.

With a collection of BaseBuilder objects, one can simply write the equation
using normal python syntax:
> A = ArgumentBuilder(name = "A")
> a = ArgumentBuilder(name = "a")
> x = ArgumentBuilder(name = "x")
> # sin is defined in this module as an OperatorBuilder
> sin = get_builder("sin")
> beq = A*sin(a*x)
> eq = beq.get_equation()

The equation builder can also handle scalar constants. Staring with the above
setup:
> beq2 = A*sin(a*x) + 3
> eq2 = beq2.get_equation()
Here, we didn't have to wrap '3' in an ArgumentBuilder. Non scalars, constant
or otherwise, must be wrapped as ArgumentBuilders in order to be used in this
way.

BaseBuilder can make use of user-defined functions.  Any callable python
object can be wrapped as an OperatorBuilder with the wrap_function function.
For example.
> _f = lambda a, b : (a-b)/(a+b)
> f = wrap_function("f", _f)
> # Using BaseBuilder
> a = ArgumentBuilder(name = "a")
> b = ArgumentBuilder(name = "b")
> c = ArgumentBuilder(name = "c")
> beq = c*f(a,b)
> eq = beq.get_equation()
"""

import inspect
import io
import numbers
import token
import tokenize

import numpy

import diffpy.srfit.equation.literals as literals
from diffpy.srfit.equation.equationmod import Equation
from diffpy.srfit.equation.literals.literal import Literal
from diffpy.utils._deprecator import build_deprecation_message, deprecated

__all__ = [
    "EquationFactory",
    "BaseBuilder",
    "ArgumentBuilder",
    "OperatorBuilder",
    "wrap_argument",
    "wrap_operator",
    "wrap_function",
    "get_builder",
    "wrapArgument",
    "wrapOperator",
    "wrapFunction",
    "getBuilder",
]

# NOTE - the builder cannot handle numpy arrays on the left of a binary
# operation because the array will automatically loop the operator of the
# right-side over its arguments. This results in an array of BaseBuilder
# instances, not an BaseBuilder that contains an array.

_builders = {}

EquationFactory_base = "diffpy.srfit.equation.builder.EquationFactory"
builder_base = "diffpy.srfit.equation.builder"
removal_version = "4.0.0"

registerFunction_dep_msg = build_deprecation_message(
    EquationFactory_base,
    "registerFunction",
    "register_function",
    removal_version,
)
makeEquation_dep_msg = build_deprecation_message(
    EquationFactory_base, "makeEquation", "make_equation", removal_version
)
registerConstant_dep_msg = build_deprecation_message(
    EquationFactory_base,
    "registerConstant",
    "register_constant",
    removal_version,
)
registerArgument_dep_msg = build_deprecation_message(
    EquationFactory_base,
    "registerArgument",
    "register_argument",
    removal_version,
)
registerOperator_dep_msg = build_deprecation_message(
    EquationFactory_base,
    "registerOperator",
    "register_operator",
    removal_version,
)
registerBuilder_dep_msg = build_deprecation_message(
    EquationFactory_base,
    "registerBuilder",
    "register_builder",
    removal_version,
)
deRegisterBuilder_dep_msg = build_deprecation_message(
    EquationFactory_base,
    "deRegisterBuilder",
    "deregister_builder",
    removal_version,
)
wrapArgument_dep_msg = build_deprecation_message(
    builder_base, "wrapArgument", "wrap_argument", removal_version
)
wrapOperator_dep_msg = build_deprecation_message(
    builder_base, "wrapOperator", "wrap_operator", removal_version
)
wrapFunction_dep_msg = build_deprecation_message(
    builder_base, "wrapFunction", "wrap_function", removal_version
)
getBuilder_dep_msg = build_deprecation_message(
    builder_base, "getBuilder", "get_builder", removal_version
)


class EquationFactory(object):
    """A Factory for equations.

    Attributes
    ----------
    builders : dict
        The BaseBuilders registered with the factory, indexed by name.
    newargs : set
        The new arguments created by the last call to ``make_equation``.
    equations : set
        The equations that have been built by the factory.
    """

    symbols = ("+", "-", "*", "/", "**", "%", "|")
    ignore = ("(", ",", ")")

    def __init__(self):
        """Initialize the factory and register "pi" and "e" as
        constants."""
        self.builders = dict(_builders)
        self.newargs = set()
        self.equations = set()
        self.register_constant("pi", numpy.pi)
        self.register_constant("e", numpy.e)
        return

    def make_equation(
        self, eqstr, buildargs=True, argclass=literals.Argument, argkw={}
    ):
        """Make an equation from an equation string.

        Parameters
        ----------
        eqstr : str
            The equation in string form using standard python syntax. The
            equation string can use any registered literal or function,
            including the numpy ufuncs that are registered automatically.
        buildargs : bool, optional
            The flag that is True (default) to let the factory create
            missing arguments. If False, undefined names in ``eqstr``
            raise a ValueError.
        argclass : type, optional
            The class used to create new arguments. Its constructor must
            accept the ``name`` keyword. Default is
            ``diffpy.srfit.equation.literals.Argument``.
        argkw : dict, optional
            The keyword arguments passed to the ``argclass`` constructor.
            Default is an empty dict.

        Returns
        -------
        Equation
            The callable equation built from ``eqstr``.

        Raises
        ------
        ValueError
            If ``buildargs`` is False and ``eqstr`` uses undefined names.
        SyntaxError
            If ``eqstr`` is not valid python syntax.
        """
        self._prepare_builders(eqstr, buildargs, argclass, argkw)
        beq = eval(eqstr, {}, self.builders)
        # handle scalar numbers or numpy arrays
        if isinstance(beq, (numbers.Number, numpy.ndarray)):
            lit = literals.Argument(value=beq, const=True)
            eq = Equation(name="", root=lit)
        else:
            eq = beq.get_equation()
            self.equations.add(eq)
        return eq

    def register_constant(self, name, value):
        """Register a named constant with the factory.

        Parameters
        ----------
        name : str
            The name of the constant.
        value : object
            The value of the constant.

        Returns
        -------
        ArgumentBuilder
            The registered builder.
        """
        arg = literals.Argument(name=name, value=value, const=True)
        return self.register_argument(name, arg)

    def register_argument(self, name, arg):
        """Register a named Argument with the factory.

        Parameters
        ----------
        name : str
            The name to register the argument under.
        arg : Argument
            The argument to register.

        Returns
        -------
        ArgumentBuilder
            The registered builder.
        """
        argbuilder = wrap_argument(name, arg)
        return self.register_builder(name, argbuilder)

    def register_operator(self, name, op):
        """Register an Operator literal with the factory.

        Operators can be used with or without arguments (or parentheses)
        in an equation string.  If used with arguments, then the
        Operator will use the passed arguments as arguments for the
        operation. If used without arguments, it is assumed that the
        operator is already populated with arguments, and those will be
        used.

        Parameters
        ----------
        name : str
            The name to register the operator under.
        op : Operator
            The operator to register.

        Returns
        -------
        OperatorBuilder
            The registered builder.
        """
        opbuilder = wrap_operator(name, op)
        return self.register_builder(name, opbuilder)

    def register_function(self, name, func, argnames):
        """Register a named function with the factory.

        This will register a builder for the function.

        Parameters
        ----------
        name : str
            The name of the function.
        func : callable
            The callable python object.
        argnames : list of str
            The argument names for func. If these names do not
            correspond to builders, then new constants with value 0
            will be created for each name.

        Returns
        -------
        registered_builder : OperatorBuilder
            The registered builder.
        """
        for n in argnames:
            if n not in self.builders:
                self.register_constant(n, 0)
        opbuilder = wrap_function(name, func, len(argnames))
        for argname in argnames:
            builder = self.builders[argname]
            argliteral = builder.literal
            opbuilder.literal.add_literal(argliteral)
        registered_builder = self.register_builder(name, opbuilder)
        return registered_builder

    def register_builder(self, name, builder):
        """Register a builder so it can be used in ``make_equation``.

        If an extant builder with the given name is already registered,
        this will replace all instances of the old builder's literal in
        the factory's equation set with the new builder's literal. Note
        that this may lead to errors if one of the replacements causes a
        self-reference.

        Parameters
        ----------
        name : str
            The name to register the builder under.
        builder : BaseBuilder
            The builder to register.

        Returns
        -------
        BaseBuilder
            The registered builder.

        Raises
        ------
        TypeError
            If ``name`` is not a string or ``builder`` is not a
            BaseBuilder.
        ValueError
            If the new builder's literal causes a self-reference in an
            existing equation.
        """
        if not isinstance(name, str):
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

    def deregister_builder(self, name):
        """Deregister a builder by name.

        This does not change the equations that use the Literal wrapped
        by the builder. Nothing happens if no builder has this name.

        Parameters
        ----------
        name : str
            The name of the builder to deregister.
        """
        if name in self.builders:
            del self.builders[name]
        return

    def wipeout(self, eq):
        """Invalidate the specified equation and remove it from the
        factory.

        This will remove the equation from the purview of the factory
        and also change its formula to return NaN.  This ensures that eq
        does not observe any object in the factory and thus prevents its
        indirect pickling with the factory because of observer callback
        function.

        Parameters
        ----------
        eq : Equation or None
            The equation to invalidate. Nothing happens if this is None.
        """
        if eq is None:
            assert eq not in self.equations
            return
        self.equations.discard(eq)
        # invalidate this equation to clean up any observer relations of
        # objects in the factory towards its literals tree.
        nan = literals.Argument("nan", value=numpy.nan, const=True)
        eq.set_root(nan)
        return

    @deprecated(makeEquation_dep_msg)
    def makeEquation(
        self, eqstr, buildargs=True, argclass=literals.Argument, argkw={}
    ):
        """This function has been deprecated and will be removed in
        version 4.0.0.

        Please use
        diffpy.srfit.equation.builder.EquationFactory.make_equation
        instead.
        """
        return self.make_equation(eqstr, buildargs, argclass, argkw)

    @deprecated(registerConstant_dep_msg)
    def registerConstant(self, name, value):
        """This function has been deprecated and will be removed in
        version 4.0.0.

        Please use
        diffpy.srfit.equation.builder.EquationFactory.register_constant
        instead.
        """
        return self.register_constant(name, value)

    @deprecated(registerArgument_dep_msg)
    def registerArgument(self, name, arg):
        """This function has been deprecated and will be removed in
        version 4.0.0.

        Please use
        diffpy.srfit.equation.builder.EquationFactory.register_argument
        instead.
        """
        return self.register_argument(name, arg)

    @deprecated(registerOperator_dep_msg)
    def registerOperator(self, name, op):
        """This function has been deprecated and will be removed in
        version 4.0.0.

        Please use
        diffpy.srfit.equation.builder.EquationFactory.register_operator
        instead.
        """
        return self.register_operator(name, op)

    @deprecated(registerFunction_dep_msg)
    def registerFunction(self, name, func, argnames):
        """This function has been deprecated and will be removed in
        version 4.0.0.

        Please use
        diffpy.srfit.equation.builder.EquationFactory.register_function
        instead.
        """
        return self.register_function(name, func, argnames)

    @deprecated(registerBuilder_dep_msg)
    def registerBuilder(self, name, builder):
        """This function has been deprecated and will be removed in
        version 4.0.0.

        Please use
        diffpy.srfit.equation.builder.EquationFactory.register_builder
        instead.
        """
        return self.register_builder(name, builder)

    @deprecated(deRegisterBuilder_dep_msg)
    def deRegisterBuilder(self, name):
        """This function has been deprecated and will be removed in
        version 4.0.0.

        Please use
        diffpy.srfit.equation.builder.EquationFactory.deregister_builder
        instead.
        """
        return self.deregister_builder(name)

    def _prepare_builders(self, eqstr, buildargs, argclass, argkw):
        """Prepare builders so that equation string can be evaluated.

        This method checks the equation string for errors and missing
        arguments, and creates new arguments if allowed. In the process it
        rebuilds the newargs attribute.

        Parameters
        ----------
        eqstr : str
            The equation string, as in ``make_equation``.
        buildargs : bool
            The flag that is True to let the factory create missing
            arguments.
        argclass : type
            The class used to create new arguments.
        argkw : dict
            The keyword arguments passed to the ``argclass`` constructor.

        Raises
        ------
        ValueError
            If new arguments must be created but ``buildargs`` is False.
        SyntaxError
            If the equation string uses invalid syntax.
        """
        eqargs = self._get_undefined_args(eqstr)

        # Raise an error if there are arguments that need to be created, but
        # this is disallowed.
        if not buildargs and eqargs:
            eqargsstr = ", ".join(eqargs)
            msg = "The equation contains undefined arguments: %s" % eqargsstr
            raise ValueError(msg)

        # Make the arguments
        newargs = set()
        for argname in eqargs:
            arg = argclass(name=argname, **argkw)
            argbuilder = ArgumentBuilder(name=argname, arg=arg)
            newargs.add(arg)
            self.register_builder(argname, argbuilder)

        self.newargs = newargs

        return

    def _get_undefined_args(self, eqstr):
        """Get the undefined arguments from eqstr.

        This tokenizes eqstr and extracts undefined arguments. An
        undefined argument is defined as any token that is not a special
        character that does not correspond to a builder.

        Raises SyntaxError if the equation string uses invalid syntax.
        """
        interface = io.StringIO(eqstr).readline
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
            m = "invalid syntax: '%s'" % eqstr
            raise SyntaxError(m)

        # Scan the tokens for names that do not correspond to registered
        # builders. These will be treated as arguments that need to be
        # generated.
        for tok in set(args):
            # Move genuine variables to the eqargs dictionary
            if (
                # Check registered builders
                tok in self.builders
                or
                # Check symbols
                tok in EquationFactory.symbols
                or
                # Check ignored characters
                tok in EquationFactory.ignore
            ):
                args.remove(tok)

        return args


# End class EquationFactory

base_basebuilder = "diffpy.srfit.equation.builder.BaseBuilder"
removal_version = "4.0.0"

getequation_dep_msg = build_deprecation_message(
    base_basebuilder,
    "getEquation",
    "get_equation",
    removal_version,
)


class BaseBuilder(object):
    """Class for building equations.

    Equation builder objects can be composed like a normal function where the
    arguments can be other BaseBuilder instances or constants.

    Attributes
    ----------
    literal : Literal or None
        The equation Literal being built by this instance.
    """

    def __init__(self):
        self.literal = None
        return

    def __call__(self, *args):
        """Raise a TypeError, because a BaseBuilder takes no
        arguments."""
        m = "%s (%s) cannot accept arguments" % (
            self.literal.name,
            self.__class__.__name__,
        )
        raise TypeError(m)

    def get_equation(self):
        """Get the equation built by this object.

        The equation will given the name "_eq_<root>" where "<root>" is
        the name of the root node.

        Returns
        -------
        Equation
            The equation whose root is this builder's literal.
        """
        # We need to make a name for this, so we name it after its root
        name = "_eq_%s" % self.literal.name
        eq = Equation(name, self.literal)
        return eq

    @deprecated(getequation_dep_msg)
    def getEquation(self):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use diffpy.srfit.equation.builder.BaseBuilder.get_equation
        instead.
        """
        return self.get_equation()

    def __eval_binary(self, other, OperatorClass, onleft=True):
        """Evaluate a binary function.

        Other can be an BaseBuilder or a constant.

        Parameters
        ----------
        onleft
            Indicates that the operator was passed on the left side
            (default True).
        """
        # Create the Operator
        op = OperatorClass()

        # onleft takes care of non-commutative operators, and assures that the
        # ordering is preserved.
        if onleft:
            # Add the literals to the operator
            op.add_literal(self.literal)

        # Deal with constants
        if isinstance(other, BaseBuilder):
            literal = other.literal
        elif isinstance(other, Literal):
            literal = other
        else:
            literal = literals.Argument(value=other, const=True)
        op.add_literal(literal)

        if not onleft:
            # Add the literals to the operator
            op.add_literal(self.literal)

        # Create a new OperatorBuilder for the Operator
        opbuilder = OperatorBuilder(op.name)
        opbuilder.literal = op
        return opbuilder

    def __eval_unary(self, OperatorClass):
        """Evaluate a unary function."""
        op = OperatorClass()
        op.add_literal(self.literal)
        opbuilder = OperatorBuilder(op.name)
        opbuilder.literal = op
        return opbuilder

    def __add__(self, other):
        return self.__eval_binary(other, literals.AdditionOperator)

    def __radd__(self, other):
        return self.__eval_binary(other, literals.AdditionOperator, False)

    def __sub__(self, other):
        return self.__eval_binary(other, literals.SubtractionOperator)

    def __rsub__(self, other):
        return self.__eval_binary(other, literals.SubtractionOperator, False)

    def __mul__(self, other):
        return self.__eval_binary(other, literals.MultiplicationOperator)

    def __rmul__(self, other):
        return self.__eval_binary(
            other, literals.MultiplicationOperator, False
        )

    def __truediv__(self, other):
        return self.__eval_binary(other, literals.DivisionOperator)

    def __rtruediv__(self, other):
        return self.__eval_binary(other, literals.DivisionOperator, False)

    def __pow__(self, other):
        return self.__eval_binary(other, literals.ExponentiationOperator)

    def __rpow__(self, other):
        return self.__eval_binary(
            other, literals.ExponentiationOperator, False
        )

    def __mod__(self, other):
        return self.__eval_binary(other, literals.RemainderOperator)

    def __rmod__(self, other):
        return self.__eval_binary(other, literals.RemainderOperator, False)

    def __neg__(self):
        return self.__eval_unary(literals.NegationOperator)


# These are used by the class.


class ArgumentBuilder(BaseBuilder):
    """BaseBuilder wrapper around an Argument literal.

    Equation builder objects can be composed like a normal function where the
    arguments can be other BaseBuilder instances or constants.

    Parameters
    ----------
    value : object, optional
        The value of the new Argument. Default is None.
    name : str, optional
        The name of the new Argument. Default is None.
    const : bool, optional
        The flag that is True if the new Argument is constant. Default is
        False.
    arg : Argument, optional
        The existing Argument to wrap. If None (default), a new Argument
        is created from ``value``, ``name`` and ``const``.

    Attributes
    ----------
    literal : Argument
        The Argument wrapped by this instance.
    """

    def __init__(self, value=None, name=None, const=False, arg=None):
        BaseBuilder.__init__(self)
        if arg is None:
            self.literal = literals.Argument(
                value=value, name=name, const=const
            )
        else:
            self.literal = arg
        return


# end class ArgumentBuilder


class OperatorBuilder(BaseBuilder):
    """BaseBuilder wrapper around an Operator literal.

    Parameters
    ----------
    name : str
        The name of the wrapped Operator.
    op : Operator, optional
        The Operator to wrap. If None (default), ``name`` is taken to be
        the name of a numpy ufunc, which defines the Operator when the
        builder is called.

    Attributes
    ----------
    literal : Operator or None
        The Operator wrapped by this instance.
    name : str
        The name of the wrapped Operator.
    """

    def __init__(self, name, op=None):
        BaseBuilder.__init__(self)
        self.name = name
        self.literal = op
        return

    def __call__(self, *args):
        """Call the operator builder.

        This creates a new builder that encapsulates the operation.

        Parameters
        ----------
        *args
            The arguments of the operation.

        Returns
        -------
        OperatorBuilder
            The new builder for the operation applied to ``args``.

        Raises
        ------
        ValueError
            If the operator takes a fixed number of arguments and
            ``args`` has a different length.
        """
        newobj = OperatorBuilder(self.name)

        # If all we have is a name, then we assume that it is the name of a
        # numpy operator, and we use the corresponding Operator.
        if self.literal is None:
            ufunc = getattr(numpy, self.name)
            self.literal = literals.UFuncOperator(ufunc)
        # Here the Operator is already specified.  We can copy its attributes
        # to a new Operator inside of the new OperatorBuilder.
        op = literals.make_operator(
            name=self.literal.name,
            symbol=self.literal.symbol,
            nin=self.literal.nin,
            nout=self.literal.nout,
            operation=self.literal.operation,
        )
        newobj.literal = op

        # Now that we have a literal, let's check our inputs
        literal = newobj.literal
        if literal.nin >= 0 and len(args) != literal.nin:
            raise ValueError(
                "%s takes %i arguments (%i given)"
                % (self.literal, self.literal.nin, len(args))
            )

        # Wrap scalar arguments
        for i, arg in enumerate(args):
            # Wrap the argument if it is not already
            if not isinstance(arg, BaseBuilder):
                name = self.name + "_%i" % i
                arg = ArgumentBuilder(value=arg, name=name, const=True)
            newobj.literal.add_literal(arg.literal)

        return newobj


# end class OperatorBuilder

# Utility functions


def wrap_argument(name, arg):
    """Wrap an Argument as a builder.

    Parameters
    ----------
    name : str
        The name of the argument. It is not used, because the builder
        takes its name from ``arg``.
    arg : Argument
        The argument to wrap.

    Returns
    -------
    ArgumentBuilder
        The builder that wraps ``arg``.
    """
    argbuilder = ArgumentBuilder(arg=arg)
    return argbuilder


def wrap_operator(name, op):
    """Wrap an Operator as a builder.

    Parameters
    ----------
    name : str
        The name of the builder.
    op : Operator
        The operator to wrap.

    Returns
    -------
    OperatorBuilder
        The builder that wraps ``op``.
    """
    opbuilder = OperatorBuilder(name, op)
    return opbuilder


def wrap_function(name, func, nin=2, nout=1):
    """Wrap a function in an OperatorBuilder instance.

    Parameters
    ----------
    name : str
        The name of the function.
    func : callable
        The callable python object.
    nin : int, optional
        The number of input arguments. Default is 2.
    nout : int, optional
        The number of return values. Default is 1.

    Returns
    -------
    OperatorBuilder
        The builder that wraps the function.
    """
    op = literals.make_operator(
        name=name, symbol=name, nin=nin, nout=nout, operation=func
    )

    # Create the OperatorBuilder
    opbuilder = OperatorBuilder(name, op)

    return opbuilder


def get_builder(name):
    """Get an operator builder from the global builders dictionary.

    Parameters
    ----------
    name : str
        The name of the builder, such as "sin" or "add".

    Returns
    -------
    OperatorBuilder
        The global builder registered under ``name``.

    Raises
    ------
    KeyError
        If no global builder has this name.
    """
    return _builders[name]


@deprecated(wrapArgument_dep_msg)
def wrapArgument(name, arg):
    """This function has been deprecated and will be removed in version
    4.0.0.

    Please use diffpy.srfit.equation.builder.wrap_argument instead.
    """
    return wrap_argument(name, arg)


@deprecated(wrapOperator_dep_msg)
def wrapOperator(name, op):
    """This function has been deprecated and will be removed in version
    4.0.0.

    Please use diffpy.srfit.equation.builder.wrap_operator instead.
    """
    return wrap_operator(name, op)


@deprecated(wrapFunction_dep_msg)
def wrapFunction(name, func, nin=2, nout=1):
    """This function has been deprecated and will be removed in version
    4.0.0.

    Please use diffpy.srfit.equation.builder.wrap_function instead.
    """
    return wrap_function(name, func, nin, nout)


@deprecated(getBuilder_dep_msg)
def getBuilder(name):
    """This function has been deprecated and will be removed in version
    4.0.0.

    Please use diffpy.srfit.equation.builder.get_builder instead.
    """
    return get_builder(name)


def __wrap_numpy_operators():
    """Export all numpy operators as OperatorBuilder instances in the
    module namespace."""
    for name in dir(numpy):
        op = getattr(numpy, name)
        if isinstance(op, numpy.ufunc):
            _builders[name] = OperatorBuilder(name)
    return


__wrap_numpy_operators()


# Register other functions as well
def __wrap_srfit_operators():
    """Export all non-base operators from the
    diffpy.srfit.equation.literals.operators module as OperatorBuilder
    instances in the module namespace."""
    opmod = literals.operators
    excluded_types = set((opmod.CustomOperator, opmod.UFuncOperator))

    # check if opmod member should be wrapped as OperatorBuilder
    def _is_exported_type(cls):
        return (
            inspect.isclass(cls)
            and issubclass(cls, opmod.Operator)
            and not inspect.isabstract(cls)
            and cls not in excluded_types
        )

    # create OperatorBuilder objects
    for nm, opclass in inspect.getmembers(opmod, _is_exported_type):
        op = opclass()
        assert op.name, "Unnamed Operator should never appear here."
        _builders[op.name] = OperatorBuilder(op.name, op)
    return


__wrap_srfit_operators()

# End of file
