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
"""Classes and utilities for creating Equations objects.

The makeEquation function turns the string-representation of an equation into
an Equation instance.  For example
> eq = makeEquation("A*sin(a*x)")
will create an Equation with Arguments A, a and x that evaluates as
"A*sin(a*x)". 

Scalar constants can be included in the equation:
> consts = {"offset" : 3 }
> eq = makeEquation("A*sin(a*x) + offset", consts=consts)
This includes a constant offset in the equation. Similarly, one could have
written
> eq = makeEquation("A*sin(a*x) + 3")
to get the same effect.

The EquationBuilder class is at the core of makeEquation, and it can be used
directly to create Equations. EquationBuilder instances overload normal
arithmetic operations so that they build an Equation object instead.
EquationBuilder is specified in the ArgumentBuilder, GeneratorBuilder,
OperatorBuilder and PartitionBuilder classes.  All of the numpy ufunc operators
are overloaded within this module as OperatorBuilder instances. You can
register new builder objects with the 'registerBuilder' method, which will
allow those builders to be used by name in the makeEquation method.

With a collection of EquationBuilder objects, one can simply write
the Equation using normal python syntax:
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
otherwise, must be wrapped as ArgumentBuilders in order to be used in
this way.

Both makeEquation and EquationBuilder can make use of user-defined functions.
Any callable python object can be wrapped as an OperatorBuilder with the
wrapFunction function. The wrapped function is registered with this module so it
is usable in other imports and by the makeEquationFunction. For example.
> _f = lambda a, b : (a-b)/(a+b)
> f = wrapFunction("f", _f)
> # The builder, "f", which is a wrapper around _f is now registered with this
> # module as a custom builder.
> # Using makeEquation
> eq = makeEquation("c*f(a,b)")
> # gives an Equation with c, a and b as arguments.
> # Using EquationBuilder
> a = ArgumentBuilder(name = "a")
> b = ArgumentBuilder(name = "b")
> c = ArgumentBuilder(name = "c")
> beq = c*f(a,b)
> eq = beq.makeEquation()

Tags can be passed to Operators using both makeEquation and EquationBuilder
(see the diffpy.srfit.equation.literals.operator module).  When an
OperatorBuilder is called with extra arguments, these arguments are interpreted
as tags. For example, in the above example, the last two lines could be
replaced by
> beq = c*f(a,b,"tag1")
> eq = beq.makeEquation()
and the 'f' function would only operate on Partitions Arguments that have a
"tag1" tag.

Also included in this module is the 'custombuilders' dictionary of wrapped
EquationBuilder objects. 
"""

# NOTE - the builder cannot handle numpy arrays on the left of a binary
# operation because the array will automatically loop the operator of the
# right-side over its arguments. This results in an array of EquationBuilder
# instances, not an EquationBuilder that contains an array.


import sys
import numpy

from .Equation import Equation
import diffpy.srfit.equation.literals as literals

custombuilders = {}

def makeEquation(eqstr, consts = {}):
    """Make an equation from an equation string.

    Arguments
    eqstr   --  An equation in string form using standard python syntax. The
                equation string can use any function that exists as a numpy
                ufunc or that has been created with the wrapFunction method. In
                the latter case, the name of a custom function used in the
                equation string must match the name used in the wrapFunction
                method.
    consts  --  A dictionary of named scalar constants used in the equation.

    Returns an Equation instance representing the equation string.
    """
    ns = _makeNamespace(eqstr, consts)
    beq = eval(eqstr, ns)
    return beq.getEquation()

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
            # Check custom builders
            tok not in custombuilders and
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
    # Get the operators, partitions and generators
    for opname in eqops.values():
        #print opname
        if opname in custombuilders:
            opbuilder = custombuilders[opname]
        else:
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

    Attributes
    literal     --  The root of the Equation being built by this instance.
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
    """EquationBuilder wrapper around an Argument literal.

    Equation builder objects can be composed like a normal function where the
    arguments can be other EquationBuilder instances or constants.

    Attributes
    literal     --  The Argument wrapped by this instance.
    """

    def __init__(self, value=None, name=None, const=False):
        """Create an ArgumentBuilder instance, containing a new Argument.

        Arguments
        value   --  The value of the wrapped Argument (float, default None)
        name    --  The name of the wrapped Argument (string, default None)
        const   --  Flag indicating whether the Argument is constant (bool,
                    default False)
        """
        EquationBuilder.__init__(self)
        self.literal = literals.Argument(value=value, name=name, const=const)
        return

# end class ArgumentBuilder

class PartitionBuilder(EquationBuilder):
    """EquationBuilder wrapper around a Partition literal.

    Equation builder objects can be composed like a normal function where the
    arguments can be other EquationBuilder instances or constants.

    Attributes
    literal     --  The Partition wrapped by this instance.
    """

    def __init__(self, part):
        """Wrap the Partition.
        
        Arguments
        part    --  The Partition instance to be wrapped.
        """
        EquationBuilder.__init__(self)
        self.literal = part
        return

# end class PartitionBuilder

class GeneratorBuilder(EquationBuilder):
    """EquationBuilder wrapper around a Generator literal.

    Equation builder objects can be composed like a normal function where the
    arguments can be other EquationBuilder instances or constants.

    Attributes
    literal     --  The Generator wrapped by this instance.
    """

    def __init__(self, gen):
        """Wrap the Generator.
        
        Arguments
        gen --  The Generator instance to be wrapped.
        """
        EquationBuilder.__init__(self)
        self.literal = gen
        return

# end class GeneratorBuilder

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
        
        args    --  Arguments of the operation. Args past nin are considered
                    tags.
        kw      --  Key-word arguments. If "combine" appears and is True, the
                    operation is allowed to combine any partitions it
                    encounters after the operation.
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
            op = literals.Operator()
            op.name = self.literal.name
            op.symbol = self.literal.name
            op.nin = self.literal.nin
            op.nout = self.literal.nout
            op.operation = self.literal.operation
            newobj.literal = op

        combine = kw.get("combine", False)
        # Wrap scalar arguments and process tags
        for arg in args[:newobj.literal.nin]:
            # Wrap the argument if it is not already
            if not isinstance(arg, EquationBuilder):
                arg = ArgumentBuilder(value=arg, const=True)
            newobj.literal.addLiteral(arg.literal)
        newobj.literal.addTags(*args[newobj.literal.nin:])
        newobj.literal.setCombine(combine)
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

    # Register it with the module
    registerBuilder(name, opbuilder)

    # Return it
    return opbuilder

def wrapPartition(name, part):
    """Wrap a Partition as a builder and register it with the module.

    Returns the wrapped Partition.
    """
    part.name = name
    partbuilder = PartitionBuilder(part)
    registerBuilder(name, partbuilder)
    return partbuilder

def wrapGenerator(name, gen):
    """Wrap a Generator as a builder and register it with the module.
    
    Returns the wrapped Generator.
    """
    gen.name = name
    genbuilder = GeneratorBuilder(gen)
    registerBuilder(name, genbuilder)
    return genbuilder

def registerBuilder(name, builder):
    """Register builder in this module so it can be used in makeEquation."""
    global custom
    custombuilders[name] = builder

# Export all numpy operators as OperatorBuilder instances in the module
# namespace.
for name in dir(numpy):
    op = getattr(numpy, name)
    if isinstance(op, numpy.ufunc):
        setattr(sys.modules[__name__], name, OperatorBuilder(name))

# version
__id__ = "$Id$"

#
# End of file
