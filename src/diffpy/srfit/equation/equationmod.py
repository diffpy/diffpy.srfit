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

"""The Equation class for holding and evaluating an equation.

Equation is a functor that holds a Literal tree that defines an equation. It's
__call__ method evaluates the equation at the most recent value of its
Arguments. The non-constant arguments are accessible as attributes of the
Equation instance and can be passed as arguments to __call__.

Example
> # make a Literal tree. Here's a simple one
> add = AdditionOperator()
> a = Argument(name="a") # Don't forget to name these!
> b = Argument(name="b")
> add.addLiteral(a)
> add.addLiteral(b)
> # make an Equation instance and pass the root
> eq = Equation(root = add)
> eq(a=3, b=4) # returns 7
> eq(a=2) # remembers b=4, returns 6
> eq.a.setValue(-3)
> eq.b.setValue(3)
> eq() # uses last assignment of a and b, returns 0

See the class documentation for more information.
"""

__all__ = ["Equation"]

from collections import OrderedDict

from diffpy.srfit.equation.visitors import validate, getArgs, swap
from diffpy.srfit.equation.literals.operators import Operator
from diffpy.srfit.equation.literals.literal import Literal

class Equation(Operator):
    """Class for holding and evaluating a Literal tree.

    Instances have attributes that are the non-const Arguments of the tree
    (accessed by name) and a __call__ method that evaluates the tree.  It is
    assumed, but not enforced that Arguments have unique names.  If this is not
    the case, then one should keep its own list of Arguments.

    The tree is scanned for errors when it is added. Thus, the tree should be
    complete before putting it inside an Equation.

    Equations can act as Operator nodes within a literal tree. In this context,
    they evaluate as the root node, but provide a calling interface that
    accepts new argument values for the literal tree.

    Attributes
    root    --  The root Literal of the equation tree
    argdict --  An OrderedDict of Arguments from the root.
    args    --  Property that gets the values of argdict.

    Operator Attributes
    args    --  List of Literal arguments, set with 'addLiteral'
    name    --  A name for this operator. e.g. "add" or "sin"
    nin     --  Number of inputs (<1 means this is variable)
    nout    --  Number of outputs
    operation   --  Function that performs the operation. e.g. numpy.add.
    symbol  --  The symbolic representation. e.g. "+" or "sin".
    _value  --  The value of the Operator.
    value   --  Property for 'getValue'.
    """

    # define abstract attributes from the Operator base.
    nin = None
    nout = 1

    def __init__(self, name = None, root = None):
        """Initialize.

        name    --  A name for this Equation.
        root    --  The root node of the Literal tree (default None). If root
                    is not passed here, you must call the 'setRoot' method to
                    set or change the root node.

        """
        # Operator stuff. We circumvent Operator.__init__ since we're using
        # args as a property. We cannot set it, as the Operator tries to do.
        if name is None and root is not None:
            name = "eq_%s"%root.name
        Literal.__init__(self, name)
        self.root = None
        self.argdict = OrderedDict()
        if root is not None:
            self.setRoot(root)
        return


    @property
    def symbol(self):
        return self.name


    def operation(self, *args, **kw):
        """Evaluate this Equation object.

        Same as the __call__ method.  This method is used via
        the Operator interface.

        Return the result of __call__(*args, **kw).
        """
        return self.__call__(*args, **kw)


    def _getArgs(self):
        return self.argdict.values()

    args = property(_getArgs)

    def __getattr__(self, name):
        """Gives access to the Arguments as attributes."""
        # Avoid infinite loop on argdict lookup.
        argdict = object.__getattribute__(self, 'argdict')
        if not name in argdict:
            raise AttributeError("No argument named '%s' here" % name)
        return argdict[name]


    # Ensure there is no __dir__ override in the base class.
    assert (getattr(Operator, '__dir__', None) is
            getattr(object, '__dir__', None))


    def __dir__(self):
        "Return sorted list of attributes for this object."
        rv = set(dir(type(self)))
        rv.update(self.__dict__, self.argdict)
        rv = sorted(rv)
        return rv


    def setRoot(self, root):
        """Set the root of the Literal tree.

        Raises:
        ValueError if errors are found in the Literal tree.

        """

        # Validate the new root
        validate(root)

        # Stop observing the old root
        if self.root is not None:
            self.root.removeObserver(self._flush)

        # Add the new root
        self.root = root
        self.root.addObserver(self._flush)
        self._flush(other=(self,))

        # Get the args
        args = getArgs(root, getconsts=False)
        self.argdict = OrderedDict( [(arg.name, arg) for arg in args] )

        # Set Operator attributes
        self.nin = len(self.args)

        return


    def __call__(self, *args, **kw):
        """Call the equation.

        New Argument values are acceped as arguments or keyword assignments (or
        both).  The order of accepted arguments is given by the args
        attribute.  The Equation will remember values set in this way.

        Raises
        ValueError when a passed argument cannot be found
        """
        # Process args
        for idx, val in enumerate(args):
            if idx >= len(self.argdict):
                raise ValueError("Too many arguments")
            arg = self.args[idx]
            arg.setValue(val)

        # Process kw
        for name, val in kw.items():
            arg = self.argdict.get(name)
            if arg is None:
                raise ValueError("No argument named '%s' here"%name)
            arg.setValue(val)

        self._value = self.root.getValue()
        return self._value

    def swap(self, oldlit, newlit):
        """Swap a literal in the equation for another.

        Note that this may change the root and the operation interface

        """
        newroot = swap(self.root, oldlit, newlit)
        self.setRoot(newroot)
        return

    # Operator methods

    def addLiteral(self, literal):
        """Cannot add a literal to an Equation."""
        raise RuntimeError("Cannot add literals to an Equation.")

    # Visitors can treat us differently than an Operator.

    def identify(self, visitor):
        """Identify self to a visitor."""
        return visitor.onEquation(self)

# End of file
