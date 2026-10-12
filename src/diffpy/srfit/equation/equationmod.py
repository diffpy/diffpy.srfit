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

Equation is a functor that holds a Literal tree that defines an
equation. Its __call__ method evaluates the equation at the most recent
value of its Arguments. The non-constant arguments are accessible as
attributes of the Equation instance and can be passed as arguments to
__call__.

Example::

    # make a Literal tree. Here's a simple one
    add = AdditionOperator()
    a = Argument(name="a")  # Don't forget to name these!
    b = Argument(name="b")
    add.add_literal(a)
    add.add_literal(b)
    # make an Equation instance and pass the root
    eq = Equation(root=add)
    eq(a=3, b=4)  # returns 7
    eq(a=2)  # remembers b=4, returns 6
    eq.a.set_value(-3)
    eq.b.set_value(3)
    eq()  # uses last assignment of a and b, returns 0

See the class documentation for more information.
"""

__all__ = ["Equation"]

from collections import OrderedDict

from diffpy.srfit.equation.literals.literal import Literal
from diffpy.srfit.equation.literals.operators import Operator
from diffpy.srfit.equation.visitors import get_args, swap, validate
from diffpy.utils._deprecator import build_deprecation_message, deprecated

equation_base = "diffpy.srfit.equation.Equation"
removal_version = "4.0.0"
setRoot_dep_msg = build_deprecation_message(
    equation_base, "setRoot", "set_root", removal_version
)


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

    Parameters
    ----------
    name : str, optional
        The name of the Equation. If None (default) and ``root`` is
        given, the name is "eq_<root name>".
    root : Literal, optional
        The root node of the Literal tree. If None (default), call
        ``set_root`` to set it.

    Attributes
    ----------
    root : Literal or None
        The root Literal of the equation tree.
    argdict : OrderedDict
        The non-constant Arguments of the tree, indexed by name.
    args : list
        The values of ``argdict``, in the same order.
    name : str
        The name of the Equation.
    nin : int
        The number of non-constant Arguments.
    nout : int
        The number of outputs, which is always 1.
    symbol : str
        The symbolic representation, which is the same as ``name``.
    value : object
        The value of the Equation, read with ``get_value``.
    """

    # define abstract attributes from the Operator base.
    nin = None
    nout = 1

    def __init__(self, name=None, root=None):
        # Operator stuff. We circumvent Operator.__init__ since we're using
        # args as a property. We cannot set it, as the Operator tries to do.
        if name is None and root is not None:
            name = "eq_%s" % root.name
        Literal.__init__(self, name)
        self.root = None
        self.argdict = OrderedDict()
        if root is not None:
            self.set_root(root)
        return

    @property
    def symbol(self):
        return self.name

    def operation(self, *args, **kw):
        """Evaluate this Equation object.

        Same as the __call__ method.  This method is used via the
        Operator interface.

        Parameters
        ----------
        *args
            The new values of the arguments, in the order of ``args``.
        **kw
            The new values of the arguments, by name.

        Returns
        -------
        object
            The value of the Equation.
        """
        return self.__call__(*args, **kw)

    def _get_args(self):
        return list(self.argdict.values())

    args = property(_get_args)

    def __getattr__(self, name):
        """Gives access to the Arguments as attributes."""
        # Avoid infinite loop on argdict lookup.
        argdict = object.__getattribute__(self, "argdict")
        if name not in argdict:
            raise AttributeError("No argument named '%s' here" % name)
        return argdict[name]

    # Ensure there is no __dir__ override in the base class.
    assert getattr(Operator, "__dir__", None) is getattr(
        object, "__dir__", None
    )

    def __dir__(self):
        """Return sorted list of attributes for this object."""
        rv = set(dir(type(self)))
        rv.update(self.__dict__, self.argdict)
        rv = sorted(rv)
        return rv

    def set_root(self, root):
        """Set the root of the Literal tree.

        Parameters
        ----------
        root : Literal
            The new root of the Literal tree.

        Raises
        ------
        ValueError
            If errors are found in the Literal tree.
        """
        # Validate the new root
        validate(root)

        # Stop observing the old root
        if self.root is not None:
            self.root.remove_observer(self._flush)

        # Add the new root
        self.root = root
        self.root.add_observer(self._flush)
        self._flush(other=(self,))

        # Get the args
        args = get_args(root, getconsts=False)
        self.argdict = OrderedDict([(arg.name, arg) for arg in args])

        # Set Operator attributes
        self.nin = len(self.args)

        return

    @deprecated(setRoot_dep_msg)
    def setRoot(self, root):
        """This function has been deprecated and will be removed in
        version 4.0.0.

        Please use diffpy.srfit.equation.Equation.set_root instead.
        """
        return self.set_root(root)

    def __call__(self, *args, **kw):
        """Call the equation.

        New Argument values are accepted as arguments or keyword
        assignments (or both).  The order of accepted arguments is given
        by the args attribute. The Equation will remember values set in
        this way.

        Parameters
        ----------
        *args
            The new values of the arguments, in the order of ``args``.
        **kw
            The new values of the arguments, by name.

        Returns
        -------
        object
            The value of the Equation.

        Raises
        ------
        ValueError
            If there are more positional values than arguments, or a
            keyword does not name an argument.
        """
        # Process args
        for idx, val in enumerate(args):
            if idx >= len(self.argdict):
                raise ValueError("Too many arguments")
            arg = self.args[idx]
            arg.set_value(val)

        # Process kw
        for name, val in kw.items():
            arg = self.argdict.get(name)
            if arg is None:
                raise ValueError("No argument named '%s' here" % name)
            arg.set_value(val)

        self._value = self.root.get_value()
        return self._value

    def swap(self, oldlit, newlit):
        """Swap a literal in the equation for another.

        Note that this may change the root and the operation interface.

        Parameters
        ----------
        oldlit : Literal
            The literal to replace.
        newlit : Literal
            The literal to put in its place.

        Raises
        ------
        ValueError
            If ``newlit`` causes a self-reference.
        """
        newroot = swap(self.root, oldlit, newlit)
        self.set_root(newroot)
        return

    # Operator methods

    def add_literal(self, literal):
        """Raise a RuntimeError, because an Equation takes its literals
        from its root.

        Parameters
        ----------
        literal : Literal
            The literal that cannot be added.

        Raises
        ------
        RuntimeError
            Always.
        """
        raise RuntimeError("Cannot add literals to an Equation.")

    # Visitors can treat us differently than an Operator.

    def identify(self, visitor):
        """Identify this Equation to a visitor.

        Parameters
        ----------
        visitor : Visitor
            The visitor to identify to.

        Returns
        -------
        object
            The result of ``visitor.on_equation(self)``.
        """
        return visitor.on_equation(self)


# End of file
