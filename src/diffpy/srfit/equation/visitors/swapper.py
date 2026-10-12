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
"""Swapper for replacing a Literal in an equation with another
Literals."""

__all__ = ["Swapper"]

from diffpy.srfit.equation.visitors.visitor import Visitor


class Swapper(Visitor):
    """Swapper for swapping out one literal for another in a literal
    tree.

    Note that this cannot swap out a root node of a literal tree. This case
    must be tested for explicitly.

    Parameters
    ----------
    oldlit : Literal
        The literal to be replaced.
    newlit : Literal
        The literal to be placed into the literal tree.

    Attributes
    ----------
    newlit : Literal
        The literal to be placed into the literal tree.
    oldlit : Literal
        The literal to be replaced.
    """

    def __init__(self, oldlit, newlit):
        self.newlit = newlit
        self.oldlit = oldlit

        self._swap = False

        return

    def on_argument(self, arg):
        """Process an Argument node.

        Tell the parent to swap the old Argument with the replacement
        Literal.

        Parameters
        ----------
        arg : Argument
            The Argument to check.
        """
        if arg is self.oldlit:
            self._swap = True

        return

    def on_operator(self, op):
        """Process an Operator node.

        Tell the parent to swap the old Operator with the replacement
        Literal, and swap any of its arguments in place.

        Parameters
        ----------
        op : Operator
            The Operator to check and traverse.

        Raises
        ------
        ValueError
            If the replacement Literal causes a self-reference.
        """
        # Check to see if we need to swap out this Operator. If so, then we
        # don't need to traverse into its arguments.
        if op is self.oldlit:
            self._swap = True
            return

        # Now traverse into the arguments.
        for literal in op.args:
            literal.identify(self)

        # If we've been told to swap out a literal, then we must do it in-place
        # because the order of op.args matters.
        if self._swap:

            oldlit = self.oldlit
            newlit = self.newlit

            while oldlit in op.args:

                # Record the index
                idx = op.args.index(oldlit)
                # Remove the literal
                del op.args[idx]
                # Remove op as an observer. A KeyError will be raised if we
                # attempt to remove the same observer more than once, which
                # might happen if the oldlit appears multiple times in op.args.
                try:
                    oldlit.remove_observer(op._flush)
                except KeyError:
                    pass

                # Validate the new literal. If it fails, we need to restore the
                # old one
                try:
                    op._loop_check(newlit)
                except ValueError:
                    # Restore the old literal
                    op.args.insert(idx, oldlit)
                    oldlit.add_observer(op._flush)
                    raise

                # If we got here, then go on with replacing the literal
                op.args.insert(idx, newlit)
                newlit.add_observer(op._flush)
                op._flush(other=())

            self._swap = False

        return

    def on_equation(self, eq):
        """Process an Equation node.

        This looks at the equation itself as well as the root.

        Parameters
        ----------
        eq : Equation
            The Equation to check and traverse.

        Raises
        ------
        ValueError
            If the replacement Literal causes a self-reference.
        """
        if eq is self.oldlit:
            self._swap = True
            return

        # If the newlit is the root, then swap that out and move on.
        if eq.root is self.oldlit:
            eq.set_root(self.newlit)
            return

        # Now move into the equation. We have to do a _loop_check to make sure
        # that we won't have any loops in the equation.
        eq._loop_check(self.newlit)
        eq.root.identify(self)

        # Reset the root in case anything changed underneath.
        eq.set_root(eq.root)

        return


# End of file
