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

"""Swapper for replacing a Literal in an equation with another Literals.
"""

__all__ = ["Swapper"]

from diffpy.srfit.equation.visitors.visitor import Visitor

class Swapper(Visitor):
    """Swapper for swapping out one literal for another in a literal tree.

    Note that this cannot swap out a root node of a literal tree. This case
    must be tested for explicitly.

    Attributes:
    newlit  --  The literal to be placed into the literal tree.
    oldlit  --  The literal to be replaced.

    """

    def __init__(self, oldlit, newlit):
        """Initialize.

        oldlit  --  The literal to be replaced.
        newlit  --  The literal to be placed into the literal tree. See the
                    class for how the replacement takes place.

        """

        self.newlit = newlit
        self.oldlit = oldlit

        self._swap = False

        return

    def onArgument(self, arg):
        """Process an Argument node.

        Tell the parent to swap the old Argument with the replacement Literal.

        """

        if arg is self.oldlit:
            self._swap = True

        return

    def onOperator(self, op):
        """Process an Operator node.

        Tell the parent to swap the old Operator with the replacement Literal.

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
                    oldlit.removeObserver(op._flush)
                except KeyError:
                    pass

                # Validate the new literal. If it fails, we need to restore the
                # old one
                try:
                    op._loopCheck(newlit)
                except ValueError:
                    # Restore the old literal
                    op.args.insert(idx, oldlit)
                    oldlit.addObserver(op._flush)
                    raise

                # If we got here, then go on with replacing the literal
                op.args.insert(idx, newlit)
                newlit.addObserver(op._flush)
                op._flush(other=())


            self._swap = False

        return

    def onEquation(self, eq):
        """Process an Equation node.

        This looks at the equation itself as well as the root.

        """
        if eq is self.oldlit:
            self._swap = True
            return

        # If the newlit is the root, then swap that out and move on.
        if eq.root is self.oldlit:
            eq.setRoot(self.newlit)
            return

        # Now move into the equation. We have to do a _loopCheck to make sure
        # that we won't have any loops in the equation.
        eq._loopCheck(self.newlit)
        eq.root.identify(self)

        # Reset the root in case anything changed underneath.
        eq.setRoot(eq.root)

        return

# End of file
