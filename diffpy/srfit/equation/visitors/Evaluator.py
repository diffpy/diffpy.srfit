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
"""Evaluator visitor for evaluating a tree of Literals.

The Evaluator walks an equation tree composed of Literals and computes the value
of the equation. It does this by computing value of each Operator node by
performing its operation on the values of its arguments.  An Evaluator has a
clicker that it compares to the clicker of each Operator node. When the
Operator's clicker is greater than the Evaluator's, the Evaluator will recompute
the value of the node. Otherwise, the most recent value of the node is used in
the evaluation. In this way, the Evaluator can quickly re-evaluate an equation
when few of the Arguments have changed since the last evaluation. To make use of
this feature the Evaluator's clicker must be clicked after it walks a tree. 
For example
> # make an equation called 'equation'...
> # make an Evaluator
> evaluator = Evaluator()
> equation.identify(evaluator)
> value = evaluator.value
> # Do something with value ...
> # Click evaluator clicker
> evaluator.clicker.click()

Evaluator does not check of the validity of the expression. See the Validator
visitor for that.
"""

from diffpy.srfit.equation.visitors.Visitor import Visitor

from diffpy.srfit.equation.literals import Clicker

class Evaluator(Visitor):
    """Evaluator visitor for computing the value of an expression.
    
    Arguments
    value   --  The value of the expression
    clicker --  A reference
    """

    def __init__(self):
        """Initialize."""
        self.value = 0
        self.clicker = Clicker()
        self.vals = []
        return

    def onArgument(self, arg):
        """Process an Argument node."""
        self.value = arg.value
        self.vals = [arg.value]
        return

    def onOperator(self, op):
        """Process an Operator node."""

        # We have to re-evaluate if the operator's clicker is greater than the
        # evaluator's. 
        if op.clicker > self.clicker:

            # Visit each literal in the operator's arg list and update its
            # value.
            localvals = []
            for literal in op.args:
                literal.identify(self)
                localvals.extend(self.vals)

            # Now evaluate the operation. This does not check for validity.
            op.value = op.operation(*localvals)

        self.value = op.value
        if(op.nout == 1):
            self.vals = [self.value]
        elif(op.nout > 1):
            self.vals = self.value

        return


# version
__id__ = "$Id$"

#
# End of file
