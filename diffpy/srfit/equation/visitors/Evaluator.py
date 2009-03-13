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

from .Visitor import Visitor
from ..literals import Partition
from ..literals import Argument

from .. import Clicker

class Evaluator(Visitor):
    """Evaluator visitor for computing the value of an expression.
    
    Arguments
    value   --  The value of the expression
    clicker --  A reference
    """

    def __init__(self):
        """Initialize."""
        self.value = 0
        self._clicker = Clicker()

        # Are we at the root?
        self._atroot = True
        return

    def click(self):
        """Click the clicker and reset non-public data."""
        self._atroot = True
        self._clicker.click()
        return

    def onArgument(self, arg):
        """Process an Argument node."""
        return

    def onPartition(self, part):
        """Process a Partition node."""
        return

    def onOperator(self, op):
        """Process an Operator node.

        This version deals with a single partition.
        """
        if op.clicker > self._clicker:

            # This evaluates the operator and creates an Argument or Partition
            # that can be used to get its value.
            self._createProxy(op)

        # Doesn't do anything right now
        # op._proxy.identify(self)
        return

    def _createProxy(self, op):
        """Create a proxy that can be called in place of an Operator.

        The type of proxy depends on the arguments of the Operator. It uses the
        Helper class for actually evaluating the operation.

        Partition Logic
        * If two or more arguments of the Operator are Partitions, then each
          will have its combine method called.
        * If there is a single Partition argument, then the proxy will be
          Partition-like
        * Otherwise it will be Argument-like.

        Returns a new Literal
        """

        # This tells us whether we can combine after the operation and
        # self._combine will be propogated to the next level of the tree. In
        # this way, combine will only be true until the first tagged operator is
        # encounterd.
        atroot = self._atroot
        self._atroot = False

        # Evaluation helper exclusive to this operation.
        helper = Evaluator.Helper()
        helper.op = op

        for arg in op.args:
            arg.identify(self)
            arg.identify(helper)

        helper.evaluate()

        # Combine if we can. This has no effect if the helper evaluation didn't
        # result in a Partition.
        if atroot or op._cancombine:
            helper.combine()

        # Make the proxy
        if helper.part is not None:
            if op._proxy is None:
                op._proxy = Evaluator.PseudoPartition()
            op._proxy.combine = helper.part.combine
            op._proxy.tags = helper.part.tags
            op._proxy.tagmap = helper.part.tagmap
            op._proxy._partvals = helper.partvals
        else:
            if op._proxy is None:
                op._proxy = Argument()
            op._proxy.setValue(helper.value)
            self.value = helper.value

        return

    class Helper(object):
        """Class for evaluating an operation that might involve Partitions.

        The purpose of the helper class is to evaluate the operation on its
        arguments. On first glance, it seems redundant to have a class for this
        when the Evaluator can do this job. However, the Evaluator's roll is to
        traverse the entire tree where as the Helper stays at one level of the
        tree. This allows the use of instance attributes for managing the data
        necessary for the operation without managing the attributes themselves.
        This makes the code much cleaner and logically encapsulated. Think of
        this class as an "operator helper" rather than an "evaluator helper".
        """

        def __init__(self):
            self.op = None
            # Literals of the operation
            self.args = []
            self.argvals = []
            self.partidx = []
            # Partition and Argument data
            self.part = None
            self.partvals = []
            self.value = None
            return

        # Methods for collecting and evaluating arguments

        def onArgument(self, arg):
            """Record an Argument."""
            self.args.append(arg)
            self.argvals.append(arg.value)
            return

        def onPartition(self, part):
            """Record a Partition."""
            part._prepare()
            self.args.append(part)
            self.argvals.append(None) # A place-holder
            # Record the location of the partition
            self.partidx.append(len(self.args)-1)
            return

        def onOperator(self, op):
            """Process an Operator through its proxy."""
            op._proxy.identify(self)
            return

        def evaluate(self):
            """Evaluate the operation."""
            # 1. call 'combine' on Partitions, if there is more than one.
            # 2. Evaluate the operator on the converted arguments.

            # 1. Transform arguments
            numparts = len(self.partidx)
            if numparts == 1:
                # Exactly one Partition
                self.evaluatePartition()
                return

            if numparts > 1:
                # More than one Partition, combine
                for idx in self.partidx:
                    part = self.args[idx]
                    self.argvals[idx] = part.combine(part._partvals)

            self.evaluateArgument()
            return

        def evaluateArgument(self):
            """Evaluate with Argument arguments."""
            self.value = self.op.operation(*self.argvals)
            return

        def evaluatePartition(self):
            """Evaluate with a partition as an argument."""

            # Get the partition values we're working with
            pidx = self.partidx[0]
            self.part = self.args[pidx]
            self.partvals = self.part._partvals[:]

            # Create the set of indices that are consistent with the operator's
            # tags
            idxlist = set()
            for tag in self.op.tags:
                idxlist.update( self.part.tagmap.get(tag, []) )
            if not idxlist:
                idxlist = range(len(self.partvals))

            # Now evaluate the function on each part of the partition with the
            # proper tags.
            for idx in idxlist:
                self.argvals[pidx] = self.partvals[idx]
                self.partvals[idx] = self.op.operation(*self.argvals)

            return

        def combine(self):
            """Combine if we have a Partition to combine."""
            if not self.part: 
                return

            self.value = self.part.combine(self.partvals)
            self.part = None
            return

    # End class Helper

    class PseudoPartition(Partition):
        """A partition-like object that can be used to hold the output of an
        operated partition.
        """

        def _prepare(self):
            """Prepare this Partition...which means do nothing.

            The _partvals list is set by the onOperator method, so this method
            does nothing.
            """
            return

    # End class PseudoPartition

# End class Evaluator

# version
__id__ = "$Id$"

#
# End of file
