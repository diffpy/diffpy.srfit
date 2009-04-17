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
"""Example of using a partition."""

import numpy
from numpy import exp
from pylab import plot, show

from diffpy.srfit.equation.literals import Partition, Argument
from diffpy.srfit.equation.builder import EquationFactory

# The data
x = numpy.arange(0, 10, 0.1)

# Make three gaussians that will be held in the partition
pi = numpy.pi
g1 = exp(-((x-2)/0.4)**2) * ((2*pi)**0.5 * 0.4)**-1
g2 = exp(-((x-4)/0.3)**2) * ((2*pi)**0.5 * 0.3)**-1
g3 = exp(-((x-6)/0.2)**2) * ((2*pi)**0.5 * 0.2)**-1

# Make arguments holding these
a1 = Argument(name="g1", value=g1, const=True)
a2 = Argument(name="g2", value=g2, const=True)
a3 = Argument(name="g3", value=g3, const=True)

# Make the partition. It holds each gaussian separately, but can compose them
# as well.
gpart = Partition()
gpart.addArgument(a1, "A", "A-A")
gpart.addArgument(a2, "A", "B", "A-B", "B-A")
gpart.addArgument(a3, "B", "B-B")

# Make this partition known to the factory so we can use it in an equation.
factory = EquationFactory()
factory.registerPartition("gpart", gpart)

#### Now the fun begins ###

# Make the equation "A0 * gpart + c".
eq = factory.makeEquation("A0 * gpart + c")

y1 = eq(A0=1, c=0)
y2 = eq(A0=2, c=0)
plot(x, y1, x, y2)
show()
raw_input("Press anything")

# Make the same equation, but only modify the amplitude of the "B-B" pair
eq = factory.makeEquation("multiply(A0, gpart, 'B-B') + c")

y3 = eq(A0=2, c=0)
plot(x, y1, x, y2, x, y3)
show()
raw_input("Press anything")

# Do this again, but add an arbitrary constant to the "A" pairs
eq = factory.makeEquation("add( multiply(A0, gpart, 'B-B'), c, 'A')")

y4 = eq(A0=2, c=1)
plot(x, y1, x, y2, x, y3, x, y4)
show()
raw_input("Press anything")

# What went wrong?
# c=1 gets added to both g1 and g2 after the multiplication operation, even the
# zero entries. When the different parts of the partition are assembled, this
# offset appears twice, so it is doubled in the final result. Here's the result
# in normal notation
# answer = (g1 + c) + (g2 + c) + (A0 * g3)
#        = g1 + g2 + A0*g3 + 2*c
# We must think of partitions as the generalization of a numpy array.
# Operations are applied to each element.



# version
__id__ = "$Id$"

#
# End of file
