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
"""The Profile class containing the physical and calculated data.

Profile holds the arrays representing an observed profile, a selected subset of
the observed profile and a calculated profile. Profiles are used by Calculators
to store a calculated signal, and by Contributions to help calculate a residual
equation.

"""

import numpy

from .parameter import Parameter
from diffpy.srfit.equation import Clicker

# This is the roundoff tolerance for selecting bounds on arrays.
epsilon = 1e-8


class Profile(object):
    """Observed and calculated profile container.

    Attributes

    xobs    --  A numpy array of the observed independent variable (default
                None)
    yobs    --  A numpy array of the observed signal (default None)
    dyobs   --  A numpy array of the uncertainty of the observed signal (default
                None, optional).
    y       --  The profile over the calculation range (default None).
    dy      --  The uncertainty in the profile over the calculation range
                (default None).
    x       --  A numpy array of the calculated independent variable (default
                None)
    ycalc   --  A numpy array of the calculated signal (default None).
    xpar    --  A Parameter that stores x (named "x").
    ypar    --  A Parameter that stores y (named "y").
    dypar   --  A Parameter that stores dy (named "dy").

    """

    def __init__(self):
        """Initialize the attributes."""
        self.xobs = None
        self.yobs = None
        self.dyobs = None
        self.y = None
        self.dy = None
        self.x = None
        self.ycalc = None
        self.xpar = Parameter("x")
        self.ypar = Parameter("y")
        self.dypar = Parameter("dy")
        return

    def setObservedProfile(self, xobs, yobs, dyobs = None):
        """Set the observed profile.

        Arguments
        xobs    --  Numpy array of the independent variable
        yobs    --  Numpy array of the observed signal.
        dyobs   --  Numpy array of the uncertainty in the observed signal. If
                    dyobs is None (default), it will be set to 1 at each
                    observed xobs.

        Raises ValueError if len(yobs) != len(xobs)
        Raises ValueError if dyobs != None and len(dyobs) != len(xobs)

        """
        if len(yobs) != len(xobs):
            raise ValueError("xobs and yobs are different lengths")
        if dyobs is not None and len(dyobs) != len(xobs):
            raise ValueError("xobs and dyobs are different lengths")

        self.xobs = numpy.asarray(xobs, dtype=float)
        self.yobs = numpy.asarray(yobs)

        if dyobs is None:
            self.dyobs = numpy.ones_like(xobs)
        else:
            self.dyobs = numpy.asarray(dyobs, dtype=float)

        # Set the default calculation points
        if self.x is None:
            self.setCalculationPoints(self.xobs)
        else:
            self.setCalculationPoints(self.x)

        return

    def setCalculationRange(self, xmin = None, xmax = None, dx = None):
        """Set the calculation range

        Arguments
        xmin    --  The minimum value of the independent variable.
                    If xmin is None (default), the minimum observed value will
                    be used. This is clipped to the minimum observed x.
        xmax    --  The maximum value of the independent variable.
                    If xmax is None (default), the maximum observed value will
                    be used. This is clipped to the maximum observed x.
        dx      --  The sample spacing in the independent variable. If dx is
                    None (default), then the spacing in the observed points
                    will be preserved.

        Note that xmin is always inclusive (unless clipped). xmax is inclusive
        if it is within the bounds of the observed data.

        raises AttributeError if there is no observed profile
        raises ValueError if xmin > xmax
        raises ValueError if dx > xmax-xmin
        raises ValueError if dx <= 0

        """
        clip = dx is None

        if self.xobs is None:
            raise AttributeError("No observed profile")

        if xmin is None and xmax is None and dx is None:
            self.x = self.xobs
            self.y = self.yobs
            self.dy = self.dyobs
            self._updatePars()
            return

        if xmin is None:
            xmin = self.xobs[0]
        else:
            xmin = float(xmin)

        if xmax is None:
            xmax = self.xobs[-1]
        else:
            xmax = float(xmax)

        if dx is None:
            dx = (self.xobs[-1] - self.xobs[0]) / len(self.xobs)
        else:
            dx = float(dx)

        if xmin > xmax:
            raise ValueError("xmax must be greater than xmin")
        if dx > xmax - xmin:
            raise ValueError("dx must be less than xmax-xmin")
        if dx <= 0:
            raise ValueError("dx must be positive")

        if clip:
            x = self.xobs
            indices = numpy.logical_and( xmin - epsilon <= x , x <= xmax +
                    epsilon )
            self.x = self.xobs[indices]
            self.y = self.yobs[indices]
            self.dy = self.dyobs[indices]
        else:
            self.setCalculationPoints(numpy.arange(xmin, xmax+0.5*dx, dx))

        self._updatePars()
        return

    def setCalculationPoints(self, x):
        """Set the calculation points.

        Arguments
        x   --  A non-empty numpy array containing the calculation points. If
                xobs exists, the bounds of x will be limited to its bounds.

        This will create y and dy on the specified grid if xobs, yobs and
        dyobs exist.

        """
        x = numpy.asarray(x)
        if self.xobs is not None:
            x = x[ x >= self.xobs[0] - epsilon ]
            x = x[ x <= self.xobs[-1] + epsilon ]
        self.x = x
        if self.yobs is not None:
            self.y = rebinArray(self.yobs, self.xobs, self.x)
        if self.dyobs is not None:
            # FIXME - This does not follow error propogation rules and it
            # introduces (more) correlation between the data points.
            self.dy = rebinArray(self.dyobs, self.xobs, self.x)

        self._updatePars()
        return

    def _updatePars(self):
        """Update the paramters holding the observed arrays."""
        self.xpar.setValue(self.x)
        self.ypar.setValue(self.y)
        self.dypar.setValue(self.dy)
        return

# End class Profile

def rebinArray(A, xold, xnew):
    """Rebin the an array by interpolating over the new x range.

    Arguments:
    A       --  Array to interpolate
    xold    --  Old sampling array
    xnew    --  New sampling array

    This uses cubic spline interpolation.
    
    Returns: A new array over the new sampling array.
    """
    if numpy.array_equal(xold, xnew):
        return A
    from scipy.interpolate import splrep, splev
    finterp = splrep(xold, A, s=0)
    return splev(xnew, finterp, der=0)

__id__ = "$Id$"
