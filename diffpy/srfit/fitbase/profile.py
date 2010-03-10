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
to store a calculated signal, and by FitContributions to help calculate a
residual equation.

"""
__all__ = ["ProfileParameter", "Profile"]

import numpy

from diffpy.srfit.util.observable import Observable
from .parameter import Parameter

# This is the roundoff tolerance for selecting bounds on arrays.
epsilon = 1e-8

class ProfileParameter(Parameter):
    """A Parameter for profiles that can have a None value."""

    def getValue(self):
        """Get the value, even if it is None."""
        return self._value

# End class ProfileParameter

class Profile(Observable):
    """Observed and calculated profile container.

    Profile is an Observable. The xpar, ypar and dypar attributes are observed
    by the Profile, which can in turn be observed by some other object.

    Attributes

    _xobs   --  A numpy array of the observed independent variable (default
                None)
    xobs    --  Read-only property of _xobs.
    _yobs   --  A numpy array of the observed signal (default None)
    yobs    --  Read-only property of _yobs.
    _dyobs  --  A numpy array of the uncertainty of the observed signal (default
                None, optional).
    dyobs   --  Read-only property of _dyobs.
    x       --  A numpy array of the calculated independent variable (default
                None, property for xpar accessors).
    y       --  The profile over the calculation range (default None, property
                for ypar accessors).
    dy      --  The uncertainty in the profile over the calculation range
                (default None, property for dypar accessors).
    ycalc   --  A numpy array of the calculated signal (default None).
    xpar    --  A ProfileParameter that stores x (named "x").
    ypar    --  A ProfileParameter that stores y (named "y").
    dypar   --  A ProfileParameter that stores dy (named "dy").
    ycpar   --  A ProfileParameter that stores ycalc (named "ycalc"). This is
                not observed by the profile, but it is present so it can be
                observed from elsewhere.
    meta    --  A dictionary of metadata. This is only set if provided by a
                parser.

    """

    def __init__(self):
        """Initialize the attributes."""
        Observable.__init__(self)
        self._xobs = None
        self._yobs = None
        self._dyobs = None
        self.xpar = ProfileParameter("x")
        self.ypar = ProfileParameter("y")
        self.dypar = ProfileParameter("dy")
        self.ycpar = ProfileParameter("ycalc")
        self.meta = {}

        # Observable
        self.xpar.addObserver(self._flush)
        self.ypar.addObserver(self._flush)
        self.dypar.addObserver(self._flush)
        return

    # We want x, y, ycalc and dy to stay in-sync with xpar, ypar and dypar
    x = property( lambda self : self.xpar.getValue(),
                  lambda self, val : self.xpar.setValue(val) )
    y = property( lambda self : self.ypar.getValue(),
                  lambda self, val : self.ypar.setValue(val) )
    dy = property( lambda self : self.dypar.getValue(),
                   lambda self, val : self.dypar.setValue(val) )
    ycalc = property( lambda self : self.ycpar.getValue(),
                  lambda self, val : self.ycpar.setValue(val) )

    # We want xobs, yobs and dyobs to be read-only
    xobs = property( lambda self: self._xobs )
    yobs = property( lambda self: self._yobs )
    dyobs = property( lambda self: self._dyobs )

    def loadParsedData(self, parser):
        """Load parsed data from a ProfileParser.

        This sets the xobs, yobs, dyobs arrays as well as the metadata.

        """
        x, y, junk, dy = parser.getData()
        self.meta = dict(parser.getMetaData())
        self.setObservedProfile(x, y, dy)
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

        self._xobs = numpy.asarray(xobs, dtype=float)
        self._yobs = numpy.asarray(yobs, dtype=float)

        if dyobs is None:
            self._dyobs = numpy.ones_like(xobs)
        else:
            self._dyobs = numpy.asarray(dyobs, dtype=float)

        # Set the default calculation points
        if self.x is None:
            self.setCalculationPoints(self._xobs)
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

        return

    def _flush(self, other):
        """Invalidate cached state.

        This will force any observer to invalidate its state.

        """
        self.ycalc = None
        self.notify()
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
