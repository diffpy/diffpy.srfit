#!/usr/bin/env python
########################################################################
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
########################################################################

"""The Profile class containing the physical and calculated data.

Profile holds the arrays representing an observed profile, a selected subset of
the observed profile and a calculated profile. Profiles are used by Calculators
to store a calculated signal, and by FitContributions to help calculate a
residual equation.
"""

__all__ = ["Parameter", "Profile"]

import numpy

from diffpy.srfit.util.observable import Observable
from diffpy.srfit.fitbase.parameter import Parameter
from diffpy.srfit.fitbase.validatable import Validatable
from diffpy.srfit.exceptions import SrFitError

# This is the roundoff tolerance for selecting bounds on arrays.
epsilon = 1e-8

class Profile(Observable, Validatable):
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
    xpar    --  A Parameter that stores x (named "x").
    ypar    --  A Parameter that stores y (named "y").
    dypar   --  A Parameter that stores dy (named "dy").
    ycpar   --  A Parameter that stores ycalc (named "ycalc"). This is
                not observed by the profile, but it is present so it can be
                constrained to.
    meta    --  A dictionary of metadata. This is only set if provided by a
                parser.

    """

    def __init__(self):
        """Initialize the attributes."""
        Observable.__init__(self)
        self._xobs = None
        self._yobs = None
        self._dyobs = None
        self.xpar = Parameter("x")
        self.ypar = Parameter("y")
        self.dypar = Parameter("dy")
        self.ycpar = Parameter("ycalc")
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
            # work around for interpolation issue making some of these non-1
            if (self.dyobs == 1).all():
                self.dy = numpy.ones_like(self.x)
            else:
            # FIXME - This does not follow error propogation rules and it
            # introduces (more) correlation between the data points.
                self.dy = rebinArray(self.dyobs, self.xobs, self.x)

        return

    def loadtxt(self, *args, **kw):
        """Use numpy.loadtxt to load data.

        Arguments are passed to numpy.loadtxt.
        unpack = True is enforced.
        The first two arrays returned by numpy.loadtxt are assumed to be x and
        y.  If there is a third array, it is assumed to by dy. Any other arrays
        are ignored. These are passed to setObservedProfile.

        Raises ValueError if the call to numpy.loadtxt returns fewer than 2
        arrays.

        Returns the x, y and dy arrays loaded from the file

        """
        if len(args) == 8 and not args[-1]:
            args = list(args)
            args[-1] = True
        else:
            kw["unpack"] = True
        cols = numpy.loadtxt(*args, **kw)

        x = y = dy = None
        # Due to using 'unpack', a single column will come out as a single
        # array, thus the second check.
        if len(cols) < 2 or not isinstance(cols[0], numpy.ndarray):
            raise ValueError("numpy.loadtxt returned fewer than 2 arrays")
        x = cols[0]
        y = cols[1]
        if len(cols) > 2:
            dy = cols[2]

        self.setObservedProfile(x, y, dy)
        return x, y, dy

    def savetxt(self, fname, fmt='%.18e', delimiter=' '):
        """Call numpy.savetxt with x, ycalc, y, dy

        Arguments are passed to numpy.savetxt.

        """
        x = self.x
        ycalc = self.ycalc
        if ycalc is None:
            raise AttributeError("ycalc is None")
        y = self.y
        dy = self.dy

        # Add the header
        if not hasattr(fname, 'write'):
            fname = file(fname, 'w')
        if fname.closed:
            raise ValueError("I/O operation on closed file")
        header = "# x           ycalc           y           dy\n"
        fname.write(header)
        numpy.savetxt(fname, zip(x, ycalc, y, dy), fmt, delimiter)
        return

    def _flush(self, other):
        """Invalidate cached state.

        This will force any observer to invalidate its state.

        """
        self.ycalc = None
        self.notify(other)
        return

    def _validate(self):
        """Validate my state.

        This validates that x, y, dy, xobx, yobs and dyobs are not None.
        This validates that x, y, and dy are the same length.

        Raises SrFitError if validation fails.

        """
        datanotset = any(v is None for v in
                [self.x, self.y, self.dy, self.xobs, self.yobs, self.dyobs])
        if datanotset:
            raise SrFitError("Missing data")
        if len(self.x) != len(self.y) or len(self.x) != len(self.dy):
            raise SrFitError("Data are different lengths")
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
    return numpy.interp(xnew, xold, A)
