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
"""The Profile class containing the physical and calculated data.

Profile holds the arrays representing an observed profile, a selected
subset of the observed profile and a calculated profile. Profiles are
used by Calculators to store a calculated signal, and by
FitContributions to help calculate a residual equation.
"""

__all__ = ["Parameter", "Profile"]

import numpy
import six

from diffpy.srfit.exceptions import SrFitError
from diffpy.srfit.fitbase.parameter import Parameter
from diffpy.srfit.fitbase.validatable import Validatable
from diffpy.srfit.util.observable import Observable

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
    _dyobs  --  A numpy array of the uncertainty of the observed signal
                (default None, optional).
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
    x = property(
        lambda self: self.xpar.getValue(),
        lambda self, val: self.xpar.setValue(val),
    )
    y = property(
        lambda self: self.ypar.getValue(),
        lambda self, val: self.ypar.setValue(val),
    )
    dy = property(
        lambda self: self.dypar.getValue(),
        lambda self, val: self.dypar.setValue(val),
    )
    ycalc = property(
        lambda self: self.ycpar.getValue(),
        lambda self, val: self.ycpar.setValue(val),
    )

    # We want xobs, yobs and dyobs to be read-only
    xobs = property(lambda self: self._xobs)
    yobs = property(lambda self: self._yobs)
    dyobs = property(lambda self: self._dyobs)

    def loadParsedData(self, parser):
        """Load parsed data from a ProfileParser.

        This sets the xobs, yobs, dyobs arrays as well as the metadata.
        """
        x, y, junk, dy = parser.getData()
        self.meta = dict(parser.getMetaData())
        self.setObservedProfile(x, y, dy)
        return

    def setObservedProfile(self, xobs, yobs, dyobs=None):
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

    def setCalculationRange(self, xmin=None, xmax=None, dx=None):
        """Set epsilon-inclusive calculation range.

        Adhere to the observed ``xobs`` points when ``dx`` is the same
        as in the data.  ``xmin`` and ``xmax`` are clipped at the bounds
        of the observed data.

        Parameters
        ----------

        xmin : float or "obs", optional
            The minimum value of the independent variable.  Keep the
            current minimum when not specified.  If specified as "obs"
            reset to the minimum observed value.
        xmax : float or "obs", optional
            The maximum value of the independent variable.  Keep the
            current maximum when not specified.  If specified as "obs"
            reset to the maximum observed value.
        dx : float or "obs", optional
            The sample spacing in the independent variable.  When different
            from the data, resample the ``x`` as anchored at ``xmin``.

        Note that xmin is always inclusive (unless clipped). xmax is inclusive
        if it is within the bounds of the observed data.

        Raises
        ------
        AttributeError
            If there is no observed data.
        ValueError
            When xmin > xmax or if dx <= 0.  Also if dx > xmax - xmin.
        """
        if self.xobs is None:
            raise AttributeError("No observed profile")

        # local helper function
        def _isobs(a):
            if not isinstance(a, six.string_types):
                return False
            if a != "obs":
                raise ValueError('Must be either float or "obs".')
            return True

        # resolve new low and high bounds for x
        lo = (
            self.x[0]
            if xmin is None
            else self.xobs[0] if _isobs(xmin) else float(xmin)
        )
        lo = max(lo, self.xobs[0])
        hi = (
            self.x[-1]
            if xmax is None
            else self.xobs[-1] if _isobs(xmax) else float(xmax)
        )
        hi = min(hi, self.xobs[-1])
        # determine if we need to clip the original grid
        clip = True
        step = None
        ncur = len(self.x)
        stepcur = 1 if ncur < 2 else (self.x[-1] - self.x[0]) / (ncur - 1.0)
        nobs = len(self.xobs)
        stepobs = (
            1 if nobs < 2 else (self.xobs[-1] - self.xobs[0]) / (nobs - 1.0)
        )
        if dx is None:
            # check if xobs overlaps with x
            i0 = numpy.fabs(self.xobs - self.x[0]).argmin()
            n0 = min(len(self.x), len(self.xobs) - i0)
            if not numpy.allclose(self.xobs[i0 : i0 + n0], self.x[:n0]):
                clip = False
                step = stepcur if ncur > 1 else stepobs
        elif _isobs(dx):
            assert clip and step is None
        elif numpy.allclose(stepobs, dx):
            assert clip and step is None
        else:
            clip = False
            step = float(dx)
        # verify that we either clip or have the step defined.
        assert clip or step is not None
        # hi, lo, step, clip all resolved here.
        # validate arguments
        if lo > hi:
            raise ValueError("xmax must be greater than xmin.")
        if not clip:
            if step > hi - lo:
                raise ValueError("dx must be less than (xmax - xmin).")
            if step <= 0:
                raise ValueError("dx must be positive.")
        # determine epsilon extensions to the lower and upper bounds.
        epslo = abs(lo) * epsilon + epsilon
        epshi = abs(hi) * epsilon + epsilon
        # process the new grid.
        if clip:
            indices = (lo - epslo <= self.xobs) & (self.xobs <= hi + epshi)
            self.x = self.xobs[indices]
            self.y = self.yobs[indices]
            self.dy = self.dyobs[indices]
        else:
            x1 = numpy.arange(lo, hi + epshi, step)
            self.setCalculationPoints(x1)
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
            x = x[x >= self.xobs[0] - epsilon]
            x = x[x <= self.xobs[-1] + epsilon]
        self.x = x
        if self.yobs is not None:
            self.y = rebinArray(self.yobs, self.xobs, self.x)
        if self.dyobs is not None:
            # work around for interpolation issue making some of these non-1
            if (self.dyobs == 1).all():
                self.dy = numpy.ones_like(self.x)
            else:
                # FIXME - This does not follow error propagation rules and it
                # introduces (more) correlation between the data points.
                self.dy = rebinArray(self.dyobs, self.xobs, self.x)

        return

    def loadtxt(self, *args, **kw):
        """Use numpy.loadtxt to load data.

        Arguments are passed to numpy.loadtxt. unpack = True is
        enforced. The first two arrays returned by numpy.loadtxt are
        assumed to be x and y. If there is a third array, it is assumed
        to by dy. Any other arrays are ignored. These are passed to
        setObservedProfile.

        Raises ValueError if the call to numpy.loadtxt returns fewer
        than 2 arrays.

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

    def savetxt(self, fname, **kwargs):
        """Call `numpy.savetxt` with x, ycalc, y, dy.

        Parameters
        ----------
        fname : filename or file handle
            This is passed to `numpy.savetxt`.
        **kwargs
            The keyword arguments that are passed to `numpy.savetxt`.
            We preset file header "x  ycalc  y  dy".  Use ``header=''``
            to save data without any header.

        Raises
        ------
        SrFitError
            When `self.ycalc` has not been set.

        See also
        --------
        numpy.savetxt
        """
        x = self.x
        ycalc = self.ycalc
        if ycalc is None:
            raise SrFitError("ycalc is None")
        y = self.y
        dy = self.dy
        kwargs.setdefault("header", "x  ycalc  y  dy")
        data = numpy.transpose([x, ycalc, y, dy])
        numpy.savetxt(fname, data, **kwargs)
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
        datanotset = any(
            v is None
            for v in [
                self.x,
                self.y,
                self.dy,
                self.xobs,
                self.yobs,
                self.dyobs,
            ]
        )
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
