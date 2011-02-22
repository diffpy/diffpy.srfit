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
import numpy
from diffpy.srfit.fit.parameters import Fixed

__all__ = ["loadProfile"]

epsilon = 1e-8

def loadProfile(filename, fmt = 'txt', *args, **kw):
    """Load a profile with a given format.

    filename    --  Name of file to be loaded into a profile.
    fmt         --  The format of the data (string). This is used to select a
                    parser for the data. The 'parserInfo' function prints
                    information about all parsers and the 'getParser' function
                    can be used to retrieve and introspect a parser based on
                    its format.

    Remaining arguments are passed to the parser.

    Returns a configured Profile instance.

    Raises ValueError if a parser for the format cannot be found.
    Raises IOError if the file cannot be read.
    Raises ParseError if the file cannot be parsed.

    """

    from diffpy.srfit.fit.profileparser import getParser
    ParserClass = getParser(fmt)
    parser = ParserClass()
    parser.parseFile(filename, *args, **kw)

    profile = Profile()
    profile.load(parser)

    return profile

class Profile(object):
    """Observed and calculated profile container.

    Attributes

    xobs    --  A numpy array of the observed independent variable (default
                None)
    yobs    --  A numpy array of the observed signal (default None)
    dyobs   --  A numpy array of the uncertainty of the observed signal (default
                None, optional).
    x       --  xobs over the fit range (default None, property giving
                value access to xpar).
    y       --  yobs over the fit range (default None, property giving
                value access to ypar).
    dy      --  dyobs over the fit range (default None, property giving
                value access to dypar).
    xpar    --  A Parameter that stores x.
    ypar    --  A Parameter that stores y.
    dypar   --  A Parameter that stores dy.
    pars    --  Tuple of (xpar, ypar, dypar) (property).
    meta    --  A dictionary of metadata. This is only set if provided by a
                parser.

    Profile is iterable and returns x, y and dy in that order.

    """

    def __init__(self):
        """Initialize the attributes."""
        self.xobs = None
        self.yobs = None
        self.dyobs = None
        self.xpar = Fixed("x", None)
        self.ypar = Fixed("y", None)
        self.dypar = Fixed("dy", None)
        self.meta = {}
        return

    def __iter__(self):
        """Provided for easy access to the parameters.

        Returns (xpar, ypar, dypar).
        """
        return iter((self.x, self.y, self.dy))

    # Parameter properties
    x = property( lambda self : self.xpar.get(),
                  lambda self, val : self.xpar.set(val) )
    y = property( lambda self : self.ypar.get(),
                  lambda self, val : self.ypar.set(val) )
    dy = property( lambda self : self.dypar.get(),
                   lambda self, val : self.dypar.set(val) )

    pars = property( lambda self: (self.xpar, self.ypar, self.dypar) )

    def load(self, parser):
        """Load parsed data from a ProfileParser.

        This sets the xobs, yobs, dyobs arrays as well as the metadata.

        """
        x, y, junk, dy = parser.getData()
        self.meta = dict(parser.getMetaData())
        self.setObserved(x, y, dy)
        return

    def setObserved(self, xobs, yobs, dyobs = None):
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
        self.yobs = numpy.asarray(yobs, dtype=float)

        if dyobs is None:
            self.dyobs = numpy.ones_like(xobs)
        else:
            self.dyobs = numpy.asarray(dyobs, dtype=float)

        # Set the default calculation points
        if self.x is None:
            self.setPoints(self.xobs)
        else:
            self.setPoints(self.x)

        return

    def setRange(self, xmin = None, xmax = None, xstep = None):
        """Set the calculation range

        Arguments
        xmin    --  The minimum value of the independent variable.
                    If xmin is None (default), the minimum observed value will
                    be used. This is clipped to the minimum observed x.
        xmax    --  The maximum value of the independent variable.
                    If xmax is None (default), the maximum observed value will
                    be used. This is clipped to the maximum observed x.
        xstep   --  The sample spacing in the independent variable. If xstep is
                    None (default), then the spacing in the observed points
                    will be preserved.

        Note that xmin is always inclusive (unless clipped). xmax is inclusive
        if it is within the bounds of the observed data.

        raises AttributeError if there is no observed profile
        raises ValueError if xmin > xmax
        raises ValueError if xstep > xmax-xmin
        raises ValueError if xstep <= 0

        """
        clip = xstep is None

        if self.xobs is None:
            raise AttributeError("No observed profile")

        if xmin is None and xmax is None and xstep is None:
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

        if xstep is None:
            xstep = (self.xobs[-1] - self.xobs[0]) / len(self.xobs)
        else:
            xstep = float(xstep)

        if xmin > xmax:
            raise ValueError("xmax must be greater than xmin")
        if xstep > xmax - xmin:
            raise ValueError("xstep must be less than xmax-xmin")
        if xstep <= 0:
            raise ValueError("xstep must be positive")

        if clip:
            x = self.xobs
            indices = numpy.logical_and( xmin - epsilon <= x , x <= xmax +
                    epsilon )
            self.x = self.xobs[indices]
            self.y = self.yobs[indices]
            self.dy = self.dyobs[indices]
        else:
            self.setPoints(numpy.arange(xmin, xmax+0.5*xstep, xstep))

        return

    def setPoints(self, x):
        """Set the calculation points.

        Arguments
        x   --  A non-empty numpy array containing the calculation points. If
                xobs exists, a valueError will be raised if x falls out of its
                bounds.

        This will create y and dy on the specified grid if xobs, yobs and
        dyobs exist. Error in dy will be properly calculated.

        """
        x = numpy.asarray(x)
        if self.xobs is not None:
            if x[0] < self.xobs[0] or x[-1] > self.xobs[-1]:
                msg = "Points lie outside the observed profile"
                raise ValueError(msg)
        self.x = x
        if self.yobs is not None:
            if self.dyobs is None:
                self.y = numpy.interp(self.x, self.xobs, self.yobs)
            elif (self.dyobs == 1).all():
                self.y = numpy.interp(self.x, self.xobs, self.yobs)
                self.dy = numpy.ones_like(self.x)
            else:
                # This will properly propagate the interpolation errors, but
                # it's slow.
                self.y, self.dy = _interpWithError(self.x, self.xobs,
                        self.yobs, self.dyobs)

        return



# End class Profile
def _interpWithError(zarr, xarr, yarr, dyarr):
    """Interpolate data and uncertainty onto a new grid."""
    if zarr is xarr:
        return yarr, dyarr

    # Checkf or interpolation bounds
    if zarr[0] < xarr[0] or zarr[-1] > xarr[-1]:
        msg = "Cannot interpolate outside of bounds"
        raise ValueError(msg)
    if len(xarr) < 2:
        msg = "xarr too short for interpolation"
        raise ValueError(msg)

    # Check for increasing z and x
    oldz = zarr[0]
    for z in zarr[1:]:
        if z <= oldz:
            msg = "zarr is non-increasing"
            raise ValueError(msg)
        oldz = z
    oldx = xarr[0]
    for x in xarr[1:]:
        if x <= oldx:
            msg = "xarr is non-increasing"
            raiseValueError(msg)
        oldx = x


    newy = []
    newdy = []

    xidx = 0
    x1 = xarr[0]
    x2 = xarr[1]

    for z in zarr:

        # Increment x-values. We need a value on either side of z.
        while x2 < z:
            xidx += 1
            x1 = xarr[xidx]
            x2 = xarr[xidx+1]

        # Now we can interpolate
        y1 = yarr[xidx]
        y2 = yarr[xidx+1]
        dy1 = dyarr[xidx]
        dy2 = dyarr[xidx+1]
        mz = (x2 - z) / (x2 - x1)
        if mz == 0:
            yp = y2
            dyp = dy2
        if mz == 1:
            yp = y1
            dyp = dy1
        else:
            yp = y2 + (y1 - y2) * mz
            dyp = (dy2**2 + (dy2**2 + dy1**2)*mz**2)**0.5
        newy.append(yp)
        newdy.append(dyp)

    newyarr = numpy.array(newy, dtype=float)
    newdyarr = numpy.array(newdy, dtype=float)

    return newyarr, newdyarr


__id__ = "$Id$"
