#!/usr/bin/env python
##############################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2010 The Trustees of Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################
"""Form factors (characteristic functions) used in PDF nanoshape fitting.

These are used to calculate the attenuation of the PDF due to a finite
size. For a crystal-like nanoparticle, one can calculate the PDF via
Gnano(r) = f(r) Gcryst(r), where f(r) is the nanoparticle characteristic
function and Gcryst(f) is the crystal PDF.

These functions are meant to be imported and added to a FitContribution
using the 'registerFunction' method of that class.
"""

__all__ = [
    "sphericalCF",
    "spheroidalCF",
    "spheroidalCF2",
    "lognormalSphericalCF",
    "sheetCF",
    "shellCF",
    "shellCF2",
    "SASCF",
]

import numpy
from numpy import arctan as atan
from numpy import arctanh as atanh
from numpy import ceil, exp, log, log2, pi, sign, sqrt
from numpy.fft import fftfreq, ifft
from scipy.special import erf

from diffpy.srfit.fitbase.calculator import Calculator


def sphericalCF(r, psize):
    """Spherical nanoparticle characteristic function.

    r       --  distance of interaction
    psize   --  The particle diameter

    From Kodama et al., Acta Cryst. A, 62, 444-453
    (converted from radius to diameter)
    """
    f = numpy.zeros(numpy.shape(r), dtype=float)
    if psize > 0:
        x = numpy.array(r, dtype=float) / psize
        inside = x < 1.0
        xin = x[inside]
        f[inside] = 1.0 - 1.5 * xin + 0.5 * xin * xin * xin
    return f


def spheroidalCF(r, erad, prad):
    """Spheroidal characteristic function specified using radii.

    Spheroid with radii (erad, erad, prad)

    prad    --  polar radius
    erad    --  equatorial radius

    erad < prad equates to a prolate spheroid
    erad > prad equates to a oblate spheroid
    erad == prad is a sphere
    """
    psize = 2.0 * erad
    pelpt = 1.0 * prad / erad
    return spheroidalCF2(r, psize, pelpt)


def spheroidalCF2(r, psize, axrat):
    """Spheroidal nanoparticle characteristic function.

    Form factor for ellipsoid with radii (psize/2, psize/2, axrat*psize/2)

    r      --  distance of interaction
    psize  --  The equatorial diameter
    axrat  --  The ratio of axis lengths

    From Lei et al., Phys. Rev. B, 80, 024118 (2009)
    """
    pelpt = 1.0 * axrat

    if psize <= 0 or pelpt <= 0:
        return numpy.zeros_like(r)

    # to simplify the equations
    v = pelpt
    d = 1.0 * psize
    d2 = d * d
    v2 = v * v

    if v == 1:
        return sphericalCF(r, psize)

    rx = r
    if v < 1:

        r = rx[rx <= v * psize]
        r2 = r * r
        f1 = (
            1
            - 3 * r / (4 * d * v) * (1 - r2 / (4 * d2) * (1 + 2.0 / (3 * v2)))
            - 3
            * r
            / (4 * d)
            * (1 - r2 / (4 * d2))
            * v
            / sqrt(1 - v2)
            * atanh(sqrt(1 - v2))
        )

        r = rx[numpy.logical_and(rx > v * psize, rx <= psize)]
        r2 = r * r
        f2 = (
            (
                3 * d / (8 * r) * (1 + r2 / (2 * d2)) * sqrt(1 - r2 / d2)
                - 3
                * r
                / (4 * d)
                * (1 - r2 / (4 * d2))
                * atanh(sqrt(1 - r2 / d2))
            )
            * v
            / sqrt(1 - v2)
        )

        r = rx[rx > psize]
        f3 = numpy.zeros_like(r)

        f = numpy.concatenate((f1, f2, f3))

    elif v > 1:

        r = rx[rx <= psize]
        r2 = r * r
        f1 = (
            1
            - 3 * r / (4 * d * v) * (1 - r2 / (4 * d2) * (1 + 2.0 / (3 * v2)))
            - 3
            * r
            / (4 * d)
            * (1 - r2 / (4 * d2))
            * v
            / sqrt(v2 - 1)
            * atan(sqrt(v2 - 1))
        )

        r = rx[numpy.logical_and(rx > psize, rx <= v * psize)]
        r2 = r * r
        f2 = (
            1
            - 3 * r / (4 * d * v) * (1 - r2 / (4 * d2) * (1 + 2.0 / (3 * v2)))
            - 3.0
            / 8
            * (1 + r2 / (2 * d2))
            * sqrt(1 - d2 / r2)
            * v
            / sqrt(v2 - 1)
            - 3
            * r
            / (4 * d)
            * (1 - r2 / (4 * d2))
            * v
            / sqrt(v2 - 1)
            * (atan(sqrt(v2 - 1)) - atan(sqrt(r2 / d2 - 1)))
        )

        r = rx[rx > v * psize]
        f3 = numpy.zeros_like(r)

        f = numpy.concatenate((f1, f2, f3))

    return f


def lognormalSphericalCF(r, psize, psig):
    """Spherical nanoparticle characteristic function with lognormal size
    distribution.

    r      --  distance of interaction
    psize  --  The mean particle diameter
    psig   --  The log-normal width of the particle diameter

    Here, r is the independent variable, mu is the mean of the distribution
    (not of the particle size), and s is the width of the distribution. This is
    the characteristic function for the lognormal distribution of particle
    diameter:

    F(r, mu, s) = 0.5*Erfc((-mu-3*s^2+Log(r))/(sqrt(2)*s))
               + 0.25*r^3*Erfc((-mu+Log(r))/(sqrt(2)*s))*exp(-3*mu-4.5*s^2)
               - 0.75*r*Erfc((-mu-2*s^2+Log(r))/(sqrt(2)*s))*exp(-mu-2.5*s^2)

    The expectation value of the distribution gives the average particle
    diameter, psize. The variance of the distribution gives psig^2. mu and s
    can be expressed in terms of these as:

    s^2 = log((psig/psize)^2 + 1)
    mu = log(psize) - s^2/2

    Source unknown
    """
    if psize <= 0:
        return numpy.zeros_like(r)
    if psig <= 0:
        return sphericalCF(r, psize)

    sqrt2 = sqrt(2.0)
    s = sqrt(log(psig * psig / (1.0 * psize * psize) + 1))
    mu = log(psize) - s * s / 2
    if mu < 0:
        return numpy.zeros_like(r)

    return (
        0.5 * erfc((-mu - 3 * s * s + log(r)) / (sqrt2 * s))
        + 0.25
        * r
        * r
        * r
        * erfc((-mu + log(r)) / (sqrt2 * s))
        * exp(-3 * mu - 4.5 * s * s)
        - 0.75
        * r
        * erfc((-mu - 2 * s * s + log(r)) / (sqrt2 * s))
        * exp(-mu - 2.5 * s * s)
    )


def sheetCF(r, sthick):
    """Nanosheet characteristic function.

    r       --  distance of interaction
    sthick  --  Thickness of nanosheet

    From Kodama et al., Acta Cryst. A, 62, 444-453
    """
    # handle zero or negative sthick.  make it work for scalars and arrays.
    if sthick <= 0:
        return 0 * sthick
    # process scalar r
    if numpy.isscalar(r):
        rv = 1 - 0.5 * r / sthick if r < sthick else 0.5 * sthick / r
        return rv
    # handle array-type r
    ra = numpy.asarray(r)
    lo = ra < sthick
    hi = ~lo
    f = numpy.empty_like(ra, dtype=float)
    f[lo] = 1 - 0.5 * ra[lo] / sthick
    f[hi] = 0.5 * sthick / ra[hi]
    return f


def shellCF(r, radius, thickness):
    """Spherical shell characteristic function.

    radius      --  Inner radius
    thickness   --  Thickness of shell

    outer radius = radius + thickness

    From Lei et al., Phys. Rev. B, 80, 024118 (2009)
    """
    d = 1.0 * thickness
    a = 1.0 * radius + d / 2.0
    return shellCF2(r, a, d)


def shellCF2(r, a, delta):
    """Spherical shell characteristic function.

    a       --  Central radius
    delta   --  Thickness of shell

    outer radius = a + thickness/2

    From Lei et al., Phys. Rev. B, 80, 024118 (2009)
    """
    a = 1.0 * a
    d = 1.0 * delta
    a2 = a**2
    d2 = d**2
    dmr = d - r
    dmr2 = dmr**2

    f = (
        r
        * (
            16 * a * a2
            + 12 * a * d * dmr
            + 36 * a2 * (2 * d - r)
            + 3 * dmr2 * (2 * d + r)
        )
        + 2 * dmr2 * (r * (2 * d + r) - 12 * a2) * sign(dmr)
        - 2 * (2 * a - r) ** 2 * (r * (4 * a + r) - 3 * d2) * sign(2 * a - r)
        + r * (4 * a - 2 * d + r) * (2 * a - d - r) ** 2 * sign(2 * a - d - r)
    )

    f[r > 2 * a + d] = 0

    den = 8.0 * r * d * (12 * a2 + d2)
    zmask = den == 0.0
    vmask = ~zmask
    f[vmask] /= den[vmask]
    f[zmask] = 1
    return f


class SASCF(Calculator):
    """Calculator class for characteristic functions from sas-models.

    This class wraps a sas.models.BaseModel to calculate I(Q) related to
    nanoparticle shape. This I(Q) is inverted to f(r) according to:
    f(r) = 1 / (4 pi r) * SINFT(I(Q)),
    where "SINFT" represents the sine Fourier transform.

    Attributes:
    _model      --  BaseModel object this adapts.

    Managed Parameters:
    These depend on the parameters of the BaseModel object held by _model. They
    are created from the 'params' attribute of the BaseModel. If a dispersion
    is set for the BaseModel, the dispersion "width" will be accessible under
    "<parname>_width", where <parname> is the name a parameter adjusted by
    dispersion.
    """

    def __init__(self, name, model):
        """Initialize the generator.

        name    --  A name for the SASCF
        model   --  SASModel object this adapts.
        """
        Calculator.__init__(self, name)

        self._model = model

        from diffpy.srfit.sas.sasparameter import SASParameter

        # Wrap normal parameters
        for parname in model.params:
            par = SASParameter(parname, model)
            self.addParameter(par)

        # Wrap dispersion parameters
        for parname in model.dispersion:
            name = parname + "_width"
            parname += ".width"
            par = SASParameter(name, model, parname)
            self.addParameter(par)

        return

    def __call__(self, r):
        """Calculate the characteristic function from the transform of the
        BaseModel."""

        # Determine q-values.
        # We want very fine r-spacing so we can properly normalize f(r). This
        # equates to having a large qmax so that the Fourier transform is
        # finely spaced. We also want the calculation to be fast, so we pick
        # qmax such that the number of q-points is a power of 2. This allows us
        # to use the fft.
        #
        # The initial dr is somewhat arbitrary, but using dr = 0.01 allows for
        # the f(r) calculated from a particle of diameter 50, over r =
        # arange(1, 60, 0.1) to agree with the sphericalCF with Rw < 1e-4%.
        #
        # We also have to make a q-spacing small enough to compute out to at
        # least the size of the signal.
        dr = min(0.01, r[1] - r[0])
        ed = 2 * self._model.calculate_ER()

        # Check for nans. If we find any, then return zeros.
        if numpy.isnan(ed).any():
            y = numpy.zeros_like(r)
            return y

        rmax = max(ed, 2 * r[-1])
        dq = pi / rmax
        qmax = pi / dr
        numpoints = int(2 ** (ceil(log2(qmax / dq))))
        qmax = dq * numpoints

        # Calculate F(q) = q * I(q) from model
        q = fftfreq(int(qmax / dq)) * qmax
        fq = q * self._model.evalDistribution(q)

        # Calculate g(r) and the effective r-points
        rp = fftfreq(numpoints) * 2 * pi / dq
        # Note sine transform = imaginary part of ifft
        gr = ifft(fq).imag

        # Calculate full-fr for normalization
        assert rp[0] == 0.0
        frp = numpy.zeros_like(gr)
        frp[1:] = gr[1:] / rp[1:]

        # Inerpolate onto requested grid, do not use data after jump in rp
        assert numpoints % 2 == 0
        nhalf = numpoints / 2
        fr = numpy.interp(r, rp[:nhalf], gr[:nhalf])
        vmask = r != 0
        fr[vmask] /= r[vmask]

        # Normalize. We approximate fr[0] by using the fact that f(r) is linear
        # at low r. By definition, fr[0] should equal 1.
        fr0 = 2 * frp[2] - frp[1]
        fr /= fr0

        # Fix potential divide-by-zero issue, fr is 1 at r == 0
        fr[~vmask] = 1

        return fr


def erfc(x):
    return 1.0 - erf(x)


# End of file
