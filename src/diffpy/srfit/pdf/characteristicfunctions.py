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
"""Form factors (characteristic functions) used in PDF nanoshape
fitting.

These are used to calculate the attenuation of the PDF due to a finite
size. For a crystal-like nanoparticle, one can calculate the PDF via
Gnano(r) = f(r) Gcryst(r), where f(r) is the nanoparticle characteristic
function and Gcryst(f) is the crystal PDF.

These functions are meant to be imported and added to a FitContribution
using the 'register_function' method of that class.
"""

__all__ = [
    "spherical_particle",
    "spheroidal_particle",
    "lognormal_spherical_particle",
    "sheet_particle",
    "shell_particle",
    "SASCF",
    "sphericalCF",
    "spheroidalCF",
    "spheroidalCF2",
    "lognormalSphericalCF",
    "sheetCF",
    "shellCF",
    "shellCF2",
]

import warnings

import numpy
from numpy import arctan as atan
from numpy import arctanh as atanh
from numpy import ceil, exp, log, log2, pi, sign, sqrt
from numpy.fft import fftfreq, ifft
from scipy.special import erf

from diffpy.srfit.fitbase.calculator import Calculator
from diffpy.utils._deprecator import build_deprecation_message, deprecated

removal_version = "4.0.0"
cf_base = "diffpy.srfit.pdf.characteristicfunctions"


def _build_dep_msg(old_name, new_name, signature_note=None):
    """Build the deprecation message for a camel case function.

    The note describing how the signature changed, when there is one, is
    appended to the standard message from `build_deprecation_message`.
    """
    message = build_deprecation_message(
        cf_base, old_name, new_name, removal_version
    )
    if signature_note is None:
        return message
    return f"{message} {signature_note}"


sphericalCF_dep_msg = _build_dep_msg(
    "sphericalCF",
    "spherical_particle",
    "Additionally, the signature has changed. Please pass the parameter "
    "'psize' as 'particle_diameter'.",
)

spheroidalCF_dep_msg = _build_dep_msg(
    "spheroidalCF",
    "spheroidal_particle",
    "Additionally, the signature has changed. Please pass the parameters "
    "'erad' and 'prad' as 'equatorial_radius' and 'polar_radius', "
    "respectively.",
)

spheroidalCF2_dep_msg = _build_dep_msg(
    "spheroidalCF2",
    "spheroidal_particle",
    "Additionally, the parameterization has changed. 'spheroidalCF2' took "
    "the equatorial diameter 'psize' and the axis ratio 'axrat', while "
    "'spheroidal_particle' takes radii. Please pass "
    "equatorial_radius = psize / 2 and polar_radius = axrat * psize / 2.",
)

lognormalSphericalCF_dep_msg = _build_dep_msg(
    "lognormalSphericalCF",
    "lognormal_spherical_particle",
    "Additionally, the signature has changed. Please pass the parameters "
    "'psize' and 'psig' as 'particle_diameter' and "
    "'particle_diameter_sigma', respectively.",
)

sheetCF_dep_msg = _build_dep_msg(
    "sheetCF",
    "sheet_particle",
    "Additionally, the signature has changed. Please pass the parameter "
    "'sthick' as 'sheet_thickness'.",
)

shellCF_dep_msg = _build_dep_msg("shellCF", "shell_particle")

shellCF2_dep_msg = _build_dep_msg(
    "shellCF2",
    "shell_particle",
    "Additionally, the parameterization has changed. 'shellCF2' took the "
    "central radius 'a' and the shell thickness 'delta', while "
    "'shell_particle' takes the inner radius. Please pass "
    "radius = a - delta / 2 and thickness = delta.",
)


# The record of non-physical input warnings already issued in this process,
# keyed by (function name, requirement). See _warn_non_physical.
_warned_non_physical = set()


def _warn_non_physical(function_name, requirement):
    """Warn that non-physical input produced a zero characteristic
    function.

    A refinement evaluates a characteristic function on every iteration,
    so the warning is issued only once per process for each distinct
    requirement. The message must not include the offending value, since
    the `warnings` registry keys on the message text and a value that
    drifts from iteration to iteration would defeat the suppression.
    """
    key = (function_name, requirement)
    if key in _warned_non_physical:
        return
    _warned_non_physical.add(key)
    warnings.warn(
        f"In '{function_name}', {requirement}. The characteristic function "
        "was set to zero for every r, which flattens the fit residual so a "
        "refinement cannot recover on its own. Please use a physically "
        "meaningful starting value, or keep the parameter in range with "
        "FitRecipe.add_soft_bounds.",
        RuntimeWarning,
        stacklevel=3,
    )


def spherical_particle(r, particle_diameter):
    """Compute the spherical nanoparticle characteristic function.

    From Kodama et al., Acta Cryst. A, 62, 444-453 (converted from
    radius to diameter).

    Parameters
    ----------
    r : array_like
        The distance of interaction.
    particle_diameter : float
        The particle diameter.

    Returns
    -------
    numpy.ndarray
        The characteristic function values evaluated at `r`.

    Warns
    -----
    RuntimeWarning
        If `particle_diameter` is not positive, in which case the
        characteristic function is zero everywhere. The warning is issued
        only once per process.
    """
    characteristic_function = numpy.zeros(numpy.shape(r), dtype=float)
    if particle_diameter <= 0:
        _warn_non_physical(
            "spherical_particle", "'particle_diameter' must be positive"
        )
        return characteristic_function
    scaled_r = numpy.array(r, dtype=float) / particle_diameter
    inside = scaled_r < 1.0
    scaled_r_inside = scaled_r[inside]
    characteristic_function[inside] = (
        1.0 - 1.5 * scaled_r_inside + 0.5 * scaled_r_inside**3
    )
    return characteristic_function


def spheroidal_particle(r, equatorial_radius, polar_radius):
    """Compute the spheroidal nanoparticle characteristic function.

    Spheroid with radii (equatorial_radius, equatorial_radius,
    polar_radius). ``equatorial_radius < polar_radius`` equates to a
    prolate spheroid, ``equatorial_radius > polar_radius`` equates to
    an oblate spheroid, and ``equatorial_radius == polar_radius`` is a
    sphere. From Lei et al., Phys. Rev. B, 80, 024118 (2009).

    Parameters
    ----------
    r : array_like
        The distance of interaction.
    equatorial_radius : float
        The equatorial radius.
    polar_radius : float
        The polar radius.

    Returns
    -------
    numpy.ndarray
        The characteristic function values evaluated at `r`.

    Warns
    -----
    RuntimeWarning
        If `equatorial_radius` or `polar_radius` is not positive, in which
        case the characteristic function is zero everywhere. The warning is
        issued only once per process.
    """
    if equatorial_radius <= 0:
        _warn_non_physical(
            "spheroidal_particle", "'equatorial_radius' must be positive"
        )
        return numpy.zeros(numpy.shape(r), dtype=float)
    if polar_radius <= 0:
        _warn_non_physical(
            "spheroidal_particle", "'polar_radius' must be positive"
        )
        return numpy.zeros(numpy.shape(r), dtype=float)

    particle_diameter = 2.0 * equatorial_radius
    axis_ratio = 1.0 * polar_radius / equatorial_radius

    # to simplify the equations
    v = axis_ratio
    d = particle_diameter
    d2 = d * d
    v2 = v * v

    if v == 1:
        return spherical_particle(r, particle_diameter)

    rx = r
    if v < 1:

        r = rx[rx <= v * d]
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

        r = rx[numpy.logical_and(rx > v * d, rx <= d)]
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

        r = rx[rx > d]
        f3 = numpy.zeros_like(r)

        f = numpy.concatenate((f1, f2, f3))

    elif v > 1:

        r = rx[rx <= d]
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

        r = rx[numpy.logical_and(rx > d, rx <= v * d)]
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

        r = rx[rx > v * d]
        f3 = numpy.zeros_like(r)

        f = numpy.concatenate((f1, f2, f3))

    return f


def lognormal_spherical_particle(
    r, particle_diameter, particle_diameter_sigma
):
    """Compute the spherical nanoparticle characteristic function with
    lognormal size distribution.

    Here, r is the independent variable, mu is the mean of the
    distribution (not of the particle size), and s is the width of
    the distribution. This is the characteristic function for the
    lognormal distribution of particle diameter:

    F(r, mu, s) = 0.5*Erfc((-mu-3*s^2+Log(r))/(sqrt(2)*s))
               + 0.25*r^3*Erfc((-mu+Log(r))/(sqrt(2)*s))*exp(-3*mu-4.5*s^2)
               - 0.75*r*Erfc((-mu-2*s^2+Log(r))/(sqrt(2)*s))*exp(-mu-2.5*s^2)

    The expectation value of the distribution gives the average particle
    diameter, particle_diameter. The variance of the distribution gives
    particle_diameter_sigma^2. mu and s can be expressed in terms of these
    as:

    s^2 = log((particle_diameter_sigma/particle_diameter)^2 + 1)
    mu = log(particle_diameter) - s^2/2

    Source unknown.

    Parameters
    ----------
    r : array_like
        The distance of interaction.
    particle_diameter : float
        The mean particle diameter.
    particle_diameter_sigma : float
        The log-normal width of the particle diameter.

    Returns
    -------
    numpy.ndarray
        The characteristic function values evaluated at `r`.

    Warns
    -----
    RuntimeWarning
        If `particle_diameter` is not positive, if `particle_diameter_sigma`
        is negative, or if `particle_diameter` is too small for
        `particle_diameter_sigma` for the closed form to hold, in which case
        the characteristic function is zero everywhere. The warning is issued
        only once per process. A `particle_diameter_sigma` of zero is the
        sphere limit rather than an error and does not warn.
    """
    if particle_diameter <= 0:
        _warn_non_physical(
            "lognormal_spherical_particle",
            "'particle_diameter' must be positive",
        )
        return numpy.zeros(numpy.shape(r), dtype=float)
    if particle_diameter_sigma < 0:
        _warn_non_physical(
            "lognormal_spherical_particle",
            "'particle_diameter_sigma' must not be negative",
        )
        return numpy.zeros(numpy.shape(r), dtype=float)
    if particle_diameter_sigma == 0:
        return spherical_particle(r, particle_diameter)

    sqrt2 = sqrt(2.0)
    s = sqrt(
        log(
            particle_diameter_sigma
            * particle_diameter_sigma
            / (1.0 * particle_diameter * particle_diameter)
            + 1
        )
    )
    mu = log(particle_diameter) - s * s / 2
    if mu < 0:
        _warn_non_physical(
            "lognormal_spherical_particle",
            "'particle_diameter' is too small for 'particle_diameter_sigma'",
        )
        return numpy.zeros(numpy.shape(r), dtype=float)

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


def sheet_particle(r, sheet_thickness):
    """Compute the nanosheet characteristic function.

    From Kodama et al., Acta Cryst. A, 62, 444-453.

    Parameters
    ----------
    r : array_like
        The distance of interaction.
    sheet_thickness : float
        The thickness of the nanosheet.

    Returns
    -------
    numpy.ndarray
        The characteristic function values evaluated at `r`.

    Warns
    -----
    RuntimeWarning
        If `sheet_thickness` is not positive, in which case the
        characteristic function is zero everywhere. The warning is issued
        only once per process.
    """
    if sheet_thickness <= 0:
        _warn_non_physical(
            "sheet_particle", "'sheet_thickness' must be positive"
        )
        if numpy.isscalar(r):
            return 0.0
        return numpy.zeros(numpy.shape(r), dtype=float)
    # process scalar r
    if numpy.isscalar(r):
        rv = (
            1 - 0.5 * r / sheet_thickness
            if r < sheet_thickness
            else 0.5 * sheet_thickness / r
        )
        return rv
    # handle array-type r
    r_array = numpy.asarray(r)
    inside = r_array < sheet_thickness
    outside = ~inside
    characteristic_function = numpy.empty_like(r_array, dtype=float)
    characteristic_function[inside] = (
        1 - 0.5 * r_array[inside] / sheet_thickness
    )
    characteristic_function[outside] = 0.5 * sheet_thickness / r_array[outside]
    return characteristic_function


def shell_particle(r, radius, thickness):
    """Compute the spherical shell characteristic function.

    The outer radius equals ``radius + thickness``. From Lei et al.,
    Phys. Rev. B, 80, 024118 (2009).

    Parameters
    ----------
    r : array_like
        The distance of interaction.
    radius : float
        The inner radius.
    thickness : float
        The thickness of the shell.

    Returns
    -------
    numpy.ndarray
        The characteristic function values evaluated at `r`.

    Warns
    -----
    RuntimeWarning
        If `thickness` is not positive or `radius` is negative, in which case
        the characteristic function is zero everywhere. The warning is issued
        only once per process. A `radius` of zero is the solid sphere limit
        rather than an error and does not warn.
    """
    if thickness <= 0:
        _warn_non_physical("shell_particle", "'thickness' must be positive")
        return numpy.zeros(numpy.shape(r), dtype=float)
    if radius < 0:
        _warn_non_physical("shell_particle", "'radius' must not be negative")
        return numpy.zeros(numpy.shape(r), dtype=float)

    d = 1.0 * thickness
    a = 1.0 * radius + d / 2.0
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


@deprecated(sphericalCF_dep_msg)
def sphericalCF(r, psize):
    """This function has been deprecated and will be removed in version
    4.0.0.

    Please use diffpy.srfit.pdf.characteristicfunctions.spherical_particle
    instead.
    """
    return spherical_particle(r, psize)


@deprecated(spheroidalCF_dep_msg)
def spheroidalCF(r, erad, prad):
    """This function has been deprecated and will be removed in version
    4.0.0.

    Please use diffpy.srfit.pdf.characteristicfunctions.spheroidal_particle
    instead.
    """
    return spheroidal_particle(r, erad, prad)


@deprecated(spheroidalCF2_dep_msg)
def spheroidalCF2(r, psize, axrat):
    """This function has been deprecated and will be removed in version
    4.0.0.

    Please use diffpy.srfit.pdf.characteristicfunctions.spheroidal_particle
    instead.
    """
    return spheroidal_particle(r, psize / 2, axrat * psize / 2)


@deprecated(lognormalSphericalCF_dep_msg)
def lognormalSphericalCF(r, psize, psig):
    """This function has been deprecated and will be removed in version
    4.0.0.

    Please use
    diffpy.srfit.pdf.characteristicfunctions.lognormal_spherical_particle
    instead.
    """
    return lognormal_spherical_particle(r, psize, psig)


@deprecated(sheetCF_dep_msg)
def sheetCF(r, sthick):
    """This function has been deprecated and will be removed in version
    4.0.0.

    Please use diffpy.srfit.pdf.characteristicfunctions.sheet_particle
    instead.
    """
    return sheet_particle(r, sthick)


@deprecated(shellCF_dep_msg)
def shellCF(r, radius, thickness):
    """This function has been deprecated and will be removed in version
    4.0.0.

    Please use diffpy.srfit.pdf.characteristicfunctions.shell_particle
    instead.
    """
    return shell_particle(r, radius, thickness)


@deprecated(shellCF2_dep_msg)
def shellCF2(r, a, delta):
    """This function has been deprecated and will be removed in version
    4.0.0.

    Please use diffpy.srfit.pdf.characteristicfunctions.shell_particle
    instead.
    """
    return shell_particle(r, a - delta / 2, delta)


class SASCF(Calculator):
    """Calculator class for characteristic functions from sas-models.

    This class wraps a sas.models.BaseModel to calculate I(Q) related to
    nanoparticle shape. This I(Q) is inverted to f(r) according to:
    f(r) = 1 / (4 pi r) * SINFT(I(Q)),
    where "SINFT" represents the sine Fourier transform.

    Attributes
    ----------
    _model
        BaseModel object this adapts.
    Managed Parameters
        These depend on the parameters of the BaseModel object held by _model.
        They are created from the 'params' attribute of the BaseModel. If a
        dispersion is set for the BaseModel, the dispersion "width" will be
        accessible under "<parname>_width", where <parname> is the name a
        parameter adjusted by dispersion.
    """

    def __init__(self, name, model):
        """Initialize the generator.

        Parameters
        ----------
        name : str
            The name for the SASCF.
        model : BaseModel
            The sas.models.BaseModel object this adapts.
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
        """Calculate the characteristic function from the transform of
        the BaseModel."""
        # Determine q-values.
        # We want very fine r-spacing so we can properly normalize f(r). This
        # equates to having a large qmax so that the Fourier transform is
        # finely spaced. We also want the calculation to be fast, so we pick
        # qmax such that the number of q-points is a power of 2. This allows us
        # to use the fft.
        #
        # The initial dr is somewhat arbitrary, but using dr = 0.01 allows for
        # the f(r) calculated from a particle of diameter 50, over r =
        # arange(1, 60, 0.1) to agree with the spherical_particle with Rw <
        # 1e-4%.
        #
        # We also have to make a q-spacing small enough to compute out to at
        # least the size of the signal.
        raise NotImplementedError(
            "As of release 3.2.0, SAS characteristic functions are not working"
            + " but we hope to have them working again in a future release."
        )

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
    """Compute the complementary error function.

    Parameters
    ----------
    x : array_like
        The input value.

    Returns
    -------
    numpy.ndarray
        The complementary error function evaluated at `x`.
    """
    return 1.0 - erf(x)


# End of file
