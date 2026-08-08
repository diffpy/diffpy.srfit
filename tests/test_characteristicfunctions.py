#!/usr/bin/env python
##############################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2010 The Trustees of Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################
"""Tests for the nanoparticle characteristic functions and the SAS
calculator wrapper."""

# FIXME: remove this line when `docformatter` fixes the blank line bug


<<<<<<< HEAD
import re
=======
>>>>>>> origin/main
import unittest
import warnings

import numpy
import numpy.testing as npt
import pytest

try:
    from sasmodels.sasview_model import find_model, load_standard_models
except ModuleNotFoundError:
    pass

import diffpy.srfit.pdf.characteristicfunctions as cf

# # Global variables to be assigned in setUp
# cf = None
# Fixme: remove this code if imports don't break on testing

# ----------------------------------------------------------------------------


def testSphere(sas_available):
    # if not sas_available:
    # pytest.skip("sas package not available")
    pytest.skip(
        "sas characteristic functions not currently working, "
        + "remove skip when our code is refactored to use the "
        + "latest sasview API"
    )

    load_standard_models()
    radius = 25
    # Calculate sphere cf from SphereModel
    SphereModel = find_model("sphere")
    model = SphereModel()
    model.setParam("radius", radius)
    ff = cf.SASCF("sphere", model)
    r = numpy.arange(1, 60, 0.1, dtype=float)
    fr1 = ff(r)

    # Calculate sphere cf analytically
    fr2 = cf.spherical_particle(r, 2 * radius)
    diff = fr1 - fr2
    res = numpy.dot(diff, diff)
    res /= numpy.dot(fr2, fr2)
    assert res == pytest.approx(0, abs=1e-4)
    return


def testSpheroid(sas_available):
    # if not sas_available:
    # pytest.skip("sas package not available")
    pytest.skip(
        "sas characteristic functions not currently working, "
        + "remove skip when our code is refactored to use the "
        + "latest sasview API"
    )

    load_standard_models()
    prad = 20.9
    erad = 33.114
    # Calculate cf from EllipsoidModel
    EllipsoidModel = find_model("ellipsoid")
    model = EllipsoidModel()
    model.setParam("radius_polar", prad)
    model.setParam("radius_equatorial", erad)
    ff = cf.SASCF("spheroid", model)
    r = numpy.arange(0, 100, 1 / numpy.pi, dtype=float)
    fr1 = ff(r)

    # Calculate cf analytically
    fr2 = cf.spheroidal_particle(r, erad, prad)
    diff = fr1 - fr2
    res = numpy.dot(diff, diff)
    res /= numpy.dot(fr2, fr2)
    assert res == pytest.approx(0, abs=1e-4)
    return


def testShell(sas_available):
    # if not sas_available:
    # pytest.skip("sas package not available")
    pytest.skip(
        "sas characteristic functions not currently working, "
        + "remove skip when our code is refactored to use the "
        + "latest sasview API"
    )

    load_standard_models()
    radius = 19.2
    thickness = 7.8
    # Calculate cf from VesicleModel
    VesicleModel = find_model("vesicle")
    model = VesicleModel()
    model.setParam("radius", radius)
    model.setParam("thickness", thickness)
    ff = cf.SASCF("vesicle", model)
    r = numpy.arange(0, 99.45, 0.1, dtype=float)
    fr1 = ff(r)

    # Calculate sphere cf analytically
    fr2 = cf.shell_particle(r, radius, thickness)
    diff = fr1 - fr2
    res = numpy.dot(diff, diff)
    res /= numpy.dot(fr2, fr2)
    assert res == pytest.approx(0, abs=1e-4)
    return


def testCylinder(sas_available):
    # if not sas_available:
    # pytest.skip("sas package not available")
    pytest.skip(
        "sas characteristic functions not currently working, "
        + "remove skip when our code is refactored to use the "
        + "latest sasview API"
    )
    load_standard_models()
    # Make sure cylinder works over different r-ranges.
    radius = 100
    length = 30

    CylinderModel = find_model("cylinder")
    model = CylinderModel()
    model.setParam("radius", radius)
    model.setParam("length", length)

    ff = cf.SASCF("cylinder", model)

    r1 = numpy.arange(0, 10, 0.1, dtype=float)
    r2 = numpy.arange(0, 50, 0.1, dtype=float)
    r3 = numpy.arange(0, 100, 0.1, dtype=float)
    r4 = numpy.arange(0, 500, 0.1, dtype=float)

    fr1 = ff(r1)
    fr2 = ff(r2)
    fr3 = ff(r3)
    fr4 = ff(r4)

    d = fr1 - numpy.interp(r1, r2, fr2)
    res12 = numpy.dot(d, d)
    res12 /= numpy.dot(fr1, fr1)
    assert res12 == pytest.approx(0, abs=1e-4)

    d = fr1 - numpy.interp(r1, r3, fr3)
    res13 = numpy.dot(d, d)
    res13 /= numpy.dot(fr1, fr1)
    assert res13 == pytest.approx(0, abs=1e-4)

    d = fr1 - numpy.interp(r1, r4, fr4)
    res14 = numpy.dot(d, d)
    res14 /= numpy.dot(fr1, fr1)
    assert res14 == pytest.approx(0, abs=1e-4)

    d = fr2 - numpy.interp(r2, r3, fr3)
    res23 = numpy.dot(d, d)
    res23 /= numpy.dot(fr2, fr2)
    assert res23 == pytest.approx(0, abs=1e-4)

    d = fr2 - numpy.interp(r2, r4, fr4)
    res24 = numpy.dot(d, d)
    res24 /= numpy.dot(fr2, fr2)
    assert res24 == pytest.approx(0, abs=1e-4)

    d = fr3 - numpy.interp(r3, r4, fr4)
    res34 = numpy.dot(d, d)
    res34 /= numpy.dot(fr3, fr3)
    assert res34 == pytest.approx(0, abs=1e-4)
    return


# End of class TestSASCF

# ----------------------------------------------------------------------------
# The characteristic functions below compute the attenuation of the PDF due
# to a finite nanoparticle shape. Each has a known value of 1 at r = 0 (the
# particle always overlaps itself at zero distance) and drops to 0 once r
# reaches the largest distance spanned by the shape.


@pytest.mark.parametrize(
    "input_r, input_particle_diameter, expected_characteristic_function",
    [
        # C1: r ranges over values inside and outside the particle.
        # Expected: the polynomial form 1 - 1.5x + 0.5x^3 (x = r /
        # particle_diameter) applies inside the particle, and the function
        # is zero at and beyond the particle diameter.
        (
            numpy.array([0.0, 5.0, 10.0, 15.0]),
            10.0,
            [1.0, 0.3125, 0.0, 0.0],
        ),
    ],
)
def test_spherical_particle(
    input_r, input_particle_diameter, expected_characteristic_function
):
    actual_characteristic_function = cf.spherical_particle(
        input_r, input_particle_diameter
    )
    npt.assert_allclose(
        actual_characteristic_function, expected_characteristic_function
    )


@pytest.mark.parametrize(
    "input_r, input_equatorial_radius, input_polar_radius, "
    "expected_characteristic_function",
    [
        # C1: equatorial_radius equals polar_radius, i.e. a sphere.
        # Expected: matches spherical_particle for the equivalent diameter.
        (
            numpy.linspace(0, 25, 6),
            10.0,
            10.0,
            cf.spherical_particle(numpy.linspace(0, 25, 6), 20.0),
        ),
        # C2: polar_radius > equatorial_radius, i.e. a prolate spheroid.
        # Expected: values from the Lei et al. characteristic function.
        (
            numpy.array([0.0, 10.0, 20.0, 30.0, 40.0, 50.0]),
            10.0,
            20.0,
            [
                1.0,
                0.40106264900760513,
                0.054200238412168256,
                0.003650768214115807,
                0.0,
                0.0,
            ],
        ),
        # C3: polar_radius < equatorial_radius, i.e. an oblate spheroid.
        # Expected: values from the Lei et al. characteristic function.
        (
            numpy.array([0.0, 5.0, 10.0, 15.0, 20.0, 25.0]),
            20.0,
            10.0,
            [
                1.0,
                0.744181556741916,
                0.5061470768546105,
                0.30368052370886184,
                0.15456586067544859,
                0.0675583474738014,
            ],
        ),
    ],
)
def test_spheroidal_particle(
    input_r,
    input_equatorial_radius,
    input_polar_radius,
    expected_characteristic_function,
):
    actual_characteristic_function = cf.spheroidal_particle(
        input_r, input_equatorial_radius, input_polar_radius
    )
    npt.assert_allclose(
        actual_characteristic_function, expected_characteristic_function
    )


@pytest.mark.parametrize(
    "input_r, input_particle_diameter, input_particle_diameter_sigma, "
    "expected_characteristic_function",
    [
        # C1: the lognormal distribution has no width.
        # Expected: matches spherical_particle for the same diameter.
        (
            numpy.array([0.0, 5.0, 10.0]),
            10.0,
            0.0,
            cf.spherical_particle(numpy.array([0.0, 5.0, 10.0]), 10.0),
        ),
        # C2: a genuine lognormal size distribution is applied.
        # Expected: values from the closed-form lognormal-averaged
        # spherical characteristic function.
        (
            numpy.array([1.0, 5.0, 10.0, 15.0]),
            10.0,
            2.0,
            [
                0.8617610662266728,
                0.362144891493816,
                0.039068792680198694,
                0.0009056141323687261,
            ],
        ),
    ],
)
def test_lognormal_spherical_particle(
    input_r,
    input_particle_diameter,
    input_particle_diameter_sigma,
    expected_characteristic_function,
):
    actual_characteristic_function = cf.lognormal_spherical_particle(
        input_r, input_particle_diameter, input_particle_diameter_sigma
    )
    npt.assert_allclose(
        actual_characteristic_function, expected_characteristic_function
    )


@pytest.mark.parametrize(
    "input_r, input_sheet_thickness, expected_characteristic_function",
    [
        # C1: r is an array with values inside and outside the sheet
        # thickness.
        # Expected: 1 - 0.5*r/thickness inside the sheet, and
        # 0.5*thickness/r beyond it.
        (
            numpy.array([0.0, 2.0, 4.0, 8.0]),
            4.0,
            [1.0, 0.75, 0.5, 0.25],
        ),
        # C2: r is a scalar inside the sheet thickness.
        # Expected: 1 - 0.5*r/thickness.
        (2.0, 4.0, 0.75),
    ],
)
def test_sheet_particle(
    input_r, input_sheet_thickness, expected_characteristic_function
):
    actual_characteristic_function = cf.sheet_particle(
        input_r, input_sheet_thickness
    )
    npt.assert_allclose(
        actual_characteristic_function, expected_characteristic_function
    )


@pytest.mark.parametrize(
    "input_r, input_radius, input_thickness, "
    "expected_characteristic_function",
    [
        # C1: r ranges from the shell center out past the outer radius.
        # Expected: values from the Lei et al. shell characteristic
        # function; the function is 1 at r = 0 and 0 once r exceeds the
        # outer radius (radius + thickness).
        (
            numpy.array([0.0, 5.0, 12.5, 15.0, 20.0, 30.0]),
            10.0,
            5.0,
            [
                1.0,
                0.4934210526315789,
                0.19736842105263158,
                0.16447368421052633,
                0.12335526315789473,
                0.0,
            ],
        ),
    ],
)
def test_shell_particle(
    input_r, input_radius, input_thickness, expected_characteristic_function
):
    actual_characteristic_function = cf.shell_particle(
        input_r, input_radius, input_thickness
    )
    npt.assert_allclose(
        actual_characteristic_function, expected_characteristic_function
    )


# ----------------------------------------------------------------------------
# A refinement can step a shape parameter into a region that has no physical
# meaning, such as a negative diameter or thickness. The characteristic
# functions must not raise there, because that would abort a refinement that
# would otherwise recover on its own, so they return zero for every r instead.
# A zero return is flat in the shape parameter and so stalls the optimizer,
# which the user cannot see from the fit output alone; each function therefore
# also warns. The warning is issued only once per process for each distinct
# problem so that a refinement loop does not bury the user in repeated
# messages.
#
# Every such warning states which requirement was violated and then closes
# with the same explanation of the consequence and of how to avoid it. Only
# the requirement varies between cases.


@pytest.mark.parametrize(
    "function_name, input_args, expected_requirement",
    [
        # C1: the sphere diameter is zero.
        # Expected: zero everywhere, warning that the diameter must be
        # positive.
        (
            "spherical_particle",
            (numpy.array([0.0, 5.0]), 0.0),
            "'particle_diameter' must be positive",
        ),
        # C2: the sphere diameter is negative.
        # Expected: zero everywhere, warning that the diameter must be
        # positive.
        (
            "spherical_particle",
            (numpy.array([1.0, 5.0, 10.0]), -10.0),
            "'particle_diameter' must be positive",
        ),
        # C3: the spheroid equatorial radius is zero, which the axis ratio
        # would otherwise be divided by.
        # Expected: zero everywhere, warning that the equatorial radius must
        # be positive, and no ZeroDivisionError.
        (
            "spheroidal_particle",
            (numpy.array([1.0, 5.0]), 0.0, 5.0),
            "'equatorial_radius' must be positive",
        ),
        # C4: the spheroid polar radius is negative.
        # Expected: zero everywhere, warning that the polar radius must be
        # positive.
        (
            "spheroidal_particle",
            (numpy.array([1.0, 5.0]), 10.0, -5.0),
            "'polar_radius' must be positive",
        ),
        # C5: the mean particle diameter of the lognormal distribution is
        # zero.
        # Expected: zero everywhere, warning that the diameter must be
        # positive.
        (
            "lognormal_spherical_particle",
            (numpy.array([1.0, 2.0]), 0.0, 1.0),
            "'particle_diameter' must be positive",
        ),
        # C6: the width of the lognormal distribution is negative. A negative
        # width is meaningless rather than a limiting case, so it must not be
        # quietly treated as a distribution of zero width.
        # Expected: zero everywhere, warning that the width must not be
        # negative.
        (
            "lognormal_spherical_particle",
            (numpy.array([1.0, 5.0, 10.0]), 10.0, -2.0),
            "'particle_diameter_sigma' must not be negative",
        ),
        # C7: the mean particle diameter is too small for the requested
        # distribution width, which puts the underlying lognormal mean below
        # zero and invalidates the closed form.
        # Expected: zero everywhere, warning that the diameter is too small
        # for the width.
        (
            "lognormal_spherical_particle",
            (numpy.array([1.0, 5.0]), 2.0, 4.0),
            "'particle_diameter' is too small for 'particle_diameter_sigma'",
        ),
        # C8: the sheet thickness is zero.
        # Expected: zero everywhere, warning that the thickness must be
        # positive.
        (
            "sheet_particle",
            (numpy.array([1.0, 2.0]), 0.0),
            "'sheet_thickness' must be positive",
        ),
        # C9: the sheet thickness is negative.
        # Expected: zero everywhere, warning that the thickness must be
        # positive.
        (
            "sheet_particle",
            (numpy.array([1.0, 2.0]), -1.0),
            "'sheet_thickness' must be positive",
        ),
        # C10: the shell thickness is zero, which would otherwise make the
        # normalizing denominator vanish and return an unattenuated one at
        # every r.
        # Expected: zero everywhere, warning that the thickness must be
        # positive.
        (
            "shell_particle",
            (numpy.array([1.0, 5.0, 10.0]), 10.0, 0.0),
            "'thickness' must be positive",
        ),
        # C11: the shell thickness is negative, which would otherwise return
        # a smoothly varying branch with values below negative one that an
        # optimizer can mistake for a real solution.
        # Expected: zero everywhere, warning that the thickness must be
        # positive.
        (
            "shell_particle",
            (numpy.array([1.0, 5.0, 10.0]), 10.0, -5.0),
            "'thickness' must be positive",
        ),
        # C12: the shell inner radius is negative.
        # Expected: zero everywhere, warning that the radius must not be
        # negative.
        (
            "shell_particle",
            (numpy.array([1.0, 5.0, 10.0]), -20.0, 5.0),
            "'radius' must not be negative",
        ),
    ],
)
def test_non_physical_input_returns_zero_and_warns(
    function_name,
    input_args,
    expected_requirement,
    reset_characteristic_function_warnings,
):
    function = getattr(cf, function_name)
    input_r = input_args[0]
    expected_remedy = (
        "The characteristic function was set to zero for every r, which "
        "flattens the fit residual so a refinement cannot recover on its own. "
        "Please use a physically meaningful starting value, or keep the "
        "parameter in range with FitRecipe.add_soft_bounds."
    )
    expected_msg = (
        f"In '{function_name}', {expected_requirement}. {expected_remedy}"
    )
    with pytest.warns(RuntimeWarning, match=re.escape(expected_msg)):
        actual_characteristic_function = function(*input_args)
    npt.assert_allclose(
        actual_characteristic_function, numpy.zeros_like(input_r)
    )


@pytest.mark.parametrize(
    "function_name, input_args, expected_characteristic_function",
    [
        # C1: a lognormal distribution of zero width is the sphere limit
        # rather than a mistake.
        # Expected: matches spherical_particle for the same diameter, with no
        # warning.
        (
            "lognormal_spherical_particle",
            (numpy.array([1.0, 5.0, 10.0]), 10.0, 0.0),
            cf.spherical_particle(numpy.array([1.0, 5.0, 10.0]), 10.0),
        ),
        # C2: a shell of zero inner radius is a solid sphere rather than a
        # mistake.
        # Expected: matches spherical_particle for twice the thickness, with
        # no warning.
        (
            "shell_particle",
            (numpy.array([1.0, 5.0, 10.0]), 0.0, 5.0),
            cf.spherical_particle(numpy.array([1.0, 5.0, 10.0]), 10.0),
        ),
    ],
)
def test_limiting_input_does_not_warn(
    function_name,
    input_args,
    expected_characteristic_function,
    reset_characteristic_function_warnings,
):
    function = getattr(cf, function_name)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        actual_characteristic_function = function(*input_args)
    npt.assert_allclose(
        actual_characteristic_function, expected_characteristic_function
    )


def test_non_physical_input_warns_only_once(
    reset_characteristic_function_warnings,
):
    # A refinement evaluates the characteristic function on every iteration,
    # so a function that warns on each call would flood the user. The
    # suppression must not rely on the warnings registry, which callers can
    # disable, so this asserts under the "always" filter.
    with warnings.catch_warnings(record=True) as raised_warnings:
        warnings.simplefilter("always")
        for _ in range(5):
            cf.spherical_particle(numpy.array([1.0, 2.0]), -10.0)
    actual_warning_count = len(raised_warnings)
    expected_warning_count = 1
    assert actual_warning_count == expected_warning_count


# ----------------------------------------------------------------------------
# sphericalCF, spheroidalCF, spheroidalCF2, lognormalSphericalCF, sheetCF,
# shellCF, and shellCF2 are deprecated in favor of the snake_case functions
# above, and forward to them directly (translating arguments where the old
# and new parameterizations differ). Consequently they also inherit the
# replacement's handling of non-physical input in full: values that the old,
# pre-3.3.0 implementation used to accept silently, or handled with its own
# quirks, now warn and return zero exactly like the replacement. Each old
# name must still work, forward the (translated) arguments, and emit a
# DeprecationWarning naming its replacement.

# The warning must also tell the caller how to translate the old arguments,
# so each case pins the exact note appended to the standard message. A case
# whose parameters did not change appends nothing.


@pytest.mark.parametrize(
    "old_name, new_name, expected_signature_note, old_args, new_args",
    [
        # C1: sphericalCF forwards its arguments unchanged, only renamed.
        # Expected: The note names the renamed parameter.
        (
            "sphericalCF",
            "spherical_particle",
            "Additionally, the signature has changed. Please pass the "
            "parameter 'psize' as 'particle_diameter'.",
            (numpy.array([0.0, 5.0, 10.0]), 10.0),
            (numpy.array([0.0, 5.0, 10.0]), 10.0),
        ),
        # C2: spheroidalCF forwards its arguments unchanged, only renamed.
        # Expected: The note names both renamed parameters.
        (
            "spheroidalCF",
            "spheroidal_particle",
            "Additionally, the signature has changed. Please pass the "
            "parameters 'erad' and 'prad' as 'equatorial_radius' and "
            "'polar_radius', respectively.",
            (numpy.array([0.0, 5.0, 10.0]), 10.0, 15.0),
            (numpy.array([0.0, 5.0, 10.0]), 10.0, 15.0),
        ),
        # C3: spheroidalCF2 used the raw (diameter, axis ratio)
        # parameterization, which must be converted to (equatorial_radius,
        # polar_radius) before forwarding to spheroidal_particle.
        # Expected: The note gives the conversion, not a rename.
        (
            "spheroidalCF2",
            "spheroidal_particle",
            "Additionally, the parameterization has changed. 'spheroidalCF2' "
            "took the equatorial diameter 'psize' and the axis ratio "
            "'axrat', while 'spheroidal_particle' takes radii. Please pass "
            "equatorial_radius = psize / 2 and "
            "polar_radius = axrat * psize / 2.",
            (numpy.array([0.0, 5.0, 10.0]), 20.0, 1.5),
            (numpy.array([0.0, 5.0, 10.0]), 10.0, 15.0),
        ),
        # C4: lognormalSphericalCF forwards its arguments unchanged, only
        # renamed.
        # Expected: The note names both renamed parameters.
        (
            "lognormalSphericalCF",
            "lognormal_spherical_particle",
            "Additionally, the signature has changed. Please pass the "
            "parameters 'psize' and 'psig' as 'particle_diameter' and "
            "'particle_diameter_sigma', respectively.",
            (numpy.array([1.0, 5.0, 10.0]), 10.0, 2.0),
            (numpy.array([1.0, 5.0, 10.0]), 10.0, 2.0),
        ),
        # C5: sheetCF forwards its arguments unchanged, only renamed.
        # Expected: The note names the renamed parameter.
        (
            "sheetCF",
            "sheet_particle",
            "Additionally, the signature has changed. Please pass the "
            "parameter 'sthick' as 'sheet_thickness'.",
            (numpy.array([0.0, 2.0, 4.0]), 4.0),
            (numpy.array([0.0, 2.0, 4.0]), 4.0),
        ),
        # C6: shellCF, whose parameters kept their names and meanings.
        # Expected: The standard message with no note appended.
        (
            "shellCF",
            "shell_particle",
            "",
            (numpy.array([0.0, 5.0, 12.5]), 10.0, 5.0),
            (numpy.array([0.0, 5.0, 12.5]), 10.0, 5.0),
        ),
        # C7: shellCF2 used the raw (central_radius, thickness)
        # parameterization, which must be converted to (radius, thickness)
        # before forwarding to shell_particle.
        # Expected: The note gives the conversion, not a rename.
        (
            "shellCF2",
            "shell_particle",
            "Additionally, the parameterization has changed. 'shellCF2' took "
            "the central radius 'a' and the shell thickness 'delta', while "
            "'shell_particle' takes the inner radius. Please pass "
            "radius = a - delta / 2 and thickness = delta.",
            (numpy.array([0.0, 5.0, 12.5]), 12.5, 5.0),
            (numpy.array([0.0, 5.0, 12.5]), 10.0, 5.0),
        ),
    ],
)
def test_deprecated_functions_warn_and_forward(
    old_name,
    new_name,
    expected_signature_note,
    old_args,
    new_args,
):
    old_function = getattr(cf, old_name)
    new_function = getattr(cf, new_name)
    module_path = "diffpy.srfit.pdf.characteristicfunctions"
    expected_characteristic_function = new_function(*new_args)
    expected_msg = (
        f"'{module_path}.{old_name}' is deprecated and will be removed "
        f"in version 4.0.0. Please use '{module_path}.{new_name}' "
        "instead."
    )
    if expected_signature_note:
        expected_msg = f"{expected_msg} {expected_signature_note}"
    # Capture every warning rather than using pytest.warns, which only
    # checks that a matching warning occurred somewhere and would not
    # notice a second, irrelevant DeprecationWarning from a call chained
    # through another deprecated function.
    with warnings.catch_warnings(record=True) as raised_warnings:
        warnings.simplefilter("always")
        actual_characteristic_function = old_function(*old_args)
    deprecation_warnings = [
        raised
        for raised in raised_warnings
        if raised.category is DeprecationWarning
    ]
    actual_deprecation_count = len(deprecation_warnings)
    expected_deprecation_count = 1
    assert actual_deprecation_count == expected_deprecation_count
    actual_msg = str(deprecation_warnings[0].message)
    assert actual_msg == expected_msg
    npt.assert_allclose(
        actual_characteristic_function, expected_characteristic_function
    )


# Because the deprecated functions forward to their replacements, they must
# also forward the replacement's non-physical-input handling: a
# RuntimeWarning alongside the DeprecationWarning, and a result of zero (or,
# for a translated case, zero reached via the translated argument).


@pytest.mark.parametrize(
    "old_name, new_name, old_args, new_args, expected_requirement",
    [
        # C1: sphericalCF forwards a non-positive diameter unchanged.
        # Expected: the RuntimeWarning names the forwarded parameter.
        (
            "sphericalCF",
            "spherical_particle",
            (numpy.array([1.0, 5.0]), -10.0),
            (numpy.array([1.0, 5.0]), -10.0),
            "'particle_diameter' must be positive",
        ),
        # C2: spheroidalCF2 translates a non-positive equatorial diameter
        # into a non-positive equatorial radius.
        # Expected: the RuntimeWarning names the translated parameter.
        (
            "spheroidalCF2",
            "spheroidal_particle",
            (numpy.array([1.0, 5.0]), -20.0, 1.5),
            (numpy.array([1.0, 5.0]), -10.0, -15.0),
            "'equatorial_radius' must be positive",
        ),
        # C3: lognormalSphericalCF forwards a negative distribution width
        # unchanged. This is the case where the old, pre-3.3.0 code and the
        # new function most disagreed: the old code treated a negative width
        # as the zero-width (sphere) limit, where the new function rejects
        # it outright.
        # Expected: the RuntimeWarning that the old code never raised.
        (
            "lognormalSphericalCF",
            "lognormal_spherical_particle",
            (numpy.array([1.0, 5.0, 10.0]), 10.0, -2.0),
            (numpy.array([1.0, 5.0, 10.0]), 10.0, -2.0),
            "'particle_diameter_sigma' must not be negative",
        ),
        # C4: sheetCF forwards a non-positive thickness unchanged. The old
        # code returned a bare scalar 0 regardless of the shape of r, where
        # the new function returns an array shaped like r.
        # Expected: an array of zeros shaped like r.
        (
            "sheetCF",
            "sheet_particle",
            (numpy.array([1.0, 5.0]), 0.0),
            (numpy.array([1.0, 5.0]), 0.0),
            "'sheet_thickness' must be positive",
        ),
        # C5: shellCF forwards a zero thickness unchanged. The old code hit
        # a vanishing normalization denominator and substituted one, where
        # the new function rejects the input outright.
        # Expected: zero everywhere rather than the old fallback of one.
        (
            "shellCF",
            "shell_particle",
            (numpy.array([1.0, 5.0]), 10.0, 0.0),
            (numpy.array([1.0, 5.0]), 10.0, 0.0),
            "'thickness' must be positive",
        ),
        # C6: shellCF2 translates a central radius smaller than half the
        # thickness into a negative inner radius.
        # Expected: the RuntimeWarning names the translated parameter.
        (
            "shellCF2",
            "shell_particle",
            (numpy.array([1.0, 5.0]), 1.0, 5.0),
            (numpy.array([1.0, 5.0]), -1.5, 5.0),
            "'radius' must not be negative",
        ),
    ],
)
def test_deprecated_functions_forward_non_physical_input(
    old_name,
    new_name,
    old_args,
    new_args,
    expected_requirement,
    reset_characteristic_function_warnings,
):
    old_function = getattr(cf, old_name)
    new_function = getattr(cf, new_name)
    expected_characteristic_function = new_function(*new_args)

    # The oracle call above already consumed the once-per-process
    # RuntimeWarning for this (function, requirement) pair; clear the
    # record so the call under test below is free to raise it again.
    cf._warned_non_physical.clear()

    with warnings.catch_warnings(record=True) as raised_warnings:
        warnings.simplefilter("always")
        actual_characteristic_function = old_function(*old_args)

    # Exactly two warnings are expected: the DeprecationWarning for the old
    # name, and the RuntimeWarning forwarded from the new function. A stray
    # third warning, such as a second DeprecationWarning from a call
    # chained through another deprecated function, must fail this check.
    actual_warning_count = len(raised_warnings)
    expected_warning_count = 2
    assert actual_warning_count == expected_warning_count

    actual_categories = {raised.category for raised in raised_warnings}
    expected_categories = {DeprecationWarning, RuntimeWarning}
    assert actual_categories == expected_categories

    # The RuntimeWarning always names the new function, since that is where
    # it is actually raised, regardless of which name the caller used.
    expected_remedy = (
        "The characteristic function was set to zero for every r, which "
        "flattens the fit residual so a refinement cannot recover on its "
        "own. Please use a physically meaningful starting value, or keep "
        "the parameter in range with FitRecipe.add_soft_bounds."
    )
    expected_runtime_msg = (
        f"In '{new_name}', {expected_requirement}. {expected_remedy}"
    )
    runtime_warning = next(
        raised
        for raised in raised_warnings
        if raised.category is RuntimeWarning
    )
    actual_runtime_msg = str(runtime_warning.message)
    assert actual_runtime_msg == expected_runtime_msg

    npt.assert_allclose(
        actual_characteristic_function, expected_characteristic_function
    )


if __name__ == "__main__":
    unittest.main()
