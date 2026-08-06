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


import re
import unittest

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
        # C2: the particle diameter is non-positive.
        # Expected: the characteristic function is zero everywhere.
        (numpy.array([0.0, 5.0]), 0.0, [0.0, 0.0]),
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
        # C4: the polar radius is non-positive.
        # Expected: the characteristic function is zero everywhere.
        (numpy.array([1.0, 2.0]), 10.0, 0.0, [0.0, 0.0]),
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
        # C3: the mean particle diameter is non-positive.
        # Expected: the characteristic function is zero everywhere.
        (numpy.array([1.0, 2.0]), 0.0, 1.0, [0.0, 0.0]),
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
    "input_sheet_thickness",
    [
        # C1: zero thickness.
        # Expected: the characteristic function is zero, regardless of r.
        0.0,
        # C2: negative thickness.
        # Expected: the characteristic function is zero, regardless of r.
        -1.0,
    ],
)
def test_sheet_particle_non_positive_thickness(input_sheet_thickness):
    actual_characteristic_function = cf.sheet_particle(
        numpy.array([1.0, 2.0]), input_sheet_thickness
    )
    expected_characteristic_function = 0
    assert actual_characteristic_function == expected_characteristic_function


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
# sphericalCF, spheroidalCF, spheroidalCF2, lognormalSphericalCF, sheetCF,
# shellCF, and shellCF2 are deprecated in favor of the snake_case functions
# above. Each old name must still work, emit a DeprecationWarning naming its
# replacement, and forward to the new implementation (translating arguments
# where the old and new parameterizations differ).

module_path = "diffpy.srfit.pdf.characteristicfunctions"


@pytest.mark.parametrize(
    "old_name, new_name, old_args, expected_characteristic_function",
    [
        # C1: sphericalCF forwards directly to spherical_particle.
        (
            "sphericalCF",
            "spherical_particle",
            (numpy.array([0.0, 5.0, 10.0]), 10.0),
            cf.spherical_particle(numpy.array([0.0, 5.0, 10.0]), 10.0),
        ),
        # C2: spheroidalCF forwards directly to spheroidal_particle.
        (
            "spheroidalCF",
            "spheroidal_particle",
            (numpy.array([0.0, 5.0, 10.0]), 10.0, 15.0),
            cf.spheroidal_particle(numpy.array([0.0, 5.0, 10.0]), 10.0, 15.0),
        ),
        # C3: spheroidalCF2 used the raw (diameter, axis ratio)
        # parameterization, which must be converted to (equatorial_radius,
        # polar_radius) before forwarding to spheroidal_particle.
        (
            "spheroidalCF2",
            "spheroidal_particle",
            (numpy.array([0.0, 5.0, 10.0]), 20.0, 1.5),
            cf.spheroidal_particle(numpy.array([0.0, 5.0, 10.0]), 10.0, 15.0),
        ),
        # C4: lognormalSphericalCF forwards directly to
        # lognormal_spherical_particle.
        (
            "lognormalSphericalCF",
            "lognormal_spherical_particle",
            (numpy.array([1.0, 5.0, 10.0]), 10.0, 2.0),
            cf.lognormal_spherical_particle(
                numpy.array([1.0, 5.0, 10.0]), 10.0, 2.0
            ),
        ),
        # C5: sheetCF forwards directly to sheet_particle.
        (
            "sheetCF",
            "sheet_particle",
            (numpy.array([0.0, 2.0, 4.0]), 4.0),
            cf.sheet_particle(numpy.array([0.0, 2.0, 4.0]), 4.0),
        ),
        # C6: shellCF forwards directly to shell_particle.
        (
            "shellCF",
            "shell_particle",
            (numpy.array([0.0, 5.0, 12.5]), 10.0, 5.0),
            cf.shell_particle(numpy.array([0.0, 5.0, 12.5]), 10.0, 5.0),
        ),
        # C7: shellCF2 used the raw (central_radius, thickness)
        # parameterization, which must be converted to (radius, thickness)
        # before forwarding to shell_particle.
        (
            "shellCF2",
            "shell_particle",
            (numpy.array([0.0, 5.0, 12.5]), 12.5, 5.0),
            cf.shell_particle(numpy.array([0.0, 5.0, 12.5]), 10.0, 5.0),
        ),
    ],
)
def test_deprecated_functions_warn_and_forward(
    old_name, new_name, old_args, expected_characteristic_function
):
    old_function = getattr(cf, old_name)
    expected_msg = (
        f"'{module_path}.{old_name}' is deprecated and will be removed "
        f"in version 4.0.0. Please use '{module_path}.{new_name}' "
        "instead."
    )
    with pytest.warns(DeprecationWarning, match=re.escape(expected_msg)):
        actual_characteristic_function = old_function(*old_args)
    npt.assert_allclose(
        actual_characteristic_function, expected_characteristic_function
    )


if __name__ == "__main__":
    unittest.main()
