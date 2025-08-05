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
"""Tests for sas package."""

import unittest

import numpy
import pytest

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

    from sasmodels.sasview_model import find_model, load_standard_models

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
    fr2 = cf.sphericalCF(r, 2 * radius)
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

    from sasmodels.sasview_model import find_model, load_standard_models

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
    fr2 = cf.spheroidalCF(r, erad, prad)
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

    from sasmodels.sasview_model import find_model, load_standard_models

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
    fr2 = cf.shellCF(r, radius, thickness)
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
    from sasmodels.sasview_model import find_model, load_standard_models

    load_standard_models()
    """Make sure cylinder works over different r-ranges."""
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

if __name__ == "__main__":
    unittest.main()
