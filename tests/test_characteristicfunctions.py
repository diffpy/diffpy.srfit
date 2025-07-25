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
from diffpy.srfit.sas.sasimport import sasimport

# # Global variables to be assigned in setUp
# cf = None
# Fixme: remove this code if imports don't break on testing

# ----------------------------------------------------------------------------


def testSphere(sas_available):
    if not sas_available:
        pytest.skip("sas package not available")
    radius = 25
    # Calculate sphere cf from SphereModel
    SphereModel = sasimport("sas.models.SphereModel").SphereModel
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
    if not sas_available:
        pytest.skip("sas package not available")
    prad = 20.9
    erad = 33.114
    # Calculate cf from EllipsoidModel
    EllipsoidModel = sasimport("sas.models.EllipsoidModel").EllipsoidModel
    model = EllipsoidModel()
    model.setParam("radius_a", prad)
    model.setParam("radius_b", erad)
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
    if not sas_available:
        pytest.skip("sas package not available")
    radius = 19.2
    thickness = 7.8
    # Calculate cf from VesicleModel
    VesicleModel = sasimport("sas.models.VesicleModel").VesicleModel
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
    if not sas_available:
        pytest.skip("sas package not available")
    """Make sure cylinder works over different r-ranges."""
    radius = 100
    length = 30

    CylinderModel = sasimport("sas.models.CylinderModel").CylinderModel
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
