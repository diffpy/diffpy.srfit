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

import numpy
import pytest

from diffpy.srfit.sas import SASGenerator, SASParser, SASProfile

# ----------------------------------------------------------------------------
# FIXME: adjust sensitivity of the pytest.approx statements when ready to test
# with sasview installed.


def testParser(sas_available, datafile):
    if not sas_available:
        pytest.skip("sas package not available")

    data = datafile("sas_ascii_test_1.txt")
    parser = SASParser()
    parser.parseFile(data)
    x, y, dx, dy = parser.getData()
    testx = numpy.array(
        [
            0.002618,
            0.007854,
            0.01309,
            0.01832,
            0.02356,
            0.02879,
            0.03402,
            0.03925,
            0.04448,
            0.0497,
        ]
    )
    diff = testx - x
    res = numpy.dot(diff, diff)
    assert 0 == pytest.approx(res)

    testy = numpy.array(
        [
            0.02198,
            0.02201,
            0.02695,
            0.02645,
            0.03024,
            0.3927,
            7.305,
            17.43,
            13.43,
            8.346,
        ]
    )
    diff = testy - y
    res = numpy.dot(diff, diff)
    assert 0 == pytest.approx(res)

    testdy = numpy.array(
        [
            0.002704,
            0.001643,
            0.002452,
            0.001769,
            0.001531,
            0.1697,
            1.006,
            0.5351,
            0.3677,
            0.191,
        ]
    )
    diff = testdy - dy
    res = numpy.dot(diff, diff)
    assert 0 == pytest.approx(res)

    testdx = numpy.array(
        [
            0.0004091,
            0.005587,
            0.005598,
            0.005624,
            0.005707,
            0.005975,
            0.006264,
            0.006344,
            0.006424,
            0.006516,
        ]
    )
    diff = testdx - dx
    res = numpy.dot(diff, diff)
    assert 0 == pytest.approx(res)
    return


# End of class TestSASParser


def test_generator(sas_available):
    if not sas_available:
        pytest.skip("sas package not available")
    from sasmodels.sasview_model import find_model, load_standard_models

    load_standard_models()
    SphereModel = find_model("sphere")
    model = SphereModel()
    gen = SASGenerator("sphere", model)
    for pname in model.params:
        defval = model.getParam(pname)
        par = gen.get(pname)
        assert defval == par.getValue()
        # Test setting values
        par.setValue(1.0)
        assert 1.0 == par.getValue()
        assert 1.0 == model.getParam(pname)
        par.setValue(defval)
        assert defval == par.getValue()
        assert defval == model.getParam(pname)

    r = numpy.arange(1, 10, 0.1, dtype=float)
    y = gen(r)
    refy = model.evalDistribution(r)
    diff = y - refy
    res = numpy.dot(diff, diff)
    assert 0 == pytest.approx(res)
    return


def testGenerator2(sas_available, datafile):
    if not sas_available:
        pytest.skip("sas package not available")
    from sasmodels.sasview_model import find_model, load_standard_models

    load_standard_models()
    EllipsoidModel = find_model("ellipsoid")
    model = EllipsoidModel()
    gen = SASGenerator("ellipsoid", model)

    # Load the data using SAS tools
    import sasdata.dataloader.loader as sas_dataloader

    Loader = sas_dataloader.Loader
    loader = Loader()
    data = datafile("sas_ellipsoid_testdata.txt")
    datainfo = loader.load(str(data))
    profile = SASProfile(datainfo)

    gen.setProfile(profile)
    gen.scale.value = 1.0
    gen.radius_polar.value = 20
    gen.radius_equatorial.value = 400
    gen.background.value = 0.01

    y = gen(profile.xobs)
    diff = profile.yobs - y
    res = numpy.dot(diff, diff)
    # FIXME: go back to default tolerance when we figure out why
    # the models are not identical
    assert 0 == pytest.approx(res, abs=1e-3)
    return
