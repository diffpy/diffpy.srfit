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

from diffpy.srfit.sas import SASGenerator, SASParser, SASProfile
from diffpy.srfit.tests.utils import TestCaseSaS, datafile
from diffpy.srfit.sas.sasimport import sasimport


class TestSASParser(TestCaseSaS):

    def testParser(self):
        data = datafile("sas_ascii_test_1.txt")
        parser = SASParser()
        parser.parseFile(data)

        x, y, dx, dy = parser.getData()

        testx = numpy.array([0.002618, 0.007854, 0.01309, 0.01832, 0.02356,
            0.02879, 0.03402, 0.03925, 0.04448, 0.0497])
        diff = testx - x
        res = numpy.dot(diff, diff)
        self.assertAlmostEqual(0, res)

        testy = numpy.array([ 0.02198, 0.02201, 0.02695, 0.02645, 0.03024,
            0.3927, 7.305, 17.43, 13.43, 8.346])
        diff = testy - y
        res = numpy.dot(diff, diff)
        self.assertAlmostEqual(0, res)

        testdy = numpy.array([ 0.002704, 0.001643, 0.002452, 0.001769,
            0.001531, 0.1697, 1.006, 0.5351, 0.3677, 0.191])
        diff = testdy - dy
        res = numpy.dot(diff, diff)
        self.assertAlmostEqual(0, res)

        testdx = numpy.array([0.0004091, 0.005587, 0.005598, 0.005624,
            0.005707, 0.005975, 0.006264, 0.006344, 0.006424, 0.006516])
        diff = testdx - dx
        res = numpy.dot(diff, diff)
        self.assertAlmostEqual(0, res)

        return


class TestSASGenerator(TestCaseSaS):

    def testGenerator(self):

        # Test generator output
        SphereModel = sasimport('sas.models.SphereModel').SphereModel
        model = SphereModel()
        gen = SASGenerator("sphere", model)

        for pname in model.params:
            defval = model.getParam(pname)
            par = gen.get(pname)
            self.assertEqual(defval, par.getValue())
            # Test setting values
            par.setValue(1.0)
            self.assertEqual(1.0, par.getValue())
            self.assertEqual(1.0, model.getParam(pname))
            par.setValue(defval)
            self.assertEqual(defval, par.getValue())
            self.assertEqual(defval, model.getParam(pname))


        r = numpy.arange(1, 10, 0.1, dtype = float)
        y = gen(r)
        refy = model.evalDistribution(r)
        diff = y - refy
        res = numpy.dot(diff, diff)
        self.assertAlmostEqual(0, res)

        return

    def testGenerator2(self):

        # Test generator with a profile
        EllipsoidModel = sasimport('sas.models.EllipsoidModel').EllipsoidModel
        model = EllipsoidModel()
        gen = SASGenerator("ellipsoid", model)

        # Load the data using SAS tools
        Loader = sasimport('sas.dataloader.loader').Loader
        loader = Loader()
        data = datafile("sas_ellipsoid_testdata.txt")
        datainfo = loader.load(data)
        profile = SASProfile(datainfo)

        gen.setProfile(profile)
        gen.scale.value = 1.0
        gen.radius_a.value = 20
        gen.radius_b.value = 400
        gen.background.value = 0.01

        y = gen(profile.xobs)
        diff = profile.yobs - y
        res = numpy.dot(diff, diff)
        self.assertAlmostEqual(0, res)
        return


if __name__ == "__main__":
    unittest.main()
