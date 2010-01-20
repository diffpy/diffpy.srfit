#!/usr/bin/env python
"""Tests for sas package."""

import unittest

import numpy

from diffpy.srfit.sas import SASGenerator, SASParser, SASProfile
from diffpy.srfit.fitbase import Profile

class TestSASParser(unittest.TestCase):

    def testParser(self):
        data = "testdata/sas_ascii_test_1.txt"
        parser = SASParser()
        parser.parseFile(data)

        meta = parser._meta

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


class TestSASGenerator(unittest.TestCase):

    def testGenerator(self):

        # Test generator output
        from sans.models.SphereModel import SphereModel
        model = SphereModel()
        gen = SASGenerator("sphere", model)

        for pname in model.params:
            defval = model.getParam(pname)
            par = gen.get(pname)
            self.assertEquals(defval, par.getValue())
            # Test setting values
            par.setValue(1.0)
            self.assertEquals(1.0, par.getValue())
            self.assertEquals(1.0, model.getParam(pname))
            par.setValue(defval)
            self.assertEquals(defval, par.getValue())
            self.assertEquals(defval, model.getParam(pname))


        r = numpy.arange(1, 10, 0.1, dtype = float)
        y = gen(r)
        refy = model.evalDistribution(r)
        diff = y - refy
        res = numpy.dot(diff, diff)
        self.assertAlmostEqual(0, res)

        return

    def testGenerator2(self):

        # Test generator with a profile
        from sans.models.EllipsoidModel import EllipsoidModel
        model = EllipsoidModel()
        gen = SASGenerator("ellipsoid", model)

        # Load the data using SAS tools
        from DataLoader.loader import Loader
        loader = Loader()
        datainfo = loader.load("testdata/sas_ellipsoid_testdata.txt")
        profile = SASProfile(datainfo)

        gen.setProfile(profile)
        gen.scale = 1.0
        gen.radius_a = 20
        gen.radius_b = 400
        gen.contrast = 3e-6
        gen.background = 0.01

        y = gen(profile.xobs)
        diff = profile.yobs - y
        res = numpy.dot(diff, diff)
        self.assertAlmostEqual(0, res)
        return


if __name__ == "__main__":
    unittest.main()
