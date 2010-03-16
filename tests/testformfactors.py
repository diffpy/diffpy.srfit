#!/usr/bin/env python
"""Tests for sas package."""

import unittest
import os.path

import numpy

import diffpy.srfit.pdf.nanoformfactors as nff

class TestSASFormFactor(unittest.TestCase):

    def testSphere(self):
        radius = 25
        # Calculate sphere ff from SphereModel
        from sans.models.SphereModel import SphereModel
        model = SphereModel()
        model.setParam("radius", radius)
        ff = nff.SASFormFactor("sphere", model)
        r = numpy.arange(1, 60, 0.1, dtype = float)
        fr1 = ff(r)

        # Calculate sphere ff analytically
        fr2 = nff.sphericalFF(r, 2*radius)
        diff = fr1 - fr2
        res = numpy.dot(diff, diff)
        res /= numpy.dot(fr2, fr2)
        self.assertAlmostEqual(0, res, 4)
        return

    def testSpheroid(self):
        prad = 20.9
        erad = 33.114
        # Calculate ff from EllipsoidModel
        from sans.models.EllipsoidModel import EllipsoidModel
        model = EllipsoidModel()
        model.setParam("radius_a", prad)
        model.setParam("radius_b", erad)
        ff = nff.SASFormFactor("spheroid", model)
        r = numpy.arange(0, 100, 1/numpy.pi, dtype = float)
        fr1 = ff(r)

        # Calculate ff analytically
        fr2 = nff.spheroidalFF(r, erad, prad)
        diff = fr1 - fr2
        res = numpy.dot(diff, diff)
        res /= numpy.dot(fr2, fr2)
        self.assertAlmostEqual(0, res, 4)
        return

    def testShell(self):
        radius = 19.2
        thickness = 7.8
        # Calculate ff from VesicleModel
        from sans.models.VesicleModel import VesicleModel
        model = VesicleModel()
        model.setParam("radius", radius)
        model.setParam("thickness", thickness)
        ff = nff.SASFormFactor("vesicle", model)
        r = numpy.arange(0, 99.45, 0.1, dtype = float)
        fr1 = ff(r)

        # Calculate sphere ff analytically
        fr2 = nff.shellFF(r, radius, thickness)
        diff = fr1 - fr2
        res = numpy.dot(diff, diff)
        res /= numpy.dot(fr2, fr2)
        self.assertAlmostEqual(0, res, 4)
        return

if __name__ == "__main__":
    unittest.main()
