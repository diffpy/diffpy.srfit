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

    def testCylinder(self):
        """Make sure cylinder works over different r-ranges"""
        radius = 100
        length = 30

        from sans.models.CylinderModel import CylinderModel
        model = CylinderModel()
        model.setParam("radius", radius)
        model.setParam("length", length)

        ff = nff.SASFormFactor("cylinder", model)

        r1 = numpy.arange(0, 10, 0.1, dtype = float)
        r2 = numpy.arange(0, 50, 0.1, dtype = float)
        r3 = numpy.arange(0, 100, 0.1, dtype = float)
        r4 = numpy.arange(0, 500, 0.1, dtype = float)

        fr1 = ff(r1)
        fr2 = ff(r2)
        fr3 = ff(r3)
        fr4 = ff(r4)

        d = fr1 - numpy.interp(r1, r2, fr2)
        res12 = numpy.dot(d,d)
        res12 /= numpy.dot(fr1, fr1)
        self.assertAlmostEqual(0, res12, 4)

        d = fr1 - numpy.interp(r1, r3, fr3)
        res13 = numpy.dot(d,d)
        res13 /= numpy.dot(fr1, fr1)
        self.assertAlmostEqual(0, res13, 4)

        d = fr1 - numpy.interp(r1, r4, fr4)
        res14 = numpy.dot(d,d)
        res14 /= numpy.dot(fr1, fr1)
        self.assertAlmostEqual(0, res14, 4)

        d = fr2 - numpy.interp(r2, r3, fr3)
        res23 = numpy.dot(d,d)
        res23 /= numpy.dot(fr2, fr2)
        self.assertAlmostEqual(0, res23, 4)

        d = fr2 - numpy.interp(r2, r4, fr4)
        res24 = numpy.dot(d,d)
        res24 /= numpy.dot(fr2, fr2)
        self.assertAlmostEqual(0, res24, 4)

        d = fr3 - numpy.interp(r3, r4, fr4)
        res34 = numpy.dot(d,d)
        res34 /= numpy.dot(fr3, fr3)
        self.assertAlmostEqual(0, res34, 4)
        return



if __name__ == "__main__":
    unittest.main()
