#!/usr/bin/env python
##############################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2010 Trustees of the Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
##############################################################################
"""Tests for sas package."""

import unittest

import numpy

from utils import testcase, TestCaseSaS, TestCasePdf

# Global variables to be assigned in setUp
cf = None


class TestSASCF(testcase(TestCaseSaS,TestCasePdf)):

    def setUp(self):
        global cf
        import diffpy.srfit.pdf.characteristicfunctions as cf

    def testSphere(self):
        radius = 25
        # Calculate sphere cf from SphereModel
        from sans.models.SphereModel import SphereModel
        model = SphereModel()
        model.setParam("radius", radius)
        ff = cf.SASCF("sphere", model)
        r = numpy.arange(1, 60, 0.1, dtype = float)
        fr1 = ff(r)

        # Calculate sphere cf analytically
        fr2 = cf.sphericalCF(r, 2*radius)
        diff = fr1 - fr2
        res = numpy.dot(diff, diff)
        res /= numpy.dot(fr2, fr2)
        self.assertAlmostEqual(0, res, 4)
        return

    def testSpheroid(self):
        prad = 20.9
        erad = 33.114
        # Calculate cf from EllipsoidModel
        from sans.models.EllipsoidModel import EllipsoidModel
        model = EllipsoidModel()
        model.setParam("radius_a", prad)
        model.setParam("radius_b", erad)
        ff = cf.SASCF("spheroid", model)
        r = numpy.arange(0, 100, 1/numpy.pi, dtype = float)
        fr1 = ff(r)

        # Calculate cf analytically
        fr2 = cf.spheroidalCF(r, erad, prad)
        diff = fr1 - fr2
        res = numpy.dot(diff, diff)
        res /= numpy.dot(fr2, fr2)
        self.assertAlmostEqual(0, res, 4)
        return

    def testShell(self):
        radius = 19.2
        thickness = 7.8
        # Calculate cf from VesicleModel
        from sans.models.VesicleModel import VesicleModel
        model = VesicleModel()
        model.setParam("radius", radius)
        model.setParam("thickness", thickness)
        ff = cf.SASCF("vesicle", model)
        r = numpy.arange(0, 99.45, 0.1, dtype = float)
        fr1 = ff(r)

        # Calculate sphere cf analytically
        fr2 = cf.shellCF(r, radius, thickness)
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

        ff = cf.SASCF("cylinder", model)

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
