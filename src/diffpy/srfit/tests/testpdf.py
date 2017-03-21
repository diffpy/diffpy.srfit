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

"""Tests for pdf package."""

import unittest

import numpy

from diffpy.srfit.tests.utils import datafile
from diffpy.srfit.tests.utils import testoptional
from diffpy.srfit.tests.utils import TestCasePDF, TestCaseStructure

# Global variables to be assigned in setUp
PDFGenerator = PDFParser = None


class TestPDFParset(testoptional(TestCasePDF)):

    def setUp(self):
        global PDFParser
        from diffpy.srfit.pdf import PDFParser

    def testParser1(self):
        data = datafile("ni-q27r100-neutron.gr")
        parser = PDFParser()
        parser.parseFile(data)

        meta = parser._meta

        self.assertEqual(data, meta['filename'])
        self.assertEqual(1, meta['nbanks'])
        self.assertEqual('N', meta['stype'])
        self.assertEqual(27, meta['qmax'])
        self.assertEqual(300, meta.get('temperature'))
        self.assertEqual(None, meta.get('qdamp'))
        self.assertEqual(None, meta.get('qbroad'))
        self.assertEqual(None, meta.get('spdiameter'))
        self.assertEqual(None, meta.get('scale'))
        self.assertEqual(None, meta.get('doping'))

        x, y, dx, dy = parser.getData()
        self.assertTrue(dx is None)
        self.assertTrue(dy is None)

        testx = numpy.linspace(0.01, 100, 10000)
        diff = testx - x
        res = numpy.dot(diff, diff)
        self.assertAlmostEqual(0, res)

        testy = numpy.array([1.144, 2.258, 3.312, 4.279, 5.135, 5.862, 6.445,
            6.875, 7.150, 7.272])
        diff = testy - y[:10]
        res = numpy.dot(diff, diff)
        self.assertAlmostEqual(0, res)

        return

    def testParser2(self):
        data = datafile("si-q27r60-xray.gr")
        parser = PDFParser()
        parser.parseFile(data)

        meta = parser._meta

        self.assertEqual(data, meta['filename'])
        self.assertEqual(1, meta['nbanks'])
        self.assertEqual('X', meta['stype'])
        self.assertEqual(27, meta['qmax'])
        self.assertEqual(300, meta.get('temperature'))
        self.assertEqual(None, meta.get('qdamp'))
        self.assertEqual(None, meta.get('qbroad'))
        self.assertEqual(None, meta.get('spdiameter'))
        self.assertEqual(None, meta.get('scale'))
        self.assertEqual(None, meta.get('doping'))

        x, y, dx, dy = parser.getData()
        testx = numpy.linspace(0.01, 60, 5999, endpoint=False)
        diff = testx - x
        res = numpy.dot(diff, diff)
        self.assertAlmostEqual(0, res)

        testy = numpy.array([0.1105784, 0.2199684, 0.3270088, 0.4305913,
            0.5296853, 0.6233606, 0.7108060, 0.7913456, 0.8644501, 0.9297440])
        diff = testy - y[:10]
        res = numpy.dot(diff, diff)
        self.assertAlmostEqual(0, res)

        testdy = numpy.array([0.001802192, 0.003521449, 0.005079115,
            0.006404892, 0.007440527, 0.008142955, 0.008486813, 0.008466340,
            0.008096858, 0.007416456])
        diff = testdy - dy[:10]
        res = numpy.dot(diff, diff)
        self.assertAlmostEqual(0, res)

        self.assertTrue(dx is None)
        return

class TestPDFGenerator(testoptional(TestCaseStructure, TestCasePDF)):

    def setUp(self):
        global PDFGenerator
        from diffpy.srfit.pdf import PDFGenerator
        self.gen = PDFGenerator()
        return


    def testGenerator(self):
        qmax = 27.0
        gen = self.gen
        gen.setScatteringType('N')
        self.assertEqual('N', gen.getScatteringType())
        gen.setQmax(qmax)
        self.assertAlmostEqual(qmax, gen.getQmax())
        from diffpy.Structure import PDFFitStructure
        stru = PDFFitStructure()
        ciffile = datafile("ni.cif")
        stru.read(ciffile)
        for i in range(4):
            stru[i].Bisoequiv = 1
        gen.setStructure(stru)

        calc = gen._calc
        # Test parameters
        for par in gen.iterPars(recurse=False):
            pname = par.name
            defval = calc._getDoubleAttr(pname)
            self.assertEqual(defval, par.getValue())
            # Test setting values
            par.setValue(1.0)
            self.assertEqual(1.0, par.getValue())
            par.setValue(defval)
            self.assertEqual(defval, par.getValue())

        r = numpy.arange(0, 10, 0.1)
        y = gen(r)

        # Now create a reference PDF. Since the calculator is testing its
        # output, we just have to make sure we can calculate from the
        # PDFGenerator interface.
        from diffpy.srreal.pdfcalculator import PDFCalculator
        calc = PDFCalculator()
        calc.rstep = r[1] - r[0]
        calc.rmin = r[0]
        calc.rmax = r[-1] + 0.5 * calc.rstep
        calc.qmax = qmax
        calc.setScatteringFactorTableByType('N')
        calc.eval(stru)
        yref = calc.pdf

        diff = y - yref
        res = numpy.dot(diff, diff)
        self.assertAlmostEqual(0, res)
        return


    def test_setQmin(self):
        """Verify qmin is propagated to the calculator object.
        """
        gen = self.gen
        self.assertEqual(0, gen.getQmin())
        self.assertEqual(0, gen._calc.qmin)
        gen.setQmin(0.93)
        self.assertEqual(0.93, gen.getQmin())
        self.assertEqual(0.93, gen._calc.qmin)
        return


if __name__ == "__main__":
    unittest.main()
