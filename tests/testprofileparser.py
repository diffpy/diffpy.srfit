#!/usr/bin/env python
"""Tests for refinableobj module."""

import unittest
import os.path

from diffpy.srfit.fit.profileparser import TextParser, ParseError
from diffpy.srfit.fit.profileparser import getParser

thisfile = locals().get('__file__', 'testprofileparser.py')
tests_dir = os.path.dirname(os.path.abspath(thisfile))
testdata_dir = os.path.join(tests_dir, 'testdata')

class TestTextParser(unittest.TestCase):

    def testGetParser(self):
        """Test getParser."""
        parser = getParser("txt")
        self.assertTrue(parser is TextParser)
        self.assertRaises(ValueError, getParser, "alksdfkj")
        return

    def testParse(self):
        """Test the parser."""
        parser = TextParser()
        filename = os.path.join(testdata_dir, "testdata.txt")

        val0 = 1e-2
        val1 = 1.105784e-1
        val2 = -7.539822e-4
        val3 = 1.802192e-3

        # Test normal load. It assumes that dy is the third array and dx is the
        # fourth.
        parser.parseFile(filename)
        x, y, dx, dy = parser.getData()
        self.assertAlmostEqual(val0, x[0])
        self.assertAlmostEqual(val1, y[0])
        self.assertAlmostEqual(val3, dx[0])
        self.assertAlmostEqual(val2, dy[0])

        # Test normal load. Try to pass unpack = False to verify that it is
        # ignored.
        parser.parseFile(filename, unpack = False)
        x, y, dx, dy = parser.getData()
        self.assertAlmostEqual(val0, x[0])
        self.assertAlmostEqual(val1, y[0])
        self.assertAlmostEqual(val3, dx[0])
        self.assertAlmostEqual(val2, dy[0])

        # Load the arrays in the desired order
        parser.parseFile(filename, usecols=(0,1,3,2))
        x, y, dx, dy = parser.getData()
        self.assertAlmostEqual(val0, x[0])
        self.assertAlmostEqual(val1, y[0])
        self.assertAlmostEqual(val2, dx[0])
        self.assertAlmostEqual(val3, dy[0])

        # Try not including dy
        parser.parseFile(filename, usecols=(0,1))
        x, y, dx, dy = parser.getData()
        self.assertAlmostEqual(val0, x[0])
        self.assertAlmostEqual(val1, y[0])
        self.assertTrue(dy is None)
        self.assertTrue(dx is None)

        # Try to include too little
        self.assertRaises(ParseError, parser.parseFile, filename, usecols=(0,))


        return

if __name__ == "__main__":

    unittest.main()

