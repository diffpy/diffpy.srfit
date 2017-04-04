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

"""Tests for fitresults module."""

import unittest

from diffpy.srfit.fitbase.fitrecipe import FitRecipe
from diffpy.srfit.fitbase.fitresults import initializeRecipe
from diffpy.srfit.tests.utils import datafile


class TestInitializeRecipe(unittest.TestCase):

    def setUp(self):
        self.recipe = recipe = FitRecipe("recipe")
        recipe.newVar("A", 0)
        recipe.newVar("sig", 0)
        recipe.newVar("x0", 0)
        self.filename = datafile("results.res")

        self.Aval = 5.77619823e-01
        self.sigval = -9.22758690e-01
        self.x0val = 6.12422115e+00
        return

    def testInitializeFromFileName(self):
        recipe = self.recipe
        self.assertEqual(0, recipe.A.value)
        self.assertEqual(0, recipe.sig.value)
        self.assertEqual(0, recipe.x0.value)
        initializeRecipe(recipe, self.filename)
        self.assertAlmostEqual(self.Aval, recipe.A.value)
        self.assertAlmostEqual(self.sigval, recipe.sig.value)
        self.assertAlmostEqual(self.x0val, recipe.x0.value)
        return

    def testInitializeFromFileObj(self):
        recipe = self.recipe
        self.assertEqual(0, recipe.A.value)
        self.assertEqual(0, recipe.sig.value)
        self.assertEqual(0, recipe.x0.value)
        infile = open(self.filename, 'r')
        initializeRecipe(recipe, infile)
        self.assertFalse(infile.closed)
        infile.close()
        self.assertAlmostEqual(self.Aval, recipe.A.value)
        self.assertAlmostEqual(self.sigval, recipe.sig.value)
        self.assertAlmostEqual(self.x0val, recipe.x0.value)
        return


    def testInitializeFromString(self):
        recipe = self.recipe
        self.assertEqual(0, recipe.A.value)
        self.assertEqual(0, recipe.sig.value)
        self.assertEqual(0, recipe.x0.value)
        infile = open(self.filename, 'r')
        resstr = infile.read()
        infile.close()
        initializeRecipe(recipe, resstr)
        self.assertAlmostEqual(self.Aval, recipe.A.value)
        self.assertAlmostEqual(self.sigval, recipe.sig.value)
        self.assertAlmostEqual(self.x0val, recipe.x0.value)
        return

if __name__ == "__main__":

    unittest.main()
