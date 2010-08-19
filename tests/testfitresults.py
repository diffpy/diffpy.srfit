#!/usr/bin/env python
"""Tests for fitresults module."""

import unittest
import os.path

from diffpy.srfit.fitbase.fitrecipe import FitRecipe
from diffpy.srfit.fitbase.fitresults import initializeRecipe

thisfile = locals().get('__file__', 'testpdf.py')
tests_dir = os.path.dirname(os.path.abspath(thisfile))
testdata_dir = os.path.join(tests_dir, 'testdata')

class TestInitializeRecipe(unittest.TestCase):

    def setUp(self):
        self.recipe = recipe = FitRecipe("recipe")
        recipe.newVar("A", 0)
        recipe.newVar("sig", 0)
        recipe.newVar("x0", 0)
        self.filename = "results.res"

        self.Aval = 5.77619823e-01
        self.sigval = -9.22758690e-01
        self.x0val = 6.12422115e+00
        return

    def testInitializeFromFileName(self):
        recipe = self.recipe
        self.assertEquals(0, recipe.A.value)
        self.assertEquals(0, recipe.sig.value)
        self.assertEquals(0, recipe.x0.value)
        initializeRecipe(recipe, os.path.join(testdata_dir, self.filename))
        self.assertAlmostEquals(self.Aval, recipe.A.value)
        self.assertAlmostEquals(self.sigval, recipe.sig.value)
        self.assertAlmostEquals(self.x0val, recipe.x0.value)
        return

    def testInitializeFromFileObj(self):
        recipe = self.recipe
        self.assertEquals(0, recipe.A.value)
        self.assertEquals(0, recipe.sig.value)
        self.assertEquals(0, recipe.x0.value)
        infile = file(os.path.join(testdata_dir, self.filename), 'r')
        initializeRecipe(recipe, infile)
        self.assertFalse(infile.closed)
        infile.close()
        self.assertAlmostEquals(self.Aval, recipe.A.value)
        self.assertAlmostEquals(self.sigval, recipe.sig.value)
        self.assertAlmostEquals(self.x0val, recipe.x0.value)
        return


    def testInitializeFromString(self):
        recipe = self.recipe
        self.assertEquals(0, recipe.A.value)
        self.assertEquals(0, recipe.sig.value)
        self.assertEquals(0, recipe.x0.value)
        infile = file(os.path.join(testdata_dir, self.filename), 'r')
        resstr = infile.read()
        infile.close()
        initializeRecipe(recipe, resstr)
        self.assertAlmostEquals(self.Aval, recipe.A.value)
        self.assertAlmostEquals(self.sigval, recipe.sig.value)
        self.assertAlmostEquals(self.x0val, recipe.x0.value)
        return

if __name__ == "__main__":

    unittest.main()
