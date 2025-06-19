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

import pytest

from diffpy.srfit.fitbase.fitrecipe import FitRecipe
from diffpy.srfit.fitbase.fitresults import initializeRecipe


def testInitializeFromFileName(datafile):
    recipe = FitRecipe("recipe")
    recipe.newVar("A", 0)
    recipe.newVar("sig", 0)
    recipe.newVar("x0", 0)
    filename = datafile("results.res")
    Aval = 5.77619823e-01
    sigval = -9.22758690e-01
    x0val = 6.12422115e00

    assert 0 == recipe.A.value
    assert 0 == recipe.sig.value
    assert 0 == recipe.x0.value
    initializeRecipe(recipe, filename)
    assert Aval == pytest.approx(recipe.A.value)
    assert sigval == pytest.approx(recipe.sig.value)
    assert x0val == pytest.approx(recipe.x0.value)
    return


def testInitializeFromFileObj(datafile):
    recipe = FitRecipe("recipe")
    recipe.newVar("A", 0)
    recipe.newVar("sig", 0)
    recipe.newVar("x0", 0)
    filename = datafile("results.res")
    Aval = 5.77619823e-01
    sigval = -9.22758690e-01
    x0val = 6.12422115e00

    assert 0 == recipe.A.value
    assert 0 == recipe.sig.value
    assert 0 == recipe.x0.value
    infile = open(filename, "r")
    initializeRecipe(recipe, infile)
    assert not infile.closed
    infile.close()
    assert Aval == pytest.approx(recipe.A.value)
    assert sigval == pytest.approx(recipe.sig.value)
    assert x0val == pytest.approx(recipe.x0.value)
    return


def testInitializeFromString(datafile):
    recipe = FitRecipe("recipe")
    recipe.newVar("A", 0)
    recipe.newVar("sig", 0)
    recipe.newVar("x0", 0)
    filename = datafile("results.res")
    Aval = 5.77619823e-01
    sigval = -9.22758690e-01
    x0val = 6.12422115e00

    assert 0 == recipe.A.value
    assert 0 == recipe.sig.value
    assert 0 == recipe.x0.value
    infile = open(filename, "r")
    resstr = infile.read()
    infile.close()
    initializeRecipe(recipe, resstr)
    assert Aval == pytest.approx(recipe.A.value)
    assert sigval == pytest.approx(recipe.sig.value)
    assert x0val == pytest.approx(recipe.x0.value)
    return


if __name__ == "__main__":

    unittest.main()
