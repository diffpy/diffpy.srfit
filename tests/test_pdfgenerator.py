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
"""Tests for pdf.pdfgenerator module."""

import numpy
import pytest

from diffpy.srfit.pdf import PDFGenerator

# ----------------------------------------------------------------------------


def testGenerator(diffpy_srreal_available, datafile):
    if not diffpy_srreal_available:
        pytest.skip("diffpy.srreal package not available")

    from diffpy.srreal.pdfcalculator import PDFCalculator
    from diffpy.structure import PDFFitStructure

    qmax = 27.0
    gen = PDFGenerator()
    gen.setScatteringType("N")
    assert "N" == gen.getScatteringType()
    gen.setQmax(qmax)
    assert qmax == pytest.approx(gen.getQmax())

    stru = PDFFitStructure()
    ciffile = datafile("ni.cif")
    cif_path = str(ciffile)
    stru.read(cif_path)
    for i in range(4):
        stru[i].Bisoequiv = 1
    gen.setStructure(stru)

    calc = gen._calc
    # Test parameters
    for par in gen.iterPars(recurse=False):
        pname = par.name
        defval = calc._getDoubleAttr(pname)
        assert defval == par.getValue()
        # Test setting values
        par.set_value(1.0)
        assert 1.0 == par.getValue()
        par.set_value(defval)
        assert defval == par.getValue()

    r = numpy.arange(0, 10, 0.1)
    y = gen(r)

    # Now create a reference PDF. Since the calculator is testing its
    # output, we just have to make sure we can calculate from the
    # PDFGenerator interface.

    calc = PDFCalculator()
    calc.rstep = r[1] - r[0]
    calc.rmin = r[0]
    calc.rmax = r[-1] + 0.5 * calc.rstep
    calc.qmax = qmax
    calc.setScatteringFactorTableByType("N")
    calc.eval(stru)
    yref = calc.pdf

    diff = y - yref
    res = numpy.dot(diff, diff)
    assert 0 == pytest.approx(res)
    return


def test_setQmin(diffpy_srreal_available):
    """Verify qmin is propagated to the calculator object."""
    if not diffpy_srreal_available:
        pytest.skip("diffpy.srreal package not available")

    gen = PDFGenerator()
    assert 0 == gen.getQmin()
    assert 0 == gen._calc.qmin
    gen.setQmin(0.93)
    assert 0.93 == gen.getQmin()
    assert 0.93 == gen._calc.qmin
    return
