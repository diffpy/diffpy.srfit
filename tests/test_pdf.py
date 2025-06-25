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

import io
import pickle
import unittest
from itertools import chain

import numpy
import pytest

from diffpy.srfit.exceptions import SrFitError
from diffpy.srfit.pdf import PDFContribution, PDFGenerator, PDFParser

# ----------------------------------------------------------------------------


def testParser1(datafile):
    data = datafile("ni-q27r100-neutron.gr")
    parser = PDFParser()
    parser.parseFile(data)

    meta = parser._meta

    assert data == meta["filename"]
    assert 1 == meta["nbanks"]
    assert "N" == meta["stype"]
    assert 27 == meta["qmax"]
    assert 300 == meta.get("temperature")
    assert meta.get("qdamp") is None
    assert meta.get("qbroad") is None
    assert meta.get("spdiameter") is None
    assert meta.get("scale") is None
    assert meta.get("doping") is None

    x, y, dx, dy = parser.getData()
    assert dx is None
    assert dy is None

    testx = numpy.linspace(0.01, 100, 10000)
    diff = testx - x
    res = numpy.dot(diff, diff)
    assert 0 == pytest.approx(res)

    testy = numpy.array(
        [
            1.144,
            2.258,
            3.312,
            4.279,
            5.135,
            5.862,
            6.445,
            6.875,
            7.150,
            7.272,
        ]
    )
    diff = testy - y[:10]
    res = numpy.dot(diff, diff)
    assert 0 == pytest.approx(res)

    return


def testParser2(datafile):
    data = datafile("si-q27r60-xray.gr")
    parser = PDFParser()
    parser.parseFile(data)

    meta = parser._meta

    assert data == meta["filename"]
    assert 1 == meta["nbanks"]
    assert "X" == meta["stype"]
    assert 27 == meta["qmax"]
    assert 300 == meta.get("temperature")
    assert meta.get("qdamp") is None
    assert meta.get("qbroad") is None
    assert meta.get("spdiameter") is None
    assert meta.get("scale") is None
    assert meta.get("doping") is None

    x, y, dx, dy = parser.getData()
    testx = numpy.linspace(0.01, 60, 5999, endpoint=False)
    diff = testx - x
    res = numpy.dot(diff, diff)
    assert 0 == pytest.approx(res)

    testy = numpy.array(
        [
            0.1105784,
            0.2199684,
            0.3270088,
            0.4305913,
            0.5296853,
            0.6233606,
            0.7108060,
            0.7913456,
            0.8644501,
            0.9297440,
        ]
    )
    diff = testy - y[:10]
    res = numpy.dot(diff, diff)
    assert 0 == pytest.approx(res)

    testdy = numpy.array(
        [
            0.001802192,
            0.003521449,
            0.005079115,
            0.006404892,
            0.007440527,
            0.008142955,
            0.008486813,
            0.008466340,
            0.008096858,
            0.007416456,
        ]
    )
    diff = testdy - dy[:10]
    res = numpy.dot(diff, diff)
    assert 0 == pytest.approx(res)

    assert dx is None
    return


def testGenerator(
    diffpy_srreal_available, diffpy_structure_available, datafile
):
    if not diffpy_structure_available:
        pytest.skip("diffpy.structure package not available")
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
    stru.read(ciffile)
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
        par.setValue(1.0)
        assert 1.0 == par.getValue()
        par.setValue(defval)
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


def test_setQmin(diffpy_structure_available, diffpy_srreal_available):
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


def test_setQmax(diffpy_structure_available, diffpy_srreal_available):
    """Check PDFContribution.setQmax()"""
    if not diffpy_structure_available:
        pytest.skip("diffpy.structure package not available")
    from diffpy.structure import Structure

    if not diffpy_srreal_available:
        pytest.skip("diffpy.structure package not available")

    pc = PDFContribution("pdf")
    pc.setQmax(21)
    pc.addStructure("empty", Structure())
    assert 21 == pc.empty.getQmax()
    pc.setQmax(22)
    assert 22 == pc.getQmax()
    assert 22 == pc.empty.getQmax()
    return


def test_getQmax(diffpy_structure_available, diffpy_srreal_available):
    """Check PDFContribution.getQmax()"""
    if not diffpy_structure_available:
        pytest.skip("diffpy.structure package not available")
    from diffpy.structure import Structure

    if not diffpy_srreal_available:
        pytest.skip("diffpy.structure package not available")

    # cover all code branches in PDFContribution._getMetaValue
    # (1) contribution metadata
    pc1 = PDFContribution("pdf")
    assert pc1.getQmax() is None
    pc1.setQmax(17)
    assert 17 == pc1.getQmax()
    # (2) contribution metadata
    pc2 = PDFContribution("pdf")
    pc2.addStructure("empty", Structure())
    pc2.empty.setQmax(18)
    assert 18 == pc2.getQmax()
    # (3) profile metadata
    pc3 = PDFContribution("pdf")
    pc3.profile.meta["qmax"] = 19
    assert 19 == pc3.getQmax()
    return


def test_savetxt(
    diffpy_structure_available, diffpy_srreal_available, datafile
):
    "check PDFContribution.savetxt()"
    if not diffpy_structure_available:
        pytest.skip("diffpy.structure package not available")
    from diffpy.structure import Structure

    if not diffpy_srreal_available:
        pytest.skip("diffpy.structure package not available")

    pc = PDFContribution("pdf")
    pc.loadData(datafile("si-q27r60-xray.gr"))
    pc.setCalculationRange(0, 10)
    pc.addStructure("empty", Structure())
    fp = io.BytesIO()
    with pytest.raises(SrFitError):
        pc.savetxt(fp)
    pc.evaluate()
    pc.savetxt(fp)
    txt = fp.getvalue().decode()
    nlines = len(txt.strip().split("\n"))
    assert 1001 == nlines
    return


def test_pickling(
    diffpy_structure_available, diffpy_srreal_available, datafile
):
    "validate PDFContribution.residual() after pickling."
    if not diffpy_structure_available:
        pytest.skip("diffpy.structure package not available")
    from diffpy.structure import loadStructure

    if not diffpy_srreal_available:
        pytest.skip("diffpy.structure package not available")

    pc = PDFContribution("pdf")
    pc.loadData(datafile("ni-q27r100-neutron.gr"))
    ni = loadStructure(datafile("ni.cif"))
    ni.Uisoequiv = 0.003
    pc.addStructure("ni", ni)
    pc.setCalculationRange(0, 10)
    pc2 = pickle.loads(pickle.dumps(pc))
    res0 = pc.residual()
    assert numpy.array_equal(res0, pc2.residual())
    for p in chain(pc.iterPars("Uiso"), pc2.iterPars("Uiso")):
        p.value = 0.004
    res1 = pc.residual()
    assert not numpy.allclose(res0, res1)
    assert numpy.array_equal(res1, pc2.residual())
    return


if __name__ == "__main__":
    unittest.main()
