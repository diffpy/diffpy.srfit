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
"""Tests for pdf.pdfcontribution module."""

import io
import pickle
from itertools import chain

import numpy
import pytest

from diffpy.srfit.exceptions import SrFitError
from diffpy.srfit.pdf import PDFContribution

# ----------------------------------------------------------------------------


def test_pdfcontribution_loadData(datafile):
    """LoadData passes the PDF metadata on to the built-in profile."""
    contribution = PDFContribution("pdf")
    contribution.loadData(datafile("si-q27r60-xray.gr"))

    expected_metadata = {
        "version": "diffpy.pdfgetx-2.4.0",
        "dataformat": "QA",
        "outputtype": "gr",
        "stype": "X",
        "composition": "Si",
        "bgscale": 1.0,
        "rpoly": 0.9,
        "qmaxinst": 29.0,
        "qmin": 0.01,
        "qmax": 27.0,
        "rmin": 0.0,
        "rmax": 60.0,
        "rstep": 0.01,
        "temperature": 300.0,
        "filename": str(datafile("si-q27r60-xray.gr")),
        "bank": 0,
        "nbanks": 1,
    }
    actual_metadata = contribution.profile.meta
    assert actual_metadata == expected_metadata
    actual_point_count = len(contribution.profile.xobs)
    expected_point_count = 5999
    assert actual_point_count == expected_point_count


def test_setQmax(diffpy_srreal_available):
    """Check PDFContribution.setQmax()"""
    from diffpy.structure import Structure

    if not diffpy_srreal_available:
        pytest.skip("diffpy.srreal package not available")

    pc = PDFContribution("pdf")
    pc.setQmax(21)
    pc.addStructure("empty", Structure())
    assert 21 == pc.empty.getQmax()
    pc.setQmax(22)
    assert 22 == pc.getQmax()
    assert 22 == pc.empty.getQmax()
    return


def test_getQmax(diffpy_srreal_available):
    """Check PDFContribution.getQmax()"""
    from diffpy.structure import Structure

    if not diffpy_srreal_available:
        pytest.skip("diffpy.srreal package not available")

    # cover all code branches in PDFContribution._get_meta_value
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


def test_savetxt(diffpy_srreal_available, datafile):
    "check PDFContribution.savetxt()"
    from diffpy.structure import Structure

    if not diffpy_srreal_available:
        pytest.skip("diffpy.srreal package not available")

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


def test_pickling(diffpy_srreal_available, datafile):
    "validate PDFContribution.residual() after pickling."
    from diffpy.structure import loadStructure

    if not diffpy_srreal_available:
        pytest.skip("diffpy.srreal package not available")

    pc = PDFContribution("pdf")
    pc.loadData(datafile("ni-q27r100-neutron.gr"))
    ciffile = datafile("ni.cif")
    cif_path = str(ciffile)
    ni = loadStructure(cif_path)
    ni.Uisoequiv = 0.003
    pc.addStructure("ni", ni)
    pc.setCalculationRange(0, 10)
    pc2 = pickle.loads(pickle.dumps(pc))
    res0 = pc.residual()
    assert numpy.array_equal(res0, pc2.residual())
    for p in chain(
        pc.iterate_over_parameters("Uiso"), pc2.iterate_over_parameters("Uiso")
    ):
        p.value = 0.004
    res1 = pc.residual()
    assert not numpy.allclose(res0, res1)
    assert numpy.array_equal(res1, pc2.residual())
    return
