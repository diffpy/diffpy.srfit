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
        assert defval == par.get_value()
        # Test setting values
        par.set_value(1.0)
        assert 1.0 == par.get_value()
        par.set_value(defval)
        assert defval == par.get_value()

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
    calc.scatteringfactortable = "N"
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


@pytest.mark.parametrize(
    "input_stype, expected_scattering_type",
    [
        # Test that the scattering type sets the calculator's scattering
        # factor table.
        # C1: X-ray scattering type given as a short code.
        # Expected: The X-ray table is used and "X" is reported.
        ("X", "X"),
        # C2: Neutron scattering type given as a short code.
        # Expected: The neutron table is used and "N" is reported.
        ("N", "N"),
        # C3: Electron scattering type given as a short code.
        # Expected: The electron table is used and "E" is reported.
        ("E", "E"),
        # C4: Neutron scattering type given as its registered name.
        # Expected: The neutron table is used and "N" is reported.
        ("neutron", "N"),
    ],
)
def test_set_scattering_type(
    diffpy_srreal_available, input_stype, expected_scattering_type
):
    if not diffpy_srreal_available:
        pytest.skip("diffpy.srreal package not available")
    generator = PDFGenerator()
    generator.setScatteringType(input_stype)
    actual_scattering_type = generator.getScatteringType()
    actual_table_type = generator._calc.scatteringfactortable.radiationType()
    actual_meta_stype = generator.meta["stype"]
    assert actual_scattering_type == expected_scattering_type
    assert actual_table_type == expected_scattering_type
    assert actual_meta_stype == expected_scattering_type


def test_set_scattering_type_bad(diffpy_srreal_available):
    # C1: Unknown scattering type.
    # Expected: ValueError is raised and the previous scattering type is
    # kept.
    if not diffpy_srreal_available:
        pytest.skip("diffpy.srreal package not available")
    generator = PDFGenerator()
    generator.setScatteringType("N")
    with pytest.raises(ValueError, match="Unknown type 'bogus'"):
        generator.setScatteringType("bogus")
    actual_scattering_type = generator.getScatteringType()
    actual_meta_stype = generator.meta["stype"]
    assert actual_scattering_type == "N"
    assert actual_meta_stype == "N"
