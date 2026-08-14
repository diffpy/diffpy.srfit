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
"""Tests for pdf.pdfparser module."""

import numpy
import pytest

from diffpy.srfit.pdf import PDFParser

# ----------------------------------------------------------------------------


def approx_or_none(expected_values):
    """Wrap expected_values in pytest.approx, unless it is None."""
    if expected_values is None:
        return None
    return pytest.approx(expected_values)


@pytest.mark.parametrize(
    "input_filename, expected_x, expected_y, expected_dy",
    [
        # C1: A neutron PDF written by PDFgetN, which has no dx or dy
        # columns.
        # Expected: x and y are read correctly, and dx and dy are None.
        (
            "ni-q27r100-neutron.gr",
            numpy.linspace(0.01, 100, 10000),
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
            ],
            None,
        ),
        # C2: An x-ray PDF written by PDFgetX2, which has a dy column
        # and a negative dx column.
        # Expected: x, y, and dy are read correctly, and the invalid
        # negative dx column is dropped.
        (
            "si-q27r60-xray.gr",
            numpy.linspace(0.01, 60, 5999, endpoint=False),
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
            ],
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
            ],
        ),
    ],
)
def test_pdfparser_data(
    datafile, as_list, input_filename, expected_x, expected_y, expected_dy
):
    """PDFParser reads the x, y, and dy arrays correctly, and always
    drops the invalid dx column."""
    parser = PDFParser()
    parser.parse_file(datafile(input_filename))

    actual_x, actual_y, actual_dx, actual_dy = parser.get_data()
    actual_dy = as_list(actual_dy)
    if actual_dy is not None:
        # Compare only the first 10 values
        actual_dy = actual_dy[:10]
    assert actual_dx is None
    assert actual_x.tolist() == pytest.approx(expected_x.tolist())
    assert actual_y[:10].tolist() == pytest.approx(expected_y)
    assert actual_dy == approx_or_none(expected_dy)


# PDFParser inherits ProfileParser's hooks unchanged: PDFgetX and
# PDFgetN headers are already plain name = value pairs, including
# stype = X or stype = N for the scattering type. The metadata below
# reaches PDFGenerator, which uses stype, qmin and qmax to set the
# scattering type and the Q range, so losing a key silently changes a
# refinement.
@pytest.mark.parametrize(
    "input_filename, expected_metadata",
    [
        # C1: An x-ray PDF written by PDFgetX2.
        # Expected: The header yields the x-ray scattering type,
        # the Q range and the rest of the diffpy.pdfgetx config.
        (
            "si-q27r60-xray.gr",
            {
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
                "bank": 0,
                "nbanks": 1,
            },
        ),
        # C2: A neutron PDF written by PDFgetN.
        # Expected: The header yields the neutron scattering type,
        # the Q range and the rest of the xPDFsuite config.
        (
            "ni-q27r100-neutron.gr",
            {
                "wavelength": 1.333,
                "dataformat": "QA",
                "inputfile": "npdf_03315.chi",
                "backgroundfile": "npdf_03001.chi",
                "stype": "N",
                "bgscale": 1.0,
                "composition": "Ni",
                "outputtype": "gr",
                "qmaxinst": 27.0,
                "qmin": 0.87,
                "qmax": 27.0,
                "temperature": 300.0,
                "rmax": 100.0,
                "rmin": 0.0,
                "rstep": 0.01,
                "rpoly": 0.9,
                "inputdir": "/data/npdf/chi",
                "savedir": "/data/npdf/gr",
                "bank": 0,
                "nbanks": 1,
            },
        ),
    ],
)
def test_pdfparser_metadata(datafile, input_filename, expected_metadata):
    """PDF specific metadata survives the load_data based parse_file."""
    parser = PDFParser()
    parser.parse_file(datafile(input_filename))
    actual_metadata = parser.get_metadata()
    # add the filename key to the expected metadata for comparison
    expected_metadata["filename"] = str(datafile(input_filename))
    assert actual_metadata == expected_metadata


# ----------------------------------------------------------------------------
# parseFile is deprecated in favor of parse_file. The old name must still
# work, emit a DeprecationWarning, and forward to the new implementation.


def test_pdfparser_deprecated_parseFile(datafile):
    """The deprecated parseFile warns and delegates to parse_file."""
    input_filename = datafile("si-q27r60-xray.gr")
    expected_parser = PDFParser()
    expected_parser.parse_file(input_filename)
    actual_parser = PDFParser()
    with pytest.warns(DeprecationWarning):
        actual_parser.parseFile(input_filename)

    actual_metadata = actual_parser.get_metadata()
    expected_metadata = expected_parser.get_metadata()
    assert actual_metadata == expected_metadata

    actual_x, actual_y, actual_dx, actual_dy = actual_parser.get_data()
    expected_x, expected_y, expected_dx, expected_dy = (
        expected_parser.get_data()
    )
    assert actual_x.tolist() == expected_x.tolist()
    assert actual_y.tolist() == expected_y.tolist()
    assert actual_dx == expected_dx
    assert actual_dy.tolist() == expected_dy.tolist()
