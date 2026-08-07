#!/usr/bin/env python
##############################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2008 The Trustees of Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################
"""This module contains parsers for PDF data.

PDFParser is suitable for parsing data generated from PDFGetN and
PDFGetX.

See the class documentation for more information.
"""

__all__ = ["PDFParser"]

from diffpy.srfit.fitbase.profileparser import ProfileParser


class PDFParser(ProfileParser):
    """Parser for PDF diffraction pattern data.

    PDFgetX and PDFgetN write their header as plain ``name = value``
    pairs, including ``stype = X`` or ``stype = N`` for the scattering
    type, so this class parses files identically to ``ProfileParser``
    and only sets ``_format`` to identify PDF data.

    Attributes
    ----------
    _format
        The name of the data format that this parses (string, default
        ``""``). The format string is a unique identifier for the data
        format handled by the parser.
    _banks
        The data from each bank. Each bank contains a
        (x, y, dx, dy) tuple:
        x
            A numpy array containing the independent
            variable read from the file.
        y
            A numpy array containing the profile
            from the file.
        dx
            A numpy array containing the uncertainty in x
            read from the file. This is None if the
            uncertainty cannot be read.
        dy
            A numpy array containing the uncertainty read
            from the file. This is None if the uncertainty
            cannot be read.
    _x
        The independent variable from the chosen bank.
    _y
        The profile from the chosen bank.
    _dx
        The uncertainty in independent variable from the chosen bank.
    _dy
        The uncertainty in profile from the chosen bank.
    _meta
        A dictionary containing metadata read from the file.

    General Metadata
    -----------------
    filename
        The name of the file from which data was parsed. This key
        will not exist if data was not read from file.
    nbanks
        The number of banks parsed.
    bank
        The chosen bank number.

    Metadata
    --------
    stype
        The scattering type ("X", "N").
    qmin
        The minimum scattering vector (float).
    qmax
        The maximum scattering vector (float).

    These, along with any other ``name = value`` pairs in the header,
    may appear in the metadata dictionary.
    """

    _format = "PDF"


# End of PDFParser
