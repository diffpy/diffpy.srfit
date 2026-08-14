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

import re
from pathlib import Path

from diffpy.srfit.fitbase.profileparser import ProfileParser

_FLOAT_RX = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"


class PDFParser(ProfileParser):
    """Parser for PDF diffraction pattern data.

    PDFgetX and PDFgetN write their header as plain ``name = value``
    pairs, including ``stype = X`` or ``stype = N`` for the scattering
    type, so this class parses files identically to ``ProfileParser``
    for those. Some facilities (e.g. NOMAD at SNS) instead prepend a
    free-text instrument comment, so this class also falls back to
    scanning that comment for the scattering type, ``qmin``, ``qmax``,
    ``qdamp``, and ``qbroad`` when they are not already present as
    ``name = value`` pairs.

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
    qdamp
        The Q-resolution damping factor (float).
    qbroad
        The Q-resolution broadening factor (float).

    These, along with any other ``name = value`` pairs in the header,
    may appear in the metadata dictionary.
    """

    _format = "PDF"

    def _parse_metadata(self, filename):
        """Return the metadata read from a PDFgetX or PDFgetN header.

        This calls ``ProfileParser``'s ``name = value`` based parsing
        first, then falls back to scanning the free-text instrument
        comments some facilities prepend to their files for the
        scattering type and Q-resolution parameters that
        such comments are not already covered by a ``name = value``
        pair.

        Parameters
        ----------
        filename : str or Path
            The name of the file to parse.

        Returns
        -------
        dict
            The metadata read from the file header.
        """
        metadata = super()._parse_metadata(filename)
        self._parse_comment_metadata(Path(filename).read_text(), metadata)
        return metadata

    @staticmethod
    def _parse_comment_metadata(header_text, metadata):
        """Fill in stype, qmin, qmax, qdamp, and qbroad from free-text
        instrument comments, without overwriting values already found
        by the ``name = value`` based parsing.

        Parameters
        ----------
        header_text : str
            The full text of the file being parsed.
        meta : dict
            The metadata dictionary to update in place.

        Returns
        -------
        dict
            The updated metadata dictionary.
        """
        if "stype" not in metadata:
            if re.search(r"(x-?ray|PDFgetX)", header_text, re.I):
                metadata["stype"] = "X"
            elif re.search(r"(neutron|PDFgetN)", header_text, re.I):
                metadata["stype"] = "N"
        regexes = {
            "qmin": r"\bqmin *= *(%s)\b" % _FLOAT_RX,
            "qmax": r"\bqmax *= *(%s)\b" % _FLOAT_RX,
            "qdamp": r"\b(?:qdamp|qsig) *= *(%s)\b" % _FLOAT_RX,
            "qbroad": r"\b(?:qbroad|qalp) *= *(%s)\b" % _FLOAT_RX,
        }
        for key, pattern in regexes.items():
            if key in metadata:
                continue
            res = re.search(pattern, header_text, re.I)
            if res:
                metadata[key] = float(res.group(1))
        return metadata


# End of PDFParser
