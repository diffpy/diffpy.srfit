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


class PDFParser(ProfileParser):
    """Class for holding a diffraction pattern.

    Attributes
    ----------
    _format
        Name of the data format that this parses (string, default
        ""). The format string is a unique identifier for the data
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
        Independent variable from the chosen bank
    _y
        Profile from the chosen bank
    _dx
        Uncertainty in independent variable from the chosen bank
    _dy
        Uncertainty in profile from the chosen bank
    _meta
        A dictionary containing metadata read from the file.

    General Metadata

    Attributes
    ----------
    filename
        The name of the file from which data was parsed. This key
        will not exist if data was not read from file.
    nbanks
        The number of banks parsed.
    bank
        The chosen bank number.

    Metadata
    ----------
    stype
        The scattering type ("X", "N")
    qmin
        Minimum scattering vector (float)
    qmax
        Maximum scattering vector (float)
    qdamp
        Resolution damping factor (float)
    qbroad
        Resolution broadening factor (float)
    spdiameter
        Nanoparticle diameter (float)
    scale
        Data scale (float)
    temperature
        Temperature (float)
    doping
        Doping (float)


    These may appear in the metadata dictionary.
    """

    _format = "PDF"

    def _parse_metadata(self, filename):
        """Return the metadata read from a PDFgetX or PDFgetN header."""
        return self._parse_header(Path(filename).read_text())

    def _parse_header(self, patstring):
        """Return the metadata read from the header part of a
        pattern."""
        rx = {"f": r"[-+]?(\d+(\.\d*)?|\d*\.\d+)([eE][-+]?\d+)?"}
        header = self._split_header(patstring, rx)

        # find where the metadata starts
        metadata = ""
        res = re.search(r"^#+\ +metadata\b\n", header, re.M)
        if res:
            metadata = header[res.end() :]
            header = header[: res.start()]

        # parse header
        meta = {}
        # stype
        if re.search("(x-?ray|PDFgetX)", header, re.I):
            meta["stype"] = "X"
        elif re.search("(neutron|PDFgetN)", header, re.I):
            meta["stype"] = "N"
        # qmin
        regexp = r"\bqmin *= *(%(f)s)\b" % rx
        res = re.search(regexp, header, re.I)
        if res:
            meta["qmin"] = float(res.groups()[0])
        # qmax
        regexp = r"\bqmax *= *(%(f)s)\b" % rx
        res = re.search(regexp, header, re.I)
        if res:
            meta["qmax"] = float(res.groups()[0])
        # qdamp
        regexp = r"\b(?:qdamp|qsig) *= *(%(f)s)\b" % rx
        res = re.search(regexp, header, re.I)
        if res:
            meta["qdamp"] = float(res.groups()[0])
        # qbroad
        regexp = r"\b(?:qbroad|qalp) *= *(%(f)s)\b" % rx
        res = re.search(regexp, header, re.I)
        if res:
            meta["qbroad"] = float(res.groups()[0])
        # spdiameter
        regexp = r"\bspdiameter *= *(%(f)s)\b" % rx
        res = re.search(regexp, header, re.I)
        if res:
            meta["spdiameter"] = float(res.groups()[0])
        # dscale
        regexp = r"\bdscale *= *(%(f)s)\b" % rx
        res = re.search(regexp, header, re.I)
        if res:
            meta["scale"] = float(res.groups()[0])
        # temperature
        regexp = r"\b(?:temp|temperature|T)\ *=\ *(%(f)s)\b" % rx
        res = re.search(regexp, header)
        if res:
            meta["temperature"] = float(res.groups()[0])
        # doping
        regexp = r"\b(?:x|doping)\ *=\ *(%(f)s)\b" % rx
        res = re.search(regexp, header)
        if res:
            meta["doping"] = float(res.groups()[0])

        # parsing general metadata
        if metadata:
            regexp = r"\b(\w+)\ *=\ *(%(f)s)\b" % rx
            while True:
                res = re.search(regexp, metadata, re.M)
                if res:
                    meta[res.groups()[0]] = float(res.groups()[1])
                    metadata = metadata[res.end() :]
                else:
                    break

        return meta

    @staticmethod
    def _split_header(patstring, rx):
        """Return the header part of a pattern."""
        res = re.search(r"^#+ start data\s*(?:#.*\s+)*", patstring, re.M)
        # start_data is position where the first data line starts
        if res:
            start_data = res.end()
        else:
            # find line that starts with a floating point number
            regexp = r"^\s*%(f)s" % rx
            res = re.search(regexp, patstring, re.M)
            if res:
                start_data = res.start()
            else:
                start_data = 0
        return patstring[:start_data]


# End of PDFParser
