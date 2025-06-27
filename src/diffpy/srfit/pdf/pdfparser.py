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

import numpy

from diffpy.srfit.exceptions import ParseError
from diffpy.srfit.fitbase.profileparser import ProfileParser


class PDFParser(ProfileParser):
    """Class for holding a diffraction pattern.

    Attributes

    _format     --  Name of the data format that this parses (string, default
                    ""). The format string is a unique identifier for the data
                    format handled by the parser.
    _banks      --  The data from each bank. Each bank contains a
                    (x, y, dx, dy) tuple:
                    x       --  A numpy array containing the independent
                                variable read from the file.
                    y       --  A numpy array containing the profile
                                from the file.
                    dx      --  A numpy array containing the uncertainty in x
                                read from the file. This is 0 if the
                                uncertainty cannot be read.
                    dy      --  A numpy array containing the uncertainty read
                                from the file. This is 0 if the uncertainty
                                cannot be read.
    _x          --  Independent variable from the chosen bank
    _y          --  Profile from the chosen bank
    _dx         --  Uncertainty in independent variable from the chosen bank
    _dy         --  Uncertainty in profile from the chosen bank
    _meta       --  A dictionary containing metadata read from the file.

    General Metadata

    filename    --  The name of the file from which data was parsed. This key
                    will not exist if data was not read from file.
    nbanks      --  The number of banks parsed.
    bank        --  The chosen bank number.

    Metadata - These may appear in the metadata dictionary

    stype       --  The scattering type ("X", "N")
    qmin        --  Minimum scattering vector (float)
    qmax        --  Maximum scattering vector (float)
    qdamp       --  Resolution damping factor (float)
    qbroad      --  Resolution broadening factor (float)
    spdiameter  --  Nanoparticle diameter (float)
    scale       --  Data scale (float)
    temperature --  Temperature (float)
    doping      --  Doping (float)
    """

    _format = "PDF"

    def parseString(self, patstring):
        """Parse a string and set the _x, _y, _dx, _dy and _meta variables.

        When _dx or _dy cannot be obtained in the data format it is set to 0.

        This wipes out the currently loaded data and selected bank number.

        Arguments
        patstring   --  A string containing the pattern

        Raises ParseError if the string cannot be parsed
        """
        # useful regex patterns:
        rx = {"f": r"[-+]?(\d+(\.\d*)?|\d*\.\d+)([eE][-+]?\d+)?"}
        # find where does the data start
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
        header = patstring[:start_data]
        databody = patstring[start_data:].strip()

        # find where the metadata starts
        metadata = ""
        res = re.search(r"^#+\ +metadata\b\n", header, re.M)
        if res:
            metadata = header[res.end() :]
            header = header[: res.start()]

        # parse header
        meta = self._meta
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

        # read actual data - robs, Gobs, drobs, dGobs
        inf_or_nan = re.compile("(?i)^[+-]?(NaN|Inf)\\b")
        has_drobs = True
        has_dGobs = True
        # raise ParseError if something goes wrong
        robs = []
        Gobs = []
        drobs = []
        dGobs = []
        try:
            for line in databody.split("\n"):
                v = line.split()
                # there should be at least 2 value in the line
                robs.append(float(v[0]))
                Gobs.append(float(v[1]))
                # drobs is valid if all values are defined and positive
                has_drobs = (
                    has_drobs and len(v) > 2 and not inf_or_nan.match(v[2])
                )
                if has_drobs:
                    v2 = float(v[2])
                    has_drobs = v2 > 0.0
                    drobs.append(v2)
                # dGobs is valid if all values are defined and positive
                has_dGobs = (
                    has_dGobs and len(v) > 3 and not inf_or_nan.match(v[3])
                )
                if has_dGobs:
                    v3 = float(v[3])
                    has_dGobs = v3 > 0.0
                    dGobs.append(v3)
        except (ValueError, IndexError) as err:
            raise ParseError(err)
        if has_drobs:
            drobs = numpy.asarray(drobs)
        else:
            drobs = None
        if has_dGobs:
            dGobs = numpy.asarray(dGobs)
        else:
            dGobs = None

        robs = numpy.asarray(robs)
        Gobs = numpy.asarray(Gobs)

        self._banks.append([robs, Gobs, drobs, dGobs])
        return


# End of PDFParser
