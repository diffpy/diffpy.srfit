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

"""This module contains classes for parsing profiles from files.

ProfileParser is a base class for parsing data. It can interact with a Profile
object to automatically set the Profile's data and metadata. Each specific file
format must be encapsulated in a ProfileParser subclass.

See the class documentation for more information.
"""


from diffpy.srfit.exceptions import ParseError


class ProfileParser(object):
    """Class for parsing data from a or string.

    Attributes

    _format     --  Name of the data format that this parses (string, default
                    ""). The format string is a unique identifier for the data
                    format handled by the parser.
    _banks      --  The data from each bank. Each bank contains a (x, y, dx, dy)
                    tuple:
                    x       --  A numpy array containing the independent
                                variable read from the file.
                    y       --  A numpy array containing the profile
                                from the file.
                    dx      --  A numpy array containing the uncertainty in x
                                read from the file. This is None if the
                                uncertainty cannot be read.
                    dy      --  A numpy array containing the uncertainty read
                                from the file. This is None if the uncertainty
                                cannot be read.
    _x          --  Indpendent variable from the chosen bank
    _y          --  Profile from the chosen bank
    _dx         --  Uncertainty in independent variable from the chosen bank
    _dy         --  Uncertainty in profile from the chosen bank
    _meta       --  A dictionary containing metadata read from the file.

    General Metadata

    filename    --  The name of the file from which data was parsed. This key
                    will not exist if data was not read from file.
    nbanks      --  The number of banks parsed.
    bank        --  The chosen bank number.

    """

    _format = ""

    def __init__(self):
        """Initialize the attributes."""
        self._banks = []
        self._meta = {}
        self._x = None
        self._y = None
        self._dx = None
        self._dy = None
        return

    def getFormat(self):
        """Get the format string."""
        return self._format

    def parseString(self, patstring):
        """Parse a string and set the _x, _y, _dx, _dy and _meta variables.

        When _dx or _dy cannot be obtained in the data format it is set to
        None.

        This wipes out the currently loaded data and selected bank number.

        Arguments
        patstring   --  A string containing the pattern

        Raises ParseError if the string cannot be parsed

        """
        raise NotImplementedError()

    def parseFile(self, filename):
        """Parse a file and set the _x, _y, _dx, _dy and _meta variables.

        This wipes out the currently loaded data and selected bank number.

        Arguments
        filename    --  The name of the file to parse

        Raises IOError if the file cannot be read
        Raises ParseError if the file cannot be parsed

        """
        infile = open(filename, 'r')
        self._banks = []
        self._meta = {}
        filestring = infile.read()
        self.parseString(filestring)
        infile.close()
        self._meta["filename"] = filename

        if len(self._banks) < 1:
            raise ParseError("There are no data in the banks")

        self.selectBank(0)
        return

    def getNumBanks(self):
        """Get the number of banks read by the parser."""
        return len(self._banks)

    def selectBank(self, index):
        """Select which bank to use.

        This method should only be called after the data has been parsed.  The
        chosen bank number is not persistent, and so must be re-selected if the
        parser is used to parse more data. This uses python list notation, so
        index -n returns the nth bank from the end.

        Arguments:
        index  --  index of bank (integer, starting at 0).

        Raises IndexError if requesting a bank that does not exist

        """
        if index is None:
            index = self._meta.get("bank", 0)

        numbanks = self.getNumBanks()
        if index > numbanks:
            raise IndexError("Bank index out of range")

        if index < 0:
            index += numbanks

        if index < 0:
            raise IndexError("Bank index out of range")

        self._meta["bank"] = index
        self._meta["nbanks"] = numbanks
        self._x, self._y, self._dx, self._dy = self._banks[index]
        return

    def getData(self, index = None):
        """Get the data.

        This method should only be called after the data has been parsed.  The
        chosen bank number is not persistent, and so must be re-selected if the
        parser is used to parse more data. This uses python list notation, so
        index -n returns the nth bank from the end.

        Arguments:
        index  --   index of bank (integer, starting at 0, default None). If
                    index is None then the currently selected bank is used.

        This returns (x, y, dx, dy) tuple for the bank. dx is 0 if it cannot
        be determined from the data format.

        """
        self.selectBank(index)

        return self._x, self._y, self._dx, self._dy

    def getMetaData(self):
        """Get the parsed metadata."""
        return self._meta

# End of ProfileParser
