#!/usr/bin/env python
########################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2008 Trustees of the Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
########################################################################
"""This module contains parsers for SAS data.

SASParser uses the sas DataLoader class to load data.

"""
__all__ = ["SASParser"]

import numpy

from diffpy.srfit.fitbase.profileparser import ProfileParser, ParseError

class SASParser(ProfileParser):
    """Class for parsing a sas profile.

    This uses a sas DataLoader to load the data. The DataInfo object it returns
    is held in the metadata under the name "datainfo".

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
                                read from the file. This is 0 if the uncertainty
                                cannot be read.
                    dy      --  A numpy array containing the uncertainty read
                                from the file. This is 0 if the uncertainty
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

    Metadata - These may appear in the metadata dictionary

    datainfo    --  The DataInfo object used to do the data parsing.

    """

    _format = "SAS"

    def parseFile(self, filename):
        """Parse a file and set the _x, _y, _dx, _dy and _meta variables.

        This wipes out the currently loaded data and selected bank number.

        Arguments
        filename    --  The name of the file to parse

        Raises IOError if the file cannot be read
        Raises ParseError if the file cannot be parsed

        """

        Loader = _import_sans_Loader()
        loader = Loader()

        try:
            data = loader.load(filename)
        except (RuntimeError, ValueError), e:
            raise ParseError(e)

        self._banks = []
        self._meta = {}
        self._meta["filename"] = filename
        self._meta["datainfo"] = data

        self._banks.append([data.x, data.y, data.dx, data.dy])
        self.selectBank(0)
        return

    def parseString(self, patstring):
        """Parse a string and set the _x, _y, _dx, _dy and _meta variables.

        When _dx or _dy cannot be obtained in the data format it is set to 0.

        This wipes out the currently loaded data and selected bank number.

        Arguments
        patstring   --  A string containing the pattern

        Raises ParseError if the string cannot be parsed

        """
        # This calls on parseFile, as that is how the sans data loader works.
        import tempfile
        fh, fn = tempfile.mkstemp()
        outfile = file(fn, 'w')
        fn.write(patstring)
        outfile.close()
        self.parseFile(fn)

        del self._metadata["filename"]

        # Close the temporary file and delete it
        import os
        os.close(fh)
        os.remove(fn)
        return


# End of class SASParser

# Local Helpers --------------------------------------------------------------

def _import_sans_Loader():
    """Return the Loader class from SANS DataLoader package.

    This is a workaround for a recent change in dataloader package name.
    """
    global _the_Loader
    if _the_Loader is not None:  return _the_Loader
    # first try to import from the latest interface:
    try:
        from sans.dataloader.loader import Loader as _the_Loader
        return _import_sans_Loader()
    except ImportError:
        pass
    # try to import from the old interface:
    try:
        from DataLoader.loader import Loader as _the_Loader
        return _import_sans_Loader()
    except ImportError:
        pass
    # finally use the latest interface again and let it raise ImportError
    from sans.dataloader.loader import Loader as _the_Loader
    return _import_sans_Loader()
_the_Loader = None
