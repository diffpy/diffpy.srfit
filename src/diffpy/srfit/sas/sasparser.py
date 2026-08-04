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
"""This module contains parsers for SAS data.

SASParser uses the sas DataLoader class to load data.
"""

__all__ = ["SASParser"]


from diffpy.srfit.exceptions import ParseError
from diffpy.srfit.fitbase.profileparser import ProfileParser


class SASParser(ProfileParser):
    """Class for parsing a sas profile.

    This uses a sas DataLoader to load the data. The DataInfo object it returns
    is held in the metadata under the name "datainfo".

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
            read from the file. This is 0 if the
            uncertainty cannot be read.
        dy
            A numpy array containing the uncertainty read
            from the file. This is 0 if the uncertainty
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
    ----------------
    filename
        The name of the file from which data was parsed. This key
        will not exist if data was not read from file.
    nbanks
        The number of banks parsed.
    bank
        The chosen bank number.


    Metadata
    --------
    datainfo
        The DataInfo object used to do the data parsing.


    These may appear in the metadata dictionary
    """

    _format = "SAS"

    def parse_file(self, filename, column_format=None, **kwargs):
        """Parse a file and set the _x, _y, _dx, _dy and _meta
        variables.

        The sas DataLoader reads the data and the metadata together, so
        this overrides parse_file itself rather than the
        `_parse_metadata` and `_parse_data` hooks.

        This wipes out the currently loaded data and selected bank number.

        Parameters
        ----------
        filename
            The name of the file to parse
        column_format
            Unused. Accepted so that the signature matches
            `ProfileParser.parse_file`.
        kwargs
            Unused. Accepted so that the signature matches
            `ProfileParser.parse_file`.

        Raises
        ----------
        IOError
            if the file cannot be read
        ParseError
            if the file cannot be parsed
        """
        import sasdata.dataloader.loader as sas_dataloader

        Loader = sas_dataloader.Loader
        loader = Loader()

        try:
            data = loader.load(str(filename))
        except RuntimeError as e:
            raise ParseError(e)
        except ValueError as e:
            raise ParseError(e)

        self._banks = []
        self._meta = {}
        self._meta["filename"] = filename
        self._meta["datainfo"] = data

        for data_obj in data:
            self._banks.append(
                [data_obj.x, data_obj.y, data_obj.dx, data_obj.dy]
            )
        # FIXME: Revisit when we refactor the SAS characteristic functions.
        # Why is a list imported but only the first element is taken?
        # Is this desired behavior?
        self.select_bank(0)
        return


# End of class SASParser
