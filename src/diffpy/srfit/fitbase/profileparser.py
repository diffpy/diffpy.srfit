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

ProfileParser is a base class for parsing data. It can interact with a
Profile object to automatically set the Profile's data and metadata.
Each specific file format must be encapsulated in a ProfileParser
subclass.

See the class documentation for more information.
"""

import numpy as np

from diffpy.srfit.exceptions import ParseError
from diffpy.utils._deprecator import build_deprecation_message, deprecated
from diffpy.utils.parsers import load_data

removal_verison = "4.0.0"
pp_base = "diffpy.srfit.fitbase.profileparser.ProfileParser"

parseFile_dep_msg = build_deprecation_message(
    pp_base,
    "parseFile",
    "parse_file",
    removal_verison,
)

getFormat_dep_msg = build_deprecation_message(
    pp_base,
    "getFormat",
    "get_format",
    removal_verison,
)

getNumBanks_dep_msg = build_deprecation_message(
    pp_base,
    "getNumBanks",
    "get_num_banks",
    removal_verison,
)

selectBank_dep_msg = build_deprecation_message(
    pp_base,
    "selectBank",
    "select_bank",
    removal_verison,
)

getData_dep_msg = build_deprecation_message(
    pp_base,
    "getData",
    "get_data",
    removal_verison,
)

getMetaData_dep_msg = build_deprecation_message(
    pp_base,
    "getMetaData",
    "get_metadata",
    removal_verison,
)


class ProfileParser(object):
    """Class for parsing data from a or string.

    Attributes
    ----------
    _format : str, optional
        The name of the data format that this parses (string, default
        `""`). The format string is a unique identifier for the data
        format handled by the parser.
    _banks : list of tuples
        The data from each bank. Each bank contains a (x, y, dx,
        dy)
        tuple:
        x : np.ndarray
            The independent variable read from the file.
        y : np.ndarray
            The dependent variable (profile) read
            from the file.
        dx : np.ndarray
            The uncertainties associated with x
            read from the file. This is None if the
            uncertainty cannot be read.
        dy : np.ndarray
            The uncertainties associated with y
            read from the file. This is None if the
            uncertainty cannot be read.
    _x : np.ndarray
        Independent variable from the chosen bank
    _y : np.ndarray
        Profile from the chosen bank
    _dx : np.ndarray
        Uncertainty in independent variable from the chosen bank
    _dy : np.ndarray
        Uncertainty in profile from the chosen bank
    _meta : dict
        A dictionary containing metadata read from the file.


    General Metadata
    ----------------
    filename : str or Path
        The name of the file from which data was parsed. This key
        will not exist if data was not read from file.
    nbanks : int
        The number of banks parsed.
    bank : int
        The chosen bank number.
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

    def get_format(self):
        """Get the format string.

        Returns
        -------
        str
            The unique identifier for the data format handled by this
            parser.
        """
        return self._format

    @deprecated(getFormat_dep_msg)
    def getFormat(self):
        """This function is deprecated and will be removed in version
        4.0.0.

        Please use diffpy.srfit.fitbase.ProfileParser.get_format
        instead.
        """
        return self.get_format()

    @deprecated(parseFile_dep_msg)
    def parseFile(self, filename):
        """This function is deprecated and will be removed in version
        4.0.0.

        Please use diffpy.srfit.fitbase.ProfileParser.parse_file
        instead.
        """
        return self.parse_file(filename)

    def parse_file(self, filename, column_format=None, **kwargs):
        """Parse a data file to extract data and metadata, with
        automatic handling of uncertainties.

        This is a template method. Subclasses customize a format by
        overriding the `_parse_metadata` and `_parse_data` hooks rather
        than this method. `PDFParser` in the
        `diffpy.srfit.pdf.pdfparser` module is a worked example: it
        overrides `_parse_metadata` to read PDFgetX/PDFgetN headers and
        inherits the data block handling unchanged.

        The default `_parse_data` reads a single bank:

        - For files with 2 columns: assumes (x, y) and sets dx, dy to None.
        - For files with 3 columns: assumes (x, y, dy) and sets dx to None.
        - For files with 4 columns: assumes (x, y, dx, dy).
        - For other cases: `column_format` must be explicitly specified.

        Uncertainty columns (dx, dy) are only considered valid if all values
        are positive and not NaN/Inf. Otherwise they are set to None.

        This wipes out the currently loaded data and selected bank number.

        Parameters
        ----------
        filename : str or Path
            The name of the file to parse.
        column_format : tuple of str, optional
            The order in which columns appear in the file.
            If None, the format is auto-detected based on the
            number of columns.

            Valid labels: `"x"`, `"y"`, `"dx"`, `"dy"`

            Examples:

            - `("x", "y")`
            - `("x", "y", "dy")`
            - `("x", "y", "dx", "dy")`
            - `("x", "dx", "y", "dy")`

        kwargs
            The keyword arguments passed on to
            `diffpy.utils.parsers.load_data`, such as `usecols`,
            `delimiter`, `comments` and `minrows`. Use `usecols` to
            select four columns out of a wider file, then label them
            with `column_format`.

        Raises
        ------
        ParseError
            If parsing fails or ambiguity detected.
        """
        self._banks = []
        self._meta = {}
        self._meta.update(self._parse_metadata(filename))
        self._parse_data(filename, column_format, **kwargs)
        self._meta["filename"] = str(filename)
        if len(self._banks) < 1:
            raise ParseError("There are no data in the banks")
        self.select_bank(0)

    def _parse_metadata(self, filename):
        """Return the metadata read from the header of a file.

        Override this hook to parse a format whose header is not a plain
        list of ``name = value`` pairs.
        """
        return load_data(filename, headers=True)

    def _parse_data(self, filename, column_format=None, **kwargs):
        """Append the banks read from the data block of a file.

        Override this hook to parse a format whose data block is not a
        plain matrix of columns, or one that holds several banks.
        """
        data = load_data(filename, **kwargs)
        if data.size == 0 or data.ndim == 1:
            raise ParseError(
                "Data block must have at least two columns (x, y)."
            )
        column_format = self._detect_column_format(data, column_format)
        columns = self._map_column_labels_to_data(data, column_format)
        self._banks.append(
            [
                columns["x"],
                columns["y"],
                self._validate_uncertainty(columns.get("dx")),
                self._validate_uncertainty(columns.get("dy")),
            ]
        )

    def _detect_column_format(self, data, column_format):
        """Auto-detect or validate column format."""
        num_cols = data.shape[1]

        if column_format is None:
            if num_cols == 2:
                column_format = ("x", "y")
            elif num_cols == 3:
                column_format = ("x", "y", "dy")
            elif num_cols == 4:
                column_format = ("x", "y", "dx", "dy")
            else:
                raise ParseError(
                    f"Expected 2 to 4 columns but found {num_cols}."
                )
        if len(column_format) != num_cols:
            raise ParseError(
                f"column_format has {len(column_format)} "
                f"labels but file contains {num_cols} columns."
            )
        if len(set(column_format)) != len(column_format):
            raise ParseError("column_format cannot contain duplicate labels.")
        for label in column_format:
            if label not in {"x", "y", "dx", "dy"}:
                raise ParseError(
                    f"column_format contains invalid label '{label}'. "
                    "Valid labels are 'x', 'y', 'dx', and 'dy'."
                )
        return column_format

    def _map_column_labels_to_data(self, data, column_format):
        """Map numeric data to columns by label."""
        columns = {}
        for i, label in enumerate(column_format):
            columns[label] = data[:, i]

        if "x" not in columns or "y" not in columns:
            raise ParseError(
                "Both 'x' and 'y' columns must be present in the data."
            )

        return columns

    @staticmethod
    def _validate_uncertainty(data):
        """Return the uncertainty data if valid, otherwise None."""
        if data is None or not np.all(np.isfinite(data)) or np.any(data <= 0):
            return None
        return data

    def get_num_banks(self):
        """Get the number of banks read by the parser.

        Returns
        -------
        int
            The number of banks read by the parser.
        """
        return len(self._banks)

    @deprecated(getNumBanks_dep_msg)
    def getNumBanks(self):
        """This function is deprecated and will be removed in version
        4.0.0.

        Please use diffpy.srfit.fitbase.ProfileParser.get_num_banks
        instead.
        """
        return self.get_num_banks()

    def select_bank(self, index):
        """Select which bank to use.

        This method should only be called after the data has been parsed.  The
        chosen bank number is not persistent, and so must be re-selected if the
        parser is used to parse more data. This uses python list notation, so
        index -n returns the nth bank from the end.

        Parameters
        ----------
        index
            index of bank (integer, starting at 0).

        Raises
        ----------
        IndexError
            if requesting a bank that does not exist
        """
        if index is None:
            index = self._meta.get("bank", 0)

        numbanks = self.get_num_banks()
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

    @deprecated(selectBank_dep_msg)
    def selectBank(self, index):
        """This function is deprecated and will be removed in version
        4.0.0.

        Please use diffpy.srfit.fitbase.ProfileParser.select_bank
        instead.
        """
        self.select_bank(index)
        return

    def get_data(self, index=None):
        """Get the data.

        This method should only be called after the data has been parsed.  The
        chosen bank number is not persistent, and so must be re-selected if the
        parser is used to parse more data. This uses python list notation, so
        index -n returns the nth bank from the end.

        Parameters
        ----------
        index
            index of bank (integer, starting at 0, default None). If
            index is None then the currently selected bank is used.

        Returns
        ----------
        This returns (x, y, dx, dy) tuple for the bank. dx is None if it
        cannot be determined from the data format.
        """
        self.select_bank(index)

        return self._x, self._y, self._dx, self._dy

    @deprecated(getData_dep_msg)
    def getData(self, index=None):
        """This function is deprecated and will be removed in version
        4.0.0.

        Please use diffpy.srfit.fitbase.ProfileParser.get_data instead.
        """
        return self.get_data(index)

    def get_metadata(self):
        """Get the parsed metadata.

        Returns
        -------
        dict
            A dictionary containing metadata read from the file.
        """
        return self._meta

    @deprecated(getMetaData_dep_msg)
    def getMetaData(self):
        """This function is deprecated and will be removed in version
        4.0.0.

        Please use diffpy.srfit.fitbase.ProfileParser.get_metadata
        instead.
        """
        return self._meta


# End of ProfileParser
