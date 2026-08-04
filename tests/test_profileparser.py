import re
from pathlib import Path

import pytest

from diffpy.srfit.exceptions import ParseError
from diffpy.srfit.fitbase.profileparser import ProfileParser

# UC1: User loads file with all x, y, dx, dy columns in that format
# expected: x, y, dx, dy, and metadata are all read correctly
# UC2: User loads file with x, y, dy columns in that format (dx is missing)
# expected: x, y, dy, and metadata are all read correctly
# UC3: User loads file with x, y columns in that format (dx and dy are missing)
# expected: x, y, and metadata are all read correctly
# UC4: User loads file with x, dx, y, dy columns in that format and specifies
# column_format
# expected: x, y, dx, dy, and metadata are all read correctly
# UC5: User loads file with dy and dx values containing NaN and inf values
# expected: x, y, and metadata are all read correctly and dx and dy are set to
# 0 for all values

# UC6: User loads file with only one column
# expected: ParseError is raised
# UC7: User loads file with 5 columns
# expected: ParseError is raised
# UC8: User loads file with x, y, and dy but specifies column_format with 4
# columns
# expected: ParseError is raised
# UC9: User loads file with x, y, dx, and dy but specifies column_format with 5
# columns
# expected: ParseError is raised
# UC10: User loads file with x, y, dx, and dy but specifies column_format with
# 3 columns
# expected: ParseError is raised
# UC11: User loads file with x, y, dx, and dy but specifies column_format with
# duplicate values
# expected: ParseError is raised
# UC13: User loads file with only dx, dy columns specified in column_format
# (neither x nor y is present)
# expected: ParseError is raised
# UC14: User loads a file with a format-specific parser (parse_string
# overridden, as PDFParser does) and also specifies column_format
# expected: ParseError is raised
# UC15: User loads a file with a format-specific parser whose parse_string
# does not populate any banks
# expected: ParseError is raised

EXPECTED_META = {
    "wavelength": 0.1,
    "dataformat": "QA",
    "inputfile": "input.iq",
    "backgroundfile": "backgroundfile.iq",
    "mode": "xray",
    "bgscale": 1.0,
    "composition": "TiSe2",
    "outputtype": "gr",
    "qmaxinst": 25.0,
    "qmin": 0.1,
    "qmax": 25.0,
    "rmax": 140.0,
    "rmin": 0.0,
    "rstep": 0.01,
    "rpoly": 0.7,
    "inputdir": "/my/data/dir",
    "savedir": "/my/save/dir",
    "backgroundfilefull": "/my/data/dir/backgroundfile.iq",
    "nbanks": 1,
    "bank": 0,
}


@pytest.mark.parametrize(
    "input_file, column_order, expected_x, "
    "expected_y, expected_dx, expected_dy",
    [
        # UC1: 4-column file (x, y, dx, dy) — all columns present
        # expected: x, y, dx, dy, and metadata are all read correctly
        (
            Path("four_col.gr"),
            None,
            [1.0, 1.1, 1.2],
            [2.0, 2.1, 2.2],
            [0.1, 0.3, 0.5],
            [0.2, 0.4, 0.6],
        ),
        # UC2: 3-column file (x, y, dy) — dx is missing
        # expected: x, y, dy, and metadata are all read correctly
        (
            Path("three_col.dat"),
            None,
            [1.0, 1.1, 1.2],
            [2.0, 2.1, 2.2],
            [0.0, 0.0, 0.0],
            [0.2, 0.4, 0.6],
        ),
        # UC3: 2-column file (x, y) — dx and dy are missing
        # expected: x, y, and metadata are all read correctly
        (
            Path("two_col.txt"),
            None,
            [1.0, 1.1, 1.2],
            [2.0, 2.1, 2.2],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ),
        # UC4: 4-column file in (x, dx, y, dy) order with explicit
        # column_format
        # expected: x, y, dx, dy, and metadata are all read correctly
        (
            Path("four_col_reordered.txt"),
            ("x", "dx", "y", "dy"),
            [1.0, 1.1, 1.2],
            [2.0, 2.1, 2.2],
            [0.1, 0.3, 0.5],
            [0.2, 0.4, 0.6],
        ),
        # UC5: 4-column file where dx/dy contain NaN and inf values
        # expected: x, y, and metadata are read correctly; dx and dy
        # are set to 0
        (
            Path("four_col_nan_inf.gr"),
            None,
            [1.0, 1.1, 1.2],
            [2.0, 2.1, 2.2],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ),
    ],
)
def test_parse_file(
    parser_datafiles,
    input_file,
    column_order,
    expected_x,
    expected_y,
    expected_dx,
    expected_dy,
):
    parser = ProfileParser()
    parser.parse_file(parser_datafiles / input_file, column_order)
    actual_x = parser._x.tolist()
    actual_y = parser._y.tolist()
    actual_dx = parser._dx.tolist()
    actual_dy = parser._dy.tolist()
    actual_metadata = parser._meta
    actual_metadata["filename"] = actual_metadata["filename"].split("/")[-1]

    EXPECTED_META["filename"] = str(input_file).split("/")[-1]
    assert actual_x == expected_x
    assert actual_y == expected_y
    assert actual_dx == expected_dx
    assert actual_dy == expected_dy
    assert actual_metadata == EXPECTED_META


@pytest.mark.parametrize(
    "input_file, column_order, msg",
    [
        # UC6: Only one column — cannot form x/y pair
        # expected: ParseError is raised
        (
            "one_col.gr",
            None,
            "Data block must have at least two columns (x, y).",
        ),
        # UC7: Five columns — ambiguous, no mapping defined
        # expected: ParseError is raised
        ("five_col.gr", None, "Expected 2 to 4 columns but found 5."),
        # UC8: 3-column file but column_format expects 4 columns
        # expected: ParseError is raised
        (
            "three_col.dat",
            ("x", "y", "dx", "dy"),
            "column_format has 4 labels but file contains 3 columns.",
        ),
        # UC9: 4-column file but column_format expects 5 columns
        # expected: ParseError is raised
        (
            "four_col.gr",
            ("x", "y", "dx", "dy", "extra"),
            "column_format has 5 labels but file contains 4 columns.",
        ),
        # UC10: 4-column file but column_format expects only 3 columns
        # expected: ParseError is raised
        (
            "four_col.gr",
            ("x", "y", "dy"),
            "column_format has 3 labels but file contains 4 columns.",
        ),
        # UC11: column_format contains duplicate column names
        # expected: ParseError is raised
        (
            "four_col.gr",
            ("x", "x", "dx", "dy"),
            "column_format cannot contain duplicate labels.",
        ),
        # UC12: column_format contains invalid column names
        (
            "four_col.gr",
            ("x", "y", "dx", "invalid"),
            "column_format contains invalid label 'invalid'. "
            "Valid labels are 'x', 'y', 'dx', and 'dy'.",
        ),
        # UC13: column_format specifies dx, dy but neither x nor y
        # expected: ParseError is raised
        (
            "two_col.txt",
            ("dx", "dy"),
            "Both 'x' and 'y' columns must be present in the data.",
        ),
    ],
)
def test_parse_file_bad(parser_datafiles, input_file, column_order, msg):
    parser = ProfileParser()
    with pytest.raises(ParseError, match=re.escape(msg)):
        parser.parse_file(parser_datafiles / input_file, column_order)


class _FormatSpecificParser(ProfileParser):
    """A parser with a format-specific parse_string override, mimicking
    PDFParser, used to test parse_file's dispatch logic."""

    def parse_string(self, patstring):
        raise NotImplementedError


class _NoBanksParser(ProfileParser):
    """A parser whose parse_string never populates any banks."""

    def parse_string(self, patstring):
        pass


@pytest.mark.parametrize(
    "parser, column_format, expected_msg",
    [
        # C1: column_format is specified for a parser that overrides
        # parse_string, so it determines the column layout from the file
        # format instead.
        # Expected: ParseError is raised.
        (
            _FormatSpecificParser(),
            ("x", "y", "dx", "dy"),
            "_FormatSpecificParser determines the column layout from the "
            "file format, so 'column_format' is not supported.",
        ),
        # C2: A format-specific parser's parse_string does not populate
        # any banks.
        # Expected: ParseError is raised.
        (
            _NoBanksParser(),
            None,
            "There are no data in the banks",
        ),
    ],
)
def test_parse_file_bad_parser(
    parser_datafiles, parser, column_format, expected_msg
):
    with pytest.raises(ParseError, match=re.escape(expected_msg)):
        parser.parse_file(parser_datafiles / "four_col.gr", column_format)


def test_parse_file_does_not_extract_pdf_metadata(datafile):
    # ProfileParser.parse_file extracts metadata by in the header,
    # using "=" to separate keys and values.
    # Expected: parser.get_metadata() parses treturns a dict
    parser = ProfileParser()
    # SAS data is used because the header contains metadata separated by '='.
    parser.parse_file(datafile("sas_ellipsoid_testdata.txt"))
    actual_metadata = parser.get_metadata()
    expected_metadata = {
        "pythonclass": "EllipsoidModel",
        "scale": 1.0,
        "radius_a": 20.0,
        "radius_b": 400.0,
        "contrast": 3e-06,
        "background": 0.01,
        "filename": str(datafile("sas_ellipsoid_testdata.txt")),
        "nbanks": 1,
        "bank": 0,
    }
    assert actual_metadata == expected_metadata
