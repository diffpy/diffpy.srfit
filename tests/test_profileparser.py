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
    ],
)
def test_parse_file_bad(parser_datafiles, input_file, column_order, msg):
    parser = ProfileParser()
    with pytest.raises(ParseError, match=re.escape(msg)):
        parser.parse_file(parser_datafiles / input_file, column_order)
