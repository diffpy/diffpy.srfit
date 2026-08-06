import re
from pathlib import Path

import pytest

from diffpy.srfit.exceptions import ParseError
from diffpy.srfit.fitbase.profileparser import ProfileParser

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
            None,
            [0.2, 0.4, 0.6],
        ),
        # UC3: 2-column file (x, y) — dx and dy are missing
        # expected: x, y, and metadata are all read correctly
        (
            Path("two_col.txt"),
            None,
            [1.0, 1.1, 1.2],
            [2.0, 2.1, 2.2],
            None,
            None,
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
        # are set to None
        (
            Path("four_col_nan_inf.gr"),
            None,
            [1.0, 1.1, 1.2],
            [2.0, 2.1, 2.2],
            None,
            None,
        ),
    ],
)
def test_parse_file(
    as_list,
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
    actual_dx = as_list(parser._dx)
    actual_dy = as_list(parser._dy)
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


# ProfileParser is an extension point. A subclass customizes a format by
# overriding the _parse_metadata and _parse_data hooks; parse_file itself
# is a template method that subclasses are not expected to touch. The
# tests below pin that contract so a future refactor cannot quietly turn
# parse_file back into a concrete parser.


class MetadataOnlyParser(ProfileParser):
    """A parser that customizes only the metadata, as PDFParser does."""

    _format = "metadata-only"

    def _parse_metadata(self, filename):
        return {"instrument": "custom"}


def test_parse_file_uses_metadata_hook(parser_datafiles):
    """Overriding _parse_metadata keeps the inherited column
    handling."""
    # Case: a subclass overrides _parse_metadata to customize the metadata,
    # but does not override _parse_data.
    # Expected: The subclass's metadata parser is used,
    # but the data columns are still read correctly.
    parser = MetadataOnlyParser()
    parser.parse_file(parser_datafiles / "four_col.gr")

    actual_format = parser.get_format()
    expected_format = "metadata-only"
    assert actual_format == expected_format

    # The custom hook replaces the generic header scan entirely, so the
    # name = value pairs the default would have collected are absent.
    actual_metadata = parser.get_metadata()
    assert actual_metadata["instrument"] == "custom"
    assert "wavelength" not in actual_metadata

    actual_x, actual_y, actual_dx, actual_dy = parser.get_data()
    expected_x = [1.0, 1.1, 1.2]
    expected_y = [2.0, 2.1, 2.2]
    expected_dx = [0.1, 0.3, 0.5]
    expected_dy = [0.2, 0.4, 0.6]
    assert actual_x.tolist() == expected_x
    assert actual_y.tolist() == expected_y
    assert actual_dx.tolist() == expected_dx
    assert actual_dy.tolist() == expected_dy


class TwoBankParser(ProfileParser):
    """A parser that customizes only the data, adding a second bank."""

    def _parse_data(self, filename, column_format=None, **kwargs):
        super()._parse_data(filename, column_format, **kwargs)
        input_x, input_y, input_dx, input_dy = self._banks[0]
        self._banks.append([input_x, 2 * input_y, input_dx, input_dy])


def test_parse_file_uses_data_hook(parser_datafiles):
    """Overriding _parse_data may contribute several banks."""
    # Case: a subclass overrides _parse_data to customize the data,
    # but does not override _parse_metadata.
    # Expected: The subclass's data parser is used, and the
    # metadata is still read correctly.
    parser = TwoBankParser()
    parser.parse_file(parser_datafiles / "four_col.gr")

    actual_bank_count = parser.get_num_banks()
    expected_bank_count = 2
    assert actual_bank_count == expected_bank_count

    actual_reported_bank_count = parser.get_metadata()["nbanks"]
    assert actual_reported_bank_count == expected_bank_count

    actual_second_bank_y = parser.get_data(1)[1].tolist()
    expected_second_bank_y = [4.0, 4.2, 4.4]
    assert actual_second_bank_y == expected_second_bank_y


def test_parse_file_deprecated_warns_about_profileparser(parser_datafiles):
    """The shared parseFile reports its own name, not a subclass one.

    parseFile is inherited by every parser, so a message naming a
    specific subclass would misdirect users of the others.
    """
    # Case: parseFile is called on a subclass of ProfileParser.
    # Expected: A DeprecationWarning is raised, and the message mentions
    parser = ProfileParser()
    expected_msg = (
        "'diffpy.srfit.fitbase.profileparser.ProfileParser.parseFile' is "
        "deprecated and will be removed in version 4.0.0. Please use "
        "'diffpy.srfit.fitbase.profileparser.ProfileParser.parse_file' "
        "instead."
    )
    with pytest.warns(
        DeprecationWarning,
        match=re.escape(expected_msg),
    ):
        parser.parseFile(parser_datafiles / "two_col.txt")

    actual_x = parser.get_data()[0].tolist()
    expected_x = [1.0, 1.1, 1.2]
    assert actual_x == expected_x


# parse_file accepts an optional `metadata` dict, merged into the parsed
# metadata after the file header is read. This lets a caller attach
# information that is not present in the file itself (e.g. sample name),
# without having to edit `parser.get_metadata()` by hand after parsing.


@pytest.mark.parametrize(
    "input_metadata, expected_added_metadata",
    [
        # UC13: Supplied keys are not present in the file header
        # Expected: the supplied keys/values are added to the parsed
        # metadata, alongside the keys read from the file
        (
            {"operator": "jdoe", "sample": "LaB6"},
            {"operator": "jdoe", "sample": "LaB6"},
        ),
        # UC14: A supplied key overlaps with one already parsed from the
        # file header
        # Expected: the supplied value overrides the value parsed from
        # the file
        (
            {"wavelength": 99.9},
            {"wavelength": 99.9},
        ),
    ],
)
def test_parse_file_appends_metadata(
    parser_datafiles, input_metadata, expected_added_metadata
):
    parser = ProfileParser()
    parser.parse_file(
        parser_datafiles / "four_col.gr", metadata=input_metadata
    )

    actual_metadata = parser.get_metadata()
    for key, expected_value in expected_added_metadata.items():
        assert actual_metadata[key] == expected_value


@pytest.mark.parametrize(
    "bad_metadata_input, expected_msg",
    [
        # UC15: The user supplies a metadata dict with a key that is not
        # a string.
        # Expected: A ParseError is raised with a message indicating that
        # all keys must be strings.
        (
            {"operator": "jdoe", 42: "LaB6"},
            "Key '42' in the metadata dictionary "
            "is not a string. All keys in the metadata "
            "dictionary must be strings.",
        ),
        # UC16: User supplies metadata not in the form of a dictionary.
        # Expected: A ParseError is raised with a message indicating that
        # the metadata must be a dictionary.
        (
            ["operator", "jdoe"],
            "The metadata argument must be a dictionary. "
            "Received type 'list' instead.",
        ),
    ],
)
def test_parse_file_appends_metadata_bad(
    parser_datafiles, bad_metadata_input, expected_msg
):
    parser = ProfileParser()
    with pytest.raises(ParseError, match=re.escape(expected_msg)):
        parser.parse_file(
            parser_datafiles / "four_col.gr", metadata=bad_metadata_input
        )


@pytest.mark.parametrize(
    "input_metadata, reserved_key, expected_value, expected_msg",
    [
        # UC17: Supplied metadata overlaps with the "filename" key, which
        # parse_file normally sets itself from the file argument.
        # Expected: A UserWarning is raised naming the reserved key, and
        # the supplied value overrides the automatically parsed one.
        (
            {"filename": "spoofed.dat"},
            "filename",
            "spoofed.dat",
            "'filename' is a reserved metadata key normally set by "
            "parse_file. The supplied value will override it.",
        ),
        # UC18: Supplied metadata overlaps with the "bank" key, which
        # parse_file normally sets from select_bank.
        # Expected: A UserWarning is raised naming the reserved key, and
        # the supplied value overrides the automatically parsed one.
        (
            {"bank": 5},
            "bank",
            5,
            "'bank' is a reserved metadata key normally set by "
            "parse_file. The supplied value will override it.",
        ),
        # UC19: Supplied metadata overlaps with the "nbanks" key, which
        # parse_file normally sets from the number of banks read.
        # Expected: A UserWarning is raised naming the reserved key, and
        # the supplied value overrides the automatically parsed one.
        (
            {"nbanks": 99},
            "nbanks",
            99,
            "'nbanks' is a reserved metadata key normally set by "
            "parse_file. The supplied value will override it.",
        ),
    ],
)
def test_parse_file_appends_metadata_warns_on_reserved_key(
    parser_datafiles,
    input_metadata,
    reserved_key,
    expected_value,
    expected_msg,
):
    parser = ProfileParser()
    with pytest.warns(UserWarning, match=re.escape(expected_msg)):
        parser.parse_file(
            parser_datafiles / "four_col.gr", metadata=input_metadata
        )

    actual_metadata = parser.get_metadata()
    actual_value = actual_metadata[reserved_key]
    assert actual_value == expected_value


def test_parse_file_appends_metadata_does_not_mutate_input(parser_datafiles):
    # Case: the caller's metadata dict is mutated after parse_file returns.
    # Expected: the parser's own metadata is unaffected, so parse_file must
    # have copied the dict rather than stored a reference to it.
    parser = ProfileParser()
    input_metadata = {"operator": "jdoe"}
    parser.parse_file(
        parser_datafiles / "four_col.gr", metadata=input_metadata
    )

    input_metadata["operator"] = "mutated"
    input_metadata["sample"] = "added after parsing"

    actual_metadata = parser.get_metadata()
    assert actual_metadata["operator"] == "jdoe"
    assert "sample" not in actual_metadata


def test_parse_file_usecols_selects_columns(parser_datafiles):
    """A file wider than four columns is read by selecting columns.

    column_format has to label every column that is loaded, so a wider
    file is narrowed first with the load_data usecols argument.
    """
    # Case: A file with five columns is loaded, but only the first two are
    # used.
    # Expected: The first two columns are read as x and y, and dx and dy are
    # None.
    parser = ProfileParser()
    parser.parse_file(
        parser_datafiles / "five_col.gr", ("x", "y"), usecols=(0, 1)
    )

    actual_x, actual_y, actual_dx, actual_dy = parser.get_data()
    expected_x = [1.0, 1.1, 1.2]
    expected_y = [2.0, 2.1, 2.2]
    assert actual_x.tolist() == expected_x
    assert actual_y.tolist() == expected_y
    assert actual_dx is None
    assert actual_dy is None
