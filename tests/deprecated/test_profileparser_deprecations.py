"""Tests for the deprecated camelCase names in
diffpy.srfit.fitbase.profileparser.

Each test calls the deprecated name, checks that it warns, and checks
that it does the same as the new name.

Remove this file in version 4.0.0, together with the deprecated names.
"""

import re

import pytest

from diffpy.srfit import ProfileParser


def expected_message(base, old_name, new_name):
    return re.escape(
        f"'{base}.{old_name}' is deprecated and will be removed in version "
        f"4.0.0. Please use '{base}.{new_name}' instead."
    )


PROFILE_PARSER = "diffpy.srfit.fitbase.profileparser.ProfileParser"


class _SubclassParser(ProfileParser):
    """A parser that inherits every method from ProfileParser."""


def _make_parser(parser_class, method_name, parser_datafiles):
    """Return a parser to call the method on.

    parse_file is called on a new parser. The other methods are called
    on a parser that has parsed a four-column file.
    """
    parser = parser_class()
    if method_name != "parse_file":
        parser.parse_file(parser_datafiles / "four_col.gr")
    return parser


def _make_call_args(method_name, parser_datafiles):
    """Return fresh arguments for a ProfileParser method."""
    if method_name == "parse_file":
        return (parser_datafiles / "two_col.txt",)
    if method_name == "select_bank":
        return (0,)
    return ()


def _describe_value(value, as_list):
    if isinstance(value, tuple):
        return tuple(as_list(array) for array in value)
    return value


def _describe_parser(parser, as_list):
    return (
        parser.get_format(),
        parser.get_num_banks(),
        _describe_value(parser.get_data(), as_list),
        dict(parser.get_metadata()),
    )


# ProfileParser: each deprecated method does the same as its new name, and
# its warning names ProfileParser even when it is called on a subclass,
# since the method is shared by every parser. See test_profileparser.py.
@pytest.mark.parametrize(
    "input_parser_class",
    [
        # C1: The parser is a ProfileParser.
        # Expected: The warning names ProfileParser.
        ProfileParser,
        # C2: The parser is a subclass that inherits the method.
        # Expected: The warning still names ProfileParser.
        _SubclassParser,
    ],
)
@pytest.mark.parametrize(
    "input_old_name, input_new_name",
    [
        # C1: parseFile parses a two-column file.
        # Expected: Warns and gives the same data and metadata as
        # parse_file.
        ("parseFile", "parse_file"),
        # C2: getFormat gets the format of the parser.
        # Expected: Warns and returns the same format as get_format.
        ("getFormat", "get_format"),
        # C3: getNumBanks gets the number of banks read.
        # Expected: Warns and returns the same number as get_num_banks.
        ("getNumBanks", "get_num_banks"),
        # C4: selectBank selects bank 0.
        # Expected: Warns and selects the same bank as select_bank.
        ("selectBank", "select_bank"),
        # C5: getData gets the data of the selected bank.
        # Expected: Warns and returns the same data as get_data.
        ("getData", "get_data"),
        # C6: getMetaData gets the parsed metadata.
        # Expected: Warns and returns the same metadata as get_metadata.
        ("getMetaData", "get_metadata"),
    ],
)
def test_profileparser_deprecations(
    as_list,
    parser_datafiles,
    input_parser_class,
    input_old_name,
    input_new_name,
):
    old_parser = _make_parser(
        input_parser_class, input_new_name, parser_datafiles
    )
    old_args = _make_call_args(input_new_name, parser_datafiles)
    with pytest.warns(
        DeprecationWarning,
        match=expected_message(PROFILE_PARSER, input_old_name, input_new_name),
    ):
        old_return_value = getattr(old_parser, input_old_name)(*old_args)
    new_parser = _make_parser(
        input_parser_class, input_new_name, parser_datafiles
    )
    new_args = _make_call_args(input_new_name, parser_datafiles)
    new_return_value = getattr(new_parser, input_new_name)(*new_args)
    actual_old_result = (
        _describe_value(old_return_value, as_list),
        _describe_parser(old_parser, as_list),
    )
    actual_new_result = (
        _describe_value(new_return_value, as_list),
        _describe_parser(new_parser, as_list),
    )
    assert actual_old_result == actual_new_result
