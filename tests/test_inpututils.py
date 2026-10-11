"""Unit tests for inpututils.py."""

import io

import pytest

from diffpy.srfit.util.inpututils import (
    convert_input_to_string,
    get_dict_from_results_file,
)

MULTILINE_TEXT = "a 1.0\nb 2.0\n"
LONG_TEXT = "x" * 80


# The input may be an open file, a file name, or the text itself.
@pytest.mark.parametrize(
    "input_kind, expected_string",
    [
        # C1: An open file-like object.
        # Expected: Its contents.
        ("file_object", MULTILINE_TEXT),
        # C2: The name of an existing file.
        # Expected: The file's contents.
        ("file_name", MULTILINE_TEXT),
        # C3: Text spanning several lines.
        # Expected: The text itself.
        ("multiline_text", MULTILINE_TEXT),
        # C4: A single line of 80 or more characters.
        # Expected: The text itself.
        ("long_text", LONG_TEXT),
    ],
)
def test_convert_input_to_string(tmp_path, input_kind, expected_string):
    file_path = tmp_path / "input.txt"
    file_path.write_text(MULTILINE_TEXT)
    inputs = {
        "file_object": io.StringIO(MULTILINE_TEXT),
        "file_name": str(file_path),
        "multiline_text": MULTILINE_TEXT,
        "long_text": LONG_TEXT,
    }
    actual_string = convert_input_to_string(inputs[input_kind])
    assert actual_string == expected_string


def test_convert_input_to_string_missing_file():
    # C1: A short single line that is not an existing file.
    # Expected: It is treated as a file name, so FileNotFoundError is raised.
    with pytest.raises(FileNotFoundError):
        convert_input_to_string("missing_file.txt")


def test_get_dict_from_results_file(tmp_path):
    # C1: A results file with parameter lines, a separator, a blank line,
    # a line without an uncertainty and a line with a non-numeric value.
    # Expected: Only well-formed "name value +/- uncertainty" lines are read.
    results_file = tmp_path / "results.res"
    results_file.write_text(
        "Variables\n"
        "--------------\n"
        "a 1.5 +/- 0.1\n"
        "\n"
        "b 2 +/- 0.3\n"
        "c 4.0\n"
        "d value +/- 0.1\n"
    )
    actual_results = get_dict_from_results_file(results_file)
    expected_results = {"a": 1.5, "b": 2.0}
    assert actual_results == expected_results
