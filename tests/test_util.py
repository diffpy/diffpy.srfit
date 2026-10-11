"""Unit tests for the functions in diffpy.srfit.util."""

import pytest

from diffpy.srfit.util import sort_key_for_numeric_string


# The key splits a string into text and integer parts so that numbers sort
# by value rather than character by character.
@pytest.mark.parametrize(
    "input_string, expected_key",
    [
        # C1: Text, a number, then more text.
        # Expected: The number becomes an integer between the text parts.
        ("a12b", ("a", 12, "b")),
        # C2: The string starts with a number.
        # Expected: The key starts with an empty text part.
        ("12a", ("", 12, "a")),
        # C3: The string has no digits.
        # Expected: The key holds the string unchanged.
        ("abc", ("abc",)),
    ],
)
def test_sort_key_for_numeric_string(input_string, expected_key):
    actual_key = sort_key_for_numeric_string(input_string)
    assert actual_key == expected_key


def test_sort_key_for_numeric_string_orders_numbers_by_value():
    # C1: Names whose numbers have different numbers of digits are sorted.
    # Expected: "a2" comes before "a10".
    actual_order = sorted(
        ["a10", "a2", "b1", "a1"], key=sort_key_for_numeric_string
    )
    expected_order = ["a1", "a2", "a10", "b1"]
    assert actual_order == expected_order
