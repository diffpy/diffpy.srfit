"""Unit tests for the functions in diffpy.srfit.util."""

import pytest

from diffpy.srfit.util import (
    format_value_with_uncertainty,
    sort_key_for_numeric_string,
)


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


# Test formatting a value with its uncertainty in parenthesis notation. The
# uncertainty is rounded to two significant digits and the value to the
# same decimal place, as in "3.1997(21)" for 3.1997 +/- 0.0021.
@pytest.mark.parametrize(
    "input_value, input_uncertainty, expected_formatted_value",
    [
        # C1: A user formats a fitted peak position of 3.19965644 with an
        # uncertainty of 0.0021228.
        # Expected: The uncertainty becomes 21 and the value is rounded to
        # four decimal places.
        (3.19965644, 0.0021228, "3.1997(21)"),
        # C2: A user formats a value whose uncertainty is greater than one.
        # Expected: The digits in parentheses are in the same decimal places
        # as the last digits of the value.
        # 1. An uncertainty of 1.234 is rounded to one decimal place.
        # 2. An uncertainty of 12.3 is rounded to a whole number.
        (12.345, 1.234, "12.3(12)"),
        (1234.6, 12.3, "1235(12)"),
        # C3: A user formats a negative slope of -0.298964 with an
        # uncertainty of 0.00203.
        # Expected: The value keeps its sign, and the uncertainty becomes
        # 20 rather than 2.
        (-0.298964, 0.00203, "-0.2990(20)"),
        # C4: A user formats a value whose uncertainty of 0.0996 rounds up
        # to the next power of ten.
        # Expected: The uncertainty becomes 10, and the value is rounded to
        # two decimal places rather than three.
        (1.2345, 0.0996, "1.23(10)"),
        # C5: A user formats a value whose uncertainty can't be used.
        # Expected: The value is given alone with six significant digits.
        # 1. The uncertainty is zero.
        # 2. The uncertainty is negative.
        # 3. The uncertainty is not a number, as when the covariance
        # matrix can't be computed.
        # 4. The uncertainty is infinite.
        (5.012345678, 0.0, "5.01235"),
        (5.012345678, -0.1, "5.01235"),
        (5.012345678, float("nan"), "5.01235"),
        (5.012345678, float("inf"), "5.01235"),
    ],
)
def test_format_value_with_uncertainty(
    input_value, input_uncertainty, expected_formatted_value
):
    actual_formatted_value = format_value_with_uncertainty(
        input_value, input_uncertainty
    )
    assert actual_formatted_value == expected_formatted_value
