"""Unit tests for nameutils.py."""

import re

import pytest

from diffpy.srfit.util.nameutils import is_identifier, validate_name


# A valid name starts with a letter or underscore and contains only
# letters, digits and underscores.
@pytest.mark.parametrize(
    "input_name, expected_is_identifier",
    [
        # C1: A single letter.
        # Expected: Valid.
        ("x", True),
        # C2: Letters, digits and underscores, starting with an underscore.
        # Expected: Valid.
        ("_scale_2", True),
        # C3: The name starts with a digit.
        # Expected: Not valid.
        ("2x", False),
        # C4: The name contains a hyphen.
        # Expected: Not valid.
        ("a-b", False),
        # C5: The name contains a space.
        # Expected: Not valid.
        ("a b", False),
        # C6: The name is empty.
        # Expected: Not valid.
        ("", False),
    ],
)
def test_is_identifier(input_name, expected_is_identifier):
    actual_is_identifier = is_identifier(input_name)
    assert actual_is_identifier == expected_is_identifier


@pytest.mark.parametrize(
    "input_name",
    [
        # C1: A valid name.
        # Expected: No error is raised.
        "scale",
        # C2: A valid name with digits and underscores.
        # Expected: No error is raised.
        "_scale_2",
    ],
)
def test_validate_name_valid(input_name):
    actual_result = validate_name(input_name)
    assert actual_result is None


@pytest.mark.parametrize(
    "input_name",
    [
        # C1: The name starts with a digit.
        # Expected: A ValueError explains the problem and the fix.
        "2x",
        # C2: The name contains a hyphen.
        # Expected: A ValueError explains the problem and the fix.
        "a-b",
    ],
)
def test_validate_name_invalid(input_name):
    expected_message = (
        f"Name '{input_name}' is not a valid identifier. "
        "Please use only letters, digits and underscores, "
        "and do not start the name with a digit."
    )
    with pytest.raises(ValueError, match=re.escape(expected_message)):
        validate_name(input_name)
