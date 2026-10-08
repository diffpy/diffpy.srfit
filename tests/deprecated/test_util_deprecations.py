"""Tests for the deprecated camelCase names in diffpy.srfit.util.

Each test repeats the cases of the matching test for the new name, but
calls the deprecated name and also checks that it warns.

Remove this file in version 4.0.0, together with the deprecated names.
"""

import io
import re

import pytest

from diffpy.srfit.util import sortKeyForNumericString
from diffpy.srfit.util.inpututils import inputToString
from diffpy.srfit.util.nameutils import isIdentifier, validateName
from diffpy.srfit.util.observable import Observable

MULTILINE_TEXT = "a 1.0\nb 2.0\n"
LONG_TEXT = "x" * 80


def expected_message(base, old_name, new_name):
    return re.escape(
        f"'{base}.{old_name}' is deprecated and will be removed in version "
        f"4.0.0. Please use '{base}.{new_name}' instead."
    )


OBSERVABLE = "diffpy.srfit.util.observable.Observable"
TAG_MANAGER = "diffpy.srfit.util.tagmanager.TagManager"
NAMEUTILS = "diffpy.srfit.util.nameutils"


class _Recorder:
    """Observer that records the arguments it is notified with."""

    def __init__(self):
        self.calls = []

    def record(self, semaphores):
        self.calls.append(semaphores)


# Observable: see test_observable.py::test_add_and_remove_observer.
def test_addObserver_hasObserver_removeObserver():
    # C1: An observer is added, then removed.
    # Expected: Each call warns. hasObserver is True after adding and False
    # after removing.
    observable = Observable()
    recorder = _Recorder()
    with pytest.warns(
        DeprecationWarning,
        match=expected_message(OBSERVABLE, "addObserver", "add_observer"),
    ):
        observable.addObserver(recorder.record)
    with pytest.warns(
        DeprecationWarning,
        match=expected_message(OBSERVABLE, "hasObserver", "has_observer"),
    ):
        actual_has_observer_after_add = observable.hasObserver(recorder.record)
    expected_has_observer_after_add = True
    assert actual_has_observer_after_add == expected_has_observer_after_add
    with pytest.warns(
        DeprecationWarning,
        match=expected_message(
            OBSERVABLE, "removeObserver", "remove_observer"
        ),
    ):
        observable.removeObserver(recorder.record)
    actual_has_observer_after_remove = observable.has_observer(recorder.record)
    expected_has_observer_after_remove = False
    assert (
        actual_has_observer_after_remove == expected_has_observer_after_remove
    )


# TagManager: see test_tagmanager.py::test_has_tags.
@pytest.mark.parametrize(
    "input_tags, input_silent, expected_has_tags",
    [
        # C1: The object has the one given tag.
        # Expected: Warns and returns True.
        (("3",), False, True),
        # C2: The object has both given tags.
        # Expected: Warns and returns True.
        (("3", "number"), False, True),
        # C3: The object has only one of the given tags.
        # Expected: Warns and returns False.
        (("3", "4"), False, False),
        # C4: A silent manager is asked about an unknown tag.
        # Expected: Warns and returns False.
        (("fail",), True, False),
    ],
)
def test_hasTags(tag_manager, input_tags, input_silent, expected_has_tags):
    tag_manager.silent = input_silent
    with pytest.warns(
        DeprecationWarning,
        match=expected_message(TAG_MANAGER, "hasTags", "has_tags"),
    ):
        actual_has_tags = tag_manager.hasTags(3, *input_tags)
    assert actual_has_tags == expected_has_tags


# TagManager: see test_tagmanager.py::test_has_tags_unknown_tag_bad.
def test_hasTags_unknown_tag_bad(tag_manager):
    # C1: hasTags(3, "fail") is called, but no object was ever tagged
    # "fail", and the manager is not silent.
    # Expected: Warns and raises a KeyError.
    with (
        pytest.warns(
            DeprecationWarning,
            match=expected_message(TAG_MANAGER, "hasTags", "has_tags"),
        ),
        pytest.raises(KeyError),
    ):
        tag_manager.hasTags(3, "fail")


# TagManager: see test_tagmanager.py::test_verify_tags.
@pytest.mark.parametrize(
    "input_tags, input_silent, expected_error",
    [
        # C1: All tags exist.
        # Expected: Warns and returns True.
        (("3", "number"), False, None),
        # C2: One tag does not exist and the manager is not silent.
        # Expected: Warns and raises a KeyError.
        (("3", "fail"), False, KeyError),
        # C3: One tag does not exist and the manager is silent.
        # Expected: Warns and still raises a KeyError.
        (("3", "fail"), True, KeyError),
    ],
)
def test_verifyTags(tag_manager, input_tags, input_silent, expected_error):
    tag_manager.silent = input_silent
    expected_warning_message = expected_message(
        TAG_MANAGER, "verifyTags", "verify_tags"
    )
    if expected_error is None:
        with pytest.warns(DeprecationWarning, match=expected_warning_message):
            actual_result = tag_manager.verifyTags(*input_tags)
        expected_result = True
        assert actual_result == expected_result
    else:
        with (
            pytest.warns(DeprecationWarning, match=expected_warning_message),
            pytest.raises(expected_error),
        ):
            tag_manager.verifyTags(*input_tags)


# nameutils: see test_nameutils.py::test_is_identifier.
@pytest.mark.parametrize(
    "input_name, expected_is_identifier",
    [
        # C1: A single letter.
        # Expected: Warns and returns True.
        ("x", True),
        # C2: Letters, digits and underscores, starting with an underscore.
        # Expected: Warns and returns True.
        ("_scale_2", True),
        # C3: The name starts with a digit.
        # Expected: Warns and returns False.
        ("2x", False),
        # C4: The name contains a hyphen.
        # Expected: Warns and returns False.
        ("a-b", False),
        # C5: The name contains a space.
        # Expected: Warns and returns False.
        ("a b", False),
        # C6: The name is empty.
        # Expected: Warns and returns False.
        ("", False),
    ],
)
def test_isIdentifier(input_name, expected_is_identifier):
    with pytest.warns(
        DeprecationWarning,
        match=expected_message(NAMEUTILS, "isIdentifier", "is_identifier"),
    ):
        actual_is_identifier = isIdentifier(input_name)
    assert actual_is_identifier == expected_is_identifier


# nameutils: see test_nameutils.py::test_validate_name_valid.
@pytest.mark.parametrize(
    "input_name",
    [
        # C1: A valid name.
        # Expected: Warns and raises no error.
        "scale",
        # C2: A valid name with digits and underscores.
        # Expected: Warns and raises no error.
        "_scale_2",
    ],
)
def test_validateName_valid(input_name):
    with pytest.warns(
        DeprecationWarning,
        match=expected_message(NAMEUTILS, "validateName", "validate_name"),
    ):
        actual_result = validateName(input_name)
    expected_result = None
    assert actual_result == expected_result


# nameutils: see test_nameutils.py::test_validate_name_invalid.
@pytest.mark.parametrize(
    "input_name",
    [
        # C1: The name starts with a digit.
        # Expected: Warns and raises a ValueError that explains the fix.
        "2x",
        # C2: The name contains a hyphen.
        # Expected: Warns and raises a ValueError that explains the fix.
        "a-b",
    ],
)
def test_validateName_invalid(input_name):
    expected_error_message = (
        f"Name '{input_name}' is not a valid identifier. "
        "Please use only letters, digits and underscores, "
        "and do not start the name with a digit."
    )
    with (
        pytest.warns(
            DeprecationWarning,
            match=expected_message(NAMEUTILS, "validateName", "validate_name"),
        ),
        pytest.raises(ValueError, match=re.escape(expected_error_message)),
    ):
        validateName(input_name)


# inpututils: see test_inpututils.py::test_convert_input_to_string.
@pytest.mark.parametrize(
    "input_kind, expected_string",
    [
        # C1: An open file-like object.
        # Expected: Warns and returns its contents.
        ("file_object", MULTILINE_TEXT),
        # C2: The name of an existing file.
        # Expected: Warns and returns the file's contents.
        ("file_name", MULTILINE_TEXT),
        # C3: Text spanning several lines.
        # Expected: Warns and returns the text itself.
        ("multiline_text", MULTILINE_TEXT),
        # C4: A single line of 80 or more characters.
        # Expected: Warns and returns the text itself.
        ("long_text", LONG_TEXT),
    ],
)
def test_inputToString(tmp_path, input_kind, expected_string):
    file_path = tmp_path / "input.txt"
    file_path.write_text(MULTILINE_TEXT)
    inputs = {
        "file_object": io.StringIO(MULTILINE_TEXT),
        "file_name": str(file_path),
        "multiline_text": MULTILINE_TEXT,
        "long_text": LONG_TEXT,
    }
    with pytest.warns(
        DeprecationWarning,
        match=expected_message(
            "diffpy.srfit.util.inpututils",
            "inputToString",
            "convert_input_to_string",
        ),
    ):
        actual_string = inputToString(inputs[input_kind])
    assert actual_string == expected_string


# util: see test_util.py::test_sort_key_for_numeric_string.
@pytest.mark.parametrize(
    "input_string, expected_key",
    [
        # C1: Text, a number, then more text.
        # Expected: Warns, and the number becomes an integer between the
        # text parts.
        ("a12b", ("a", 12, "b")),
        # C2: The string starts with a number.
        # Expected: Warns, and the key starts with an empty text part.
        ("12a", ("", 12, "a")),
        # C3: The string has no digits.
        # Expected: Warns, and the key holds the string unchanged.
        ("abc", ("abc",)),
    ],
)
def test_sortKeyForNumericString(input_string, expected_key):
    with pytest.warns(
        DeprecationWarning,
        match=expected_message(
            "diffpy.srfit.util",
            "sortKeyForNumericString",
            "sort_key_for_numeric_string",
        ),
    ):
        actual_key = sortKeyForNumericString(input_string)
    assert actual_key == expected_key
