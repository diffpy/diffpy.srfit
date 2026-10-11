#!/usr/bin/env python
##############################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2010 The Trustees of Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################
"""Unit tests for tagmanager.py."""

import pytest

# The tag_manager fixture (conftest.py) tags 3 with "3" and "number", and 4
# with "4" and "number". It is not silent unless a test makes it so.


# Tagging stores every tag, as a string, on the object.
@pytest.mark.parametrize(
    "input_tag_calls, expected_tags",
    [
        # C1: A new object is tagged once with one tag.
        # Expected: The object has that tag.
        ([("three",)], {"three"}),
        # C2: A new object is tagged in two calls with several tags.
        # Expected: The object has every tag from both calls.
        ([("three", "tri"), ("trois",)], {"three", "tri", "trois"}),
        # C3: A new object is tagged with a non-string tag.
        # Expected: The tag is stored as a string.
        ([(3,)], {"3"}),
    ],
)
def test_tag(tag_manager, input_tag_calls, expected_tags):
    for tags in input_tag_calls:
        tag_manager.tag("obj", *tags)
    actual_tags = set(tag_manager.tags("obj"))
    assert actual_tags == expected_tags


def test_tag_unhashable_object_bad(tag_manager):
    # C1: An unhashable object (a set) is tagged.
    # Expected: A TypeError is raised.
    with pytest.raises(TypeError):
        tag_manager.tag(set(), "unhashable")


# Untagging removes the given tags, or all tags when none are given.
@pytest.mark.parametrize(
    "input_tags, expected_tags",
    [
        # C1: One tag is removed.
        # Expected: The other tags remain.
        (("3",), {"three", "tri", "tres", "trois"}),
        # C2: Several tags are removed at once.
        # Expected: Only the tags that were not removed remain.
        (("3", "three", "tri"), {"tres", "trois"}),
        # C3: No tags are given.
        # Expected: Every tag is removed from the object.
        ((), set()),
    ],
)
def test_untag(tag_manager, input_tags, expected_tags):
    tag_manager.tag("obj", "3", "three", "tri", "tres", "trois")
    tag_manager.untag("obj", *input_tags)
    actual_tags = set(tag_manager.tags("obj"))
    assert actual_tags == expected_tags


# When not silent, untagging something that is not tagged raises KeyError.
@pytest.mark.parametrize(
    "input_object, input_tag",
    [
        # C1: The tag does not exist.
        # Expected: A KeyError is raised.
        (3, "5"),
        # C2: The tag exists but does not apply to the object.
        # Expected: A KeyError is raised.
        (4, "3"),
        # C3: Neither the tag nor the object is known.
        # Expected: A KeyError is raised.
        (5, "5"),
    ],
)
def test_untag_missing_tag_bad(tag_manager, input_object, input_tag):
    with pytest.raises(KeyError):
        tag_manager.untag(input_object, input_tag)


# Union returns objects with any of the tags; intersection returns objects
# with all of them. A silent manager treats unknown tags as empty.
@pytest.mark.parametrize(
    "input_method, input_tags, input_silent, expected_objects",
    [
        # C1: Union of no tags.
        # Expected: An empty set.
        ("union", (), False, set()),
        # C2: Union of a tag shared by both objects.
        # Expected: Both objects.
        ("union", ("number",), False, {3, 4}),
        # C3: Union of a tag on one object.
        # Expected: That object.
        ("union", ("3",), False, {3}),
        # C4: Union of tags on different objects.
        # Expected: Both objects.
        ("union", ("3", "4"), False, {3, 4}),
        # C5: Silent union of an unknown tag.
        # Expected: An empty set.
        ("union", ("fail",), True, set()),
        # C6: Silent union of an unknown tag and a known tag.
        # Expected: The object with the known tag.
        ("union", ("fail", "3"), True, {3}),
        # C7: Intersection of no tags.
        # Expected: An empty set.
        ("intersection", (), False, set()),
        # C8: Intersection of a tag shared by both objects.
        # Expected: Both objects.
        ("intersection", ("number",), False, {3, 4}),
        # C9: Intersection of a tag on one object.
        # Expected: That object.
        ("intersection", ("3",), False, {3}),
        # C10: Intersection of tags on different objects.
        # Expected: An empty set.
        ("intersection", ("3", "4"), False, set()),
        # C11: Silent intersection of an unknown tag.
        # Expected: An empty set.
        ("intersection", ("fail",), True, set()),
    ],
)
def test_union_and_intersection(
    tag_manager, input_method, input_tags, input_silent, expected_objects
):
    tag_manager.silent = input_silent
    actual_objects = getattr(tag_manager, input_method)(*input_tags)
    assert actual_objects == expected_objects


@pytest.mark.parametrize("input_method", ["union", "intersection"])
def test_union_and_intersection_unknown_tag_bad(tag_manager, input_method):
    # C1: union("fail") or intersection("fail") is called, but no object
    # was ever tagged "fail", and the manager is not silent.
    # Expected: A KeyError is raised.
    with pytest.raises(KeyError):
        getattr(tag_manager, input_method)("fail")


# has_tags is True only when the object has every given tag.
@pytest.mark.parametrize(
    "input_tags, input_silent, expected_has_tags",
    [
        # C1: The object has the one given tag.
        # Expected: True.
        (("3",), False, True),
        # C2: The object has both given tags.
        # Expected: True.
        (("3", "number"), False, True),
        # C3: The object has only one of the given tags.
        # Expected: False.
        (("3", "4"), False, False),
        # C4: A silent manager is asked about an unknown tag.
        # Expected: False.
        (("fail",), True, False),
    ],
)
def test_has_tags(tag_manager, input_tags, input_silent, expected_has_tags):
    tag_manager.silent = input_silent
    actual_has_tags = tag_manager.has_tags(3, *input_tags)
    assert actual_has_tags == expected_has_tags


def test_has_tags_unknown_tag_bad(tag_manager):
    # C1: has_tags(3, "fail") is called, but no object was ever tagged
    # "fail", and the manager is not silent.
    # Expected: A KeyError is raised.
    with pytest.raises(KeyError):
        tag_manager.has_tags(3, "fail")


# verify_tags passes for existing tags and raises KeyError for unknown ones,
# whether or not the manager is silent.
@pytest.mark.parametrize(
    "input_tags, input_silent, expected_error",
    [
        # C1: All tags exist.
        # Expected: Returns True.
        (("3", "number"), False, None),
        # C2: One tag does not exist and the manager is not silent.
        # Expected: A KeyError is raised.
        (("3", "fail"), False, KeyError),
        # C3: One tag does not exist and the manager is silent.
        # Expected: A KeyError is still raised.
        (("3", "fail"), True, KeyError),
    ],
)
def test_verify_tags(tag_manager, input_tags, input_silent, expected_error):
    tag_manager.silent = input_silent
    if expected_error is None:
        actual_result = tag_manager.verify_tags(*input_tags)
        assert actual_result is True
    else:
        with pytest.raises(expected_error):
            tag_manager.verify_tags(*input_tags)
