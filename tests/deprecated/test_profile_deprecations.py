"""Tests for the deprecated camelCase names in
diffpy.srfit.fitbase.profile.

Each test calls the deprecated name, checks that it warns, and checks
that it does the same as the new name.

Remove this file in version 4.0.0, together with the deprecated names.
"""

import re

import pytest
from numpy import arange

from diffpy.srfit import Profile, ProfileParser


def expected_message(base, old_name, new_name):
    return re.escape(
        f"'{base}.{old_name}' is deprecated and will be removed in version "
        f"4.0.0. Please use '{base}.{new_name}' instead."
    )


PROFILE = "diffpy.srfit.fitbase.profile.Profile"


def _make_profile(method_name):
    """Return a profile to call the method on.

    set_calculation_range needs a profile observed at x = 2, 2.5, ...,
    9.5 with y = dy = x. The other methods start from an empty profile.
    """
    profile = Profile()
    if method_name == "set_calculation_range":
        xobs = arange(2, 9.6, 0.5)
        profile.set_observed_profile(xobs, xobs.copy(), xobs.copy())
    return profile


def _make_call_args(method_name, parser_datafiles):
    """Return fresh arguments for a Profile method."""
    if method_name == "set_observed_profile":
        return ([1.0, 2.0, 3.0], [2.0, 4.0, 6.0])
    if method_name == "set_calculation_range":
        return (4, 7)
    if method_name == "set_calculation_points":
        return (arange(3, 12.2, 0.2),)
    parser = ProfileParser()
    parser.parse_file(parser_datafiles / "four_col.gr")
    return (parser,)


def _describe_profile(profile, as_list):
    return (
        as_list(profile.xobs),
        as_list(profile.yobs),
        as_list(profile.dyobs),
        as_list(profile.x),
        as_list(profile.y),
        as_list(profile.dy),
        profile.meta,
    )


# Profile: each deprecated method does the same as its new name. See
# test_profile.py.
@pytest.mark.parametrize(
    "input_old_name, input_new_name",
    [
        # C1: setObservedProfile sets x = [1, 2, 3] and y = 2*x without dy.
        # Expected: Warns and gives the same profile as
        # set_observed_profile.
        ("setObservedProfile", "set_observed_profile"),
        # C2: setCalculationRange sets the range [4, 7] on observed data.
        # Expected: Warns and gives the same profile as
        # set_calculation_range.
        ("setCalculationRange", "set_calculation_range"),
        # C3: setCalculationPoints sets x = 3, 3.2, ..., 12 without data.
        # Expected: Warns and gives the same profile as
        # set_calculation_points.
        ("setCalculationPoints", "set_calculation_points"),
        # C4: loadParsedData loads a parsed four-column file.
        # Expected: Warns and gives the same profile as load_parsed_data.
        ("loadParsedData", "load_parsed_data"),
    ],
)
def test_profile_deprecations(
    as_list, parser_datafiles, input_old_name, input_new_name
):
    old_profile = _make_profile(input_new_name)
    old_args = _make_call_args(input_new_name, parser_datafiles)
    with pytest.warns(
        DeprecationWarning,
        match=expected_message(PROFILE, input_old_name, input_new_name),
    ):
        old_return_value = getattr(old_profile, input_old_name)(*old_args)
    new_profile = _make_profile(input_new_name)
    new_args = _make_call_args(input_new_name, parser_datafiles)
    new_return_value = getattr(new_profile, input_new_name)(*new_args)
    actual_old_result = (
        old_return_value,
        _describe_profile(old_profile, as_list),
    )
    actual_new_result = (
        new_return_value,
        _describe_profile(new_profile, as_list),
    )
    assert actual_old_result == actual_new_result
