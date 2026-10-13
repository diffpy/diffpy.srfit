"""Tests for the deprecated camelCase names in
diffpy.srfit.fitbase.fitcontribution.

Each test calls the deprecated name, checks that it warns, and checks
that it does the same as the new name.

Remove this file in version 4.0.0, together with the deprecated names.
"""

import re

import pytest

from diffpy.srfit import FitContribution, Profile, ProfileGenerator


def expected_message(base, old_name, new_name):
    return re.escape(
        f"'{base}.{old_name}' is deprecated and will be removed in version "
        f"4.0.0. Please use '{base}.{new_name}' instead."
    )


FIT_CONTRIBUTION = "diffpy.srfit.fitbase.FitContribution"


def _make_contribution(method_name):
    """Return a contribution to call the method on.

    set_profile is called on a contribution that has the equation
    "A * x". The other methods are called on a contribution that has a
    profile, and the getters and set_residual_equation also on one that
    has the equation "A * sin(x + 5)".
    """
    contribution = FitContribution("test")
    if method_name == "set_profile":
        contribution.set_equation("A * x")
        return contribution
    contribution.set_profile(Profile())
    if method_name in (
        "get_equation",
        "set_residual_equation",
        "get_residual_equation",
    ):
        contribution.set_equation("A * sin(x + 5)")
    return contribution


def _make_call_args(method_name):
    """Return fresh arguments for a FitContribution method."""
    if method_name == "set_profile":
        return (Profile(),)
    if method_name == "add_profile_generator":
        return (ProfileGenerator("gen"), "I")
    if method_name == "set_equation":
        return ("A * x + B",)
    if method_name == "set_residual_equation":
        return ("2 * (eq - y)",)
    return ()


def _describe_contribution(contribution):
    return (
        sorted(contribution.get_names()),
        sorted(contribution._generators),
        contribution.get_equation(),
        contribution.get_residual_equation(),
    )


# FitContribution: each deprecated method does the same as its new name.
# See test_contribution.py.
@pytest.mark.parametrize(
    "input_old_name, input_new_name",
    [
        # C1: setProfile sets a profile on a contribution with an equation.
        # Expected: Warns and gives the same contribution as set_profile.
        ("setProfile", "set_profile"),
        # C2: addProfileGenerator adds a generator under the name "I".
        # Expected: Warns and gives the same contribution as
        # add_profile_generator.
        ("addProfileGenerator", "add_profile_generator"),
        # C3: setEquation sets the equation "A * x + B".
        # Expected: Warns and gives the same contribution as set_equation.
        ("setEquation", "set_equation"),
        # C4: getEquation gets the equation "A * sin(x + 5)".
        # Expected: Warns and returns the same string as get_equation.
        ("getEquation", "get_equation"),
        # C5: setResidualEquation sets the residual "2 * (eq - y)".
        # Expected: Warns and gives the same contribution as
        # set_residual_equation.
        ("setResidualEquation", "set_residual_equation"),
        # C6: getResidualEquation gets the default chiv residual.
        # Expected: Warns and returns the same string as
        # get_residual_equation.
        ("getResidualEquation", "get_residual_equation"),
    ],
)
def test_fitcontribution_deprecations(
    find_builders_with_observers, input_old_name, input_new_name
):
    old_contribution = _make_contribution(input_new_name)
    old_args = _make_call_args(input_new_name)
    with pytest.warns(
        DeprecationWarning,
        match=expected_message(
            FIT_CONTRIBUTION, input_old_name, input_new_name
        ),
    ):
        old_return_value = getattr(old_contribution, input_old_name)(*old_args)
    new_contribution = _make_contribution(input_new_name)
    new_args = _make_call_args(input_new_name)
    new_return_value = getattr(new_contribution, input_new_name)(*new_args)
    actual_old_result = (
        old_return_value,
        _describe_contribution(old_contribution),
    )
    actual_new_result = (
        new_return_value,
        _describe_contribution(new_contribution),
    )
    assert actual_old_result == actual_new_result
    actual_builders = find_builders_with_observers()
    expected_builders = []
    assert actual_builders == expected_builders
