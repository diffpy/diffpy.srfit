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
"""Tests for the fitcontribution module."""

import warnings

import pytest
from numpy import allclose, arange, array_equal, dot, nan, sin

from diffpy.srfit import FitContribution, Profile, ProfileGenerator
from diffpy.srfit.exceptions import SrFitError
from diffpy.srfit.fitbase.parameter import Parameter

XOBS = arange(0, 10, 0.5)
YOBS_SINE = 9 * sin(XOBS)


def _fsquare(x):
    return x**2


# set_profile links the profile's x, y and dy to the contribution, and sets
# the default residual equation when a profile equation already exists.
@pytest.mark.parametrize(
    "input_equation, expected_has_equation, expected_has_residual_equation",
    [
        # C1: The profile is set before any equation.
        # Expected: x, y and dy are linked, and neither the equation nor the
        # residual equation is set.
        (None, False, False),
        # C2: The profile is set after the equation "A * x".
        # Expected: x, y and dy are linked, and the residual equation is
        # set.
        ("A * x", True, True),
    ],
)
def test_set_profile(
    input_equation, expected_has_equation, expected_has_residual_equation
):
    profile = Profile()
    contribution = FitContribution("test")
    if input_equation is not None:
        contribution.set_equation(input_equation)
    contribution.set_profile(profile)
    actual_profile = contribution.profile
    expected_profile = profile
    assert actual_profile == expected_profile
    actual_parameters = [
        contribution.x.par,
        contribution.y.par,
        contribution.dy.par,
    ]
    expected_parameters = [profile.xpar, profile.ypar, profile.dypar]
    assert actual_parameters == expected_parameters
    actual_has_equation = contribution._eq is not None
    assert actual_has_equation == expected_has_equation
    actual_has_residual_equation = contribution._reseq is not None
    assert actual_has_residual_equation == expected_has_residual_equation


def test_set_profile_bad():
    # C1: A string is given instead of a Profile.
    # Expected: Raises a TypeError.
    contribution = FitContribution("test")
    with pytest.raises(TypeError):
        contribution.set_profile("invalid")


# add_profile_generator registers the generator, makes it the profile
# equation, and passes it the profile when there is one.
@pytest.mark.parametrize(
    "input_has_profile, expected_has_residual_equation",
    [
        # C1: The contribution has no profile.
        # Expected: The generator is registered and has no profile, an
        # equation is set, and no residual equation is set.
        (False, False),
        # C2: The contribution has a profile.
        # Expected: The generator is registered and gets the profile, and
        # both an equation and the residual equation are set.
        (True, True),
    ],
)
def test_add_profile_generator(
    input_has_profile, expected_has_residual_equation
):
    generator = ProfileGenerator("test")
    profile = Profile()
    contribution = FitContribution("test")
    if input_has_profile:
        contribution.set_profile(profile)
    contribution.add_profile_generator(generator, "gen")
    actual_generator_names = list(contribution._generators)
    expected_generator_names = ["test"]
    assert actual_generator_names == expected_generator_names
    actual_has_equation = contribution._eq is not None
    expected_has_equation = True
    assert actual_has_equation == expected_has_equation
    actual_generator_profile = generator.profile
    expected_generator_profile = profile if input_has_profile else None
    assert actual_generator_profile == expected_generator_profile
    actual_has_residual_equation = contribution._reseq is not None
    assert actual_has_residual_equation == expected_has_residual_equation
    actual_generator_output = generator(XOBS)
    expected_generator_output = XOBS
    assert array_equal(actual_generator_output, expected_generator_output)


def test_profile_generator_interaction():
    # C1: A profile and a generator are added, then data with y = x is
    # observed.
    # Expected: The generator sees the observed data, the profile
    # equation equals y, and the residual is zero.
    generator = ProfileGenerator("test")
    profile = Profile()
    contribution = FitContribution("test")
    contribution.set_profile(profile)
    contribution.add_profile_generator(generator, "I")
    profile.set_observed_profile(XOBS, XOBS)
    actual_generator_profile_arrays = [
        generator.profile.xobs,
        generator.profile.yobs,
        generator.profile.x,
    ]
    expected_generator_profile_arrays = [XOBS, XOBS, XOBS]
    for actual_array, expected_array in zip(
        actual_generator_profile_arrays, expected_generator_profile_arrays
    ):
        assert array_equal(actual_array, expected_array)
    actual_contribution_x = contribution.x.value
    expected_contribution_x = XOBS
    assert array_equal(actual_contribution_x, expected_contribution_x)
    actual_generator_output = generator(XOBS)
    expected_generator_output = XOBS
    assert array_equal(actual_generator_output, expected_generator_output)
    actual_profile_equation = contribution._eq()
    expected_profile_equation = XOBS
    assert array_equal(actual_profile_equation, expected_profile_equation)
    actual_residual_sum = sum(contribution._reseq())
    expected_residual_sum = 0
    assert actual_residual_sum == expected_residual_sum


def test_replace_profile():
    # C1: A contribution with a generator is evaluated on one profile, then
    # set to a second profile with different data.
    # Expected: The generator gets the second profile, the same profile and
    # residual equations are kept but their cached values are cleared, and
    # they are evaluated on the second profile.
    generator = ProfileGenerator("test")
    profile = Profile()
    profile.set_observed_profile(XOBS, XOBS)
    new_xobs = arange(0, 10, 0.8)
    new_profile = Profile()
    new_profile.set_observed_profile(new_xobs, 0.5 * new_xobs)
    contribution = FitContribution("test")
    contribution.set_profile(profile)
    contribution.add_profile_generator(generator, "I")
    actual_initial_generator_value = generator.value
    expected_initial_generator_value = XOBS
    assert array_equal(
        actual_initial_generator_value, expected_initial_generator_value
    )
    actual_initial_profile_equation = contribution._eq()
    expected_initial_profile_equation = XOBS
    assert array_equal(
        actual_initial_profile_equation, expected_initial_profile_equation
    )
    actual_initial_residual_sum = sum(contribution._reseq())
    expected_initial_residual_sum = 0
    assert actual_initial_residual_sum == expected_initial_residual_sum
    actual_initial_residual_length = len(contribution.residual())
    expected_initial_residual_length = len(XOBS)
    assert actual_initial_residual_length == expected_initial_residual_length
    expected_equations = [contribution._eq, contribution._reseq]
    contribution.set_profile(new_profile)
    actual_profiles = [contribution.profile, generator.profile]
    expected_profiles = [new_profile, new_profile]
    assert actual_profiles == expected_profiles
    actual_equations = [contribution._eq, contribution._reseq]
    assert actual_equations == expected_equations
    actual_cached_values = [
        contribution._eq._value,
        contribution._reseq._value,
    ]
    expected_cached_values = [None, None]
    assert actual_cached_values == expected_cached_values
    actual_generator_value = generator.value
    expected_generator_value = new_xobs
    assert array_equal(actual_generator_value, expected_generator_value)
    actual_profile_equation = contribution._eq()
    expected_profile_equation = new_xobs
    assert array_equal(actual_profile_equation, expected_profile_equation)
    actual_residual_length = len(contribution.residual())
    expected_residual_length = len(new_xobs)
    assert actual_residual_length == expected_residual_length


# The residual of a contribution with a generator "I" (which returns x),
# observed at x = 0, 0.5, ..., 9.5. A parameter c is added to the
# contribution with value input_c_value. The squared residual,
# dot(residual, residual), is compared.
@pytest.mark.parametrize(
    "input_equation, input_residual_equation, input_c_value, input_yobs, "
    "expected_chi2",
    [
        # C1: The equation is the generator, as set by
        # add_profile_generator, and y = x.
        # Expected: The default chiv residual is zero.
        (None, None, 2, XOBS, 0),
        # C2: The equation is changed so it no longer matches y = x,
        # 1. to "2*I", 2. to "c*I" with c = 2.
        # Expected: The chiv residual is (2*x - x)/1, so chi2 = dot(x, x).
        ("2*I", None, 2, XOBS, dot(XOBS, XOBS)),
        ("c*I", None, 2, XOBS, dot(XOBS, XOBS)),
        # C3: The equation is "c**2*sin(I)" with c = 3, and y = 9*sin(x).
        # Expected: The chiv residual is zero.
        ("c**2*sin(I)", None, 3, YOBS_SINE, 0),
        # C4: The equation is "2*I", y = 9*sin(x), and the residual is the
        # preset "resv".
        # Expected: chi2 = sum((2*x - y)**2) / sum(y**2).
        (
            "2*I",
            "resv",
            2,
            YOBS_SINE,
            sum((2 * XOBS - YOBS_SINE) ** 2) / sum(YOBS_SINE**2),
        ),
        # C5: The equation is "2*I", y = 9*sin(x), and the residual is the
        # custom "abs(eq-y)**0.5".
        # Expected: chi2 = sum(abs(2*x - y)).
        (
            "2*I",
            "abs(eq-y)**0.5",
            2,
            YOBS_SINE,
            sum(abs(2 * XOBS - YOBS_SINE)),
        ),
    ],
)
def test_residual(
    find_builders_with_observers,
    input_equation,
    input_residual_equation,
    input_c_value,
    input_yobs,
    expected_chi2,
):
    generator = ProfileGenerator("test")
    profile = Profile()
    contribution = FitContribution("test")
    contribution.set_profile(profile)
    contribution.add_profile_generator(generator, "I")
    contribution._add_parameter(Parameter("c", input_c_value))
    if input_equation is not None:
        contribution.set_equation(input_equation)
    if input_residual_equation is not None:
        contribution.set_residual_equation(input_residual_equation)
    profile.set_observed_profile(XOBS, input_yobs)
    residual = contribution.residual()
    actual_chi2 = dot(residual, residual)
    assert actual_chi2 == pytest.approx(expected_chi2)
    actual_builders = find_builders_with_observers()
    expected_builders = []
    assert actual_builders == expected_builders


# A point whose uncertainty is zero or not a number has no weight, so its
# chiv residual is zero. The profile is observed at x = 1, 2, 3, 4, 5 with
# y = x, and the equation is "A*x" with A = 2, so the residual of every
# point with an uncertainty of one is x.
@pytest.mark.parametrize(
    "input_dyobs, expected_residual",
    [
        # C1: A user loads data whose last two points have uncertainties of
        # zero, as at the end of a detector bank.
        # Expected: The residual is zero at those two points, with no
        # warning.
        ([1.0, 1.0, 1.0, 0.0, 0.0], [1.0, 2.0, 3.0, 0.0, 0.0]),
        # C2: A user loads data with one uncertainty that is not a number,
        # at x = 3.
        # Expected: The residual is zero at x = 3, with no warning.
        ([1.0, 1.0, nan, 1.0, 1.0], [1.0, 2.0, 0.0, 4.0, 5.0]),
    ],
)
def test_residual_with_invalid_uncertainties(input_dyobs, expected_residual):
    input_xobs = [1.0, 2.0, 3.0, 4.0, 5.0]
    profile = Profile()
    profile.set_observed_profile(input_xobs, input_xobs, input_dyobs)
    contribution = FitContribution("test")
    contribution.set_profile(profile)
    contribution.set_equation("A*x")
    contribution.A.set_value(2)
    # Compute the residual, recording every warning to check that none is
    # raised
    with warnings.catch_warnings(record=True) as caught_warnings:
        warnings.simplefilter("always")
        actual_residual = contribution.residual()
    actual_warning_messages = []
    for caught_warning in caught_warnings:
        actual_warning_messages.append(str(caught_warning.message))
    expected_warning_messages = []
    assert allclose(actual_residual, expected_residual)
    assert actual_warning_messages == expected_warning_messages


# The cached values of the profile and residual equations are cleared
# whenever their result may change. The contribution has a profile observed
# at x = 0, 0.5, ..., 9.5 with y = x, and a generator "I".
@pytest.mark.parametrize(
    "input_change",
    [
        # C1: Nothing has been evaluated yet.
        # Expected: Both cached values are None.
        None,
        # C2: The residual is computed, then the equation is changed,
        # 1. to "2*I", 2. to "c*I" with a new parameter c.
        # Expected: Both cached values are cleared.
        "2*I",
        "c*I",
        # C3: The residual is computed, then new data y = 9*sin(x) is
        # observed.
        # Expected: Both cached values are cleared.
        "set_observed_profile",
    ],
)
def test_residual_flushes_cached_values(input_change):
    generator = ProfileGenerator("test")
    profile = Profile()
    contribution = FitContribution("test")
    contribution.set_profile(profile)
    contribution.add_profile_generator(generator, "I")
    contribution._add_parameter(Parameter("c", 2))
    if input_change is not None:
        profile.set_observed_profile(XOBS, XOBS)
        contribution.residual()
        if input_change == "set_observed_profile":
            profile.set_observed_profile(XOBS, YOBS_SINE)
        else:
            contribution.set_equation(input_change)
    actual_cached_values = [
        contribution._eq._value,
        contribution._reseq._value,
    ]
    expected_cached_values = [None, None]
    assert actual_cached_values == expected_cached_values


# set_residual_equation needs both a profile and a profile equation.
@pytest.mark.parametrize(
    "input_has_profile",
    [
        # C1: The contribution has neither a profile nor an equation.
        # Expected: Raises an SrFitError.
        False,
        # C2: The contribution has a profile but no equation.
        # Expected: Raises an SrFitError.
        True,
    ],
)
def test_set_residual_equation_bad(input_has_profile):
    contribution = FitContribution("test")
    if input_has_profile:
        contribution.set_profile(Profile())
    with pytest.raises(SrFitError):
        contribution.set_residual_equation("chiv")


# get_residual_equation returns the residual equation as a string, with the
# profile equation written as "eq".
@pytest.mark.parametrize(
    "input_equation, input_residual_equation, expected_residual_equation",
    [
        # C1: Neither a profile nor an equation is set.
        # Expected: Returns an empty string.
        (None, None, ""),
        # C2: A profile and the equation "A * x + B" are set.
        # Expected: Returns the default chiv residual.
        ("A * x + B", None, "((eq - y) / dy)"),
        # C3: The residual equation is then set to "2 * (eq - y)".
        # Expected: Returns the new residual equation.
        ("A * x + B", "2 * (eq - y)", "(2 * (eq - y))"),
    ],
)
def test_get_residual_equation(
    input_equation, input_residual_equation, expected_residual_equation
):
    contribution = FitContribution("test")
    if input_equation is not None:
        contribution.set_profile(Profile())
        contribution.set_equation(input_equation)
    if input_residual_equation is not None:
        contribution.set_residual_equation(input_residual_equation)
    actual_residual_equation = contribution.get_residual_equation()
    assert actual_residual_equation == expected_residual_equation


# The equation factory does not hold on to equations that were replaced.
@pytest.mark.parametrize(
    "input_equation_repeats, input_residual_repeats, expected_equation_count",
    [
        # C1: Nothing is set.
        # Expected: The factory holds no equations.
        (0, 0, 0),
        # C2: set_equation is called five times.
        # Expected: The factory holds one equation.
        (5, 0, 1),
        # C3: set_equation is called five times, then a profile is set and
        # set_residual_equation is called five times.
        # Expected: The factory holds two equations.
        (5, 5, 2),
    ],
)
def test_release_old_equations(
    input_equation_repeats, input_residual_repeats, expected_equation_count
):
    contribution = FitContribution("test")
    for _ in range(input_equation_repeats):
        contribution.set_equation("A * x + B")
    if input_residual_repeats > 0:
        contribution.set_profile(Profile())
    for _ in range(input_residual_repeats):
        contribution.set_residual_equation("chiv")
    actual_equation_count = len(contribution._eqfactory.equations)
    assert actual_equation_count == expected_equation_count


# A registered function still works after set_equation is called a second
# time. The equation is first set to "fsquare" and evaluated at x = 5.
@pytest.mark.parametrize(
    "input_equation, input_x, expected_value",
    [
        # C1: The equation is set to "fsquare" again.
        # Expected: Returns 5**2 = 25.
        ("fsquare", 5, 25),
        # C2: The equation is set to "fsquare" again, and x is set to 6.
        # Expected: Returns 6**2 = 36.
        ("fsquare", 6, 36),
        # C3: The equation is set to "fsquare + 5", 1. with x = 6,
        # 2. with x = -1.
        # Expected: Returns x**2 + 5, so 1. 41, 2. 6.
        ("fsquare + 5", 6, 41),
        ("fsquare + 5", -1, 6),
    ],
)
def test_register_function(input_equation, input_x, expected_value):
    contribution = FitContribution("test")
    contribution.register_function(_fsquare, name="fsquare")
    contribution.set_equation("fsquare")
    contribution.x.set_value(5)
    contribution.evaluate()
    contribution.set_equation(input_equation)
    contribution.x << input_x
    actual_value = contribution.evaluate()
    assert actual_value == expected_value


# set_equation picks up a parameter that replaces a removed one.
@pytest.mark.parametrize(
    "input_x_replacement, expected_value",
    [
        # C1: The equation "x + 5" is evaluated with x = 2.
        # Expected: Returns 7.
        (None, 7),
        # C2: x is then removed and replaced by a new parameter x holding
        # an array.
        # Expected: Returns 5 + the array.
        (XOBS, 5 + XOBS),
    ],
)
def test_set_equation(
    find_builders_with_observers, input_x_replacement, expected_value
):
    contribution = FitContribution("test")
    contribution.set_equation("x + 5")
    contribution.x.set_value(2)
    actual_value = contribution.evaluate()
    if input_x_replacement is not None:
        contribution.remove_parameter(contribution.x)
        contribution.new_parameter("x", input_x_replacement)
        actual_value = contribution.evaluate()
    assert array_equal(actual_value, expected_value)
    actual_builders = find_builders_with_observers()
    expected_builders = []
    assert actual_builders == expected_builders


# get_equation returns the profile equation as a string.
@pytest.mark.parametrize(
    "input_equation, expected_equation",
    [
        # C1: No equation is set.
        # Expected: Returns an empty string.
        (None, ""),
        # C2: The equation "A * sin(x + 5)" is set.
        # Expected: Returns it fully parenthesized.
        ("A * sin(x + 5)", "(A * sin((x + 5)))"),
    ],
)
def test_get_equation(
    find_builders_with_observers, input_equation, expected_equation
):
    contribution = FitContribution("test")
    if input_equation is not None:
        contribution.set_equation(input_equation)
    actual_equation = contribution.get_equation()
    assert actual_equation == expected_equation
    actual_builders = find_builders_with_observers()
    expected_builders = []
    assert actual_builders == expected_builders
