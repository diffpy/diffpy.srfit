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
"""Tests for the profile module."""

import io
import re

import pytest
from numpy import allclose, arange, array, array_equal

from diffpy.srfit import FitContribution, Profile, ProfileParser
from diffpy.srfit.exceptions import SrFitError


# A new profile holds no data and no metadata.
@pytest.mark.parametrize(
    "input_attribute, expected_value",
    [
        # C1: An observed array is read, 1. xobs, 2. yobs, 3. dyobs.
        # Expected: Returns None.
        ("xobs", None),
        ("yobs", None),
        ("dyobs", None),
        # C2: A calculated array is read, 1. x, 2. y, 3. dy, 4. ycalc.
        # Expected: Returns None.
        ("x", None),
        ("y", None),
        ("dy", None),
        ("ycalc", None),
        # C3: The metadata is read.
        # Expected: Returns an empty dict.
        ("meta", {}),
    ],
)
def test_init(input_attribute, expected_value):
    profile = Profile()
    actual_value = getattr(profile, input_attribute)
    assert actual_value == expected_value


# set_observed_profile stores the observed arrays and puts the calculation
# points x, y and dy on the observed grid.
@pytest.mark.parametrize(
    "input_dyobs, expected_dyobs, expected_dy",
    [
        # C1: Uncertainties dyobs are given.
        # Expected: dyobs is stored, and dy equals dyobs.
        ([0.1, 0.2, 0.3], [0.1, 0.2, 0.3], [0.1, 0.2, 0.3]),
        # C2: dyobs is None.
        # Expected: dyobs stays None to mark the profile as unweighted,
        # and dy defaults to ones.
        (None, None, [1.0, 1.0, 1.0]),
    ],
)
def test_set_observed_profile(
    as_list, input_dyobs, expected_dyobs, expected_dy
):
    input_xobs = [1.0, 2.0, 3.0]
    input_yobs = [2.0, 4.0, 6.0]
    profile = Profile()
    profile.set_observed_profile(input_xobs, input_yobs, input_dyobs)
    actual_observed = (
        profile.xobs.tolist(),
        profile.yobs.tolist(),
        as_list(profile.dyobs),
    )
    expected_observed = (input_xobs, input_yobs, expected_dyobs)
    assert actual_observed == expected_observed
    actual_calculated = (profile.x.tolist(), profile.y.tolist())
    expected_calculated = (input_xobs, input_yobs)
    assert actual_calculated == expected_calculated
    actual_dy = profile.dy.tolist()
    assert actual_dy == expected_dy


# set_calculation_range needs an observed profile to work on.
@pytest.mark.parametrize(
    "input_args",
    [
        # C1: No observed profile is set, and set_calculation_range is called
        # 1. with no arguments, 2. with xmin, 3. with xmin and xmax,
        # 4. with xmin, xmax and dx.
        # Expected: Raises an AttributeError.
        (),
        (0,),
        (0, 5),
        (0, 5, 0.2),
    ],
)
def test_set_calculation_range_without_data_bad(input_args):
    profile = Profile()
    with pytest.raises(AttributeError):
        profile.set_calculation_range(*input_args)


# set_calculation_range keeps the current bound when one is not given,
# clips the bounds to the observed data, adheres to the observed points when
# dx matches the data, and resamples anchored at xmin otherwise. The profile
# is observed at x = 2, 2.5, ..., 9.5 with y = dy = x, and is first set to
# input_initial_range.
@pytest.mark.parametrize(
    "input_initial_range, input_range, expected_x",
    [
        # C1: No range is changed, 1. no arguments are given,
        # 2. xmin is below the observed minimum, 3. xmax is above the
        # observed maximum.
        # Expected: The calculation points stay on the full observed grid.
        ({}, {}, arange(2, 9.6, 0.5)),
        ({}, {"xmin": 0}, arange(2, 9.6, 0.5)),
        ({}, {"xmax": 100}, arange(2, 9.6, 0.5)),
        # C2: The range is narrowed by xmin and xmax on the observed grid.
        # Expected: The observed points within [4, 7] are kept.
        ({}, {"xmin": 4, "xmax": 7}, arange(4, 7.1, 0.5)),
        # C3: The range is [3, 5], then xmin="obs", xmax=7 and dx=0.001.
        # Expected: xmin is reset to the observed minimum, and the points
        # are resampled from 2 to 7 in steps of 0.001, so 5001 points.
        (
            {"xmin": 3, "xmax": 5},
            {"xmin": "obs", "xmax": 7, "dx": 0.001},
            arange(2, 7.0005, 0.001),
        ),
        # C4: The points are resampled, then xmin=4, xmax=7 and dx="obs".
        # Expected: The points return to the observed grid within [4, 7].
        (
            {"xmin": "obs", "xmax": 7, "dx": 0.001},
            {"xmin": 4, "xmax": 7, "dx": "obs"},
            arange(4, 7.1, 0.5),
        ),
        # C5: Only one bound is changed on the observed grid, 1. xmin=3
        # after [4, 7], 2. xmax=5.1 after [3, 7].
        # Expected: The other bound is kept.
        ({"xmin": 4, "xmax": 7}, {"xmin": 3}, arange(3, 7.1, 0.5)),
        ({"xmin": 3, "xmax": 7}, {"xmax": 5.1}, arange(3, 5.1, 0.5)),
        # C6: The range is [3, 5], then only dx=1 is given.
        # Expected: The points are resampled from 3 to 5 in steps of 1.
        ({"xmin": 3, "xmax": 5}, {"dx": 1}, arange(3, 5.1)),
        # C7: The points are on [3, 5] in steps of 1, then xmin=4.2, xmax=7
        # and dx=0.3, off the observed grid.
        # Expected: The points are resampled from 4.2 in steps of 0.3, up
        # to the largest point not above 7.
        (
            {"xmin": 3, "xmax": 5, "dx": 1},
            {"xmin": 4.2, "xmax": 7, "dx": 0.3},
            arange(4.2, 6.901, 0.3),
        ),
        # C8: The points are resampled from 4.2 in steps of 0.3, then only
        # xmin=4.2 and xmax=6.001 are given.
        # Expected: The current step of 0.3 is kept.
        (
            {"xmin": 4.2, "xmax": 7, "dx": 0.3},
            {"xmin": 4.2, "xmax": 6.001},
            arange(4.2, 6.001, 0.3),
        ),
        # C9: The points are on [4.2, 6] in steps of 0.3, then only dx=0.5
        # is given, which matches the observed step.
        # Expected: The observed points within [4.2, 6] are kept.
        (
            {"xmin": 4.2, "xmax": 6.001, "dx": 0.3},
            {"dx": 0.5},
            arange(4.5, 6.1, 0.5),
        ),
    ],
)
def test_set_calculation_range(
    observed_profile, input_initial_range, input_range, expected_x
):
    observed_profile.set_calculation_range(**input_initial_range)
    observed_profile.set_calculation_range(**input_range)
    actual_x = observed_profile.x
    assert array_equal(actual_x, expected_x)
    # y = dy = x in the observed data, so this holds on any grid.
    actual_y = observed_profile.y
    actual_dy = observed_profile.dy
    assert allclose(actual_y, expected_x)
    assert allclose(actual_dy, expected_x)


# set_calculation_range rejects a range it cannot build. The profile is
# observed at x = 2, 2.5, ..., 9.5.
@pytest.mark.parametrize(
    "input_range",
    [
        # C1: xmin is greater than xmax.
        # Expected: Raises a ValueError.
        {"xmin": 10, "xmax": 3},
        # C2: dx is greater than xmax - xmin.
        # Expected: Raises a ValueError.
        {"xmin": 3, "xmax": 3.9, "dx": 1.0},
        # C3: dx is not positive, 1. dx=0, 2. dx=-0.000001.
        # Expected: Raises a ValueError.
        {"dx": 0},
        {"dx": -0.000001},
        # C4: A string other than "obs" is given, 1. as xmin, 2. as xmax,
        # 3. as dx.
        # Expected: Raises a ValueError.
        {"xmin": "oobs"},
        {"xmax": "oobs"},
        {"dx": "oobs"},
    ],
)
def test_set_calculation_range_bad(observed_profile, input_range):
    with pytest.raises(ValueError):
        observed_profile.set_calculation_range(**input_range)


# set_calculation_points takes the points as given, but limits them to the
# bounds of the observed data once there is any.
@pytest.mark.parametrize(
    "input_xobs, expected_x",
    [
        # C1: No observed profile is set.
        # Expected: x is the given points, from 3 to 12 in steps of 0.2.
        (None, arange(3, 12.2, 0.2)),
        # C2: An observed profile on [2, 10] is set after the points.
        # Expected: x is limited to the observed bounds, from 3 to 10.
        (arange(2, 10.5, 0.5), arange(3, 10.1, 0.2)),
    ],
)
def test_set_calculation_points(input_xobs, expected_x):
    profile = Profile()
    profile.set_calculation_points(arange(3, 12.2, 0.2))
    if input_xobs is not None:
        profile.set_observed_profile(
            input_xobs, array(input_xobs), array(input_xobs)
        )
    actual_x = profile.x
    assert array_equal(actual_x, expected_x)


def test_savetxt():
    # C1: A profile with observed data and ycalc is saved to a file.
    # Expected: The file has a header naming the x, ycalc, y and dy
    # columns, followed by one line for each of the 21 points.
    profile = Profile()
    xobs = arange(-2, 3.01, 0.25)
    yobs = xobs**2
    profile.set_observed_profile(xobs, yobs)
    profile.ycalc = yobs.copy()
    file_object = io.BytesIO()
    profile.savetxt(file_object)
    actual_text = file_object.getvalue().decode()
    assert re.match(r"^# x +ycalc +y +dy\b", actual_text)
    actual_line_count = len(actual_text.strip().split("\n"))
    expected_line_count = 22
    assert actual_line_count == expected_line_count


def test_savetxt_bad():
    # C1: A profile without ycalc is saved.
    # Expected: Raises an SrFitError.
    profile = Profile()
    with pytest.raises(SrFitError):
        profile.savetxt("foo")


# loadtxt reads x, y and dy with numpy.loadtxt and always unpacks the
# columns. The first point of testdata.txt is x=1e-2, y=1.105784e-1 and
# dy=1.802192e-3.
@pytest.mark.parametrize(
    "input_args, input_kwargs, expected_first_point",
    [
        # C1: The x, y and dy columns are selected with usecols.
        # Expected: The first point is loaded with its dy.
        ((), {"usecols": (0, 1, 3)}, (1e-2, 1.105784e-1, 1.802192e-3)),
        # C2: unpack=False is given, 1. as a keyword, 2. as the last of
        # all eight positional arguments.
        # Expected: unpack is enforced, and the first point is loaded with
        # its dy.
        (
            (),
            {"usecols": (0, 1, 3), "unpack": False},
            (1e-2, 1.105784e-1, 1.802192e-3),
        ),
        (
            (float, "#", None, None, 0, (0, 1, 3), False),
            {},
            (1e-2, 1.105784e-1, 1.802192e-3),
        ),
        # C3: Only the x and y columns are selected.
        # Expected: The first point is loaded, and dy defaults to 1.
        ((), {"usecols": (0, 1)}, (1e-2, 1.105784e-1, 1)),
    ],
)
def test_loadtxt(datafile, input_args, input_kwargs, expected_first_point):
    profile = Profile()
    profile.loadtxt(datafile("testdata.txt"), *input_args, **input_kwargs)
    actual_first_point = (profile.x[0], profile.y[0], profile.dy[0])
    assert actual_first_point == pytest.approx(expected_first_point)


def test_loadtxt_bad(datafile):
    # C1: Only the x column is selected.
    # Expected: Raises a ValueError, since y is required.
    profile = Profile()
    with pytest.raises(ValueError):
        profile.loadtxt(datafile("testdata.txt"), usecols=(0,))


# The parsed x, y and dy arrays are copied onto the observed profile.
# Uncertainties on x are dropped, since srfit treats the independent
# variable as having no uncertainty.
@pytest.mark.parametrize(
    "input_filename, expected_xobs, expected_yobs, expected_dyobs",
    [
        # C1: File has four columns, so uncertainties are present.
        # Expected: xobs, yobs and dyobs are loaded and dx is dropped.
        (
            "four_col.gr",
            [1.0, 1.1, 1.2],
            [2.0, 2.1, 2.2],
            [0.2, 0.4, 0.6],
        ),
        # C2: File has three columns, so uncertainties are present.
        # Expected: xobs, yobs and dyobs are loaded.
        (
            "three_col.dat",
            [1.0, 1.1, 1.2],
            [2.0, 2.1, 2.2],
            [0.2, 0.4, 0.6],
        ),
        # C3: File has two columns, so no uncertainties are present.
        # Expected: xobs and yobs are loaded and dyobs stays None,
        # marking the profile as unweighted.
        (
            "two_col.txt",
            [1.0, 1.1, 1.2],
            [2.0, 2.1, 2.2],
            None,
        ),
    ],
)
def test_load_parsed_data(
    as_list,
    parser_datafiles,
    input_filename,
    expected_xobs,
    expected_yobs,
    expected_dyobs,
):
    parser = ProfileParser()
    parser.parse_file(parser_datafiles / input_filename)
    profile = Profile()
    profile.load_parsed_data(parser)
    actual_xobs = profile.xobs.tolist()
    actual_yobs = profile.yobs.tolist()
    # Unavailable uncertainties are None rather than an array.
    actual_dyobs = as_list(profile.dyobs)
    assert actual_xobs == expected_xobs
    assert actual_yobs == expected_yobs
    assert actual_dyobs == expected_dyobs


def test_residual_without_uncertainties():
    # C1: A profile is set without uncertainties, so dyobs is None, and a
    # residual is computed from the equation y = A*x with A = 2.
    # Expected: The residual is (A*x - y) with unit weights, since missing
    # uncertainties default to ones rather than being used as divisors.
    profile = Profile()
    profile.set_observed_profile(
        array([1.0, 2.0, 3.0]), array([1.0, 2.0, 3.0])
    )
    contribution = FitContribution("test")
    contribution.set_profile(profile)
    contribution.set_equation("A*x")
    contribution.A.set_value(2)
    actual_residual = contribution.residual().tolist()
    expected_residual = [1.0, 2.0, 3.0]
    assert actual_residual == expected_residual
