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
"""Tests for the restraint module."""

import numpy as np
import pytest

from diffpy.srfit.fitbase.recipeorganizer import get_equation_from_string
from diffpy.srfit.fitbase.restraint import Restraint


# Test the penalty of the restraint 1 < p1 + p2 < upper_bound, which is the
# squared distance of p1 + p2 outside the bounds, multiplied by the weight
# only when the restraint is scaled.
@pytest.mark.parametrize(
    "input_p1, input_p2, input_upper_bound, input_scaled, input_weight, "
    "expected_penalty",
    [
        # C1: p1 + p2 = 2 lies inside [1, 5].
        # Expected: No penalty.
        (1, 1, 5, False, 1.0, 0),
        # C2: p1 + p2 = 0 lies below the lower bound 1.
        # Expected: A penalty of (1 - 0)**2 = 1.
        (-1, 1, 5, False, 1.0, 1),
        # C3: p1 + p2 = 8 lies above the upper bound 5.
        # Expected: A penalty of (8 - 5)**2 = 9.
        (4, 4, 5, False, 1.0, 9),
        # C4: p1 + p2 = 8 lies above the upper bound 5, the restraint is
        # scaled and the weight is 1.5.
        # Expected: A penalty of 1.5 * (8 - 5)**2 = 13.5.
        (4, 4, 5, True, 1.5, 13.5),
        # C5: p1 + p2 = 8 lies above the upper bound 5, the restraint is not
        # scaled and the weight is 1.5.
        # Expected: The weight is ignored, a penalty of (8 - 5)**2 = 9.
        (4, 4, 5, False, 1.5, 9),
        # C6: p1 + p2 is about 1e100 and the upper bound is inf.
        # Expected: No penalty.
        (1e100, 4, np.inf, False, 1.0, 0),
    ],
)
def test_penalty(
    make_parameters_with_factory,
    input_p1,
    input_p2,
    input_upper_bound,
    input_scaled,
    input_weight,
    expected_penalty,
):
    p1, p2, factory = make_parameters_with_factory()
    equation = get_equation_from_string("p1 + p2", factory)
    restraint = Restraint(
        equation,
        lower_bound=1,
        upper_bound=input_upper_bound,
        scaled=input_scaled,
    )
    p1.set_value(input_p1)
    p2.set_value(input_p2)
    actual_penalty = restraint.penalty(input_weight)
    assert actual_penalty == expected_penalty
