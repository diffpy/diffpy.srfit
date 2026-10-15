#!/usr/bin/env python
##############################################################################
#
# (c) 2026-present The DiffPy Team. All rights reserved.
#
# File coded by: Caden Myers
#
# See GitHub contributions for a more detailed list of contributors.
# https://github.com/diffpy/diffpy.srfit/graphs/contributors
#
# See LICENSE.rst for license information.
#
##############################################################################
"""Tests for the interface module."""

from diffpy.srfit.fitbase import FitRecipe


def _make_recipe_with_variables():
    """Return a recipe with the variables p1 = 1 and p2 = 2."""
    recipe = FitRecipe("recipe")
    p1 = recipe.create_new_variable("p1", 1)
    p2 = recipe.create_new_variable("p2", 2)
    return recipe, p1, p2


# Test that "*=" constrains a Parameter to an equation.
# C1: p1 is constrained to "2*p2" with "*=", and p2 is then set to 3.
# Expected: p1 is constrained and its value is 6.
def test_constrain_operator():
    recipe, p1, p2 = _make_recipe_with_variables()
    recipe *= "p1", "2*p2"
    p2.set_value(3)
    recipe._constraints[p1].update()
    assert recipe.is_constrained(p1)
    actual_p1_value = p1.get_value()
    expected_p1_value = 6
    assert actual_p1_value == expected_p1_value


# Test that "%=" restrains an equation between soft bounds.
# C1: "p1+p2" is restrained between 0 and 10 with "%=", and p2 is then set
# to 10.
# Expected: The recipe has one restraint with a penalty of 1.
def test_restrain_operator():
    recipe, _, p2 = _make_recipe_with_variables()
    recipe %= "p1+p2", 0, 10
    p2.set_value(10)
    actual_number_of_restraints = len(recipe._restraints)
    expected_number_of_restraints = 1
    assert actual_number_of_restraints == expected_number_of_restraints
    (restraint,) = recipe._restraints
    actual_penalty = restraint.penalty()
    expected_penalty = 1
    assert actual_penalty == expected_penalty
