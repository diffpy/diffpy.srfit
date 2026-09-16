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
"""Tests for refinableobj module."""

import re
import unittest

import pytest

from diffpy.srfit.fitbase.parameter import Parameter
from diffpy.srfit.fitbase.parameterset import ParameterSet


class TestParameterSet(unittest.TestCase):

    def setUp(self):
        self.parset = ParameterSet("test")
        return

    def test_add_parameter_set(self):
        """Test the add_parameter_set method."""
        parset2 = ParameterSet("parset2")
        p1 = Parameter("parset2", 1)

        self.parset.add_parameter_set(parset2)
        self.assertTrue(self.parset.parset2 is parset2)

        self.assertRaises(ValueError, self.parset.add_parameter_set, p1)

        p1.name = "p1"
        parset2.add_parameter(p1)

        self.assertTrue(self.parset.parset2.p1 is p1)

        return


# ----------------------------------------------------------------------------
# addParameterSet is deprecated in favor of add_parameter_set. The old name
# must still work and forward to the new implementation.


class TestParameterSetDeprecated(unittest.TestCase):

    def setUp(self):
        self.parset = ParameterSet("test")
        return

    def test_add_parameter_set_deprecated(self):
        """Test the deprecated addParameterSet method.

        Remove this test after the addParameterSet is removed in version
        4.0.0.
        """
        parset2 = ParameterSet("parset2")
        p1 = Parameter("parset2", 1)

        self.parset.addParameterSet(parset2)
        self.assertTrue(self.parset.parset2 is parset2)

        self.assertRaises(ValueError, self.parset.add_parameter_set, p1)

        p1.name = "p1"
        parset2.add_parameter(p1)

        self.assertTrue(self.parset.parset2.p1 is p1)

        return


# ----------------------------------------------------------------------------
# The camelCase Parameter accessors on ParameterSet are deprecated in favour of
# their snake_case spellings. Each old name must still work, warn with a
# message that names its replacement, and have the same effect on the set.


@pytest.mark.parametrize(
    "deprecated_name, replacement_name",
    [
        # C1: A Parameter is stored under the set.
        # Expected: addParameter warns and stores it as add_parameter does.
        ("addParameter", "add_parameter"),
        # C2: A Parameter is created and stored in one call.
        # Expected: newParameter warns and creates it as new_parameter does.
        ("newParameter", "new_parameter"),
        # C3: A stored Parameter is removed from the set.
        # Expected: removeParameter warns and removes it as
        # remove_parameter does.
        ("removeParameter", "remove_parameter"),
    ],
)
def test_parameter_accessors_warn_and_forward(
    deprecated_name, replacement_name
):
    base = "diffpy.srfit.fitbase.parameterset.ParameterSet"
    expected_msg = (
        f"'{base}.{deprecated_name}' is deprecated and will be removed in "
        f"version 4.0.0. Please use '{base}.{replacement_name}' instead."
    )
    parset = ParameterSet("test")
    if deprecated_name == "newParameter":
        arguments = ("p1", 1)
    elif deprecated_name == "removeParameter":
        arguments = (parset.new_parameter("p1", 1),)
    else:
        arguments = (Parameter("p1", 1),)

    with pytest.warns(DeprecationWarning, match=re.escape(expected_msg)):
        getattr(parset, deprecated_name)(*arguments)

    expected_names = [] if deprecated_name == "removeParameter" else ["p1"]
    actual_names = [par.name for par in parset._parameters.values()]
    assert actual_names == expected_names


if __name__ == "__main__":
    unittest.main()
