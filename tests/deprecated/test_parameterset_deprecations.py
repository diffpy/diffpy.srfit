"""Tests for the deprecated camelCase names in
diffpy.srfit.fitbase.parameterset.

Each test calls the deprecated name, checks that it warns, and checks
that it gives the same result as the new name on an identical object.

Remove this file in version 4.0.0, together with the deprecated names.
"""

import re

import pytest

from diffpy.srfit.fitbase.parameter import Parameter
from diffpy.srfit.fitbase.parameterset import ParameterSet


def expected_message(base, old_name, new_name):
    return re.escape(
        f"'{base}.{old_name}' is deprecated and will be removed in version "
        f"4.0.0. Please use '{base}.{new_name}' instead."
    )


PARAMETER_SET = "diffpy.srfit.fitbase.parameterset.ParameterSet"


def _make_call_args(parset, method_name):
    """Return fresh arguments for a ParameterSet method."""
    if method_name == "add_parameter":
        return (Parameter("p2", 2),)
    if method_name == "new_parameter":
        return ("p2", 2)
    if method_name == "remove_parameter":
        return (parset.p1,)
    if method_name == "add_parameter_set":
        return (ParameterSet("other"),)
    if method_name == "remove_parameter_set":
        return (parset.child,)
    return ()


def _call_parameter_set_method(method_name, args_method_name):
    """Call a method of the set "test" that holds the Parameter p1 = 1
    and the ParameterSet "child".

    Return what the method returned, and the names, values and constant
    flags of the Parameters and the names of the ParameterSets the set
    holds afterward.
    """
    parset = ParameterSet("test")
    parset.new_parameter("p1", 1)
    parset.add_parameter_set(ParameterSet("child"))
    args = _make_call_args(parset, args_method_name)
    result = getattr(parset, method_name)(*args)
    returned = None if result is None else (result.name, result.value)
    return (
        returned,
        parset.get_names(),
        parset.get_values(),
        [par.const for par in parset.iterate_over_parameters()],
        list(parset._parsets),
    )


# ParameterSet: each deprecated method does the same as its new name, on a
# set that holds the Parameter p1 = 1 and the ParameterSet "child". See
# test_parameterset.py.
@pytest.mark.parametrize(
    "input_old_name, input_new_name",
    [
        # C1: addParameter stores the Parameter p2 = 2.
        # Expected: Warns and gives the same result as add_parameter.
        ("addParameter", "add_parameter"),
        # C2: newParameter creates the Parameter p2 = 2.
        # Expected: Warns and gives the same result as new_parameter.
        ("newParameter", "new_parameter"),
        # C3: removeParameter removes p1.
        # Expected: Warns and gives the same result as remove_parameter.
        ("removeParameter", "remove_parameter"),
        # C4: addParameterSet stores the ParameterSet "other".
        # Expected: Warns and gives the same result as add_parameter_set.
        ("addParameterSet", "add_parameter_set"),
        # C5: removeParameterSet removes "child".
        # Expected: Warns and gives the same result as remove_parameter_set.
        ("removeParameterSet", "remove_parameter_set"),
        # C6: setConst marks every Parameter in the set as constant.
        # Expected: Warns and gives the same result as set_constant.
        ("setConst", "set_constant"),
    ],
)
def test_parameter_set_deprecations(input_old_name, input_new_name):
    with pytest.warns(
        DeprecationWarning,
        match=expected_message(PARAMETER_SET, input_old_name, input_new_name),
    ):
        actual_old_result = _call_parameter_set_method(
            input_old_name, input_new_name
        )
    actual_new_result = _call_parameter_set_method(
        input_new_name, input_new_name
    )
    assert actual_old_result == actual_new_result
