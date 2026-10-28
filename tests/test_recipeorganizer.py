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
"""Tests for the recipeorganizer module."""

import warnings
from operator import attrgetter

import numpy
import pytest

from diffpy.srfit import Calculator
from diffpy.srfit.equation import Equation
from diffpy.srfit.equation.builder import EquationFactory
from diffpy.srfit.fitbase.parameter import Parameter
from diffpy.srfit.fitbase.recipeorganizer import (
    RecipeContainer,
    RecipeOrganizer,
    get_equation_from_string,
)
from diffpy.srfit.fitbase.restraint import Restraint
from diffpy.utils._deprecator import deprecated


def _make_factory_with_parameters():
    """Return a factory with p1 = 1 and p2 = 2 registered, and the
    Parameters p1 to p4, of which p3 = 3 and p4 = 4 are not
    registered."""
    parameters = {
        name: Parameter(name, value)
        for name, value in (("p1", 1), ("p2", 2), ("p3", 3), ("p4", 4))
    }
    factory = EquationFactory()
    factory.register_argument("p1", parameters["p1"])
    factory.register_argument("p2", parameters["p2"])
    return factory, parameters


# Test building an Equation from a string with the factory's Parameters and
# with Parameters passed in the namespace dictionary.
@pytest.mark.parametrize(
    "input_equation_str, input_namespace_names, expected_arg_names, "
    "expected_value",
    [
        # C1: All Parameters in the string are registered with the factory.
        # Expected: The equation uses p1 and p2 and evaluates to 3.
        ("p1+p2", [], ["p1", "p2"], 3),
        # C2: p3 is not registered with the factory but is passed in the
        # namespace dictionary.
        # Expected: The equation uses p1, p2 and p3, evaluates to 6, and p3
        # is not left registered with the factory.
        ("p1+p2+p3", ["p3"], ["p1", "p2", "p3"], 6),
    ],
)
def test_get_equation_from_string(
    input_equation_str,
    input_namespace_names,
    expected_arg_names,
    expected_value,
):
    factory, parameters = _make_factory_with_parameters()
    namespace = {name: parameters[name] for name in input_namespace_names}
    equation = get_equation_from_string(input_equation_str, factory, namespace)
    actual_args = set(equation.args)
    expected_args = {parameters[name] for name in expected_arg_names}
    assert actual_args == expected_args
    actual_value = equation()
    assert actual_value == expected_value
    actual_registered_names = {"p1", "p2", "p3"}.intersection(factory.builders)
    expected_registered_names = {"p1", "p2"}
    assert actual_registered_names == expected_registered_names


# Test that building an Equation from a string fails when a Parameter is
# unknown or the namespace dictionary redefines a registered Parameter.
@pytest.mark.parametrize(
    "input_equation_str, input_namespace",
    [
        # C1: p3 is neither registered nor in the namespace dictionary.
        # Expected: Raises a ValueError.
        ("p1+p2+p3", {}),
        # C2: p3 is in the namespace dictionary, but p4 is neither
        # registered nor in the namespace dictionary.
        # Expected: Raises a ValueError.
        ("p1+p2+p3+p4", {"p3": "p3"}),
        # C3: The namespace dictionary redefines the registered p2.
        # Expected: Raises a ValueError.
        ("p1+p2", {"p2": "p4"}),
    ],
)
def test_get_equation_from_string_bad(input_equation_str, input_namespace):
    factory, parameters = _make_factory_with_parameters()
    namespace = {
        name: parameters[parameter_name]
        for name, parameter_name in input_namespace.items()
    }
    with pytest.raises(ValueError):
        get_equation_from_string(input_equation_str, factory, namespace)


def _make_container_tree():
    """Return the container m1 holding the Parameter p1 and the
    container m2, which holds the Parameter p2, along with the unmanaged
    Parameter p3."""
    m1 = RecipeContainer("m1")
    m1._containers = {}
    m1._manage(m1._containers)
    p1 = Parameter("p1", 1)
    m1._add_object(p1, m1._parameters)
    m2 = RecipeContainer("m2")
    p2 = Parameter("p2", 2)
    m2._add_object(p2, m2._parameters)
    m1._add_object(m2, m1._containers)
    p3 = Parameter("p3", 3)
    return {"m1": m1, "m2": m2, "p1": p1, "p2": p2, "p3": p3}


# Test access to the managed objects of the container m1 as attributes.
@pytest.mark.parametrize(
    "input_attribute_path, expected_object_name",
    [
        # C1: A managed Parameter is accessed as an attribute, m1.p1.
        # Expected: Returns the Parameter.
        ("p1", "p1"),
        # C2: A managed container is accessed as an attribute, m1.m2.
        # Expected: Returns the container.
        ("m2", "m2"),
        # C3: A Parameter of a managed container is accessed through the
        # attributes of both containers, m1.m2.p2.
        # Expected: Returns the Parameter.
        ("m2.p2", "p2"),
    ],
)
def test_container_attribute_access(
    input_attribute_path, expected_object_name
):
    tree = _make_container_tree()
    actual_object = attrgetter(input_attribute_path)(tree["m1"])
    expected_object = tree[expected_object_name]
    assert actual_object == expected_object


# Test access to the Parameters of a container by index.
@pytest.mark.parametrize(
    "input_container_name, expected_parameter_name",
    [
        # C1: The first Parameter of m1, which also holds a container, is
        # accessed as m1[0].
        # Expected: Returns p1.
        ("m1", "p1"),
        # C2: The first Parameter of m2 is accessed as m2[0].
        # Expected: Returns p2.
        ("m2", "p2"),
    ],
)
def test_container_index_access(input_container_name, expected_parameter_name):
    tree = _make_container_tree()
    actual_parameter = tree[input_container_name][0]
    expected_parameter = tree[expected_parameter_name]
    assert actual_parameter == expected_parameter


# Test that slicing and len consider only the Parameters of a container.
@pytest.mark.parametrize(
    "input_container_name, expected_parameter_names",
    [
        # C1: A container holding one Parameter and one container.
        # Expected: Only the Parameter is counted and sliced.
        ("m1", ["p1"]),
        # C2: A container holding one Parameter only.
        # Expected: The Parameter is counted and sliced.
        ("m2", ["p2"]),
    ],
)
def test_container_slice_and_length(
    input_container_name, expected_parameter_names
):
    tree = _make_container_tree()
    container = tree[input_container_name]
    actual_parameters = container[0:]
    expected_parameters = [tree[name] for name in expected_parameter_names]
    assert actual_parameters == expected_parameters
    actual_length = len(container)
    expected_length = len(expected_parameter_names)
    assert actual_length == expected_length


# Test finding the path from a container to a managed object.
@pytest.mark.parametrize(
    "input_container_name, input_object_name, expected_path_names",
    [
        # C1: m1 locates its own Parameter p1.
        # Expected: The path is m1, p1.
        ("m1", "p1", ["m1", "p1"]),
        # C2: m1 locates its own container m2.
        # Expected: The path is m1, m2.
        ("m1", "m2", ["m1", "m2"]),
        # C3: m2 locates its own Parameter p2.
        # Expected: The path is m2, p2.
        ("m2", "p2", ["m2", "p2"]),
        # C4: m1 locates p2, which is held by its container m2.
        # Expected: The path is m1, m2, p2.
        ("m1", "p2", ["m1", "m2", "p2"]),
        # C5: m1 or m2 locates p3, which neither of them manages.
        # Expected: The path is empty.
        ("m1", "p3", []),
        ("m2", "p3", []),
    ],
)
def test_locate_managed_object(
    input_container_name, input_object_name, expected_path_names
):
    tree = _make_container_tree()
    container = tree[input_container_name]
    actual_path = container._locate_managed_object(tree[input_object_name])
    expected_path = [tree[name] for name in expected_path_names]
    assert actual_path == expected_path


def _make_organizer():
    """Return an organizer with a managed dictionary of containers."""
    organizer = RecipeOrganizer("test")
    organizer._containers = {}
    organizer._manage(organizer._containers)
    return organizer


def test_new_parameter():
    # C1: A Parameter is created with a name not used in the organizer.
    # Expected: The new Parameter is returned and is managed under its
    # name.
    organizer = _make_organizer()
    new_parameter = organizer._new_parameter("p2", 0)
    actual_parameter = organizer.p2
    expected_parameter = new_parameter
    assert actual_parameter == expected_parameter


def test_new_parameter_bad():
    # C1: A Parameter is created with the name of a managed Parameter.
    # Expected: Raises a ValueError.
    organizer = _make_organizer()
    organizer._add_parameter(Parameter("p1", 1))
    with pytest.raises(ValueError):
        organizer._new_parameter("p1", 0)


# Test storing a Parameter and registering it with the equation factory.
@pytest.mark.parametrize(
    "input_existing_names, input_check",
    [
        # C1: The name of the Parameter is not used in the organizer.
        # Expected: The Parameter is managed and registered under its name.
        ([], True),
        # C2: A Parameter of the same name exists and the check is off.
        # Expected: The new Parameter replaces the existing one and stays
        # registered under its name.
        (["p1"], False),
    ],
)
def test_add_parameter(input_existing_names, input_check):
    organizer = _make_organizer()
    for name in input_existing_names:
        organizer._add_parameter(Parameter(name, 1))
    new_parameter = Parameter("p1", 2)
    organizer._add_parameter(new_parameter, check=input_check)
    actual_parameter = organizer.p1
    expected_parameter = new_parameter
    assert actual_parameter == expected_parameter
    actual_registered = "p1" in organizer._eqfactory.builders
    expected_registered = True
    assert actual_registered == expected_registered


# Test that a Parameter cannot be stored under a name already in use.
@pytest.mark.parametrize(
    "input_name",
    [
        # C1: A Parameter of the same name is managed.
        # Expected: Raises a ValueError.
        "p1",
        # C2: A container of the same name is managed.
        # Expected: Raises a ValueError.
        "test",
    ],
)
def test_add_parameter_bad(input_name):
    organizer = _make_organizer()
    organizer._add_parameter(Parameter("p1", 1))
    organizer._add_object(RecipeContainer("test"), organizer._containers)
    with pytest.raises(ValueError):
        organizer._add_parameter(Parameter(input_name, 0))


def test_add_parameter_keeps_order():
    # C1: p1 is added, replaced by another p1 with the check off, then xyz
    # is added.
    # Expected: The replacement keeps the position of p1, and xyz follows.
    organizer = _make_organizer()
    organizer._add_parameter(Parameter("p1", 1))
    replacement = Parameter("p1", 2)
    organizer._add_parameter(replacement, check=False)
    xyz = Parameter("xyz", 0)
    organizer._add_parameter(xyz)
    actual_names = list(organizer._parameters.keys())
    expected_names = ["p1", "xyz"]
    assert actual_names == expected_names
    actual_parameters = list(organizer._parameters.values())
    expected_parameters = [replacement, xyz]
    assert actual_parameters == expected_parameters


def test_remove_parameter():
    # C1: A managed Parameter is removed.
    # Expected: Its name is no longer registered with the equation factory.
    organizer = _make_organizer()
    parameter = Parameter("p1", 1)
    organizer._add_parameter(parameter)
    organizer._remove_parameter(parameter)
    actual_registered = "p1" in organizer._eqfactory.builders
    expected_registered = False
    assert actual_registered == expected_registered


# Test that only managed Parameters can be removed.
@pytest.mark.parametrize(
    "input_kind",
    [
        # C1: A different Parameter with the name of a managed Parameter.
        # Expected: Raises a ValueError.
        "same_name_parameter",
        # C2: A Parameter that has already been removed.
        # Expected: Raises a ValueError.
        "removed_parameter",
        # C3: A container rather than a Parameter.
        # Expected: Raises a ValueError.
        "container",
    ],
)
def test_remove_parameter_bad(input_kind):
    organizer = _make_organizer()
    parameter = Parameter("p1", 1)
    organizer._add_parameter(parameter)
    if input_kind == "removed_parameter":
        organizer._remove_parameter(parameter)
    inputs = {
        "same_name_parameter": Parameter("p1", 2),
        "removed_parameter": parameter,
        "container": RecipeContainer("test"),
    }
    with pytest.raises(ValueError):
        organizer._remove_parameter(inputs[input_kind])


def _get_constraint_state(organizer, parameter):
    """Return whether the Parameter and the organizer say the Parameter
    is constrained, the names of the constrained Parameters and the
    number of constraints."""
    return (
        parameter.constrained,
        organizer.is_constrained(parameter),
        [
            constrained.name
            for constrained in organizer.get_constrained_parameters()
        ],
        len(organizer._constraints),
    )


def test_constraint_state():
    # Test the constraint state of p1 as a constraint is added and removed.
    # Each step changes the organizer left by the step before it.
    organizer = _make_organizer()
    p1 = organizer._new_parameter("p1", 1)
    organizer._new_parameter("p2", 2)
    # C1: No constraint is added.
    # Expected: p1 is not constrained.
    actual_states = [_get_constraint_state(organizer, p1)]
    # C2: p1 is constrained to "2*p2".
    # Expected: p1 is constrained and is the only constrained Parameter.
    organizer.add_constraint(p1, "2*p2")
    actual_states.append(_get_constraint_state(organizer, p1))
    # C3: The constraint on p1 is removed.
    # Expected: p1 is not constrained.
    organizer.remove_constraint(p1)
    actual_states.append(_get_constraint_state(organizer, p1))
    # C4: p1 is constrained again, then all constraints are cleared.
    # Expected: p1 is not constrained.
    organizer.add_constraint(p1, "2*p2")
    organizer.clear_all_constraints()
    actual_states.append(_get_constraint_state(organizer, p1))
    expected_states = [
        (False, False, [], 0),
        (True, True, ["p1"], 1),
        (False, False, [], 0),
        (False, False, [], 0),
    ]
    assert actual_states == expected_states


# Test that a constrained Parameter follows its constraint.
@pytest.mark.parametrize(
    "input_constraint_kind, input_p2_value, expected_p1_value",
    [
        # C1: p1 is constrained to the equation string "2*p2" and p2 is set
        # to 10.
        # Expected: p1 becomes 20.
        ("equation_string", 10, 20),
        # C2: p1 is constrained to the Parameter p2 and p2 is set to 7.
        # Expected: p1 becomes 7.
        ("parameter", 7, 7),
    ],
)
def test_add_constraint(
    input_constraint_kind, input_p2_value, expected_p1_value
):
    organizer = _make_organizer()
    p1 = organizer._new_parameter("p1", 1)
    p2 = organizer._new_parameter("p2", 2)
    constraints = {"equation_string": "2*p2", "parameter": p2}
    organizer.add_constraint(p1, constraints[input_constraint_kind])
    p2.set_value(input_p2_value)
    organizer._constraints[p1].update()
    actual_p1_value = p1.get_value()
    assert actual_p1_value == expected_p1_value


# Test that a constraint fails on unknown or redefined Parameters.
@pytest.mark.parametrize(
    "input_constraint_eq, input_params_names",
    [
        # C1: The constraint uses p3, which is not managed.
        # Expected: Raises a ValueError.
        ("2*p3", []),
        # C2: The params dictionary redefines the managed p2.
        # Expected: Raises a ValueError.
        ("2*p2", ["p2"]),
    ],
)
def test_add_constraint_bad(input_constraint_eq, input_params_names):
    organizer = _make_organizer()
    p1 = organizer._new_parameter("p1", 1)
    organizer._new_parameter("p2", 2)
    p3 = Parameter("p3", 3)
    params = {name: p3 for name in input_params_names}
    with pytest.raises(ValueError):
        organizer.add_constraint(p1, input_constraint_eq, params)


def _make_organizer_with_registered_arguments():
    """Return an organizer whose equation factory knows p1 = 1 and
    p2 = 2."""
    organizer = _make_organizer()
    p1 = Parameter("p1", 1)
    p2 = Parameter("p2", 2)
    organizer._eqfactory.register_argument("p1", p1)
    organizer._eqfactory.register_argument("p2", p2)
    return organizer, p1, p2


# Test the penalty of a soft bound on an equation string or a Parameter.
@pytest.mark.parametrize(
    "input_restrained_kind, input_changed_name, input_changed_value",
    [
        # C1: "p1+p2" has an upper bound of 10, and p2 is set to 10.
        # Expected: One restraint with a penalty of 1.
        ("equation_string", "p2", 10),
        # C2: The Parameter p1 has an upper bound of 10, and p1 is set
        # to 11.
        # Expected: One restraint with a penalty of 1.
        ("parameter", "p1", 11),
    ],
)
def test_add_soft_bounds(
    input_restrained_kind, input_changed_name, input_changed_value
):
    organizer, p1, p2 = _make_organizer_with_registered_arguments()
    restrained = {"equation_string": "p1+p2", "parameter": p1}
    restraint = organizer.add_soft_bounds(
        restrained[input_restrained_kind], upper_bound=10
    )
    {"p1": p1, "p2": p2}[input_changed_name].set_value(input_changed_value)
    actual_number_of_restraints = len(organizer._restraints)
    expected_number_of_restraints = 1
    assert actual_number_of_restraints == expected_number_of_restraints
    actual_penalty = restraint.penalty()
    expected_penalty = 1
    assert actual_penalty == expected_penalty


def test_remove_soft_bounds():
    # C1: A soft bound is added, then removed.
    # Expected: The organizer has no restraints.
    organizer, _, _ = _make_organizer_with_registered_arguments()
    restraint = organizer.add_soft_bounds("p1+p2", upper_bound=10)
    organizer.remove_soft_bounds(restraint)
    actual_number_of_restraints = len(organizer._restraints)
    expected_number_of_restraints = 0
    assert actual_number_of_restraints == expected_number_of_restraints


# Test that a soft bound fails on unknown or redefined Parameters.
@pytest.mark.parametrize(
    "input_equation_str, input_params_names",
    [
        # C1: The equation uses p3, which is not registered.
        # Expected: Raises a ValueError.
        ("2*p3", []),
        # C2: The params dictionary redefines the registered p2.
        # Expected: Raises a ValueError.
        ("2*p2", ["p2"]),
    ],
)
def test_add_soft_bounds_bad(input_equation_str, input_params_names):
    organizer, _, _ = _make_organizer_with_registered_arguments()
    p3 = Parameter("p3", 3)
    params = {name: p3 for name in input_params_names}
    with pytest.raises(ValueError):
        organizer.add_soft_bounds(input_equation_str, params=params)


def _make_nested_organizers():
    """Return an organizer holding p1 and p2 and the organizer m2, which
    holds p3 and p4."""
    organizer = _make_organizer()
    m2 = RecipeOrganizer("m2")
    organizer._organizers = {}
    organizer._manage(organizer._organizers)
    organizer._add_object(m2, organizer._organizers)
    for name, value in (("p1", 1), ("p2", 2)):
        organizer._add_parameter(Parameter(name, value))
    for name, value in (("p3", 3), ("p4", 4)):
        m2._add_parameter(Parameter(name, value))
    return organizer, m2


def test_get_constraints():
    # C1: p1 is constrained in the organizer and p3 in the organizer m2
    # that it holds.
    # Expected: The constraints of both organizers are returned.
    organizer, m2 = _make_nested_organizers()
    organizer.add_constraint(organizer.p1, "p2")
    m2.add_constraint(m2.p3, "p4")
    actual_constrained_parameters = set(organizer._get_constraints())
    expected_constrained_parameters = {organizer.p1, m2.p3}
    assert actual_constrained_parameters == expected_constrained_parameters


def test_get_restraints():
    # C1: A restraint is added to the organizer and another to the
    # organizer m2 that it holds.
    # Expected: The restraints of both organizers are returned.
    organizer, m2 = _make_nested_organizers()
    restraint1 = organizer.add_soft_bounds("p1 + p2")
    restraint2 = m2.add_soft_bounds("2*p3 + p4")
    actual_restraints = organizer._get_restraints()
    expected_restraints = {restraint1, restraint2}
    assert actual_restraints == expected_restraints


class _GaussianCalculator(Calculator):
    """Calculator of a Gaussian with Parameters A, center and width."""

    def __init__(self, name):
        Calculator.__init__(self, name)
        self.new_parameter("A", 1.0)
        self.new_parameter("center", 0.0)
        self.new_parameter("width", 0.1)

    def __call__(self, x):
        amplitude = self.A.get_value()
        center = self.center.get_value()
        width = self.width.get_value()
        return amplitude * numpy.exp(-0.5 * ((x - center) / width) ** 2)


# Test that a registered Calculator is managed and usable in equation
# strings.
@pytest.mark.parametrize(
    "input_center",
    [
        # C1: The center of the Gaussian is set to 3 through the organizer.
        # Expected: The calculator and an equation string that uses it
        # give the Gaussian centered at 3.
        3.0,
        # C2: The center of the Gaussian is set to 5 through the organizer.
        # Expected: The calculator and an equation string that uses it
        # give the Gaussian centered at 5.
        5.0,
    ],
)
def test_register_calculator(input_center):
    organizer = _make_organizer()
    calculator = _GaussianCalculator("g")
    organizer.register_calculator(calculator)
    x = numpy.arange(0.5, 10, 0.5)
    organizer.x.set_value(x)
    organizer.g.center.set_value(input_center)
    actual_values = calculator(x)
    expected_values = numpy.exp(-0.5 * ((x - input_center) / 0.1) ** 2)
    assert numpy.array_equal(actual_values, expected_values)
    equation = organizer.register_string_function("g/x - 1", "pdf")
    actual_equation_values = equation()
    expected_equation_values = expected_values / x - 1
    assert numpy.array_equal(actual_equation_values, expected_equation_values)


def _gaussian(A, c, w, x):
    return A * numpy.exp(-0.5 * ((x - c) / w) ** 2)


def _add_one(A):
    return A + 1


def test_register_function():
    # C1: A Python function of A, c, w and x is registered as "g".
    # Expected: Its arguments become managed Parameters, the equation
    # gives the function's value, and an equation string can use it.
    organizer = _make_organizer()
    equation = organizer.register_function(_gaussian, "g")
    actual_args_are_managed = all(
        arg is organizer.get(arg.name) for arg in equation.args
    )
    expected_args_are_managed = True
    assert actual_args_are_managed == expected_args_are_managed
    x = numpy.arange(0.5, 10, 0.5)
    organizer.x.set_value(x)
    organizer.A.set_value(1.0)
    organizer.c.set_value(3.0)
    organizer.w.set_value(0.1)
    actual_values = equation()
    expected_values = numpy.exp(-0.5 * ((x - 3.0) / 0.1) ** 2)
    assert numpy.array_equal(actual_values, expected_values)
    string_equation = organizer.register_string_function("g/x - 1", "pdf")
    actual_string_values = string_equation()
    expected_string_values = expected_values / x - 1
    assert numpy.array_equal(actual_string_values, expected_string_values)


def test_register_function_replaces_existing():
    # C1: A function is registered as "g", then A + 1 is registered as "g"
    # with A = 1.
    # Expected: The equation made for the first "g" evaluates to 2.
    organizer = _make_organizer()
    equation = organizer.register_function(_gaussian, "g")
    organizer.A.set_value(1.0)
    organizer.register_function(_add_one, "g")
    actual_value = equation()
    expected_value = 2.0
    assert actual_value == pytest.approx(expected_value)


class _ValueHolder:
    def eval(self):
        return 1.23

    def __call__(self):
        return 4.56


# Test registering callables other than plain functions.
@pytest.mark.parametrize(
    "input_callable_kind, expected_value",
    [
        # C1: A bound method.
        # Expected: The equation gives the method's value.
        ("bound_method", 1.23),
        # C2: A callable object.
        # Expected: The equation gives the object's value.
        ("callable_object", 4.56),
    ],
)
def test_register_function_callables(input_callable_kind, expected_value):
    organizer = _make_organizer()
    holder = _ValueHolder()
    callables = {"bound_method": holder.eval, "callable_object": holder}
    equation = organizer.register_function(
        callables[input_callable_kind], input_callable_kind
    )
    actual_value = equation()
    assert actual_value == pytest.approx(expected_value)


# `register_function` extracts the name and the argument names from the code
# object of the function it is given. A decorator such as `deprecated`
# replaces that code object with the one of its (*args, **kwargs) wrapper, so
# a decorated function must be introspected through the wrapper to give the
# same result as the undecorated function.
@pytest.mark.parametrize(
    "input_function",
    [
        # C1: Plain undecorated function.
        # Expected: The name and the argument names come from the function.
        _gaussian,
        # C2: The same function behind a deprecation decorator.
        # Expected: The name and the argument names of the wrapped function,
        # not those of the wrapper.
        deprecated("'_gaussian' is deprecated.")(_gaussian),
    ],
)
def test_register_function_introspects_through_a_decorator(input_function):
    organizer = RecipeOrganizer("organizer")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        equation = organizer.register_function(input_function)
    actual_name = equation.name
    expected_name = "_eq__gaussian"
    actual_argnames = list(equation.argdict.keys())
    expected_argnames = ["A", "c", "w", "x"]
    assert actual_name == expected_name
    assert actual_argnames == expected_argnames


def _make_organizer_with_eq1():
    """Return an organizer with eq1 = x**2 + 3 registered, x = 0 and
    y = 3."""
    organizer = _make_organizer()
    eq1 = organizer.register_string_function("x**2 + 3", "eq1")
    eq1.x.set_value(0)
    organizer._new_parameter("y", 3.0)
    return organizer, eq1


def test_register_string_function_manages_parameters():
    # C1: "x**2 + 3" is registered as eq1 and x is set to 0 through eq1.
    # Expected: x is a managed Parameter, so the organizer sees x = 0.
    organizer, eq1 = _make_organizer_with_eq1()
    actual_args_are_managed = all(
        arg is organizer.get(arg.name) for arg in eq1.args
    )
    expected_args_are_managed = True
    assert actual_args_are_managed == expected_args_are_managed
    actual_x_value = organizer.x.get_value()
    expected_x_value = 0
    assert actual_x_value == expected_x_value


# Test registering string functions that use eq1 = x**2 + 3, with x = 0 and
# y = 3.
@pytest.mark.parametrize(
    "input_functions, expected_value",
    [
        # C1: eq2 = eq1 - 3 uses eq1 at x = 0.
        # Expected: eq2 evaluates to 0.
        ([("eq2", "eq1 - 3")], 0),
        # C2: eq3 = eq1(y) - 3 passes y = 3 as the argument of eq1.
        # Expected: eq3 evaluates to y**2 = 9.
        ([("eq3", "eq1(y) - 3")], 9.0),
        # C3: eq4 = 2*eq3 embeds eq1 two levels deep.
        # Expected: eq4 evaluates to 2*y**2 = 18.
        ([("eq3", "eq1(y) - 3"), ("eq4", "2*eq3")], 18.0),
    ],
)
def test_register_string_function(input_functions, expected_value):
    organizer, _ = _make_organizer_with_eq1()
    for name, function_str in input_functions:
        equation = organizer.register_string_function(function_str, name)
    actual_value = equation()
    assert actual_value == expected_value


def test_register_string_function_injects_arguments():
    # C1: eq2 = eq1 - 3 and eq3 = eq1(y) - 3 are registered, and eq3 is
    # evaluated before eq2.
    # Expected: eq3 sets x, the argument of eq1, to y = 3, so eq2 then
    # evaluates to 9.
    organizer, _ = _make_organizer_with_eq1()
    eq2 = organizer.register_string_function("eq1 - 3", "eq2")
    eq3 = organizer.register_string_function("eq1(y) - 3", "eq3")
    eq3()
    actual_value = eq2()
    expected_value = 9.0
    assert actual_value == expected_value


def test_evaluate_equation():
    # C1: "3 * x" is evaluated with x = 12.
    # Expected: Returns 36, and the equation factory does not keep the
    # temporary equation.
    organizer = _make_organizer()
    organizer._new_parameter("x", 12)
    actual_value = organizer.evaluate_equation("3 * x")
    expected_value = 36
    assert actual_value == expected_value
    actual_number_of_equations = len(organizer._eqfactory.equations)
    expected_number_of_equations = 0
    assert actual_number_of_equations == expected_number_of_equations


_SHOW_STAGES = ["parameters", "constraints", "restraints", "child"]


def _make_show_organizer(stage):
    """Return an organizer built up to the given stage of _SHOW_STAGES.

    Each stage adds to the previous ones: the Parameters x = 1 and
    y = 2, then z = 7 with the constraint y = 3 * z, then a soft bound
    on z, then a child organizer foo holding bar = 13.
    """
    organizer = _make_organizer()
    stages = _SHOW_STAGES[: _SHOW_STAGES.index(stage) + 1]
    if "parameters" in stages:
        organizer._new_parameter("x", 1)
        organizer._new_parameter("y", 2)
    if "constraints" in stages:
        organizer._new_parameter("z", 7)
        organizer.add_constraint("y", "3 * z")
    if "restraints" in stages:
        organizer.add_soft_bounds("z", lower_bound=2, upper_bound=3, sig=0.001)
    if "child" in stages:
        organizer._add_object(RecipeOrganizer("foo"), organizer._containers)
        organizer.foo._new_parameter("bar", 13)
    return organizer


def test_show_empty(capturestdout):
    # C1: The organizer holds nothing.
    # Expected: Nothing is printed.
    organizer = _make_organizer()
    actual_output = capturestdout(organizer.show)
    expected_output = ""
    assert actual_output == expected_output


# Test which sections show prints as the organizer is built up.
@pytest.mark.parametrize(
    "input_stage, expected_line_count, expected_sections",
    [
        # C1: The organizer holds Parameters only.
        # Expected: 4 lines with a Parameters section only.
        ("parameters", 4, {"Parameters"}),
        # C2: A constraint is added.
        # Expected: 9 lines with Parameters and Constraints sections.
        ("constraints", 9, {"Parameters", "Constraints"}),
        # C3: A restraint is added.
        # Expected: 13 lines with all three sections.
        ("restraints", 13, {"Parameters", "Constraints", "Restraints"}),
    ],
)
def test_show(
    capturestdout, input_stage, expected_line_count, expected_sections
):
    organizer = _make_show_organizer(input_stage)
    lines = capturestdout(organizer.show).strip().split("\n")
    actual_line_count = len(lines)
    assert actual_line_count == expected_line_count
    actual_sections = {
        section
        for section in ("Parameters", "Constraints", "Restraints")
        if section in lines
    }
    assert actual_sections == expected_sections


def test_show_pattern_filters_parameters(capturestdout):
    # C1: The organizer with Parameters, a constraint and a restraint is
    # shown with the pattern "x".
    # Expected: 9 lines, since only the Parameter x is listed.
    organizer = _make_show_organizer("restraints")
    lines = capturestdout(organizer.show, pattern="x").strip().split("\n")
    actual_line_count = len(lines)
    expected_line_count = 9
    assert actual_line_count == expected_line_count


def test_show_child_parameters(capturestdout):
    # C1: The organizer holds a child organizer foo with the Parameter bar.
    # Expected: bar is listed by its full name foo.bar.
    organizer = _make_show_organizer("child")
    actual_output = capturestdout(organizer.show)
    assert "foo.bar" in actual_output


# Test patterns that keep every line of the organizer without a child.
@pytest.mark.parametrize(
    "input_stage, input_pattern",
    [
        # C1: The pattern "^" matches every name.
        # Expected: The same output as without a pattern.
        ("restraints", "^"),
        # C2: The organizer also holds foo.bar, and the pattern excludes
        # names starting with foo.
        # Expected: The same output as for the organizer without foo.
        ("child", "^(?!foo)."),
    ],
)
def test_show_pattern(capturestdout, input_stage, input_pattern):
    organizer = _make_show_organizer(input_stage)
    actual_output = capturestdout(organizer.show, input_pattern)
    expected_output = capturestdout(_make_show_organizer("restraints").show)
    assert actual_output == expected_output


def test_show_added_soft_bounds(capturestdout):
    # C1: The soft bound 0 < p < 2 is added with add_soft_bounds.
    # Expected: The Restraints section lists it by the Parameter name p.
    organizer = _make_organizer()
    organizer._new_parameter("p", 1)
    organizer.add_soft_bounds("p", lower_bound=0, upper_bound=2)
    lines = capturestdout(organizer.show, textwidth=0).splitlines()
    actual_restraint_lines = lines[lines.index("Restraints") + 2 :]
    expected_restraint_lines = [
        "p: lower_bound = 0.000000, upper_bound = 2.000000, "
        "sig = 1.000000, scaled = False"
    ]
    assert actual_restraint_lines == expected_restraint_lines


def test_show_registered_soft_bounds(capturestdout):
    # C1: A Restraint 0 < p < 2 made directly from an Equation on p is
    # added with register_soft_bounds.
    # Expected: The Restraints section lists it by its equation name eq_p.
    organizer = _make_organizer()
    parameter = organizer._new_parameter("p", 1)
    restraint = Restraint(Equation(root=parameter), 0, 2)
    organizer.register_soft_bounds(restraint)
    lines = capturestdout(organizer.show, textwidth=0).splitlines()
    actual_restraint_lines = lines[lines.index("Restraints") + 2 :]
    expected_restraint_lines = [
        "eq_p: lower_bound = 0.000000, upper_bound = 2.000000, "
        "sig = 1.000000, scaled = False"
    ]
    assert actual_restraint_lines == expected_restraint_lines


def _make_iterpars_tree():
    """Build a small hierarchy for iterate_over_parameters tests."""
    root = RecipeContainer("root")
    root._containers = {}
    root._manage(root._containers)

    root_biso = Parameter("Biso", 10)
    root._add_object(root_biso, root._parameters)

    ni0 = RecipeContainer("Ni0")
    ni0_biso = Parameter("Biso", 20)
    ni0_uiso = Parameter("Uiso", 30)
    ni0._add_object(ni0_biso, ni0._parameters)
    ni0._add_object(ni0_uiso, ni0._parameters)

    ni1 = RecipeContainer("Ni1")
    ni1_biso = Parameter("Biso", 40)
    ni1._add_object(ni1_biso, ni1._parameters)

    o0 = RecipeContainer("O0")
    o0_biso = Parameter("Biso", 50)
    o0._add_object(o0_biso, o0._parameters)

    root._add_object(ni0, root._containers)
    root._add_object(ni1, root._containers)
    root._add_object(o0, root._containers)

    return {
        "root": root,
        "root_biso": root_biso,
        "ni0": ni0,
        "ni0_biso": ni0_biso,
        "ni0_uiso": ni0_uiso,
        "ni1": ni1,
        "ni1_biso": ni1_biso,
        "o0": o0,
        "o0_biso": o0_biso,
    }


# Test leaf-name and fullname matching in iterate_over_parameters.
@pytest.mark.parametrize(
    "input_pattern, input_kwargs, expected_values",
    [
        # C1: Match leaf parameter names without fullnames.
        # Expected: All Biso parameters in the hierarchy are returned.
        (r"^Biso$", {}, [10, 20, 40, 50]),
        # C2: Match hierarchical names without fullnames.
        # Expected: No leaf names match the hierarchical pattern.
        (r"^Ni\d+\.Biso$", {}, []),
        # C3: Match hierarchical names with fullnames enabled.
        # Expected: Matching Ni Biso parameters are returned.
        (r"^Ni\d+\.Biso$", {"fullnames": True}, [20, 40]),
        # C4: Match one hierarchical Uiso name.
        # Expected: Only Ni0.Uiso is returned.
        (r"^Ni0\.Uiso$", {"fullnames": True}, [30]),
        # C5: Match one hierarchical Biso name outside Ni containers.
        # Expected: Only O0.Biso is returned.
        (r"^O0\.Biso$", {"fullnames": True}, [50]),
        # C6: Disable recursion while matching child fullnames.
        # Expected: No child parameters are returned.
        (r"^Ni\d+\.Biso$", {"fullnames": True, "recurse": False}, []),
        # C7: Disable recursion while matching root fullname.
        # Expected: Only the root-level Biso parameter is returned.
        (r"^Biso$", {"fullnames": True, "recurse": False}, [10]),
    ],
)
def test_iterate_over_parameters_fullname_matching(
    input_pattern, input_kwargs, expected_values
):
    root = _make_iterpars_tree()["root"]
    actual_values = [
        parameter.value
        for parameter in root.iterate_over_parameters(
            input_pattern, **input_kwargs
        )
    ]
    assert actual_values == expected_values


# Test that fullnames are relative to the container that is iterated over.
@pytest.mark.parametrize(
    "input_pattern, expected_names",
    [
        # C1: Match Biso relative to the called Ni0 container.
        # Expected: Ni0.Biso is returned without the Ni0 prefix.
        (r"^Biso$", ["Biso"]),
        # C2: Match Uiso relative to the called Ni0 container.
        # Expected: Ni0.Uiso is returned without the Ni0 prefix.
        (r"^Uiso$", ["Uiso"]),
        # C3: Match with the parent container prefix from inside Ni0.
        # Expected: No parameter is returned because fullnames are relative
        # to the container on which iterate_over_parameters is called.
        (r"^Ni0\.Biso$", []),
    ],
)
def test_iterate_over_parameters_fullnames_are_relative_to_called_container(
    input_pattern, expected_names
):
    ni0 = _make_iterpars_tree()["ni0"]
    actual_names = [
        parameter.name
        for parameter in ni0.iterate_over_parameters(
            input_pattern, fullnames=True
        )
    ]
    assert actual_names == expected_names
