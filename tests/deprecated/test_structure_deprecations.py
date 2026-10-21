"""Tests for the deprecation of diffpy.srfit.structure, which moved to
diffpy.cmistructure.

Each test calls a deprecated name, checks that it warns with a message
pointing to diffpy.cmistructure, and checks that it gives the same
result as the new name in srfit. The last tests check that every name a
message points to exists in diffpy.cmistructure.

Remove this file in version 4.0.0, together with diffpy.srfit.structure.
"""

import importlib
import re

import pytest

import diffpy.srfit.structure
from diffpy.srfit.structure import (
    structure_to_parameter_set,
    struToParameterSet,
    struToParameterSet_dep_msg,
)
from diffpy.srfit.structure.basestructureparset import stru_dep_msg
from diffpy.srfit.structure.diffpyparset import DiffpyStructureParSet
from diffpy.srfit.structure.srrealparset import restrainBVS_dep_msg
from diffpy.structure import Atom, Lattice, Structure


def build_expected_message(old_base, old_name, new_base, new_name):
    return re.escape(
        f"'{old_base}.{old_name}' is deprecated and will be removed in "
        f"version 4.0.0. Please use '{new_base}.{new_name}' instead."
    )


def get_new_name(message):
    """Return the dotted name that a deprecation message points to."""
    new_name = re.search(r"Please use '(.+)' instead\.", message).group(1)
    return new_name


def resolve_dotted_name(dotted_name):
    """Import the longest importable module prefix of a dotted name and
    look up the rest as attributes.

    Return None if any part is missing.
    """
    parts = dotted_name.split(".")
    for split in range(len(parts), 0, -1):
        try:
            resolved = importlib.import_module(".".join(parts[:split]))
        except ImportError:
            continue
        for attribute in parts[split:]:
            resolved = getattr(resolved, attribute, None)
            if resolved is None:
                return None
        return resolved
    return None


def make_carbon_structure():
    structure = Structure([Atom("C", [0, 0.2, 0.5])])
    return structure


def make_nacl_parameter_set():
    lattice = Lattice(5.6, 5.6, 5.6, 90, 90, 90)
    atoms = [Atom("Na", [0, 0, 0]), Atom("Cl", [0.5, 0.5, 0.5])]
    parameter_set = DiffpyStructureParSet(
        "nacl", Structure(atoms, lattice=lattice)
    )
    return parameter_set


def get_bvs_restraint_state(method_name):
    """Restrain the BVS of a fresh NaCl parameter set with the named
    method.

    Return the restraint's type, sig and scaled flag, and how many
    restraints the parameter set holds afterward.
    """
    parameter_set = make_nacl_parameter_set()
    restraint = getattr(parameter_set, method_name)(sig=0.5, scaled=True)
    return (
        type(restraint).__name__,
        restraint.sig,
        restraint.scaled,
        len(parameter_set._restraints),
    )


def test_structure_package_import_warns():
    # Test the warning given when diffpy.srfit.structure is imported.
    # C1: The diffpy.srfit.structure package is imported.
    # Expected: A DeprecationWarning says it moved to diffpy.cmistructure.
    expected_msg = (
        "'diffpy.srfit.structure' is deprecated and will be removed in "
        "version 4.0.0. Please use 'diffpy.cmistructure' instead."
    )
    with pytest.warns(DeprecationWarning, match=re.escape(expected_msg)):
        importlib.reload(diffpy.srfit.structure)


def test_stru_to_parameter_set():
    # Test the deprecated diffpy.srfit.structure.struToParameterSet.
    # C1: struToParameterSet adapts a diffpy Structure.
    # Expected: Warns and gives the same result as
    # structure_to_parameter_set.
    structure = make_carbon_structure()
    expected_msg = build_expected_message(
        "diffpy.srfit.structure",
        "struToParameterSet",
        "diffpy.cmistructure",
        "structure_to_parameter_set",
    )
    with pytest.warns(DeprecationWarning, match=expected_msg):
        old_parameter_set = struToParameterSet("carbon", structure)
    new_parameter_set = structure_to_parameter_set("carbon", structure)
    actual_old_result = (
        type(old_parameter_set),
        old_parameter_set.name,
        old_parameter_set.structure,
    )
    actual_new_result = (
        type(new_parameter_set),
        new_parameter_set.name,
        new_parameter_set.structure,
    )
    assert actual_old_result == actual_new_result


def test_restrain_bvs(diffpy_srreal_available):
    # Test the deprecated SrRealParSet.restrainBVS.
    # C1: restrainBVS adds a scaled BVS restraint with sig = 0.5 to a
    # NaCl parameter set.
    # Expected: Warns and gives the same result as restrain_bvs.
    if not diffpy_srreal_available:
        pytest.skip("diffpy.srreal package not available")
    expected_msg = build_expected_message(
        "diffpy.srfit.structure.srrealparset.SrRealParSet",
        "restrainBVS",
        "diffpy.cmistructure.srrealparset.SrRealParSet",
        "restrain_bvs",
    )
    with pytest.warns(DeprecationWarning, match=expected_msg):
        actual_old_result = get_bvs_restraint_state("restrainBVS")
    actual_new_result = get_bvs_restraint_state("restrain_bvs")
    assert actual_old_result == actual_new_result


def test_stru_attribute_get():
    # Test reading the deprecated BaseStructureParSet.stru attribute.
    # C1: .stru is read from a DiffpyStructureParSet.
    # Expected: Warns and gives the same object as .structure.
    structure = make_carbon_structure()
    parameter_set = DiffpyStructureParSet("carbon", structure)
    expected_msg = build_expected_message(
        "diffpy.srfit.structure.basestructureparset.BaseStructureParSet",
        "stru",
        "diffpy.cmistructure.basestructureparset.BaseStructureParSet",
        "structure",
    )
    with pytest.warns(DeprecationWarning, match=expected_msg):
        actual_old_result = parameter_set.stru
    actual_new_result = parameter_set.structure
    assert actual_old_result == actual_new_result


def test_stru_attribute_set():
    # Test setting the deprecated BaseStructureParSet.stru attribute.
    # C1: .stru is set on a DiffpyStructureParSet.
    # Expected: Warns and sets .structure.
    parameter_set = DiffpyStructureParSet("carbon", make_carbon_structure())
    new_structure = make_carbon_structure()
    expected_msg = build_expected_message(
        "diffpy.srfit.structure.basestructureparset.BaseStructureParSet",
        "stru",
        "diffpy.cmistructure.basestructureparset.BaseStructureParSet",
        "structure",
    )
    with pytest.warns(DeprecationWarning, match=expected_msg):
        parameter_set.stru = new_structure
    actual_structure = parameter_set.structure
    expected_structure = new_structure
    assert actual_structure == expected_structure


@pytest.mark.parametrize(
    "input_message",
    [
        # Test that each deprecation message points to a function or
        # method that exists in diffpy.cmistructure.
        # C1: The message of struToParameterSet.
        # Expected: diffpy.cmistructure.structure_to_parameter_set exists.
        struToParameterSet_dep_msg,
        # C2: The message of SrRealParSet.restrainBVS.
        # Expected: diffpy.cmistructure.srrealparset.SrRealParSet.restrain_bvs
        # exists.
        restrainBVS_dep_msg,
    ],
)
def test_new_callable_exists_in_cmistructure(
    diffpy_cmistructure_available, input_message
):
    if not diffpy_cmistructure_available:
        pytest.skip("diffpy.cmistructure package not available")
    new_name = get_new_name(input_message)
    actual_is_callable = callable(resolve_dotted_name(new_name))
    expected_is_callable = True
    assert actual_is_callable == expected_is_callable


def test_new_attribute_exists_in_cmistructure(diffpy_cmistructure_available):
    # Test that the .stru message points to an attribute that exists in
    # diffpy.cmistructure.
    # C1: A diffpy.cmistructure DiffpyStructureParSet adapts a structure.
    # Expected: It derives from the class the message names, and that
    # class's attribute holds the adapted structure.
    if not diffpy_cmistructure_available:
        pytest.skip("diffpy.cmistructure package not available")
    from diffpy.cmistructure.diffpyparset import (
        DiffpyStructureParSet as CmiDiffpyStructureParSet,
    )

    new_class_name, new_attribute_name = get_new_name(stru_dep_msg).rsplit(
        ".", 1
    )
    structure = make_carbon_structure()
    parameter_set = CmiDiffpyStructureParSet("carbon", structure)
    actual_result = (
        issubclass(type(parameter_set), resolve_dotted_name(new_class_name)),
        getattr(parameter_set, new_attribute_name),
    )
    expected_result = (True, structure)
    assert actual_result == expected_result
