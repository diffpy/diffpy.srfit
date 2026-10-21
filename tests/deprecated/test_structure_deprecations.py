"""Tests for the deprecation of diffpy.srfit.structure, which moved to
diffpy.cmistructure.

Remove this file in version 4.0.0, together with diffpy.srfit.structure.
"""

import importlib

import pytest

import diffpy.srfit.structure
from diffpy.srfit.structure import struToParameterSet
from diffpy.srfit.structure.diffpyparset import DiffpyStructureParSet
from diffpy.structure import Atom, Lattice, Structure


@pytest.fixture
def nacl_structure():
    lattice = Lattice(5.6, 5.6, 5.6, 90, 90, 90)
    atoms = [Atom("Na", [0, 0, 0]), Atom("Cl", [0.5, 0.5, 0.5])]
    return Structure(atoms, lattice=lattice)


def test_structure_package_import_warns():
    # Case: diffpy.srfit.structure is imported.
    # Expected: A DeprecationWarning points to diffpy.cmistructure.
    expected_msg = (
        "'diffpy.srfit.structure' is deprecated and will be removed in "
        "version 4.0.0. Please use 'diffpy.cmistructure' instead."
    )
    with pytest.warns(DeprecationWarning) as record:
        importlib.reload(diffpy.srfit.structure)
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg


def test_struToParameterSet(nacl_structure):
    # Case: struToParameterSet adapts a diffpy Structure.
    # Expected: A DeprecationWarning points to
    # diffpy.cmistructure.structure_to_parameter_set, and a
    # DiffpyStructureParSet is still returned.
    expected_msg = (
        "'diffpy.srfit.structure.struToParameterSet' is deprecated and will "
        "be removed in version 4.0.0. Please use "
        "'diffpy.cmistructure.structure_to_parameter_set' instead."
    )
    with pytest.warns(DeprecationWarning) as record:
        parameter_set = struToParameterSet("nacl", nacl_structure)
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg
    actual_type = type(parameter_set)
    expected_type = DiffpyStructureParSet
    assert actual_type == expected_type
    actual_structure = parameter_set.structure
    expected_structure = nacl_structure
    assert actual_structure == expected_structure


def test_restrainBVS(diffpy_srreal_available, nacl_structure):
    # Case: restrainBVS is called on a DiffpyStructureParSet.
    # Expected: A DeprecationWarning points to
    # diffpy.cmistructure's restrain_bvs, and the restraint is still added.
    if not diffpy_srreal_available:
        pytest.skip("diffpy.srreal package not available")
    parameter_set = DiffpyStructureParSet("nacl", nacl_structure)
    expected_msg = (
        "'diffpy.srfit.structure.srrealparset.SrRealParSet.restrainBVS' is "
        "deprecated and will be removed in version 4.0.0. Please use "
        "'diffpy.cmistructure.srrealparset.SrRealParSet.restrain_bvs' "
        "instead."
    )
    with pytest.warns(DeprecationWarning) as record:
        restraint = parameter_set.restrainBVS(sig=0.5)
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg
    actual_sig = restraint.sig
    expected_sig = 0.5
    assert actual_sig == expected_sig


def test_stru_get(nacl_structure):
    # Case: .stru is read from a DiffpyStructureParSet.
    # Expected: A DeprecationWarning points to .structure in
    # diffpy.cmistructure, and the adapted structure is returned.
    parameter_set = DiffpyStructureParSet("nacl", nacl_structure)
    expected_msg = (
        "'diffpy.srfit.structure.basestructureparset.BaseStructureParSet."
        "stru' is deprecated and will be removed in version 4.0.0. Please "
        "use 'diffpy.cmistructure.basestructureparset.BaseStructureParSet."
        "structure' instead."
    )
    with pytest.warns(DeprecationWarning) as record:
        actual_structure = parameter_set.stru
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg
    expected_structure = nacl_structure
    assert actual_structure == expected_structure


def test_stru_set(nacl_structure):
    # Case: .stru is set on a DiffpyStructureParSet.
    # Expected: A DeprecationWarning points to .structure in
    # diffpy.cmistructure, and .structure is updated.
    parameter_set = DiffpyStructureParSet("nacl", Structure())
    expected_msg = (
        "'diffpy.srfit.structure.basestructureparset.BaseStructureParSet."
        "stru' is deprecated and will be removed in version 4.0.0. Please "
        "use 'diffpy.cmistructure.basestructureparset.BaseStructureParSet."
        "structure' instead."
    )
    with pytest.warns(DeprecationWarning) as record:
        parameter_set.stru = nacl_structure
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg
    actual_structure = parameter_set.structure
    expected_structure = nacl_structure
    assert actual_structure == expected_structure
