"""Tests for the deprecation of diffpy.srfit.structure, which moved to
diffpy.cmistructure.

Remove this file in version 4.0.0, together with diffpy.srfit.structure.
"""

import importlib
from types import SimpleNamespace

import pytest

import diffpy.srfit.structure
from diffpy.srfit.fitbase.parameter import Parameter
from diffpy.srfit.structure import constrain_as_space_group, struToParameterSet
from diffpy.srfit.structure.cctbxparset import CCTBXScattererParSet
from diffpy.srfit.structure.diffpyparset import DiffpyStructureParSet
from diffpy.srfit.structure.sgconstraints import BaseSpaceGroupParameters
from diffpy.structure import Atom, Lattice, Structure


@pytest.fixture
def nacl_structure():
    lattice = Lattice(5.6, 5.6, 5.6, 90, 90, 90)
    atoms = [Atom("Na", [0, 0, 0]), Atom("Cl", [0.5, 0.5, 0.5])]
    return Structure(atoms, lattice=lattice)


class CCTBXStructureStandIn:
    """Stand in for a cctbx structure that holds one scatterer."""

    def __init__(self):
        self._scatterers = [
            SimpleNamespace(site=(0.0, 0.0, 0.0), occupancy=1.0, u_iso=0.01)
        ]

    def scatterers(self):
        return self._scatterers


def make_crystal_parameter_set():
    """Make an ObjCrystCrystalParSet holding a four-atom molecule."""
    from pyobjcryst.crystal import Crystal
    from pyobjcryst.molecule import Molecule
    from pyobjcryst.scatteringpower import ScatteringPowerAtom

    from diffpy.srfit.structure.objcrystparset import ObjCrystCrystalParSet

    crystal = Crystal(10, 10, 10, "P1")
    molecule = Molecule(crystal, "molecule")
    crystal.AddScatterer(molecule)
    scattering_power = ScatteringPowerAtom("C", "C")
    atom_positions = [(0, 0, 0), (1.5, 0, 0), (1.5, 1.5, 0), (1.5, 1.5, 1.5)]
    for index, (x, y, z) in enumerate(atom_positions):
        molecule.AddAtom(x, y, z, scattering_power, f"C{index}")
    crystal_parameter_set = ObjCrystCrystalParSet("crystal", crystal)
    return crystal_parameter_set


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


@pytest.mark.parametrize(
    "old_name, new_name, expected_msg",
    [
        # Test the deprecated getters of the structure parameter sets.
        # C1: getLattice is called on a DiffpyStructureParSet.
        # Expected: Warns that get_lattice in diffpy.cmistructure replaces
        # it, and returns the same lattice as get_lattice.
        (
            "getLattice",
            "get_lattice",
            "'diffpy.srfit.structure.basestructureparset.BaseStructureParSet."
            "getLattice' is deprecated and will be removed in version 4.0.0. "
            "Please use 'diffpy.cmistructure.basestructureparset."
            "BaseStructureParSet.get_lattice' instead.",
        ),
        # C2: getScatterers is called on a DiffpyStructureParSet.
        # Expected: Warns that get_scatterers in diffpy.cmistructure
        # replaces it, and returns the same scatterers as get_scatterers.
        (
            "getScatterers",
            "get_scatterers",
            "'diffpy.srfit.structure.basestructureparset.BaseStructureParSet."
            "getScatterers' is deprecated and will be removed in version "
            "4.0.0. Please use 'diffpy.cmistructure.basestructureparset."
            "BaseStructureParSet.get_scatterers' instead.",
        ),
        # C3: usingSymmetry is called on a DiffpyStructureParSet.
        # Expected: Warns that using_symmetry in diffpy.cmistructure
        # replaces it, and returns the same flag as using_symmetry.
        (
            "usingSymmetry",
            "using_symmetry",
            "'diffpy.srfit.structure.srrealparset.SrRealParSet.usingSymmetry' "
            "is deprecated and will be removed in version 4.0.0. Please use "
            "'diffpy.cmistructure.srrealparset.SrRealParSet.using_symmetry' "
            "instead.",
        ),
    ],
)
def test_structure_parameter_set_getters(
    nacl_structure, old_name, new_name, expected_msg
):
    parameter_set = DiffpyStructureParSet("nacl", nacl_structure)
    with pytest.warns(DeprecationWarning) as record:
        actual_result = getattr(parameter_set, old_name)()
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg
    expected_result = getattr(parameter_set, new_name)()
    assert actual_result == expected_result


def test_canAdapt(nacl_structure):
    # Case: canAdapt is called on DiffpyStructureParSet with a diffpy
    # Structure.
    # Expected: A DeprecationWarning points to can_adapt in
    # diffpy.cmistructure, and True is returned.
    expected_msg = (
        "'diffpy.srfit.structure.basestructureparset.BaseStructureParSet."
        "canAdapt' is deprecated and will be removed in version 4.0.0. "
        "Please use 'diffpy.cmistructure.basestructureparset."
        "BaseStructureParSet.can_adapt' instead."
    )
    with pytest.warns(DeprecationWarning) as record:
        actual_result = DiffpyStructureParSet.canAdapt(nacl_structure)
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg
    expected_result = True
    assert actual_result == expected_result


def test_useSymmetry(nacl_structure):
    # Case: useSymmetry(False) is called on a DiffpyStructureParSet.
    # Expected: A DeprecationWarning points to use_symmetry in
    # diffpy.cmistructure, and using_symmetry returns False.
    parameter_set = DiffpyStructureParSet("nacl", nacl_structure)
    expected_msg = (
        "'diffpy.srfit.structure.srrealparset.SrRealParSet.useSymmetry' is "
        "deprecated and will be removed in version 4.0.0. Please use "
        "'diffpy.cmistructure.srrealparset.SrRealParSet.use_symmetry' "
        "instead."
    )
    with pytest.warns(DeprecationWarning) as record:
        parameter_set.useSymmetry(False)
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg
    actual_using_symmetry = parameter_set.using_symmetry()
    expected_using_symmetry = False
    assert actual_using_symmetry == expected_using_symmetry


def test_angunits_diffpy_lattice(nacl_structure):
    # Case: .angunits is read and then set on a DiffpyLatticeParSet.
    # Expected: Both warn that .angle_units in diffpy.cmistructure
    # replaces it. The read returns "deg", and the set updates
    # .angle_units.
    lattice = DiffpyStructureParSet("nacl", nacl_structure).lattice
    expected_msg = (
        "'diffpy.srfit.structure.diffpyparset.DiffpyLatticeParSet.angunits' "
        "is deprecated and will be removed in version 4.0.0. Please use "
        "'diffpy.cmistructure.diffpyparset.DiffpyLatticeParSet.angle_units' "
        "instead."
    )
    with pytest.warns(DeprecationWarning) as record:
        actual_angle_units = lattice.angunits
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg
    expected_angle_units = "deg"
    assert actual_angle_units == expected_angle_units
    with pytest.warns(DeprecationWarning) as record:
        lattice.angunits = "rad"
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg
    actual_angle_units = lattice.angle_units
    expected_angle_units = "rad"
    assert actual_angle_units == expected_angle_units


def test_sg(nacl_structure):
    # Case: .sg is read and then set on SpaceGroupParameters.
    # Expected: Both warn that .space_group in diffpy.cmistructure
    # replaces it. The read returns the space group the parameters were
    # made with, and the set updates .space_group.
    parameter_set = DiffpyStructureParSet("nacl", nacl_structure)
    space_group_parameters = constrain_as_space_group(parameter_set, "Fm-3m")
    expected_msg = (
        "'diffpy.srfit.structure.sgconstraints.SpaceGroupParameters.sg' is "
        "deprecated and will be removed in version 4.0.0. Please use "
        "'diffpy.cmistructure.sgconstraints.SpaceGroupParameters."
        "space_group' instead."
    )
    with pytest.warns(DeprecationWarning) as record:
        actual_space_group = space_group_parameters.sg
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg
    expected_space_group = space_group_parameters.space_group
    assert actual_space_group == expected_space_group
    with pytest.warns(DeprecationWarning) as record:
        space_group_parameters.sg = "P1"
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg
    actual_space_group = space_group_parameters.space_group
    expected_space_group = "P1"
    assert actual_space_group == expected_space_group


def test_addParameter():
    # Case: addParameter stores a Parameter in BaseSpaceGroupParameters.
    # Expected: A DeprecationWarning points to add_parameter in
    # diffpy.cmistructure, and the Parameter is stored.
    space_group_parameters = BaseSpaceGroupParameters()
    parameter = Parameter("a", 4.0)
    expected_msg = (
        "'diffpy.srfit.structure.sgconstraints.BaseSpaceGroupParameters."
        "addParameter' is deprecated and will be removed in version 4.0.0. "
        "Please use 'diffpy.cmistructure.sgconstraints."
        "BaseSpaceGroupParameters.add_parameter' instead."
    )
    with pytest.warns(DeprecationWarning) as record:
        space_group_parameters.addParameter(parameter)
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg
    actual_parameter = space_group_parameters.a
    expected_parameter = parameter
    assert actual_parameter == expected_parameter


@pytest.mark.parametrize(
    "old_name, new_name, input_value, expected_msg",
    [
        # Test the deprecated attributes of CCTBXScattererParSet.
        # C1: .strups is read and then set.
        # Expected: Both warn that .structure_parameter_set in
        # diffpy.cmistructure replaces it, and they read and set the same
        # value as .structure_parameter_set.
        (
            "strups",
            "structure_parameter_set",
            "new crystal",
            "'diffpy.srfit.structure.cctbxparset.CCTBXScattererParSet."
            "strups' is deprecated and will be removed in version 4.0.0. "
            "Please use 'diffpy.cmistructure.cctbxparset."
            "CCTBXScattererParSet.structure_parameter_set' instead.",
        ),
        # C2: .idx is read and then set.
        # Expected: Both warn that .index in diffpy.cmistructure replaces
        # it, and they read and set the same value as .index.
        (
            "idx",
            "index",
            3,
            "'diffpy.srfit.structure.cctbxparset.CCTBXScattererParSet.idx' "
            "is deprecated and will be removed in version 4.0.0. Please use "
            "'diffpy.cmistructure.cctbxparset.CCTBXScattererParSet.index' "
            "instead.",
        ),
    ],
)
def test_cctbx_scatterer_attributes(
    old_name, new_name, input_value, expected_msg
):
    crystal = SimpleNamespace(structure=CCTBXStructureStandIn())
    scatterer = CCTBXScattererParSet("atom", crystal, 0)
    with pytest.warns(DeprecationWarning) as record:
        actual_value = getattr(scatterer, old_name)
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg
    expected_value = getattr(scatterer, new_name)
    assert actual_value == expected_value
    with pytest.warns(DeprecationWarning) as record:
        setattr(scatterer, old_name, input_value)
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg
    actual_value = getattr(scatterer, new_name)
    expected_value = input_value
    assert actual_value == expected_value


@pytest.mark.parametrize(
    "old_name, new_name, input_args, expected_msg",
    [
        # Test the deprecated methods of ObjCrystMoleculeParSet that
        # return a value or change nothing.
        # C1: getScatterers is called.
        # Expected: Warns that get_scatterers in diffpy.cmistructure
        # replaces it, and returns the same atoms as get_scatterers.
        (
            "getScatterers",
            "get_scatterers",
            (),
            "'diffpy.srfit.structure.objcrystparset.ObjCrystMoleculeParSet."
            "getScatterers' is deprecated and will be removed in version "
            "4.0.0. Please use 'diffpy.cmistructure.objcrystparset."
            "ObjCrystMoleculeParSet.get_scatterers' instead.",
        ),
        # C2: usingSymmetry is called.
        # Expected: Warns that using_symmetry in diffpy.cmistructure
        # replaces it, and returns False, as using_symmetry does.
        (
            "usingSymmetry",
            "using_symmetry",
            (),
            "'diffpy.srfit.structure.objcrystparset.ObjCrystMoleculeParSet."
            "usingSymmetry' is deprecated and will be removed in version "
            "4.0.0. Please use 'diffpy.cmistructure.objcrystparset."
            "ObjCrystMoleculeParSet.using_symmetry' instead.",
        ),
        # C3: useSymmetry(False) is called.
        # Expected: Warns that use_symmetry in diffpy.cmistructure
        # replaces it, and returns None, as use_symmetry does.
        (
            "useSymmetry",
            "use_symmetry",
            (False,),
            "'diffpy.srfit.structure.objcrystparset.ObjCrystMoleculeParSet."
            "useSymmetry' is deprecated and will be removed in version "
            "4.0.0. Please use 'diffpy.cmistructure.objcrystparset."
            "ObjCrystMoleculeParSet.use_symmetry' instead.",
        ),
        # C4: hasScatterers is called.
        # Expected: Warns that has_scatterers in diffpy.cmistructure
        # replaces it, and returns True, as has_scatterers does.
        (
            "hasScatterers",
            "has_scatterers",
            (),
            "'diffpy.srfit.structure.objcrystparset.ObjCrystScattererParSet."
            "hasScatterers' is deprecated and will be removed in version "
            "4.0.0. Please use 'diffpy.cmistructure.objcrystparset."
            "ObjCrystScattererParSet.has_scatterers' instead.",
        ),
        # C5: isDummy is called.
        # Expected: Warns that is_dummy in diffpy.cmistructure replaces it,
        # and returns False, as is_dummy does.
        (
            "isDummy",
            "is_dummy",
            (),
            "'diffpy.srfit.structure.objcrystparset.ObjCrystScattererParSet."
            "isDummy' is deprecated and will be removed in version 4.0.0. "
            "Please use 'diffpy.cmistructure.objcrystparset."
            "ObjCrystScattererParSet.is_dummy' instead.",
        ),
        # C6: wrapRestraints is called on a molecule without restraints.
        # Expected: Warns that wrap_restraints in diffpy.cmistructure
        # replaces it, and returns None, as wrap_restraints does.
        (
            "wrapRestraints",
            "wrap_restraints",
            (),
            "'diffpy.srfit.structure.objcrystparset.ObjCrystMoleculeParSet."
            "wrapRestraints' is deprecated and will be removed in version "
            "4.0.0. Please use 'diffpy.cmistructure.objcrystparset."
            "ObjCrystMoleculeParSet.wrap_restraints' instead.",
        ),
        # C7: wrapStretchModeParameters is called on a molecule without
        # stretch modes.
        # Expected: Warns that wrap_stretch_mode_parameters in
        # diffpy.cmistructure replaces it, and returns None, as
        # wrap_stretch_mode_parameters does.
        (
            "wrapStretchModeParameters",
            "wrap_stretch_mode_parameters",
            (),
            "'diffpy.srfit.structure.objcrystparset.ObjCrystMoleculeParSet."
            "wrapStretchModeParameters' is deprecated and will be removed in "
            "version 4.0.0. Please use 'diffpy.cmistructure.objcrystparset."
            "ObjCrystMoleculeParSet.wrap_stretch_mode_parameters' instead.",
        ),
    ],
)
def test_molecule_methods(
    pyobjcryst_available, old_name, new_name, input_args, expected_msg
):
    if not pyobjcryst_available:
        pytest.skip("pyobjcryst package not available")
    molecule = make_crystal_parameter_set().molecule
    with pytest.warns(DeprecationWarning) as record:
        actual_result = getattr(molecule, old_name)(*input_args)
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg
    expected_result = getattr(molecule, new_name)(*input_args)
    assert actual_result == expected_result


def test_molecule_getLattice(pyobjcryst_available):
    # Case: getLattice is called on an ObjCrystMoleculeParSet.
    # Expected: A DeprecationWarning points to get_lattice in
    # diffpy.cmistructure, and the lattice has the same Parameters as the
    # one get_lattice returns.
    if not pyobjcryst_available:
        pytest.skip("pyobjcryst package not available")
    molecule = make_crystal_parameter_set().molecule
    expected_msg = (
        "'diffpy.srfit.structure.objcrystparset.ObjCrystMoleculeParSet."
        "getLattice' is deprecated and will be removed in version 4.0.0. "
        "Please use 'diffpy.cmistructure.objcrystparset."
        "ObjCrystMoleculeParSet.get_lattice' instead."
    )
    with pytest.warns(DeprecationWarning) as record:
        actual_lattice = molecule.getLattice()
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg
    expected_lattice = molecule.get_lattice()
    actual_names = actual_lattice.get_names()
    expected_names = expected_lattice.get_names()
    assert actual_names == expected_names
    actual_values = actual_lattice.get_values()
    expected_values = expected_lattice.get_values()
    assert actual_values == expected_values


def test_molecule_canAdapt(pyobjcryst_available):
    # Case: canAdapt is called on ObjCrystMoleculeParSet with a pyobjcryst
    # Molecule.
    # Expected: A DeprecationWarning points to can_adapt in
    # diffpy.cmistructure, and True is returned.
    if not pyobjcryst_available:
        pytest.skip("pyobjcryst package not available")
    from diffpy.srfit.structure.objcrystparset import ObjCrystMoleculeParSet

    molecule = make_crystal_parameter_set().molecule
    expected_msg = (
        "'diffpy.srfit.structure.objcrystparset.ObjCrystMoleculeParSet."
        "canAdapt' is deprecated and will be removed in version 4.0.0. "
        "Please use 'diffpy.cmistructure.objcrystparset."
        "ObjCrystMoleculeParSet.can_adapt' instead."
    )
    with pytest.warns(DeprecationWarning) as record:
        actual_result = ObjCrystMoleculeParSet.canAdapt(molecule.scatterer)
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg
    expected_result = True
    assert actual_result == expected_result


@pytest.mark.parametrize(
    "old_name, new_name, atom_count, expected_msg",
    [
        # Test the deprecated methods that add stretch-mode Parameters to
        # an ObjCrystMoleculeParSet.
        # C1: addBondLengthParameter is called with two atoms.
        # Expected: Warns that add_bond_length_parameter in
        # diffpy.cmistructure replaces it, and the Parameter has the same
        # value as one add_bond_length_parameter makes.
        (
            "addBondLengthParameter",
            "add_bond_length_parameter",
            2,
            "'diffpy.srfit.structure.objcrystparset.ObjCrystMoleculeParSet."
            "addBondLengthParameter' is deprecated and will be removed in "
            "version 4.0.0. Please use 'diffpy.cmistructure.objcrystparset."
            "ObjCrystMoleculeParSet.add_bond_length_parameter' instead.",
        ),
        # C2: addBondAngleParameter is called with three atoms.
        # Expected: Warns that add_bond_angle_parameter in
        # diffpy.cmistructure replaces it, and the Parameter has the same
        # value as one add_bond_angle_parameter makes.
        (
            "addBondAngleParameter",
            "add_bond_angle_parameter",
            3,
            "'diffpy.srfit.structure.objcrystparset.ObjCrystMoleculeParSet."
            "addBondAngleParameter' is deprecated and will be removed in "
            "version 4.0.0. Please use 'diffpy.cmistructure.objcrystparset."
            "ObjCrystMoleculeParSet.add_bond_angle_parameter' instead.",
        ),
        # C3: addDihedralAngleParameter is called with four atoms.
        # Expected: Warns that add_dihedral_angle_parameter in
        # diffpy.cmistructure replaces it, and the Parameter has the same
        # value as one add_dihedral_angle_parameter makes.
        (
            "addDihedralAngleParameter",
            "add_dihedral_angle_parameter",
            4,
            "'diffpy.srfit.structure.objcrystparset.ObjCrystMoleculeParSet."
            "addDihedralAngleParameter' is deprecated and will be removed in "
            "version 4.0.0. Please use 'diffpy.cmistructure.objcrystparset."
            "ObjCrystMoleculeParSet.add_dihedral_angle_parameter' instead.",
        ),
    ],
)
def test_molecule_add_parameter_methods(
    pyobjcryst_available, old_name, new_name, atom_count, expected_msg
):
    if not pyobjcryst_available:
        pytest.skip("pyobjcryst package not available")
    old_molecule = make_crystal_parameter_set().molecule
    new_molecule = make_crystal_parameter_set().molecule
    with pytest.warns(DeprecationWarning) as record:
        actual_parameter = getattr(old_molecule, old_name)(
            "stretch", *old_molecule.atoms[:atom_count]
        )
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg
    expected_parameter = getattr(new_molecule, new_name)(
        "stretch", *new_molecule.atoms[:atom_count]
    )
    actual_value = actual_parameter.get_value()
    expected_value = expected_parameter.get_value()
    assert actual_value == expected_value


@pytest.mark.parametrize(
    "old_name, new_name, atom_count, target_value, expected_msg",
    [
        # Test the deprecated methods that restrain atoms of an
        # ObjCrystMoleculeParSet.
        # C1: restrainBondLength is called with two atoms.
        # Expected: Warns that restrain_bond_length in diffpy.cmistructure
        # replaces it, and the restraint has the same penalty as one
        # restrain_bond_length makes.
        (
            "restrainBondLength",
            "restrain_bond_length",
            2,
            1.4,
            "'diffpy.srfit.structure.objcrystparset.ObjCrystMoleculeParSet."
            "restrainBondLength' is deprecated and will be removed in "
            "version 4.0.0. Please use 'diffpy.cmistructure.objcrystparset."
            "ObjCrystMoleculeParSet.restrain_bond_length' instead.",
        ),
        # C2: restrainBondAngle is called with three atoms.
        # Expected: Warns that restrain_bond_angle in diffpy.cmistructure
        # replaces it, and the restraint has the same penalty as one
        # restrain_bond_angle makes.
        (
            "restrainBondAngle",
            "restrain_bond_angle",
            3,
            1.5,
            "'diffpy.srfit.structure.objcrystparset.ObjCrystMoleculeParSet."
            "restrainBondAngle' is deprecated and will be removed in version "
            "4.0.0. Please use 'diffpy.cmistructure.objcrystparset."
            "ObjCrystMoleculeParSet.restrain_bond_angle' instead.",
        ),
        # C3: restrainDihedralAngle is called with four atoms.
        # Expected: Warns that restrain_dihedral_angle in
        # diffpy.cmistructure replaces it, and the restraint has the same
        # penalty as one restrain_dihedral_angle makes.
        (
            "restrainDihedralAngle",
            "restrain_dihedral_angle",
            4,
            1.5,
            "'diffpy.srfit.structure.objcrystparset.ObjCrystMoleculeParSet."
            "restrainDihedralAngle' is deprecated and will be removed in "
            "version 4.0.0. Please use 'diffpy.cmistructure.objcrystparset."
            "ObjCrystMoleculeParSet.restrain_dihedral_angle' instead.",
        ),
    ],
)
def test_molecule_restrain_methods(
    pyobjcryst_available,
    old_name,
    new_name,
    atom_count,
    target_value,
    expected_msg,
):
    if not pyobjcryst_available:
        pytest.skip("pyobjcryst package not available")
    old_molecule = make_crystal_parameter_set().molecule
    new_molecule = make_crystal_parameter_set().molecule
    with pytest.warns(DeprecationWarning) as record:
        actual_restraint = getattr(old_molecule, old_name)(
            *old_molecule.atoms[:atom_count], target_value, 0.1, 0.1
        )
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg
    expected_restraint = getattr(new_molecule, new_name)(
        *new_molecule.atoms[:atom_count], target_value, 0.1, 0.1
    )
    actual_penalty = actual_restraint.penalty()
    expected_penalty = expected_restraint.penalty()
    assert actual_penalty == expected_penalty


@pytest.mark.parametrize(
    "parameter_adder, old_name, new_name, atom_count, target_value, "
    "expected_msg",
    [
        # Test the deprecated methods that restrain stretch-mode
        # Parameters of an ObjCrystMoleculeParSet.
        # C1: restrainBondLengthParameter is called with a bond length
        # Parameter.
        # Expected: Warns that restrain_bond_length_parameter in
        # diffpy.cmistructure replaces it, and the restraint has the same
        # penalty as one restrain_bond_length_parameter makes.
        (
            "add_bond_length_parameter",
            "restrainBondLengthParameter",
            "restrain_bond_length_parameter",
            2,
            1.4,
            "'diffpy.srfit.structure.objcrystparset.ObjCrystMoleculeParSet."
            "restrainBondLengthParameter' is deprecated and will be removed "
            "in version 4.0.0. Please use 'diffpy.cmistructure."
            "objcrystparset.ObjCrystMoleculeParSet."
            "restrain_bond_length_parameter' instead.",
        ),
        # C2: restrainBondAngleParameter is called with a bond angle
        # Parameter.
        # Expected: Warns that restrain_bond_angle_parameter in
        # diffpy.cmistructure replaces it, and the restraint has the same
        # penalty as one restrain_bond_angle_parameter makes.
        (
            "add_bond_angle_parameter",
            "restrainBondAngleParameter",
            "restrain_bond_angle_parameter",
            3,
            1.5,
            "'diffpy.srfit.structure.objcrystparset.ObjCrystMoleculeParSet."
            "restrainBondAngleParameter' is deprecated and will be removed "
            "in version 4.0.0. Please use 'diffpy.cmistructure."
            "objcrystparset.ObjCrystMoleculeParSet."
            "restrain_bond_angle_parameter' instead.",
        ),
        # C3: restrainDihedralAngleParameter is called with a dihedral
        # angle Parameter.
        # Expected: Warns that restrain_dihedral_angle_parameter in
        # diffpy.cmistructure replaces it, and the restraint has the same
        # penalty as one restrain_dihedral_angle_parameter makes.
        (
            "add_dihedral_angle_parameter",
            "restrainDihedralAngleParameter",
            "restrain_dihedral_angle_parameter",
            4,
            1.5,
            "'diffpy.srfit.structure.objcrystparset.ObjCrystMoleculeParSet."
            "restrainDihedralAngleParameter' is deprecated and will be "
            "removed in version 4.0.0. Please use 'diffpy.cmistructure."
            "objcrystparset.ObjCrystMoleculeParSet."
            "restrain_dihedral_angle_parameter' instead.",
        ),
    ],
)
def test_molecule_restrain_parameter_methods(
    pyobjcryst_available,
    parameter_adder,
    old_name,
    new_name,
    atom_count,
    target_value,
    expected_msg,
):
    if not pyobjcryst_available:
        pytest.skip("pyobjcryst package not available")
    old_molecule = make_crystal_parameter_set().molecule
    new_molecule = make_crystal_parameter_set().molecule
    old_parameter = getattr(old_molecule, parameter_adder)(
        "stretch", *old_molecule.atoms[:atom_count]
    )
    new_parameter = getattr(new_molecule, parameter_adder)(
        "stretch", *new_molecule.atoms[:atom_count]
    )
    with pytest.warns(DeprecationWarning) as record:
        actual_restraint = getattr(old_molecule, old_name)(
            old_parameter, target_value, 0.1, 0.1
        )
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg
    expected_restraint = getattr(new_molecule, new_name)(
        new_parameter, target_value, 0.1, 0.1
    )
    actual_penalty = actual_restraint.penalty()
    expected_penalty = expected_restraint.penalty()
    assert actual_penalty == expected_penalty


def test_addAtoms(pyobjcryst_available):
    # Case: addAtoms adds a third atom to a bond length Parameter.
    # Expected: A DeprecationWarning points to add_atoms in
    # diffpy.cmistructure, and the Parameter mutates the second and third
    # atoms.
    if not pyobjcryst_available:
        pytest.skip("pyobjcryst package not available")
    molecule = make_crystal_parameter_set().molecule
    parameter = molecule.add_bond_length_parameter(
        "bond", molecule.C0, molecule.C1
    )
    expected_msg = (
        "'diffpy.srfit.structure.objcrystparset.StretchModeParameter."
        "addAtoms' is deprecated and will be removed in version 4.0.0. "
        "Please use 'diffpy.cmistructure.objcrystparset."
        "StretchModeParameter.add_atoms' instead."
    )
    with pytest.warns(DeprecationWarning) as record:
        parameter.addAtoms([molecule.C2])
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg
    actual_mutated_atoms = parameter.mutated_atoms
    expected_mutated_atoms = {molecule.C1, molecule.C2}
    assert actual_mutated_atoms == expected_mutated_atoms


def test_angunits_objcryst_crystal(pyobjcryst_available):
    # Case: .angunits is read and then set on an ObjCrystCrystalParSet.
    # Expected: Both warn that .angle_units in diffpy.cmistructure
    # replaces it. The read returns "rad", and the set updates
    # .angle_units.
    if not pyobjcryst_available:
        pytest.skip("pyobjcryst package not available")
    crystal = make_crystal_parameter_set()
    expected_msg = (
        "'diffpy.srfit.structure.objcrystparset.ObjCrystCrystalParSet."
        "angunits' is deprecated and will be removed in version 4.0.0. "
        "Please use 'diffpy.cmistructure.objcrystparset."
        "ObjCrystCrystalParSet.angle_units' instead."
    )
    with pytest.warns(DeprecationWarning) as record:
        actual_angle_units = crystal.angunits
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg
    expected_angle_units = "rad"
    assert actual_angle_units == expected_angle_units
    with pytest.warns(DeprecationWarning) as record:
        crystal.angunits = "deg"
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg
    actual_angle_units = crystal.angle_units
    expected_angle_units = "deg"
    assert actual_angle_units == expected_angle_units


def test_scat(pyobjcryst_available):
    # Case: .scat is read and then set on an ObjCrystMoleculeParSet.
    # Expected: Both warn that .scatterer in diffpy.cmistructure replaces
    # it. The read returns the adapted pyobjcryst Molecule, and the set
    # updates .scatterer.
    if not pyobjcryst_available:
        pytest.skip("pyobjcryst package not available")
    molecule = make_crystal_parameter_set().molecule
    expected_msg = (
        "'diffpy.srfit.structure.objcrystparset.ObjCrystScattererParSet."
        "scat' is deprecated and will be removed in version 4.0.0. Please "
        "use 'diffpy.cmistructure.objcrystparset.ObjCrystScattererParSet."
        "scatterer' instead."
    )
    with pytest.warns(DeprecationWarning) as record:
        actual_scatterer = molecule.scat
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg
    expected_scatterer = molecule.scatterer
    assert actual_scatterer == expected_scatterer
    other_scatterer = make_crystal_parameter_set().molecule.scatterer
    with pytest.warns(DeprecationWarning) as record:
        molecule.scat = other_scatterer
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg
    actual_scatterer = molecule.scatterer
    expected_scatterer = other_scatterer
    assert actual_scatterer == expected_scatterer


def test_res(pyobjcryst_available):
    # Case: .res is read and then set on a bond length restraint.
    # Expected: Both warn that .restraint in diffpy.cmistructure replaces
    # it. The read returns the wrapped pyobjcryst restraint, and the set
    # updates .restraint.
    if not pyobjcryst_available:
        pytest.skip("pyobjcryst package not available")
    molecule = make_crystal_parameter_set().molecule
    restraint = molecule.restrain_bond_length(
        molecule.C0, molecule.C1, 1.4, 0.1, 0.1
    )
    other_restraint = molecule.restrain_bond_length(
        molecule.C1, molecule.C2, 1.4, 0.1, 0.1
    )
    expected_msg = (
        "'diffpy.srfit.structure.objcrystparset.ObjCrystMoleculeRestraint."
        "res' is deprecated and will be removed in version 4.0.0. Please use "
        "'diffpy.cmistructure.objcrystparset.ObjCrystMoleculeRestraint."
        "restraint' instead."
    )
    with pytest.warns(DeprecationWarning) as record:
        actual_restraint = restraint.res
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg
    expected_restraint = restraint.restraint
    assert actual_restraint == expected_restraint
    with pytest.warns(DeprecationWarning) as record:
        restraint.res = other_restraint.restraint
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg
    actual_restraint = restraint.restraint
    expected_restraint = other_restraint.restraint
    assert actual_restraint == expected_restraint


def test_matoms(pyobjcryst_available):
    # Case: .matoms is read and then set on a bond length Parameter.
    # Expected: Both warn that .mutated_atoms in diffpy.cmistructure
    # replaces it. The read returns the mutated atoms, and the set updates
    # .mutated_atoms.
    if not pyobjcryst_available:
        pytest.skip("pyobjcryst package not available")
    molecule = make_crystal_parameter_set().molecule
    parameter = molecule.add_bond_length_parameter(
        "bond", molecule.C0, molecule.C1
    )
    expected_msg = (
        "'diffpy.srfit.structure.objcrystparset.StretchModeParameter."
        "matoms' is deprecated and will be removed in version 4.0.0. Please "
        "use 'diffpy.cmistructure.objcrystparset.StretchModeParameter."
        "mutated_atoms' instead."
    )
    with pytest.warns(DeprecationWarning) as record:
        actual_mutated_atoms = parameter.matoms
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg
    expected_mutated_atoms = {molecule.C1}
    assert actual_mutated_atoms == expected_mutated_atoms
    with pytest.warns(DeprecationWarning) as record:
        parameter.matoms = {molecule.C1, molecule.C2}
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg
    actual_mutated_atoms = parameter.mutated_atoms
    expected_mutated_atoms = {molecule.C1, molecule.C2}
    assert actual_mutated_atoms == expected_mutated_atoms
