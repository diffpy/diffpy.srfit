#!/usr/bin/env python
##############################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2009 The Trustees of Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################
"""Modules and classes that adapt structure representations to the
ParameterSet interface and automatic structure constraint generation
from space group information."""

import warnings

from diffpy.srfit.structure.sgconstraints import (
    constrain_as_space_group,
    constrainAsSpaceGroup,
)
from diffpy.utils._deprecator import build_deprecation_message, deprecated

removal_version = "4.0.0"

warnings.warn(
    build_deprecation_message(
        "diffpy.srfit",
        "structure",
        "cmistructure",
        removal_version,
        new_base="diffpy",
    ),
    DeprecationWarning,
    stacklevel=2,
)

struToParameterSet_dep_msg = build_deprecation_message(
    "diffpy.srfit.structure",
    "struToParameterSet",
    "structure_to_parameter_set",
    removal_version,
    new_base="diffpy.cmistructure",
)


def structure_to_parameter_set(name, structure):
    """Create a ParameterSet adapted to a structure object.

    The adapter is chosen from the type of `structure`. Supported types
    are diffpy.structure.Structure, pyobjcryst.crystal.Crystal,
    pyobjcryst.molecule.Molecule and
    cctbx.crystal.special_position_settings.

    Parameters
    ----------
    name : str
        The name to give the structure.
    structure : object
        The structure object to adapt.

    Returns
    -------
    BaseStructureParSet
        The ParameterSet adapting `structure`.

    Raises
    ------
    TypeError
        If `structure` is not one of the supported structure types.
    """
    from diffpy.srfit.structure.diffpyparset import DiffpyStructureParSet

    if DiffpyStructureParSet.canAdapt(structure):
        return DiffpyStructureParSet(name, structure)

    from diffpy.srfit.structure.objcrystparset import ObjCrystCrystalParSet

    if ObjCrystCrystalParSet.canAdapt(structure):
        return ObjCrystCrystalParSet(name, structure)

    from diffpy.srfit.structure.objcrystparset import ObjCrystMoleculeParSet

    if ObjCrystMoleculeParSet.canAdapt(structure):
        return ObjCrystMoleculeParSet(name, structure)

    from diffpy.srfit.structure.cctbxparset import CCTBXCrystalParSet

    if CCTBXCrystalParSet.canAdapt(structure):
        return CCTBXCrystalParSet(name, structure)

    raise TypeError("Unadaptable structure format")


@deprecated(struToParameterSet_dep_msg)
def struToParameterSet(name, stru):
    """This function has been deprecated and will be removed in version
    4.0.0.

    Please use diffpy.cmistructure.structure_to_parameter_set instead.
    """
    return structure_to_parameter_set(name, stru)


# silence pyflakes checker
assert constrain_as_space_group
assert constrainAsSpaceGroup


# End of file
