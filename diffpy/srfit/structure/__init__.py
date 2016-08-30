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

"""Modules and classes that adapt structure representations to the ParameterSet
interface and automatic structure constraint generation from space group
information.
"""

def struToParameterSet(name, stru):
    """Creates a ParameterSet from an structure.

    This returns a ParameterSet adapted for the structure depending on its
    type.

    stru    --  a structure object known by this module
    name    --  A name to give the structure.

    Raises TypeError if stru cannot be adapted

    """
    from diffpy.srfit.structure.diffpyparset import DiffpyStructureParSet
    if DiffpyStructureParSet.canAdapt(stru):
        return DiffpyStructureParSet(name, stru)

    from diffpy.srfit.structure.objcrystparset import ObjCrystCrystalParSet
    if ObjCrystCrystalParSet.canAdapt(stru):
        return ObjCrystCrystalParSet(name, stru)

    from diffpy.srfit.structure.objcrystparset import ObjCrystMoleculeParSet
    if ObjCrystMoleculeParSet.canAdapt(stru):
        return ObjCrystMoleculeParSet(name, stru)

    from diffpy.srfit.structure.cctbxparset import CCTBXCrystalParSet
    if CCTBXCrystalParSet.canAdapt(stru):
        return CCTBXCrystalParSet(name, stru)

    raise TypeError("Unadaptable structure format")


from diffpy.srfit.structure.sgconstraints import constrainAsSpaceGroup

# silence pyflakes checker
assert constrainAsSpaceGroup


# End of file
