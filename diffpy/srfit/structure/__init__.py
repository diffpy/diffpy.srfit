#!/usr/bin/env python
########################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2009 Trustees of the Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
########################################################################
"""Classes that wrap structure representations as Parameters and ParameterSets.

"""

# package version
from diffpy.srfit.version import __version__

def struToParameterSet(name, stru):
    """Creates a ParameterSet from an structure.

    This returns a ParameterSet adapted for the structure depending on its
    type.

    stru    --  a structure object known by this module
    name    --  A name to give the structure.

    """
    from diffpy.srfit.structure.diffpystructure import StructureParSet
    if StructureParSet.canAdapt(stru):
        return StructureParSet(name, stru)

    from diffpy.srfit.structure.objcryststructure import ObjCrystParSet
    if ObjCrystParSet.canAdapt(stru):
        return ObjCrystParSet(name, stru)

    from diffpy.srfit.structure.objcryststructure import MoleculeParSet
    if MoleculeParSet.canAdapt(stru):
        return MoleculeParSet(name, stru)

    from diffpy.srfit.structure.cctbxstructure import CCTBXStructureParSet
    if CCTBXStructureParSet.canAdapt(stru):
        return CCTBXStructureParSet(name, stru)

    raise TypeError("Unadaptable structure format")

from sgconstraints import constrainAsSpaceGroup

# End of file
