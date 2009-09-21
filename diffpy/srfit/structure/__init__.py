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
"""Classes that wrap structure containers with Parameters and ParameterSets.

Modules:
    cctbxstructure      --  Home to adapters for cctbx structure objects. See
                            the module documentation for more details.
    diffpystructure     --  Home to adapters for the diffpy.Structure objects.
                            See the module documentation for more details.
    objcryststructure   --  Home to adapters for pyobjcryst structure objects.
                            See the module documentation for more details.

"""

# package version
from diffpy.srfit.version import __version__

def struToParameterSet(stru, name):
    """Creates a ParameterSet from an structure.

    This returns a ParameterSet adapted for the structure depending on its
    type.

    stru    --  a structure object known by this module
    name    --  A name to give the structure.

    """
    from diffpy.srfit.structure.diffpystructure import StructureParSet
    if StructureParSet.canAdapt(stru):
        return StructureParSet(stru, name)

    from diffpy.srfit.structure.objcryststructure import ObjCrystParSet
    if ObjCrystParSet.canAdapt(stru):
        return ObjCrystParSet(stru, name)

    from diffpy.srfit.structure.cctbxstructure import CCTBXStructureParSet
    if CCTBXStructureParSet.canAdapt(stru):
        return CCTBXStructureParSet(stru, name)

    raise TypeError("Unadapatable structure format")


# End of file
