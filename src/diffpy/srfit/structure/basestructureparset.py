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
"""Base class for adapting structures to a ParameterSet interface.

The BaseStructureParSet is a ParameterSet with functionality required by
all structure adapters.
"""

__all__ = ["BaseStructureParSet"]

from diffpy.srfit.fitbase.parameterset import ParameterSet
from diffpy.utils._deprecator import build_deprecation_message, deprecated

removal_version = "4.0.0"

stru_dep_msg = build_deprecation_message(
    "diffpy.srfit.structure.basestructureparset.BaseStructureParSet",
    "stru",
    "structure",
    removal_version,
    new_base="diffpy.cmistructure.basestructureparset.BaseStructureParSet",
)


parset_base = "diffpy.srfit.structure.basestructureparset.BaseStructureParSet"
parset_new_base = "diffpy.cmistructure.basestructureparset.BaseStructureParSet"

parset_canAdapt_dep_msg = build_deprecation_message(
    parset_base,
    "canAdapt",
    "can_adapt",
    removal_version,
    new_base=parset_new_base,
)

parset_getLattice_dep_msg = build_deprecation_message(
    parset_base,
    "getLattice",
    "get_lattice",
    removal_version,
    new_base=parset_new_base,
)

parset_getScatterers_dep_msg = build_deprecation_message(
    parset_base,
    "getScatterers",
    "get_scatterers",
    removal_version,
    new_base=parset_new_base,
)


class BaseStructureParSet(ParameterSet):
    """Base class for structure adapters.

    BaseStructureParSet derives from ParameterSet and provides methods that
    help interface the ParameterSet with the space group constraint methods in
    the sgconstraints module and to ProfileGenerators.

    Attributes
    ----------
    structure
        The adapted object
    """

    @property
    @deprecated(stru_dep_msg)
    def stru(self):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use
        diffpy.cmistructure.basestructureparset.BaseStructureParSet.structure
        instead.
        """
        return self.structure

    @stru.setter
    @deprecated(stru_dep_msg)
    def stru(self, value):
        self.structure = value

    @classmethod
    def can_adapt(self, stru):
        """Return whether the structure can be adapted by this class."""
        return False

    @classmethod
    @deprecated(parset_canAdapt_dep_msg)
    def canAdapt(self, stru):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use
        diffpy.cmistructure.basestructureparset.BaseStructureParSet.can_adapt
        instead.
        """
        return self.can_adapt(stru)

    def get_lattice(self):
        """Get a ParameterSet containing the lattice Parameters.

        The returned ParameterSet may contain other Parameters than the
        lattice Parameters. It is assumed that the lattice parameters
        are named "a", "b", "c", "alpha", "beta", "gamma".

        Lattice must also have the "angle_units" attribute, which is
        either "deg" or "rad", to signify degrees or radians.
        """
        raise NotImplementedError("The must be overloaded")

    @deprecated(parset_getLattice_dep_msg)
    def getLattice(self):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use
        diffpy.cmistructure.basestructureparset.BaseStructureParSet.get_lattice
        instead.
        """
        return self.get_lattice()

    def get_scatterers(self):
        """Get a list of ParameterSets that represents the scatterers.

        The site positions must be accessible from the list entries via
        the names "x", "y", and "z". The ADPs must be accessible as
        well, but the name and nature of the ADPs (U-factors, B-factors,
        isotropic, anisotropic) depends on the adapted structure.
        """
        raise NotImplementedError("The must be overloaded")

    @deprecated(parset_getScatterers_dep_msg)
    def getScatterers(self):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use
        diffpy.cmistructure.basestructureparset.BaseStructureParSet.get_scatterers
        instead.
        """
        return self.get_scatterers()
