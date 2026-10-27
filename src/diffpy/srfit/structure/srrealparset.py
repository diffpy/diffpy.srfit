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
"""Structure wrapper class for structures compatible with SrReal."""

__all__ = ["SrRealParSet"]

from diffpy.srfit.structure.basestructureparset import BaseStructureParSet
from diffpy.srfit.structure.bvsrestraint import BVSRestraint
from diffpy.utils._deprecator import build_deprecation_message, deprecated

removal_version = "4.0.0"

restrainBVS_dep_msg = build_deprecation_message(
    "diffpy.srfit.structure.srrealparset.SrRealParSet",
    "restrainBVS",
    "restrain_bvs",
    removal_version,
    new_base="diffpy.cmistructure.srrealparset.SrRealParSet",
)


srreal_base = "diffpy.srfit.structure.srrealparset.SrRealParSet"
srreal_new_base = "diffpy.cmistructure.srrealparset.SrRealParSet"

srreal_useSymmetry_dep_msg = build_deprecation_message(
    srreal_base,
    "useSymmetry",
    "use_symmetry",
    removal_version,
    new_base=srreal_new_base,
)

srreal_usingSymmetry_dep_msg = build_deprecation_message(
    srreal_base,
    "usingSymmetry",
    "using_symmetry",
    removal_version,
    new_base=srreal_new_base,
)


class SrRealParSet(BaseStructureParSet):
    """Base class for SrReal-compatible structure adapters.

    This derives from BaseStructureParSet and provides some extended
    functionality provided by SrReal.

    Attributes
    ----------
    structure
        The adapted object
    _usesymmetry
        A flag indicating if SrReal calculators that operate on
        this object should use symmetry. By default this is
        True.
    """

    def __init__(self, *args, **kw):
        BaseStructureParSet.__init__(self, *args, **kw)
        self._usesymmetry = True
        self.structure = None
        return

    def restrain_bvs(self, sig=1, scaled=False):
        """Restrain the bond-valence sum to zero.

        This adds a penalty to the cost function equal to
        bvmsdiff / sig**2
        where bvmsdiff is the mean-squared difference between the calculated
        and expected bond valence sums for the structure. If scaled is True,
        this is also scaled by the current point-averaged chi^2 value so the
        restraint is roughly equally weighted in the fit.

        Parameters
        ----------
        sig : float, optional
            The uncertainty on the BVS (default 1).
        scaled : bool, optional
            The flag indicating if the restraint is scaled
            (multiplied) by the unrestrained point-average chi^2
            (chi^2/numpoints) (default False).

        Returns
        -------
        BVSRestraint
            The restraint, for use with the 'unrestrain' method.
        """
        # Create the Restraint object
        res = BVSRestraint(self, sig, scaled)
        # Add it to the _restraints set
        self._restraints.add(res)
        # Our configuration changed. Notify observers.
        self._update_configuration()
        # Return the Restraint object
        return res

    @deprecated(restrainBVS_dep_msg)
    def restrainBVS(self, sig=1, scaled=False):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use diffpy.cmistructure.srrealparset.SrRealParSet.restrain_bvs
        instead.
        """
        return self.restrain_bvs(sig, scaled)

    def use_symmetry(self, use=True):
        """Set this structure to use symmetry.

        This determines how the structure is treated by SrReal
        calculators.
        """
        self._usesymmetry = bool(use)
        return

    @deprecated(srreal_useSymmetry_dep_msg)
    def useSymmetry(self, use=True):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use diffpy.cmistructure.srrealparset.SrRealParSet.use_symmetry
        instead.
        """
        return self.use_symmetry(use)

    def using_symmetry(self):
        """Check if symmetry is being used."""
        return self._usesymmetry

    @deprecated(srreal_usingSymmetry_dep_msg)
    def usingSymmetry(self):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use diffpy.cmistructure.srrealparset.SrRealParSet.using_symmetry
        instead.
        """
        return self.using_symmetry()

    def _get_srreal_structure(self):
        """Get the structure object for use with SrReal calculators.

        If this is periodic, then return the structure, otherwise, pass
        it inside of a nosymmetry wrapper.
        """
        from diffpy.srreal.structureadapter import nosymmetry

        if self._usesymmetry:
            return self.structure
        return nosymmetry(self.structure)
