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
"""Wrappers for interfacing cctbx crystal with SrFit.

This wraps a cctbx.crystal as a ParameterSet with a similar hierarchy, which
can then be used within a FitRecipe. Note that all manipulations to the
cctbx.crystal should be done before wrapping. Changes made to the cctbx.crystal
object after wrapping may not be reflected within the wrapper, which can have
unpredictable results during a structure refinement.

The following classes are adapted:

- `CCTBXCrystalParSet`: wrapper for `cctbx.crystal`.
- `CCTBXUnitCellParSet`: wrapper for the unit cell of `cctbx.crystal`.
- `CCTBXScattererParSet`: wrapper for `cctbx.xray.scatterer`.
"""

from diffpy.srfit.fitbase.parameter import ParameterAdapter
from diffpy.srfit.fitbase.parameterset import ParameterSet
from diffpy.srfit.structure.basestructureparset import BaseStructureParSet
from diffpy.utils._deprecator import build_deprecation_message, deprecated

__all__ = ["CCTBXScattererParSet", "CCTBXUnitCellParSet", "CCTBXCrystalParSet"]


removal_version = "4.0.0"

scatterer_base = "diffpy.srfit.structure.cctbxparset.CCTBXScattererParSet"
scatterer_new_base = "diffpy.cmistructure.cctbxparset.CCTBXScattererParSet"

scatterer_strups_dep_msg = build_deprecation_message(
    scatterer_base,
    "strups",
    "structure_parameter_set",
    removal_version,
    new_base=scatterer_new_base,
)

scatterer_idx_dep_msg = build_deprecation_message(
    scatterer_base,
    "idx",
    "index",
    removal_version,
    new_base=scatterer_new_base,
)

unit_cell_base = "diffpy.srfit.structure.cctbxparset.CCTBXUnitCellParSet"
unit_cell_new_base = "diffpy.cmistructure.cctbxparset.CCTBXUnitCellParSet"

unit_cell_strups_dep_msg = build_deprecation_message(
    unit_cell_base,
    "strups",
    "structure_parameter_set",
    removal_version,
    new_base=unit_cell_new_base,
)

crystal_base = "diffpy.srfit.structure.cctbxparset.CCTBXCrystalParSet"
crystal_new_base = "diffpy.cmistructure.cctbxparset.CCTBXCrystalParSet"

crystal_getSpaceGroup_dep_msg = build_deprecation_message(
    crystal_base,
    "getSpaceGroup",
    "get_space_group",
    removal_version,
    new_base=crystal_new_base,
)


class CCTBXScattererParSet(ParameterSet):
    """A wrapper for cctbx.xray.scatterer.

    This class derives from ParameterSet.

    Attributes
    ----------
    name
        Name of the scatterer. The name is always of the form
        "%s%i" % (element, number), where the number is the running
        index of that element type (starting at 0).
        x (y, z)    --  Atom position in crystal coordinates (ParameterAdapter)
    occupancy
        Occupancy of the atom on its crystal location
        (ParameterAdapter)
    Uiso
        Isotropic scattering factor (ParameterAdapter).
    """

    @property
    @deprecated(scatterer_idx_dep_msg)
    def idx(self):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use diffpy.cmistructure.cctbxparset.CCTBXScattererParSet.index
        instead.
        """
        return self.index

    @idx.setter
    @deprecated(scatterer_idx_dep_msg)
    def idx(self, value):
        self.index = value

    @property
    @deprecated(scatterer_strups_dep_msg)
    def strups(self):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use
        diffpy.cmistructure.cctbxparset.CCTBXScattererParSet.structure_parameter_set
        instead.
        """
        return self.structure_parameter_set

    @strups.setter
    @deprecated(scatterer_strups_dep_msg)
    def strups(self, value):
        self.structure_parameter_set = value

    def __init__(self, name, strups, idx):
        """Initialize.

        Parameters
        ----------
        name
            The name of this scatterer.
        strups
            The CCTBXCrystalParSet that contains the cctbx structure
        idx
            The index of the scatterer in the structure.
        """
        ParameterSet.__init__(self, name)
        self.structure_parameter_set = strups
        self.index = idx

        # x, y, z, occupancy
        self.add_parameter(
            ParameterAdapter("x", None, self._xyzgetter(0), self._xyzsetter(0))
        )
        self.add_parameter(
            ParameterAdapter("y", None, self._xyzgetter(1), self._xyzsetter(1))
        )
        self.add_parameter(
            ParameterAdapter("z", None, self._xyzgetter(2), self._xyzsetter(2))
        )
        self.add_parameter(
            ParameterAdapter("occupancy", None, self._getocc, self._setocc)
        )
        self.add_parameter(
            ParameterAdapter("Uiso", None, self._getuiso, self._setuiso)
        )
        return

    # Getters and setters

    def _xyzgetter(self, i):

        def f(dummy):
            return self.structure_parameter_set.structure.scatterers()[
                self.index
            ].site[i]

        return f

    def _xyzsetter(self, i):

        def f(dummy, value):
            xyz = list(
                self.structure_parameter_set.structure.scatterers()[
                    self.index
                ].site
            )
            xyz[i] = value
            self.structure_parameter_set.structure.scatterers()[
                self.index
            ].site = tuple(xyz)
            return

        return f

    def _getocc(self, dummy):
        return self.structure_parameter_set.structure.scatterers()[
            self.index
        ].occupancy

    def _setocc(self, dummy, value):
        self.structure_parameter_set.structure.scatterers()[
            self.index
        ].occupancy = value
        return

    def _getuiso(self, dummy):
        return self.structure_parameter_set.structure.scatterers()[
            self.index
        ].u_iso

    def _setuiso(self, dummy, value):
        self.structure_parameter_set.structure.scatterers()[
            self.index
        ].u_iso = value
        return

    def _getelem(self):
        return self.stru.element_symbol()

    element = property(_getelem)


# End class CCTBXScattererParSet


class CCTBXUnitCellParSet(ParameterSet):
    """A wrapper for cctbx unit_cell object.

    Attributes
    ----------
    name
        Always "unitcell".
    a
        Unit cell parameters (ParameterAdapter).
    b
        Unit cell parameters (ParameterAdapter).
    c
        Unit cell parameters (ParameterAdapter).
    alpha
        Unit cell parameters (ParameterAdapter).
    beta
        Unit cell parameters (ParameterAdapter).
    gamma
        Unit cell parameters (ParameterAdapter).
    """

    @property
    @deprecated(unit_cell_strups_dep_msg)
    def strups(self):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use
        diffpy.cmistructure.cctbxparset.CCTBXUnitCellParSet.structure_parameter_set
        instead.
        """
        return self.structure_parameter_set

    @strups.setter
    @deprecated(unit_cell_strups_dep_msg)
    def strups(self, value):
        self.structure_parameter_set = value

    def __init__(self, strups):
        """Initialize.

        Parameters
        ----------
        strups
            The CCTBXCrystalParSet that contains the cctbx structure
            and the unit cell we're wrapper.
        """
        ParameterSet.__init__(self, "unitcell")
        self.structure_parameter_set = strups
        self._latpars = list(
            self.structure_parameter_set.structure.unit_cell().parameters()
        )

        self.add_parameter(
            ParameterAdapter("a", None, self._latgetter(0), self._latsetter(0))
        )
        self.add_parameter(
            ParameterAdapter("b", None, self._latgetter(1), self._latsetter(1))
        )
        self.add_parameter(
            ParameterAdapter("c", None, self._latgetter(2), self._latsetter(2))
        )
        self.add_parameter(
            ParameterAdapter(
                "alpha", None, self._latgetter(3), self._latsetter(3)
            )
        )
        self.add_parameter(
            ParameterAdapter(
                "beta", None, self._latgetter(4), self._latsetter(4)
            )
        )
        self.add_parameter(
            ParameterAdapter(
                "gamma", None, self._latgetter(5), self._latsetter(5)
            )
        )

        return

    def _latgetter(self, i):

        def f(dummy):
            return self._latpars[i]

        return f

    def _latsetter(self, i):

        def f(dummy, value):
            self._latpars[i] = value
            self.structure_parameter_set._update = True
            return

        return f


# End class CCTBXUnitCellParSet

# FIXME - Special positions should be constant.


class CCTBXCrystalParSet(BaseStructureParSet):
    """A wrapper for CCTBX structure.

    Attributes
    ----------
    stru
        The adapted cctbx structure object.
    scatterers
        The list of ScattererParSets.
    unitcell
        The CCTBXUnitCellParSet for the structure.
    """

    def __init__(self, name, stru):
        """Initialize.

        Parameters
        ----------
        name
            A name for this
        stru
            A CCTBX structure instance.
        """
        ParameterSet.__init__(self, name)
        self.structure = stru
        self.add_parameter_set(CCTBXUnitCellParSet(self))
        self.scatterers = []

        self._update = False

        cdict = {}
        for s in stru.scatterers():
            el = s.element_symbol()
            i = cdict.get(el, 0)
            sname = "%s%i" % (el, i)
            cdict[el] = i + 1
            scatterer = CCTBXScattererParSet(sname, self, i)
            self.add_parameter_set(scatterer)
            self.scatterers.append(scatterer)

        # Constrain the lattice
        from diffpy.srfit.structure.sgconstraints import _constrain_space_group

        symbol = self.get_space_group()
        _constrain_space_group(self, symbol)

        return

    def update(self):
        """Update the unit_cell to a change in lattice parameters.

        This remakes the unit cell according to a change in the lattice
        parameters. Call this function before using the
        CCTBXCrystalParSet. The unit_cell will only be remade if
        necessary.
        """
        if not self._update:
            return

        self._update = False
        stru = self.structure
        sgn = stru.space_group().match_tabulated_settings().number()

        # Create the symmetry object
        from cctbx.crystal import symmetry

        symm = symmetry(
            unit_cell=self.unitcell._latpars, space_group_symbol=sgn
        )

        # Now the new structure
        newstru = stru.__class__(
            crystal_symmetry=symm, scatterers=stru.scatterers()
        )

        self.unitcell._latpars = list(newstru.unit_cell().parameters())

        self.structure = newstru
        return

    @classmethod
    def can_adapt(self, stru):
        """Return whether the structure can be adapted by this class."""
        try:
            from cctbx.crystal import special_position_settings
        except ImportError:
            return False
        return isinstance(stru, special_position_settings)

    def get_lattice(self):
        """Get the ParameterSet containing the lattice Parameters."""
        return self.unitcell

    def get_scatterers(self):
        """Get a list of ParameterSets that represents the scatterers.

        The site positions must be accessible from the list entries via
        the names "x", "y", and "z". The ADPs must be accessible as
        well, but the name and nature of the ADPs (U-factors, B-factors,
        isotropic, anisotropic) depends on the adapted structure.
        """
        return self.scatterers

    def get_space_group(self):
        """Get the HM space group symbol for the structure."""
        sg = self.structure.space_group()
        t = sg.type()
        return t.lookup_symbol()

    @deprecated(crystal_getSpaceGroup_dep_msg)
    def getSpaceGroup(self):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use
        diffpy.cmistructure.cctbxparset.CCTBXCrystalParSet.get_space_group
        instead.
        """
        return self.get_space_group()


# End class CCTBXCrystalParSet
