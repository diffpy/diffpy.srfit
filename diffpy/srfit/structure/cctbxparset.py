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

Classes:

CCTBXCrystalParSet  --  Wrapper for cctbx.crystal
CCTBXUnitCellParSet --  Wrapper for the unit cell of cctbx.crystal
CCTBXScattererParSet --  Wrapper for cctbx.xray.scatterer
"""

from diffpy.srfit.fitbase.parameter import ParameterAdapter
from diffpy.srfit.fitbase.parameterset import ParameterSet
from diffpy.srfit.structure.basestructureparset import BaseStructureParSet

__all__ = ["CCTBXScattererParSet", "CCTBXUnitCellParSet",
           "CCTBXCrystalParSet"]

class CCTBXScattererParSet(ParameterSet):
    """A wrapper for cctbx.xray.scatterer

    This class derives from ParameterSet.

    Attributes:
    name        --  Name of the scatterer. The name is always of the form
                    "%s%i" % (element, number), where the number is the running
                    index of that element type (starting at 0).
    x (y, z)    --  Atom position in crystal coordinates (ParameterAdapter)
    occupancy   --  Occupancy of the atom on its crystal location
                    (ParameterAdapter)
    Uiso        --  Isotropic scattering factor (ParameterAdapter).

    """

    def __init__(self, name, strups, idx):
        """Initialize

        name    --  The name of this scatterer.
        strups  --  The CCTBXCrystalParSet that contains the cctbx structure
        idx     --  The index of the scatterer in the structure.

        """
        ParameterSet.__init__(self, name)
        self.strups = strups
        self.idx = idx

        # x, y, z, occupancy
        self.addParameter(ParameterAdapter("x", None, self._xyzgetter(0),
            self._xyzsetter(0)))
        self.addParameter(ParameterAdapter("y", None, self._xyzgetter(1),
            self._xyzsetter(1)))
        self.addParameter(ParameterAdapter("z", None, self._xyzgetter(2),
            self._xyzsetter(2)))
        self.addParameter(ParameterAdapter("occupancy", None, self._getocc,
            self._setocc))
        self.addParameter(ParameterAdapter("Uiso", None, self._getuiso,
            self._setuiso))
        return

    # Getters and setters

    def _xyzgetter(self, i):

        def f(dummy):
            return self.strups.stru.scatterers()[self.idx].site[i]

        return f

    def _xyzsetter(self, i):

        def f(dummy, value):
            xyz = list(self.strups.stru.scatterers()[self.idx].site)
            xyz[i] = value
            self.strups.stru.scatterers()[self.idx].site = tuple(xyz)
            return

        return f

    def _getocc(self, dummy):
        return self.strups.stru.scatterers()[self.idx].occupancy

    def _setocc(self, dummy, value):
        self.strups.stru.scatterers()[self.idx].occupancy = value
        return

    def _getuiso(self, dummy):
        return self.strups.stru.scatterers()[self.idx].u_iso

    def _setuiso(self, dummy, value):
        self.strups.stru.scatterers()[self.idx].u_iso = value
        return

    def _getElem(self):
        return self.stru.element_symbol()

    element = property(_getElem)

# End class CCTBXScattererParSet

class CCTBXUnitCellParSet(ParameterSet):
    """A wrapper for cctbx unit_cell object.

    Attributes:
    name    --  Always "unitcell".
    a, b, c, alpha, beta, gamma --  Unit cell parameters (ParameterAdapter).

    """

    def __init__(self, strups):
        """Initialize

        strups  --  The CCTBXCrystalParSet that contains the cctbx structure
                    and the unit cell we're wrapper.

        """
        ParameterSet.__init__(self, "unitcell")
        self.strups = strups
        self._latpars = list(self.strups.stru.unit_cell().parameters())

        self.addParameter(ParameterAdapter("a", None, self._latgetter(0),
            self._latsetter(0)))
        self.addParameter(ParameterAdapter("b", None, self._latgetter(1),
            self._latsetter(1)))
        self.addParameter(ParameterAdapter("c", None, self._latgetter(2),
            self._latsetter(2)))
        self.addParameter(ParameterAdapter("alpha", None, self._latgetter(3),
            self._latsetter(3)))
        self.addParameter(ParameterAdapter("beta", None, self._latgetter(4),
            self._latsetter(4)))
        self.addParameter(ParameterAdapter("gamma", None, self._latgetter(5),
            self._latsetter(5)))

        return

    def _latgetter(self, i):

        def f(dummy):
            return self._latpars[i]

        return f

    def _latsetter(self, i):

        def f(dummy, value):
            self._latpars[i] = value
            self.strups._update = True
            return

        return f


# End class CCTBXUnitCellParSet

# FIXME - Special positions should be constant.

class CCTBXCrystalParSet(BaseStructureParSet):
    """A wrapper for CCTBX structure.

    Attributes:
    stru        --  The adapted cctbx structure object.
    scatterers  --  The list of ScattererParSets.
    unitcell    --  The CCTBXUnitCellParSet for the structure.

    """

    def __init__(self, name, stru):
        """Initialize

        name    --  A name for this
        stru    --  A CCTBX structure instance.

        """
        ParameterSet.__init__(self, name)
        self.stru = stru
        self.addParameterSet(CCTBXUnitCellParSet(self))
        self.scatterers = []

        self._update = False

        cdict = {}
        for s in stru.scatterers():
            el = s.element_symbol()
            i = cdict.get(el, 0)
            sname = "%s%i"%(el,i)
            cdict[el] = i+1
            scatterer = CCTBXScattererParSet(sname, self, i)
            self.addParameterSet(scatterer)
            self.scatterers.append(scatterer)

        # Constrain the lattice
        from diffpy.srfit.structure.sgconstraints import _constrainSpaceGroup
        symbol = self.getSpaceGroup()
        _constrainSpaceGroup(self, symbol)

        return

    def update(self):
        """Update the unit_cell to a change in lattice parameters.

        This remakes the unit cell according to a change in the lattice
        parameters. Call this function before using the CCTBXCrystalParSet.
        The unit_cell will only be remade if necessary.

        """
        if not self._update:
            return

        self._update = False
        stru = self.stru
        sgn = stru.space_group().match_tabulated_settings().number()

        # Create the symmetry object
        from cctbx.crystal import symmetry
        symm = symmetry(
                unit_cell = self.unitcell._latpars,
                space_group_symbol = sgn
                )

        # Now the new structure
        newstru = stru.__class__(
                crystal_symmetry = symm,
                scatterers = stru.scatterers()
                )

        self.unitcell._latpars = list(newstru.unit_cell().parameters())

        self.stru = newstru
        return

    @classmethod
    def canAdapt(self, stru):
        """Return whether the structure can be adapted by this class."""
        try:
            from cctbx.crystal import special_position_settings
        except ImportError:
            return False
        return isinstance(stru, special_position_settings)

    def getLattice(self):
        """Get the ParameterSet containing the lattice Parameters."""
        return self.unitcell

    def getScatterers(self):
        """Get a list of ParameterSets that represents the scatterers.

        The site positions must be accessible from the list entries via the
        names "x", "y", and "z". The ADPs must be accessible as well, but the
        name and nature of the ADPs (U-factors, B-factors, isotropic,
        anisotropic) depends on the adapted structure.

        """
        return self.scatterers

    def getSpaceGroup(self):
        """Get the HM space group symbol for the structure."""
        sg = self.stru.space_group()
        t = sg.type()
        return t.lookup_symbol()



# End class CCTBXCrystalParSet
