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
"""Wrappers for interfacing cctbx crystal with SrFit.

This wraps a cctbx.crystal as a ParameterSet with a similar hierarchy, which
can then be used within a FitRecipe. Note that all manipulations to the
cctbx.crystal should be done before wrapping. Changes made to the cctbx.crystal
object after wrapping may not be reflected within the wrapper, which can have
unpredictable results during a structure refinement.

Classes:

CCTBXStructureParSet    --  Wrapper for cctbx.crystal
UnitCellParSet  --  Wrapper for the unit cell of cctbx.crystal
ScattererParSet --  Wrapper for cctbx.xray.scatterer

"""
__id__ = "$Id$"

from diffpy.srfit.fitbase.parameter import Parameter, ParameterWrapper
from diffpy.srfit.fitbase.parameterset import ParameterSet
from diffpy.srfit.structure.basestructure import BaseStructure

from cctbx import crystal

class ScattererParSet(ParameterSet):
    """A wrapper for cctbx.xray.scatterer

    This class derives from ParameterSet.

    Attributes:
    name        --  Name of the scatterer. The name is always of the form
                    "%s%i" % (element, number), where the number is the running
                    index of that element type (starting at 0).
    x (y, z)    --  Atom position in crystal coordinates (ParameterWrapper)
    occupancy   --  Occupancy of the atom on its crystal location
                    (ParameterWrapper)
    Uiso        --  Isotropic scattering factor (ParameterWrapper).
    
    """

    def __init__(self, name, strups, idx):
        """Initialize

        name    --  The name of this scatterer.
        strups  --  The CCTBXStructureParSet that contains the cctbx structure
        idx     --  The index of the scatterer in the structure.

        """
        ParameterSet.__init__(self, name)
        self.strups = strups
        self.idx = idx

        # x, y, z, occupancy
        self.addParameter(ParameterWrapper("x", None, self._xyzgetter(0),
            self._xyzsetter(0)))
        self.addParameter(ParameterWrapper("y", None, self._xyzgetter(1),
            self._xyzsetter(1)))
        self.addParameter(ParameterWrapper("z", None, self._xyzgetter(2),
            self._xyzsetter(2)))
        self.addParameter(ParameterWrapper("occupancy", None, self._getocc,
            self._setocc))
        self.addParameter(ParameterWrapper("Uiso", None, self._getuiso,
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

# End class ScattererParSet

class UnitCellParSet(ParameterSet):
    """A wrapper for cctbx unit_cell object.
    
    Attributes:
    name    --  Always "unitcell".
    a, b, c, alpha, beta, gamma --  Unit cell parameters (ParameterWrapper).
    
    """

    def __init__(self, strups):
        """Initialize

        strups  --  The CCTBXStructureParSet that contains the cctbx structure
                    and the unit cell we're wrapper.

        """
        ParameterSet.__init__(self, "unitcell")
        self.strups = strups
        self._latpars = list(self.strups.stru.unit_cell().parameters())

        self.addParameter(ParameterWrapper("a", None, self._latgetter(0),
            self._latsetter(0)))
        self.addParameter(ParameterWrapper("b", None, self._latgetter(1),
            self._latsetter(1)))
        self.addParameter(ParameterWrapper("c", None, self._latgetter(2),
            self._latsetter(2)))
        self.addParameter(ParameterWrapper("alpha", None, self._latgetter(3),
            self._latsetter(3)))
        self.addParameter(ParameterWrapper("beta", None, self._latgetter(4),
            self._latsetter(4)))
        self.addParameter(ParameterWrapper("gamma", None, self._latgetter(5),
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


# End class UnitCellParSet

# FIXME - Special positions should be constant.

class CCTBXStructureParSet(BaseStructure):
    """A wrapper for CCTBX structure.

    Attributes:
    stru        --  The adapted cctbx structure object.
    scatterers  --  The list of ScattererParSets.
    unitcell    --  The UnitCellParSet for the structure.
    
    """

    def __init__(self, stru, name):
        """Initialize

        stru    --  A CCTBX structure instance.
        name    --  A name for this

        """
        ParameterSet.__init__(self, name)
        self.stru = stru
        self.addParameterSet(UnitCellParSet(self))
        self.scatterers = []

        self._update = False

        cdict = {}
        for s in stru.scatterers():
            el = s.element_symbol()
            i = cdict.get(el, 0)
            sname = "%s%i"%(el,i)
            cdict[el] = i+1
            scatterer = ScattererParSet(sname, self, i)
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
        parameters. Call this function before using the CCTBXStructureParSet.
        The unit_cell will only be remade if necessary.

        """
        if not self._update:
            return

        self._update = False
        stru = self.stru
        sgn = stru.space_group().match_tabulated_settings().number()

        # Create the symmetry object
        symm = crystal.symmetry(
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
        return isinstance(stru, crystal)

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



# End class StructureParSet

