#!/usr/bin/env python
"""Wrappers for interfacing a diffpy.Structure.Structure as a ParameterSet
with the same hierarchy.

A diffpy.Structure.Structure object is meant to be passed to a Strucure object
from this module, which can then be used as a ParameterSet. Any change to the
lattice or existing atoms will be registered with the Structure. Changes in the
number of atoms will not be recognized. Thus, the diffpy.Structure.Structure
object should be fully configured before passing it to Structure.

CCTBXStructureParSet    --  Name required. Contains a  UnitParameterSet and
                    several ScattererParSet parameter sets.
UnitCellParSet  --  Named "unitcell". Contains Parameters "a", "b", "c",
                    "alpha", "beta", "gamma".
ScattererParSet --  Named "%s%i" % (element, number). Contains Parameters "x",
                    "y", "z", "occupancy", "uiso".

"""
__id__ = "$Id$"

from diffpy.srfit.fitbase.parameter import Parameter, ParameterProxy
from diffpy.srfit.fitbase.parameter import ParameterWrapper
from diffpy.srfit.fitbase.parameterset import ParameterSet

from cctbx import crystal

class ScattererParSet(ParameterSet):
    """A wrapper for cctbx.xray.scatterer

    This class derives from ParameterSet.

    Attributes:
    x (y, z)    --  Atom position in crystal coordinates (Parameter)
    occupancy   --  Occupancy of the atom on its crystal location (Parameter)
    uiso        --  Isotropic scattering factor.
    
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
        self.addParameter(ParameterWrapper(None, "x", self._xyzgetter(0),
            self._xyzsetter(0)))
        self.addParameter(ParameterWrapper(None, "y", self._xyzgetter(1),
            self._xyzsetter(1)))
        self.addParameter(ParameterWrapper(None, "z", self._xyzgetter(2),
            self._xyzsetter(2)))
        self.addParameter(ParameterWrapper(None, "occupancy", self._getocc,
            self._setocc))
        self.addParameter(ParameterWrapper(None, "uiso", self._getuiso,
            self._setuiso))
        return

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
    """A wrapper for cctbx unit_cell object."""

    def __init__(self, strups):
        """Initialize

        strups  --  The CCTBXStructureParSet that contains the cctbx structure

        """
        ParameterSet.__init__(self, "unitcell")
        self.strups = strups
        self.addParameter(ParameterWrapper(None, "a", self._latgetter(0),
            self._latsetter(0)))
        self.addParameter(ParameterWrapper(None, "b", self._latgetter(1),
            self._latsetter(1)))
        self.addParameter(ParameterWrapper(None, "c", self._latgetter(2),
            self._latsetter(2)))
        self.addParameter(ParameterWrapper(None, "alpha", self._latgetter(3),
            self._latsetter(3)))
        self.addParameter(ParameterWrapper(None, "beta", self._latgetter(4),
            self._latsetter(4)))
        self.addParameter(ParameterWrapper(None, "gamma", self._latgetter(5),
            self._latsetter(5)))

        return

    def _latgetter(self, i):

        def f(dummy):
            return self.strups.stru.unit_cell().parameters()[i]

        return f

    def _latsetter(self, i):

        def f(dummy, value):
            pars = list(self.strups.stru.unit_cell().parameters())
            pars[i] = value
            self._remakeStructure(tuple(pars))
            return

        return f

    def _remakeStructure(self, uc):
        """Remake the structure with new unit cell parameters."""
        stru = self.strups.stru
        sgn = stru.space_group().match_tabulated_settings().number()

        # Create the symmetry object
        symm = crystal.symmetry(
                unit_cell = uc,
                space_group_symbol = sgn
                )

        # Now the new structure
        newstru = stru.__class__(
                crystal_symmetry = symm,
                scatterers = stru.scatterers()
                )

        self.strups.stru = newstru
        return


# End class UnitCellParSet

class CCTBXStructureParSet(ParameterSet):
    """A wrapper for CCTBX structure.

    scatterers   --  The list of ScattererParSets.
    
    """

    def __init__(self, stru, name):
        """Initialize

        stru    --  A CCTBX structure instance.
        name    --  A name for this
        """
        ParameterSet.__init__(self, name)
        self.stru = stru
        self.addParameter(UnitCellParSet(self))
        self.scatterers = []

        cdict = {}
        for i, s in enumerate(stru.scatterers()):
            el = s.element_symbol()
            i = cdict.get(el, 0)
            sname = "%s%i"%(el,i)
            cdict[el] = i+1
            scatterer = ScattererParSet(sname, self, i)
            self.addParameterSet(scatterer)
            self.scatterers.append(scatterer)

        return

# End class StructureParSet

