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
        self._latpars = list(self.strups.stru.unit_cell().parameters())

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

        self.constrainSpaceGroup()

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

    def constrainSpaceGroup(self):
        """Constrain the lattice parameters according to the space group.
        
        This forces the lattice parameters to conform to the space group
        symmetry. The protocol this follows is listed under Crystal Systems
        below.

        Crystal Systems:
        Triclinic       --  No constraints.
        Monoclinic      --  alpha and beta are constrained to 90 unless alpha !=
                            beta and alpha == gamma, in which case alpha and
                            gamma are constrained to 90.
        Orthorhombic    --  alpha, beta and gamma are constrained to 90
        Tetragonal      --  b is constrained to a and alpha, beta and gamma are
                            constrained to 90.
        Trigonal        --  If gamma == 120, then b is constrained to a, alpha
                            and beta are constrained to 90 and gamma is
                            constrained to 120.  Otherwise, b and c are
                            constrained to a, beta and gamma are constrained to
                            alpha.
        Hexagonal       --  b is constrained to a, alpha and beta are
                            constrained to 90 and gamma is constrained to 120.
        Cubic           --  b and c are constrained to a, and alpha, beta and
                            gamma are constrained to 90.

        """
        # First clear any constraints or constant variables in the lattice
        self._constraints = {}
        for par in self._parameters:
            par.setConst(False)

        # Now get the crystal system
        system = self.strups.stru.space_group().crystal_system()
        if system == "Undefined":
            system = "Triclinic"

        constraintMap = {
          "Triclinic"  : self._constrainTriclinic,
          "Monoclinic" : self._constrainMonoclinic,
          "Orthorhombic" : self._constrainOrthorhombic, 
          "Tetragonal" : self._constrainTetragonal,
          "Trigonal"   : self._constrainTrigonal,
          "Hexagonal"  : self._constrainHexagonal,
          "Cubic"      : self._constrainCubic
        }


        # Note that we don't actuall constrain parameters, as that is
        # redundant. We will simply mark some parameters as constant if they
        # are constrained by cctbx.

        constraintMap[system]()
        return

    def _constrainTriclinic(self):
        """Make constraints for Triclinic systems.

        This frees the current value of all parameters.
        """
        return

    def _constrainMonoclinic(self):
        """Make constraints for Monoclinic systems.
        
        alpha and beta are constrained to 90 unless alpha != beta and alpha ==
        gamma, in which case alpha and gamma are constrained to 90.
        """
        self.alpha.setConst(True, 90.0)
        beta = self.alpha.getValue()
        gamma = self.alpha.getValue()

        if 90 != beta and beta == gamma:
            self.gamma.setConst(True, 90)
        else:
            self.beta.setConst(True, 90)
        return

    def _constrainOrthorhombic(self):
        """Make constraints for Orthorhombic systems.
        
        alpha, beta and gamma are constrained to 90
        """
        self.alpha.setConst(True, 90.0)
        self.beta.setConst(True, 90.0)
        self.gamma.setConst(True, 90.0)
        return

    def _constrainTetragonal(self):
        """Make constraints for Tetragonal systems.

        b is constrained to a and alpha, beta and gamma are constrained to 90.
        """
        self.alpha.setConst(True, 90.0)
        self.beta.setConst(True, 90.0)
        self.gamma.setConst(True, 90.0)
        self.constrain(self.b, self.a)
        return

    def _constrainTrigonal(self):
        """Make constraints for Trigonal systems.

        If gamma == 120, then b is constrained to a, alpha and beta are
        constrained to 90 and gamma is constrained to 120. Otherwise, b and c
        are constrained to a, beta and gamma are constrained to alpha.
        """
        if self.gamma.getValue() == 120:
            self.alpha.setConst(True, 90.0)
            self.beta.setConst(True, 90.0)
            self.gamma.setConst(True, 120)
            self.constrain(self.b, self.a)
        else:
            self.constrain(self.b, self.a)
            self.constrain(self.c, self.a)
            self.constrain(self.beta, self.alpha)
            self.constrain(self.gamma, self.alpha)
        return

    def _constrainHexagonal(self):
        """Make constraints for Hexagonal systems.

        b is constrained to a, alpha and beta are constrained to 90 and gamma is
        constrained to 120.
        """
        self.constrain(self.b, self.a)
        self.alpha.setConst(True, 90.0)
        self.beta.setConst(True, 90.0)
        self.gamma.setConst(True, 120.0)
        return

    def _constrainCubic(self):
        """Make constraints for Cubic systems.

        b and c are constrained to a, alpha, beta and gamma are constrained to
        90.
        """
        self.constrain(self.b, self.a)
        self.constrain(self.c, self.a)
        self.alpha.setConst(True, 90.0)
        self.beta.setConst(True, 90.0)
        self.gamma.setConst(True, 90.0)
        return


# End class UnitCellParSet

# FIXME - Special positions should be constant.

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
        self.addParameterSet(UnitCellParSet(self))
        self.scatterers = []

        self._update = False

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



# End class StructureParSet

