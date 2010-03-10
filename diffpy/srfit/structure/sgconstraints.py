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
"""Code to set space group constraints for a crystal structure."""

__id__ = "$Id$"


import numpy

from diffpy.Structure import SpaceGroups
from diffpy.Structure.SymmetryUtilities import GeneratorSite, stdUsymbols
from diffpy.Structure.SymmetryUtilities import SymmetryConstraints
from diffpy.srfit.fitbase.recipeorganizer import RecipeContainer
from diffpy.srfit.fitbase.parameter import Parameter

__all__ = [ "constrainAsSpaceGroup", "SpaceGroupParameters" ]

def constrainAsSpaceGroup(phase, sgsymbol, scatterers = None, 
        sgoffset = [0, 0, 0], constrainlat = True, constrainadps = True,
        adpsymbols = stdUsymbols):
    """Constrain the structure to the space group.

    This applies space group constraints to a StructureParSet with P1
    symmetry.  Passed scatterers are explicitly constrained to the
    specified space group. The ADPs and lattice may be constrained as well.

    Arguments:
    phase       --  A BaseStructure object.
    sgsymbol    --  The space group number or symbol (compatible with
                    diffpy.Structure.SpaceGroups.GetSpaceGroup.
    sgoffset    --  Optional offset for sg origin (default [0, 0, 0]).
    scatterers  --  The scatterer ParameterSets to constrain. If scatterers
                    is None (default), then all scatterers accessible from
                    phase.getScatterers will be constrained.
    constrainlat    --  Flag indicating whether to constrain the lattice
                    (default True).
    constrainadps   --  Flag indicating whether to constrain the ADPs
                    (default True).
    adpsymbols  --  A list of the ADP names. By default this is equal to
                    diffpy.Structure.SymmetryUtilities.stdUsymbols (U11,
                    U22, etc.). The names must be given in the same order
                    as stdUsymbols.

    New Parameters that are used in constraints are created within a
    SpaceGroupParameters object, which is returned from this function.
    Constraints are created in ParameterSet that contains the constrained
    Parameter.  This will erase any constraints or constant flags on the
    scatterers, lattice or ADPs if they are to be constrained.

    The lattice constraints are applied as following.
    
    Crystal System:
    Triclinic       --  No constraints.
    Monoclinic      --  alpha and beta are fixed to 90 unless alpha != beta and
                        alpha == gamma, in which case alpha and gamma are fixed
                        to 90.
    Orthorhombic    --  alpha, beta and gamma are fixed to 90.
    Tetragonal      --  b is constrained to a and alpha, beta and gamma are
                        fixed to 90.
    Trigonal        --  If gamma == 120, then b is constrained to a, alpha
                        and beta are fixed to 90 and gamma is fixed to 120.
                        Otherwise, b and c are constrained to a, beta and gamma
                        are fixed to alpha.
    Hexagonal       --  b is constrained to a, alpha and beta are fixed to 90
                        and gamma is fixed to 120.
    Cubic           --  b and c are constrained to a, and alpha, beta and
                        gamma are fixed to 90.

    Raises ValueError if phase is not in P1 symmetry.

    """

    phasesg = SpaceGroups.GetSpaceGroup( phase.getSpaceGroup() )
    if phasesg != SpaceGroups.sg1:
        raise ValueError("Structure is not in 'P1' symmetry")

    sg = SpaceGroups.GetSpaceGroup(sgsymbol)

    if scatterers is None:
        scatterers = phase.getScatterers()

    sgp = SpaceGroupParameters(phase, sg, scatterers, sgoffset,
            constrainlat, constrainadps, adpsymbols)

    return sgp

class SpaceGroupParameters(RecipeContainer):
    """This is a class for holding space group Parameters.

    This class is used by the constrainAsSpaceGroup method to store the
    Parameters needed to contain the space group. It has the same Parameter
    attribute access of a ParameterSet.

    Attributes
    phase   --  The constrained BaseStructure object.
    sg      --  The diffpy.Structure.SpaceGroups.SpaceGroup object
                corresponding to the space group.
    sgoffset    --  Optional offset for the space group origin.
    scatterers  --  The constrained scatterer ParameterSets to constrain.
    constrainlat    --  Flag indicating whether the lattice is constrained.
    constrainadps   --  Flag indicating whether the ADPs are constrained.
    adpsymbols  --  A list of the ADP names.
    xyzpars --  List of free xyz Parameters that are constrained to.
    latpars --  List of free lattice Parameters that are constrained to.
    adppars --  List of free ADPs that are constrained to.

    """

    _idxtoij = [(0,0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]

    def __init__(self, phase, sg, scatterers, sgoffset, constrainlat,
            constrainadps, adpsymbols):
        """Create the SpaceGroupParameters object.

        Arguments:
        phase   --  A BaseStructure object.
        sg      --  The space group number or symbol (compatible with
                    diffpy.Structure.SpaceGroups.GetSpaceGroup.
        sgoffset    --  Optional offset for sg origin.
        scatterers  --  The scatterer ParameterSets to constrain. If scatterers
                    is None, then all scatterers accessible from
                    phase.getScatterers will be constrained.
        constrainlat    --  Flag indicating whether to constrain the lattice.
        constrainadps   --  Flag indicating whether to constrain the ADPs.
        adpsymbols  --  A list of the ADP names. The names must be given in the
                    same order as
                    diffpy.Structure.SymmetryUtilities.stdUsymbols.

        """
        RecipeContainer.__init__(self, "sgpars")

        self.phase = phase
        self.sg = sg
        self.sgoffset = sgoffset
        self.scatterers = scatterers
        self.constrainlat = constrainlat
        self.constrainadps = constrainadps
        self.adpsymbols = adpsymbols

        self.xyzpars = []
        self.latpars = []
        self.adppars = []

        self.__makeConstraints()
        return

    def __makeConstraints(self):
        """Constrain the structure to the space group.
        
        This works as described by the constrainAsSpaceGroup method.

        """
        phase = self.phase
        sg = self.sg
        sgoffset = self.sgoffset
        scatterers = self.scatterers
        adpsymbols = self.adpsymbols

        ## Constrain the lattice
        if self.constrainlat:

            # First clear any constraints or constant variables in the lattice
            lattice = phase.getLattice()
            latpars =  [lattice.a, lattice.b, lattice.c, lattice.alpha,
                    lattice.beta, lattice.gamma]
            for par in latpars:
                lattice.unconstrain(par)
                par.setConst(False)

            system = sg.crystal_system
            if not system:
                system = "Triclinic"
            system = system.title()
            f = _constraintMap[system]
            f(lattice)

            # Now get the unconstrained, non-constant lattice pars and create
            # Parameters for them here. We will constrain those Parameters so
            # that there is no confusion about what this function does.
            self.latpars = []
            for par in latpars:
                if not par.const and not par.constrained:
                    newpar = Parameter("sg_" + par.name, par.getValue())
                    self._addObject(newpar, self._parameters)
                    lattice.constrain(par, newpar)
                    self.latpars.append(newpar)

        ## Constrain x, y, z
        # Remove any prior constraints or constants. We do this explicitly in
        # case the scatterer ParameterSet contains more than just the site
        # information.
        positions =  []
        for scatterer in scatterers:

            for par in [scatterer.x, scatterer.y, scatterer.z]:
                scatterer.unconstrain(par)
                par.setConst(False)

            positions.append([scatterer.x.getValue(), scatterer.y.getValue(),
                scatterer.z.getValue()])

        Uijs = None
        if self.constrainadps:
            Uijs = []
            for scatterer in scatterers:
                Uij = numpy.zeros((3,3), dtype=float)
                for idx, pname in enumerate(adpsymbols):
                    par = scatterer.get(pname)
                    scatterer.unconstrain(par)
                    par.setConst(False)
                    i, j = self.__class__._idxtoij[idx]
                    Uij[i,j] = Uij[j,i] = par.getValue()

                Uijs.append(Uij)

        g = SymmetryConstraints(sg, positions, Uijs, sgoffset=sgoffset)

        # Start by creating the position Parameters
        self.xyzpars = []
        xyznames = []
        for name, val in g.pospars:
            xyznames.append(name)
            par = Parameter(name, val)
            self._addObject(par, self._parameters)
            self.xyzpars.append(par)

        formulas = g.positionFormulas(xyznames)

        for scatterer, f in zip(scatterers, formulas):

            # Extract the constraint equation from the formula
            for parname, formula in f.items():

                # Check to see if this parameter is free
                if parname == formula:
                    continue

                par = scatterer.get(parname)

                # Check to see if it is a constant
                fval = _getFloat(formula)
                if fval is not None:
                    par.setConst()
                    continue

                # If we got here, then we have a constraint equation
                scatterer.constrain(par, formula, ns = self._parameters )

        # Now constrain ADPs
        self.adppars = []
        if not self.constrainadps:
            return

        adpnames = []
        for name, val in g.Upars:
            name = name[:3] + "_" + name[3:]
            adpnames.append(name)
            par = Parameter(name, val)
            self._addObject(par, self._parameters)
            self.adppars.append(par)

        formulas = g.UFormulas(adpnames)

        for scatterer, f in zip(scatterers, formulas):

            # Extract the constraint equation from the formula
            for parname, formula in f.items():

                # Check to see if this parameter is free
                if parname == formula:
                    continue

                par = scatterer.get(parname)
                # Check to see if it is a constant
                fval = _getFloat(formula)
                if fval is not None:
                    par.setConst()
                    continue

                # If we got here, then we have a constraint equation
                scatterer.constrain(par, formula, ns = self._parameters )

        return

# End SpaceGroupParameters

def _constrainSpaceGroup(phase, sg):
    """Constrain structure Parameters according to its space group.

    This is meant to constrain a StructureParSet that has internal space group
    symmetry.  This forces the lattice parameters to conform to the space group
    symmetry.  The protocol this follows is listed under Crystal Systems below.
    It also forces related symmetry positions to be constrained or held
    constant.

    Arguments:
    phase    --  A BaseStructure object.
    sg      --  The space group number or symbol (compatible with
                diffpy.Structure.SpaceGroups.GetSpaceGroup.

    The lattice constraints are applied as following.
    
    Crystal System:
    Triclinic       --  No constraints.
    Monoclinic      --  alpha and beta are fixed to 90 unless alpha != beta and
                        alpha == gamma, in which case alpha and gamma are fixed
                        to 90.
    Orthorhombic    --  alpha, beta and gamma are fixed to 90.
    Tetragonal      --  b is constrained to a and alpha, beta and gamma are
                        fixed to 90.
    Trigonal        --  If gamma == 120, then b is constrained to a, alpha
                        and beta are fixed to 90 and gamma is fixed to 120.
                        Otherwise, b and c are constrained to a, beta and gamma
                        are fixed to alpha.
    Hexagonal       --  b is constrained to a, alpha and beta are fixed to 90
                        and gamma is fixed to 120.
    Cubic           --  b and c are constrained to a, and alpha, beta and
                        gamma are fixed to 90.

    Note that lattice constraints are applied at the level of the lattice
    ParameterSet. The scatterer constraints are applied at the level of each
    scatterer ParameterSet.

    """
    sg = SpaceGroups.GetSpaceGroup(sg)

    ## Constrain the lattice
    # First clear any constraints or constant variables in the lattice
    lattice = phase.getLattice()
    for par in [lattice.a, lattice.b, lattice.c, lattice.alpha, lattice.beta,
            lattice.gamma]:
        lattice.unconstrain(par)
        par.setConst(False)

    system = sg.crystal_system
    if not system:
        system = "Triclinic"
    system = system.title()

    f = _constraintMap[system]
    f(lattice)

    ## Constrain related positions

    # Remove any prior constraints or constants. We do this explicitly in case
    # the scatterer ParameterSet contains more than just the scatterer
    # information.
    for scatterer in phase.getScatterers():

        for par in [scatterer.x, scatterer.y, scatterer.z]:
            scatterer.unconstrain(par)
            par.setConst(False)

    # We can go now if we're in P1 symmetry
    if sg == SpaceGroups.sg1:
        return

    # Now make a list of the positions and check for constraints
    for scatterer in phase.getScatterers():

        # Get the postion
        xyz = [scatterer.x.getValue(), scatterer.y.getValue(),
                scatterer.z.getValue()]

        # Get a formula for this scatterer
        g = GeneratorSite(sg, xyz)
        f = g.positionFormula(xyz, xyzsymbols=("x","y","z"))

        # Extract the constraint equation from the formula
        for parname, formula in f.items():

            # Check to see if this parameter is free
            if parname == formula:
                continue

            par = scatterer.get(parname)

            # Check to see if it is a constant
            fval = _getFloat(formula)
            if fval is not None:
                par.setConst()
                continue

            # If we got here, then we have a constraint equation
            scatterer.constrain(par, formula)

    return

def _constrainTriclinic(lattice):
    """Make constraints for Triclinic systems.

    """
    return

def _constrainMonoclinic(lattice):
    """Make constraints for Monoclinic systems.
    
    alpha and beta are fixed to 90 unless alpha != beta and alpha == gamma, in
    which case alpha and gamma are constrained to 90.

    """
    lattice.alpha.setConst(True, 90.0)
    beta = lattice.beta.getValue()
    gamma = lattice.gamma.getValue()

    if 90 != beta and beta == gamma:
        lattice.gamma.setConst(True, 90)
    else:
        lattice.beta.setConst(True, 90)
    return

def _constrainOrthorhombic(lattice):
    """Make constraints for Orthorhombic systems.
    
    alpha, beta and gamma are constrained to 90

    """
    lattice.alpha.setConst(True, 90.0)
    lattice.beta.setConst(True, 90.0)
    lattice.gamma.setConst(True, 90.0)
    return

def _constrainTetragonal(lattice):
    """Make constraints for Tetragonal systems.

    b is constrained to a and alpha, beta and gamma are constrained to 90.

    """
    lattice.alpha.setConst(True, 90.0)
    lattice.beta.setConst(True, 90.0)
    lattice.gamma.setConst(True, 90.0)
    lattice.constrain(lattice.b, lattice.a)
    return

def _constrainTrigonal(lattice):
    """Make constraints for Trigonal systems.

    If gamma == 120, then b is constrained to a, alpha and beta are
    constrained to 90 and gamma is constrained to 120. Otherwise, b and c
    are constrained to a, beta and gamma are constrained to alpha.

    """
    if lattice.gamma.getValue() == 120:
        lattice.constrain(lattice.b, lattice.a)
        lattice.alpha.setConst(True, 90.0)
        lattice.beta.setConst(True, 90.0)
        lattice.gamma.setConst(True, 120)
    else:
        lattice.constrain(lattice.b, lattice.a)
        lattice.constrain(lattice.c, lattice.a)
        lattice.constrain(lattice.beta, lattice.alpha)
        lattice.constrain(lattice.gamma, lattice.alpha)
    return

def _constrainHexagonal(lattice):
    """Make constraints for Hexagonal systems.

    b is constrained to a, alpha and beta are constrained to 90 and gamma is
    constrained to 120.

    """
    lattice.constrain(lattice.b, lattice.a)
    lattice.alpha.setConst(True, 90.0)
    lattice.beta.setConst(True, 90.0)
    lattice.gamma.setConst(True, 120.0)
    return

def _constrainCubic(lattice):
    """Make constraints for Cubic systems.

    b and c are constrained to a, alpha, beta and gamma are constrained to 90.

    """
    lattice.constrain(lattice.b, lattice.a)
    lattice.constrain(lattice.c, lattice.a)
    lattice.alpha.setConst(True, 90.0)
    lattice.beta.setConst(True, 90.0)
    lattice.gamma.setConst(True, 90.0)
    return

# This is used to map the correct crystal system to the proper constraint
# function.
_constraintMap = {
  "Triclinic"  : _constrainTriclinic,
  "Monoclinic" : _constrainMonoclinic,
  "Orthorhombic" : _constrainOrthorhombic, 
  "Tetragonal" : _constrainTetragonal,
  "Trigonal"   : _constrainTrigonal,
  "Hexagonal"  : _constrainHexagonal,
  "Cubic"      : _constrainCubic
}

def _getFloat(formula):
    """Get a float from a formula string, or None if this is not possible."""
    try:
        return float(formula)
    except ValueError:
        return None


