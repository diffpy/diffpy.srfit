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

from __future__ import division

import numpy

from diffpy.Structure import SpaceGroups
from diffpy.Structure.SymmetryUtilities import GeneratorSite, stdUsymbols
from diffpy.Structure.SymmetryUtilities import SymmetryConstraints
from diffpy.srfit.fitbase.recipeorganizer import RecipeContainer
from diffpy.srfit.fitbase.parameter import Parameter, ParameterProxy

__all__ = ["constrainAsSpaceGroup"]

def constrainAsSpaceGroup(phase, sgsymbol, scatterers = None, 
        sgoffset = [0, 0, 0], constrainlat = True, constrainadps = True,
        adpsymbols = stdUsymbols, isosymbol = "Uiso"):
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
    isosymbol   --  Symbol for isotropic ADP (default "Uiso"). If None,
                isotropic ADPs will be constrainted via the anisotropic ADPs.

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
            constrainlat, constrainadps, adpsymbols, isosymbol)

    return sgp

# End constrainAsSpaceGroup

class BaseSpaceGroupParameters(RecipeContainer):
    """Base class for holding space group Parameters.

    This class is used to store the variable Parameters of a structure, leaving
    out those that constrained or fixed due to space group.  This class has the
    same Parameter attribute access of a ParameterSet. The purpose of this
    class is to make it easy to access the free variables of a structure for
    scripting purposes.

    Attributes
    name    --  "sgpars"
    xyzpars --  List of free xyz Parameters.
    latpars --  List of free lattice Parameters.
    adppars --  List of free ADPs. This does not include isotropic ADPs.

    """

    def __init__(self):
        """Create the BaseSpaceGroupParameters object.

        This initializes the attributes.

        """
        RecipeContainer.__init__(self, "sgpars")
        self.xyzpars = []
        self.latpars = []
        self.adppars = []
        return

    def addParameter(self, par, check = True):
        """Store a Parameter.

        par     --  The Parameter to be stored.
        check   --  If True (default), a ValueError is raised a Parameter of
                    the specified name has already been inserted.

        Raises ValueError if the Parameter has no name.

        """
        # Store the Parameter
        RecipeContainer._addObject(self, par, self._parameters, check)
        return

# End class BaseSpaceGroupParameters

class SpaceGroupParameters(BaseSpaceGroupParameters):
    """Class for holding and creating space group Parameters.

    This class is used to store the variable Parameters of a structure, leaving
    out those that constrained or fixed due to space group.  This does the work
    of the constrainAsSpaceGroup method.  This class has the same Parameter
    attribute access of a ParameterSet.

    Attributes
    name    --  "sgpars"
    phase   --  The constrained BaseStructure object.
    sg      --  The diffpy.Structure.SpaceGroups.SpaceGroup object
                corresponding to the space group.
    sgoffset    --  Optional offset for the space group origin.
    scatterers  --  The constrained scatterer ParameterSets.
    constrainlat    --  Flag indicating whether the lattice is constrained.
    constrainadps   --  Flag indicating whether the ADPs are constrained.
    adpsymbols  --  A list of the ADP names.
    xyzpars --  List of free xyz Parameters that are constrained to.
    latpars --  List of free lattice Parameters that are constrained to.
    adppars --  List of free ADPs that are constrained to.

    """

    def __init__(self, phase, sg, scatterers, sgoffset, constrainlat,
            constrainadps, adpsymbols, isosymbol):
        """Create the SpaceGroupParameters object.

        Arguments:
        phase   --  A BaseStructure object to be constrained.
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
        isosymbol   --  Symbol for isotropic ADP (default "Uiso"). If None,
                    isotropic ADPs will be constrainted via the anisotropic
                    ADPs.

        """
        BaseSpaceGroupParameters.__init__(self)

        self.phase = phase
        self.sg = sg
        self.sgoffset = sgoffset
        self.scatterers = scatterers
        self.constrainlat = constrainlat
        self.constrainadps = constrainadps
        self.adpsymbols = adpsymbols
        self.isosymbol = isosymbol

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
        isosymbol = self.isosymbol

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

            # Now get the unconstrained, non-constant lattice pars and add
            # them.
            self.latpars = []
            for par in latpars:
                if not par.const and not par.constrained:
                    # We don't have to make a proxy, but we do so for
                    # consistency.
                    newpar = ParameterProxy(par.name, par)
                    self.addParameter(newpar)
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
                    i, j = _idxtoij[idx]
                    Uij[i,j] = Uij[j,i] = par.getValue()

                if isosymbol:
                    par = scatterer.get(isosymbol)
                    scatterer.unconstrain(par)
                    par.setConst(False)

                Uijs.append(Uij)

        g = SymmetryConstraints(sg, positions, Uijs, sgoffset=sgoffset)

        # Make proxies to the free xyz parameters
        xyznames = [name[:1]+"_"+name[1:] for name, val in g.pospars]
        self.xyzpars = []
        for pname in xyznames:
            name, idx = pname.rsplit('_', 1)
            idx = int(idx)
            par = scatterers[idx].get(name)
            newpar = _addPar(par, idx, self)
            self.xyzpars.append(newpar)

        # Constrain non-free xyz parameters
        fpos = g.positionFormulas(xyznames)
        for idx, tmp in enumerate(zip(scatterers, fpos)):
            scatterer, fp = tmp

            # Extract the constraint equation from the formula
            for parname, formula in fp.items():
                _fixorfree(parname, formula, scatterer, idx, self._parameters)

        # Get out if we don't constrain adps
        if not self.constrainadps: return

        # Make proxies to the free adp parameters
        adpnames = [name[:3]+"_"+name[3:] for name, val in g.Upars]
        self.adppars = []

        # Make new parameters.
        isoidx = []
        isonames = []
        for pname in adpnames:
            name, idx = pname.rsplit('_', 1)
            idx = int(idx)
            # Check for isotropic ADPs
            scatterer = scatterers[idx]
            if isosymbol and g.Uisotropy[idx] and idx not in isoidx:
                isoidx.append(idx)
                par = scatterer.get(isosymbol)
                newpar = _addPar(par, idx, self)
                self.adppars.append(newpar)
                isonames.append(newpar.name)
            else:
                par = scatterers[idx].get(name)
                newpar = _addPar(par, idx, self)
                self.adppars.append(newpar)

        # Constrain dependent isotropics
        for idx, isoname in zip(isoidx[:], isonames):
            for j in g.coremap[idx]:
                if j == idx: continue
                isoidx.append(j)
                scatterer = scatterers[j]
                scatterer.constrain(isosymbol, isoname, ns = self._parameters)

        fadp = g.UFormulas(adpnames)
        adpmap = dict(zip(stdUsymbols, adpsymbols))

        # Constrain dependent anisotropics
        for idx, tmp in enumerate(zip(scatterers, fadp)):
            if idx in isoidx: continue
            scatterer, fa = tmp
            # Extract the constraint equation from the formula
            for stdparname, formula in fa.items():
                parname = adpmap[stdparname]
                _fixorfree(parname, formula, scatterer, idx, self._parameters)

        return

# End class SpaceGroupParameters

def _constrainSpaceGroup(phase, sg, adpsymbols = stdUsymbols, isosymbol =
        "Uiso"):
    """Constrain structure Parameters according to its space group.

    This constrains a StructureParSet that has internal space group symmetry.
    This forces the lattice parameters to conform to the space group symmetry
    according to the protocol listed under Crystal Systems below.  It also
    forces related symmetry positions and atomic displacement parameters to be
    constrained or held constant.

    Arguments:
    phase   --  A BaseStructure object.
    sg      --  The space group number or symbol (compatible with
                diffpy.Structure.SpaceGroups.GetSpaceGroup.
    adpsymbols  --  A list of the ADP names. The names must be given in the
                same order as diffpy.Structure.SymmetryUtilities.stdUsymbols.
    isosymbol   --  Symbol for isotropic ADP (default "Uiso"). If None,
                isotropic ADPs will be constrainted via the anisotropic ADPs.

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

    Returns a BaseSpaceGroupParameters object containing the variable
    parameters of the phase.  
    
    Note that lattice constraints are applied at the level of the lattice
    ParameterSet. The scatterer constraints are applied at the level of each
    scatterer ParameterSet.

    """
    sg = SpaceGroups.GetSpaceGroup(sg)
    adpmap = dict(zip(stdUsymbols, adpsymbols))

    ## Constrain the lattice
    # First clear any constraints or constant variables in the lattice
    lattice = phase.getLattice()
    latpars = [lattice.a, lattice.b, lattice.c, lattice.alpha, lattice.beta,
            lattice.gamma]
    for par in latpars:
        lattice.unconstrain(par)
        par.setConst(False)

    system = sg.crystal_system
    if not system:
        system = "Triclinic"
    system = system.title()

    # Create sgpars
    sgpars = BaseSpaceGroupParameters()

    # Constrain the lattice
    f = _constraintMap[system]
    f(lattice)

    # Grab free lattice parameters and add them to sgpars
    for par in latpars:
        if not par.const and not par.constrained:
            # We don't have to make a proxy, but we do so for consistency.
            newpar = ParameterProxy(par.name, par)
            sgpars.addParameter(newpar)
            sgpars.latpars.append(newpar)

    ## Constrain related positions

    # Now make a list of the positions and check for constraints
    for idx, scatterer in enumerate(phase.getScatterers()):

        # Remove prior constraints, get the postion
        xyz = []
        for pname in _xyzsymbols:
            par = scatterer.get(pname)
            scatterer.unconstrain(par)
            par.setConst(False)
            xyz.append(par.getValue())

        # Remove prior constraints, get the ADP. We must check whether the
        # scatterer has ADP parameters.
        Uij = numpy.zeros((3,3), dtype=float)
        hasadps = True
        for k, pname in enumerate(adpsymbols):
            par = scatterer.get(pname)
            if par is None: 
                hasadps = False
                break
            scatterer.unconstrain(par)
            par.setConst(False)
            i, j = _idxtoij[k]
            Uij[i,j] = Uij[j,i] = par.getValue()

        if isosymbol:
            par = scatterer.get(isosymbol)
            if par is not None:
                scatterer.unconstrain(par)
                par.setConst(False)

        # Get xyz and adp formulae for this scatterer
        g = GeneratorSite(sg, xyz, Uij)

        # Extract the xyz constraint equation from the formula
        fpos = g.positionFormula(xyz, xyzsymbols=_xyzsymbols)
        for parname, formula in fpos.items():
            par = _fixorfree(parname, formula, scatterer)
            if par is not None:
                newpar = _addPar(par, idx, sgpars)
                sgpars.xyzpars.append(newpar)

        # Extract the adp constraint equation from the formula
        # Check for isotropy
        if isosymbol and g.Uisotropy:
            # Constrain the isotropic parameter
            par = scatterer.get(isosymbol)
            if par is not None:
                newpar = _addPar(par, idx, sgpars)
                sgpars.adppars.append(newpar)
        elif hasadps:
            fadp = g.UFormula(xyz, Usymbols=adpsymbols)
            for stdparname, formula in fadp.items():
                parname = adpmap[stdparname]
                par = _fixorfree(parname, formula, scatterer)
                if par is not None:
                    newpar = _addPar(par, idx, sgpars)
                    sgpars.adppars.append(newpar)

    return sgpars

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

    if 90 != beta and 90 == gamma:
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

def _fixorfree(parname, formula, scatterer, idx = None, ns = {}):
    """Code used to fix or free a parameter.

    parname     --  Name of parameter
    formula     --  Constraint formula
    scatterer   --  scatterer containing par of parname
    idx         --  Index to identify scatterer from which par comes (default
                    None)
    ns          --  namespace to draw extra names from (default {})

    Returns the parameter if it is free.
    
    """
    par = scatterer.get(parname)

    compname = parname
    if idx is not None:
        compname += "_%i"%idx

    # Check to see if this parameter is free
    if compname == formula:
        return par

    # Check to see if it is a constant
    fval = _getFloat(formula)
    if fval is not None:
        par.setConst()
        return

    # If we got here, then we have a constraint equation
    scatterer.constrain(par, formula, ns = ns)
    return

def _addPar(par, idx, sgpars):
    """Constrain a parameter to sgpars via proxy with a specified name

    par     --  Parameter to constrain
    idx     --  Index to identify scatterer from which par comes
    sgpars  --  BaseSpaceGroupParameters object to add proxy to
    
    """
    name = "%s_%i"%(par.name, idx)
    newpar = ParameterProxy(name, par)
    sgpars.addParameter(newpar)
    return newpar

_idxtoij = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
_xyzsymbols = ('x', 'y', 'z')

def _getFloat(formula):
    """Get a float from a formula string, or None if this is not possible."""
    try:
        return eval(formula)
    except NameError:
        return None

# End of file

__id__ = "$Id$"
