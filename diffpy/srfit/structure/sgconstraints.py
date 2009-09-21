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

def constrainAsSpaceGroup(stru, sg, sites = None, sgoffset = [0, 0, 0],
        constrainadps = False, adpsymbols = stdUsymbols):
    """Constrain structure parameters according to a given space group.

    This works as constrainSpaceGroup, but it is assumed that the structure has
    been expanded to a P1 cell. Passed sites are explicitly constrained to the
    specified space group.

    Arguments:
    stru    --  A BaseStructure object.
    sg      --  The space group number or symbol (compatible with
                diffpy.Structure.SpaceGroups.GetSpaceGroup.
    sgoffset--  Optional offset for sg origin (default [0, 0, 0]).
    sites   --  The site ParameterSets to constrain. If sites is None
                (default), then all sites from stru will be constrained. It is
                assumed that the sites have Parameters for the ADPs named
                according to the adpsymbols.
    constrainadps   --  Flag indicating whether to constrain the ADPs.
    adpsymbols  --  A list of the ADP names. By default this is equal to
                diffpy.Structure.SymmetryUtilities.stdUsymbols. The names must
                be given in the same order as stdUsymbols.

    Parameters for the free sites are created at the structure ParameterSet
    level.  The site constraints are applied at the level of each site
    ParameterSet.  This will erase any constraints or constants on the passed
    sites before setting the space group constraints. ADP constraints will not
    be wiped unless constrainadps is True.

    Returns xyznames and uijnames, where these are the names of the new x, y, z
    parameters are in xyznames and the names of the new Uij parameters are in
    uijnames.

    Raises ValueError if stru is not in P1 symmetry.
    Raises ValueError if constrainadps is True, but the adpsymbols are not
    found in sites.

    """
    # FIXME - default name of constrained parameters are kind of crappy.
    sg = SpaceGroups.GetSpaceGroup(sg)

    if sites is None:
        sites = stru.getSites()

    strusg = SpaceGroups.GetSpaceGroup( stru.getSpaceGroup() )
    if strusg != SpaceGroups.sg1:
        raise ValueError("Structure is not in 'P1' symmetry")

    ## Constrain sites

    # Remove any prior constraints or constants. We do this explicitly in case
    # the site ParameterSet contains more than just the site information.
    positions =  []
    for site in sites:

        for par in [site.x, site.y, site.z]:
            site.unconstrain(par)
            par.setConst(False)

        positions.append([site.x.getValue(), site.y.getValue(),
            site.z.getValue()])

    Uijs = None
    if constrainadps:
        idxtoij = [(0,0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
        Uijs = []
        for site in sites:
            Uij = numpy.matrix(numpy.zeros((3,3)), dtype=float)
            for idx, pname in enumerate(adpsymbols):
                par = site.get(pname)
                site.unconstrain(par)
                par.setConst(False)
                i, j = idxtoij[idx]
                Uij[i,j] = Uij[j,i] = par.getValue()

            Uijs.append(Uij)

    g = SymmetryConstraints(sg, positions, Uijs, sgoffset)

    # Start by creating the position Parameters
    xyznames = []
    for name, val in g.pospars:
        # We want new names
        name = name[:1] + "_sg" + name[1:]
        xyznames.append(name)
        par = stru.get(name)
        if par is None:
            par = stru.newParameter(name, val)
        else:
            par.setValue(val)

    formulas = g.positionFormulas(xyznames)

    for site, f in zip(sites, formulas):

        # Extract the constraint equation from the formula
        for parname, formula in f.items():

            # Check to see if this parameter is free
            if parname == formula:
                continue

            par = site.get(parname)

            # Check to see if it is a constant
            fval = _getFloat(formula)
            if fval is not None:
                par.setConst()
                continue

            # If we got here, then we have a constraint equation
            site.constrain(par, formula, ns = stru._parameters )

    symbolmap = dict(zip(stdUsymbols, adpsymbols))
    # Now constrain ADPs
    uijnames = []
    for name, val in g.Upars:
        base = name[:3]
        pre = symbolmap[base]
        name = pre + "_sg" + name[3:]
        uijnames.append(name)
        par = stru.get(name)
        if par is None:
            par = stru.newParameter(name, val)
        else:
            par.setValue(val)

    formulas = g.UFormulas(uijnames)

    for site, f in zip(sites, formulas):

        # Extract the constraint equation from the formula
        for parname, formula in f.items():

            # Check to see if this parameter is free
            if parname == formula:
                continue

            par = site.get(parname)

            # Check to see if it is a constant
            fval = _getFloat(formula)
            if fval is not None:
                par.setConst()
                continue

            # If we got here, then we have a constraint equation
            site.constrain(par, formula, ns = stru._parameters )

    return (xyznames, uijnames)

def constrainSpaceGroup(stru, sg):
    """Constrain structure Parameters according to its space group.

    This forces the lattice parameters to conform to the space group symmetry.
    The protocol this follows is listed under Crystal Systems below. It also
    forces related symmetry positions to be constrained or held constant.

    Arguments:
    stru    --  A BaseStructure object.
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
    ParameterSet. The site constraints are applied at the level of each site
    ParameterSet.

    """
    sg = SpaceGroups.GetSpaceGroup(sg)

    ## Constrain the lattice
    # First clear any constraints or constant variables in the lattice
    lattice = stru.getLattice()
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
    # the site ParameterSet contains more than just the site information.
    for site in stru.getSites():

        for par in [site.x, site.y, site.z]:
            site.unconstrain(par)
            par.setConst(False)

    # We can go now if we're in P1 symmetry
    if sg == SpaceGroups.sg1:
        return

    # Now make a list of the positions and check for constraints
    for site in stru.getSites():

        # Get the postion
        xyz = [site.x.getValue(), site.y.getValue(), site.z.getValue()]

        # Get a formula for this site
        g = GeneratorSite(sg, xyz)
        f = g.positionFormula(xyz, xyzsymbols=("x","y","z"))

        # Extract the constraint equation from the formula
        for parname, formula in f.items():

            # Check to see if this parameter is free
            if parname == formula:
                continue

            par = site.get(parname)

            # Check to see if it is a constant
            fval = _getFloat(formula)
            if fval is not None:
                par.setConst()
                continue

            # If we got here, then we have a constraint equation
            site.constrain(par, formula)

    return

def _constrainTriclinic(lattice):
    """Make constraints for Triclinic systems.

    """
    return

def _constrainMonoclinic(lattice):
    """Make constraints for Monoclinic systems.
    
    alpha and beta are constrained to 90 unless alpha != beta and alpha ==
    gamma, in which case alpha and gamma are constrained to 90.

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

    b and c are constrained to a, alpha, beta and gamma are constrained to
    90.

    """
    lattice.constrain(lattice.b, lattice.a)
    lattice.constrain(lattice.c, lattice.a)
    lattice.alpha.setConst(True, 90.0)
    lattice.beta.setConst(True, 90.0)
    lattice.gamma.setConst(True, 90.0)
    return

# This is used to map the correct crystal system to the proper constraint
# function. It is used by constrainSpaceGroup.
_constraintMap = {
  "Triclinic"  : _constrainTriclinic,
  "Monoclinic" : _constrainMonoclinic,
  "Orthorhombic" : _constrainOrthorhombic, 
  "Tetragonal" : _constrainTetragonal,
  "Trigonal"   : _constrainTrigonal,
  "Hexagonal"  : _constrainHexagonal,
  "Cubic"      : _constrainCubic
}


def applySymmetryConstraints(spacegroup, indices, posflag, Uijflag,
        sgoffset=[0,0,0]):
    """Generate symmetry constraints for positions and thermal factors.
    Both positions and thermal factors may get corrected to reflect
    space group symmetry.  Old positional and thermal constraints get
    erased.  New parameter indices start at fist decade after the last
    used parameter.

    spacegroup  -- instance of Structure.SpaceGroup
    indices     -- list of integer indices of atoms to be expanded
    posflag     -- required bool flag for constraining positions
    Uijflag     -- required bool flag for Uij constrainment
    sgoffset    -- optional offset of space group origin [0,0,0]
    """
    if not posflag and not Uijflag:     return
    # need to do something
    from diffpy.Structure.SymmetryUtilities import SymmetryConstraints
    # get unique sorted indices
    tobeconstrained = dict.fromkeys(indices)
    uindices = tobeconstrained.keys()
    uindices.sort()
    # remove old constraints
    pospat = re.compile(r'^([xyz])\((\d+)\)')
    Uijpat = re.compile(r'^(u11|u22|u33|u12|u13|u23)\((\d+)\)')
    for var in self.constraints.keys():
        mpos = posflag and pospat.match(var)
        mUij = Uijflag and Uijpat.match(var)
        if mpos and (int(mpos.group(2)) - 1) in tobeconstrained:
            del self.constraints[var]
        elif mUij and (int(mUij.group(2)) - 1) in tobeconstrained:
            del self.constraints[var]
    # find the largest used parameter index; pidxused must have an element
    pidxused = [i for i in self.owner.updateParameters()] + [0]
    # new parameters will start at the next decade
    firstpospar = firstUijpar = 10*(max(pidxused)/10) + 11
    # dictionary of parameter indices and their values
    newparvalues = {}
    selatoms = [self.initial[i] for i in uindices]
    selpos = [a.xyz for a in selatoms]
    selUijs = [a.U for a in selatoms]
    symcon = SymmetryConstraints(spacegroup, selpos, selUijs,
            sgoffset=sgoffset, eps=self.symposeps)
    # deal with positions
    if posflag:
        # fix positions:
        for a, xyz in zip(selatoms, symcon.positions):  a.xyz = xyz
        numpospars = len(symcon.pospars)
        posparindices = [i + firstpospar for i in range(numpospars)]
        posparvalues = symcon.posparValues()
        newparvalues.update( dict(zip(posparindices, posparvalues)) )
        possymbols = [ "@%i" % i for i in posparindices ]
        eqns = symcon.positionFormulasPruned(possymbols)
        for aidx, eq in zip(uindices, eqns):
            siteidx = aidx + 1
            sortedsmbls = eq.keys()
            sortedsmbls.sort()
            for barevar, formula in eq.items():
                var = barevar + "(%i)" % siteidx
                self.constraints[var] = Constraint(formula)
        # adjust firstUijpar
        if numpospars:
            firstUijpar += numpospars
            firstUijpar = 10*(firstUijpar/10) + 11
    # deal with temperature factors
    if Uijflag:
        # fix thermals
        for a, Uij in zip(selatoms, symcon.Uijs):  a.U = Uij
        numUpars = len(symcon.Upars)
        Uparindices = [i + firstUijpar for i in range(numUpars)]
        Uparvalues = symcon.UparValues()
        newparvalues.update( dict(zip(Uparindices, Uparvalues)) )
        Usymbols = [ "@%i" % i for i in Uparindices ]
        eqns = symcon.UFormulasPruned(Usymbols)
        for aidx, eq in zip(uindices, eqns):
            siteidx = aidx + 1
            sortedsmbls = eq.keys()
            sortedsmbls.sort()
            for barevar, formula in eq.items():
                # keys in formula dictionary are uppercase
                var = barevar.lower() + "(%i)" % siteidx
                self.constraints[var] = Constraint(formula)
    # update parameter values in parent Fitting
    self.owner.updateParameters()
    for pidx, pvalue in newparvalues.iteritems():
        parobj = self.owner.parameters[pidx]
        parobj.setInitial(pvalue)
    # and finally remember this space group
    self.initial.pdffit["spcgr"] = spacegroup.short_name
    self.initial.pdffit["sgoffset"] = list(sgoffset)
    return


def _getFloat(formula):
    """Get a float from a formula string, or None if this is not possible."""
    try:
        return float(formula)
    except ValueError:
        return None


