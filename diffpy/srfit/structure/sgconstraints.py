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


from diffpy.Structure import SpaceGroups
from diffpy.Structure.SymmetryUtilities import GeneratorSite

def constrainSpaceGroup(stru, sg):
    """Constrain the Parameters according to its space group.

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

    def _getFloat(formula):
        try:
            return float(formula)
        except ValueError:
            return None

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

            # Check to see if it is a constant
            fval = _getFloat(formula)
            if fval is not None:
                site.get(parname).setConst()
                continue

            # If we got here, then we have a constraint equation
            site.constrain(parname, formula)

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

