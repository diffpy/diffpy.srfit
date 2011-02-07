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


import numpy

from diffpy.Structure import SpaceGroups
from diffpy.Structure.SymmetryUtilities import GeneratorSite, stdUsymbols
from diffpy.Structure.SymmetryUtilities import SymmetryConstraints
from diffpy.srfit.util import isVaried

__all__ = ["constrainAsSpaceGroup"]

def constrainAsSpaceGroup(phase, sgsymbol, scatterers = None, 
        sgoffset = [0, 0, 0], constrainlat = True, constrainadps = True):
    """Constrain the structure to the space group.

    This applies space group constraints to a StructureParSet with P1
    symmetry.  Passed scatterers are explicitly constrained to the
    specified space group. The ADPs and lattice may be constrained as well. All
    parameters allowed to vary by the space group will be set as such. This
    will remove existing constraints and variable status on the passed
    scatterers - the lattice parameters and ADPs will only be altered if these
    entities are to be constrained.

    Arguments:
    phase       --  An object providing the BaseStructure interface.
    sgsymbol    --  The space group number or symbol (compatible with
                    diffpy.Structure.SpaceGroups.GetSpaceGroup.
    sgoffset    --  Optional offset for sg origin (default [0, 0, 0]).
    scatterers  --  The adapted scatterers to constrain. If scatterers is None
                    (default), then all scatterers accessible from
                    phase.getScatterers will be constrained.
    constrainlat    --  Flag indicating whether to constrain the lattice
                    (default True).
    constrainadps   --  Flag indicating whether to constrain the ADPs
                    (default True).

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

    Returns a list of the variable parameters.

    """

    sg = SpaceGroups.GetSpaceGroup(sgsymbol)
    sgp = _constrainAsSpaceGroup(phase, sg, scatterers, sgoffset,
            constrainlat, constrainadps)

    return sgp

def _constrainAsSpaceGroup(phase, sg, scatterers = None, sgoffset = [0, 0, 0],
        constrainlat = True, constrainadps = True):
    """Restricted interface to constrainAsSpaceGroup.

    Arguments: As constrainAsSpaceGroup, except
    sg          --  diffpy.Structure.SpaceGroups.GetSpaceGroup instance

    """

    if scatterers is None:
        scatterers = phase.getScatterers()

    sgp = ConstrainSpaceGroup(phase, sg, scatterers, sgoffset,
            constrainlat, constrainadps)
    sgp.constrain()

    return sgp.variables

# End constrainAsSpaceGroup

class ConstrainSpaceGroup(object):
    """Class for creating space group parameters.

    This does the work of the constrainAsSpaceGroup method.

    Attributes
    phase   --  The constrained BaseStructure object.
    sg      --  The diffpy.Structure.SpaceGroups.SpaceGroup object
                corresponding to the space group.
    sgoffset    --  Optional offset for the space group origin.
    scatterers  --  The constrained scatterer ParameterSets.
    constrainlat    --  Flag indicating whether the lattice is constrained.
    constrainadps   --  Flag indicating whether the ADPs are constrained.
    variables   --  The list of varied parameters. Created by the 'constrain'
                    method.

    """

    def __init__(self, phase, sg, scatterers, sgoffset, constrainlat,
            constrainadps):
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

        """
        self.phase = phase
        self.sg = sg
        self.sgoffset = sgoffset
        self.scatterers = scatterers
        if self.scatterers is None:
            self.scatterers = phase.getScatterers()
        self.constrainlat = constrainlat
        self.constrainadps = constrainadps
        self.variables = []
        return

    def constrain(self):
        """Constrain the structure to the space group.
        
        This works as described by the constrainAsSpaceGroup method.

        """
        self.variables = []

        # Prepare positions
        scatterers = self.scatterers
        positions =  []
        for scatterer in scatterers:
            xyz = scatterer.getXYZ()
            positions.append([p.value for p in xyz])

        self._clearConstraints()
        self._constrainLattice()
        self._constrainXYZs(positions)
        self._constrainADPs(positions)
        return

    def _clearConstraints(self):
        """Clear old constraints.

        This only clears constraints on the lattice and ADPs if constraints for
        those entities has been requested.

        """
        phase = self.phase
        scatterers = self.scatterers

        # Clear xyz
        for scatterer in scatterers:
            for par in scatterer.getXYZ():
                if par.isConstrained():
                    par.unconstrain()
                par.fix()

        # Clear the lattice
        lattice = phase.getLattice()
        if self.constrainlat and lattice is not None:
            latpars = lattice.getLatPars()
            for par in latpars:
                if par.isConstrained():
                    par.unconstrain()
                par.fix()

        ## Clear ADPs
        if self.constrainadps:
            for scatterer in scatterers:
                adps = scatterer.getADPs()
                if adps is None: continue
                for par in adps:
                    if par is None: continue
                    if par.isConstrained():
                        par.unconstrain()
                    par.fix()

        return

    def _constrainLattice(self):
        """Constrain the lattice parameters."""

        if not self.constrainlat: return

        phase = self.phase
        sg = self.sg
        lattice = phase.getLattice()
        if lattice is None: return


        system = sg.crystal_system
        if not system:
            system = "Triclinic"
        system = system.title()
        # This makes the constraints
        f = _constraintMap[system]
        f(lattice)

        pars = filter(isVaried, lattice.getLatPars())
        self.variables.extend(pars)
        return

    def _constrainXYZs(self, positions):
        """Constrain the positions.

        positions   --  The coordinates of the scatterers.
        
        """

        sg = self.sg
        sgoffset = self.sgoffset
        xyzpars = {}

        # We do this without ADPs here so we can skip much complication. See
        # the _constrainADPs method for details.
        g = SymmetryConstraints(sg, positions, sgoffset=sgoffset)

        scatterers = self.scatterers

        # Store the free xyz parameters
        xyznames = [name for name, val in g.pospars]
        for forname in xyznames:
            pname = forname[:1]
            idx = int(forname[1:])
            scat = scatterers[idx]
            xyz = scat.getXYZ()
            j = _xyzidx[pname]
            par = xyz[j]
            xyzpars[forname] = par
            par.vary()

        self.variables.extend(xyzpars.values())

        # Constrain non-free xyz parameters
        fpos = g.positionFormulas(xyznames)
        for idx, tmp in enumerate(zip(scatterers, fpos)):
            scat, fp = tmp
            xyz = scat.getXYZ()

            # Extract the constraint equation from the formula
            for pname, formula in fp.items():
                j = _xyzidx[pname]
                par = xyz[j]
                # Get the formula name for this parameter.
                forname = "%s%i" % (pname, idx)
                _makeconstraint(par, forname, formula, xyzpars)

        return

    def _constrainADPs(self, positions):
        """Constrain the ADPs.

        positions   --  The coordinates of the scatterers.
        
        """
        if not self.constrainadps: return

        sg = self.sg
        sgoffset = self.sgoffset
        scatterers = self.scatterers
        adppars = {}

        # Prepare ADPs. Note that not all scatterers have constrainable ADPs.
        # For example, MoleculeParSet from objcryststructure does not. We
        # discard those.
        nonadps = []
        Uijs = []
        for sidx, scat in enumerate(scatterers):

            adps = scat.getADPs()

            if adps is None:
                nonadps.append(sidx)
                continue

            Uij = numpy.zeros((3,3), dtype=float)
            for idx, par in enumerate(adps[:-1]):
                i, j = _idxtoij[idx]
                Uij[i,j] = Uij[j,i] = par.value

            Uijs.append(Uij)

        # Discard any positions for the nonadps so we can create symmetry
        # constraints without having to worry about the nonadps.
        positions = list(positions)
        [positions.pop(idx) for idx in nonadps]

        g = SymmetryConstraints(sg, positions, Uijs, sgoffset=sgoffset)

        adpnames = [name for name, val in g.Upars]

        # Pull out the free adp parameters. We start by filtering out the
        # isotropic ones so we can use the isotropic parameter.
        isoidx = set([])
        for forname in adpnames:
            pname = forname[:3]
            idx = int(forname[3:])
            # Check for isotropic ADPs
            scat = scatterers[idx]
            adps = scat.getADPs()
            if adps is None: continue
            if idx not in isoidx:
                if g.Uisotropy[idx]:
                    par = adps[-1]
                    isoidx.add(idx)
                else:
                    paridx = _adpidx[pname]
                    par = adps[paridx]
                    adppars[forname] = par
                    par.vary()

        self.variables.extend(adppars.values())

        # Constrain dependent isotropics.
        for idx in list(isoidx):
            scat = scatterers[idx]
            adps = scat.getADPs()
            ipar = adps[-1]
            for j in g.coremap[idx]:
                scat = scatterers[j]
                adps = scat.getADPs()
                jpar = adps[-1]
                if j == idx: 
                    # This is independent. Let if vary and add it to the
                    # variables list.
                    jpar.vary()
                    self.variables.append(jpar)
                else:
                    jpar.constrain(ipar)
                    isoidx.add(j)

        fadp = g.UFormulas(adpnames)

        # Constrain dependent anisotropics. We use the fact that an
        # anisotropic cannot be dependent on an isotropic.
        for idx, tmp in enumerate(zip(scatterers, fadp)):
            if idx in isoidx: continue
            scat, fa = tmp
            adps = scat.getADPs()
            if adps is None: continue

            # Extract the constraint equation from the formula
            for pname, formula in fa.items():
                j = _adpidx[pname]
                par = adps[j]
                # Get the formula name for this parameter.
                forname = "%s%i" % (pname, idx)
                _makeconstraint(par, forname, formula, adppars)

        return

# End class SpaceGroupParameters

# crystal system rules
# ref: Benjamin, W. A., Introduction to crystallography,
# New York (1969), p.60

def _constrainTriclinic(lattice):
    """Make constraints for Triclinic systems.

    All lattice parameters are varied

    """
    a, b, c, alpha, beta, gamma = lattice.getLatPars()
    a.vary()
    b.vary()
    c.vary()
    alpha.vary()
    beta.vary()
    gamma.vary()
    return

def _constrainMonoclinic(lattice):
    """Make constraints for Monoclinic systems.
    
    alpha and beta are fixed to 90 unless alpha != beta and alpha == gamma, in
    which case alpha and gamma are constrained to 90
    a, b, and c are varied

    """
    a, b, c, alpha, beta, gamma = lattice.getLatPars()
    afactor = 1
    if lattice.angunits == "rad": afactor = deg2rad
    ang90 = 90.0 * afactor
    alpha.fix(ang90)
    beta = lattice.beta.getValue()
    gamma = lattice.gamma.getValue()

    if ang90 != beta.value and ang90 == gamma.value:
        gamma.fix(ang90)
    else:
        beta.fix(ang90)

    a.vary()
    b.vary()
    c.vary()
    return

def _constrainOrthorhombic(lattice):
    """Make constraints for Orthorhombic systems.
    
    alpha, beta and gamma are constrained to 90
    a, b, and c are varied

    """
    a, b, c, alpha, beta, gamma = lattice.getLatPars()
    afactor = 1
    if lattice.angunits == "rad": afactor = deg2rad
    ang90 = 90.0 * afactor
    alpha.fix(ang90)
    beta.fix(ang90)
    gamma.fix(ang90)
    a.vary()
    b.vary()
    c.vary()
    return

def _constrainTetragonal(lattice):
    """Make constraints for Tetragonal systems.

    alpha, beta and gamma are fixed to 90
    b is constrained to a and alpha 
    a and c are varied.

    """
    a, b, c, alpha, beta, gamma = lattice.getLatPars()
    afactor = 1
    if lattice.angunits == "rad": afactor = deg2rad
    ang90 = 90.0 * afactor
    alpha.fix(ang90)
    beta.fix(ang90)
    gamma.fix(ang90)
    b.constrain(a)
    a.vary()
    c.vary()
    return

def _constrainTrigonal(lattice):
    """Make constraints for Trigonal systems.

    if gamma == 120:
        alpha and beta are fixed to 90 
        gamma is fixed to 120 
        b is constrained to a
        a and c are varied
    else:
        beta and gamma are constrained to alpha
        b and c are constrained to a 
        a and alpha are varied

    """
    a, b, c, alpha, beta, gamma = lattice.getLatPars()
    afactor = 1
    if lattice.angunits == "rad": afactor = deg2rad
    ang90 = 90.0 * afactor
    ang120 = 120.0 * afactor
    if gamma.value == ang120:
        alpha.fix(ang90)
        beta.fix(ang90)
        gamma.fix(ang120)
        b.constrain(a)
        a.vary()
        c.vary()
    else:
        beta.constrain(alpha)
        gamma.constrain(alpha)
        b.constrain(a)
        c.constrain(a)
        a.vary()
        alpha.vary()
    return

def _constrainHexagonal(lattice):
    """Make constraints for Hexagonal systems.

    alpha and beta are fixed to 90 
    gamma is fixed to 120
    b is constrained to a 
    a and c ar varied

    """
    a, b, c, alpha, beta, gamma = lattice.getLatPars()
    afactor = 1
    if lattice.angunits == "rad": afactor = deg2rad
    ang90 = 90.0 * afactor
    ang120 = 120.0 * afactor
    alpha.fix(ang90)
    beta.fix(ang90)
    gamma.fix(ang120)
    b.constrain(a)
    a.vary()
    c.vary()
    return

def _constrainCubic(lattice):
    """Make constraints for Cubic systems.

    alpha, beta and gamma are fixed to 90
    b and c are constrained to a
    a is varied

    """
    a, b, c, alpha, beta, gamma = lattice.getLatPars()
    afactor = 1
    if lattice.angunits == "rad": afactor = deg2rad
    ang90 = 90.0 * afactor
    alpha.fix(ang90)
    beta.fix(ang90)
    gamma.fix(ang90)
    b.constrain(a)
    c.constrain(a)
    a.vary()
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

def _makeconstraint(par, forname, formula, ns):
    """Constrain a parameter according to a formula.

    par         --  Parameter to be constrained
    forname     --  Name of parameter as it appears in the constraint equation.
                    This is used to see if the parameter is constrained to
                    itself.
    formula     --  Constraint formula
    ns          --  Namespace from which to draw extra names. This consists of
                    parameters that appear in the constraint equations indexed
                    by the names used in those equations.

    Returns the parameter if it is free.
    
    """

    # Check to see if this parameter is free. We do this by seeing if it is
    # constrained to itself.
    import re
    pat = '%s *((\+|-) *\d+)?$'%forname
    if re.match(pat, formula):
        par.vary()
        return par

    # Check to see if it is a constant. We do this by trying to get a floating
    # point number out of the formula.
    fval = _getFloat(formula)
    if fval is not None:
        par.fix()
        return

    # If we got here, then we have a constraint equation. We will evaluate the
    # formula within the namespace, and then constrain the parameter to this.

    # First, fix potential division issues.
    formula = formula.replace("/", "*1.0/")
    # Evaluate the formula within the namespace.
    eq = eval(formula, ns)
    # Constrain the parameter.
    par.constrain(eq)
    return

def _getFloat(formula):
    """Get a float from a formula string, or None if this is not possible."""
    try:
        return eval(formula)
    except NameError:
        return None

# Constants needed above
_xyzidx = {"x" : 0, "y" : 1, "z" : 2}
_adpidx = {"U11" : 0, "U22" : 1, "U33" : 2, "U12" : 3, "U13" : 4, "U23" : 5,
"Uiso" : 6}
_idxtoij = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
deg2rad = numpy.pi / 180
rad2deg = 1.0 / deg2rad


# End of file

__id__ = "$Id$"
