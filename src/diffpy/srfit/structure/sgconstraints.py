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
"""Code to set space group constraints for a crystal structure."""


import re

import numpy

from diffpy.srfit.fitbase.parameter import ParameterProxy
from diffpy.srfit.fitbase.recipeorganizer import RecipeContainer

__all__ = ["constrainAsSpaceGroup"]


def constrainAsSpaceGroup(
    phase,
    spacegroup,
    scatterers=None,
    sgoffset=[0, 0, 0],
    constrainlat=True,
    constrainadps=True,
    adpsymbols=None,
    isosymbol="Uiso",
):
    """Constrain the structure to the space group.

    This applies space group constraints to a StructureParSet with P1
    symmetry.  Passed scatterers are explicitly constrained to the
    specified space group. The ADPs and lattice may be constrained as well.

    Arguments:
    phase       --  A BaseStructure object.
    spacegroup  --  The space group number, symbol or an instance of
                    SpaceGroup class from diffpy.structure package.
    sgoffset    --  Optional offset for sg origin (default [0, 0, 0]).
    scatterers  --  The scatterer ParameterSets to constrain. If scatterers
                    is None (default), then all scatterers accessible from
                    phase.getScatterers will be constrained.
    constrainlat    --  Flag indicating whether to constrain the lattice
                    (default True).
    constrainadps   --  Flag indicating whether to constrain the ADPs
                    (default True).
    adpsymbols  --  A list of the ADP names. By default this is equal to
                    diffpy.structure.symmetryutilities.stdUsymbols (U11,
                    U22, etc.). The names must be given in the same order
                    as stdUsymbols.
    isosymbol   --  Symbol for isotropic ADP (default "Uiso"). If None,
                isotropic ADPs will be constrained via the anisotropic ADPs.

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
    """

    from diffpy.structure.spacegroups import GetSpaceGroup, SpaceGroup

    sg = spacegroup
    if not isinstance(spacegroup, SpaceGroup):
        sg = GetSpaceGroup(spacegroup)
    sgp = _constrainAsSpaceGroup(
        phase,
        sg,
        scatterers,
        sgoffset,
        constrainlat,
        constrainadps,
        adpsymbols,
        isosymbol,
    )

    return sgp


def _constrainAsSpaceGroup(
    phase,
    sg,
    scatterers=None,
    sgoffset=[0, 0, 0],
    constrainlat=True,
    constrainadps=True,
    adpsymbols=None,
    isosymbol="Uiso",
):
    """Restricted interface to constrainAsSpaceGroup.

    Arguments: As constrainAsSpaceGroup, except
    sg          --  diffpy.structure.spacegroups.SpaceGroup instance
    """

    from diffpy.structure.symmetryutilities import stdUsymbols

    if scatterers is None:
        scatterers = phase.getScatterers()
    if adpsymbols is None:
        adpsymbols = stdUsymbols

    sgp = SpaceGroupParameters(
        phase,
        sg,
        scatterers,
        sgoffset,
        constrainlat,
        constrainadps,
        adpsymbols,
        isosymbol,
    )

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
    """

    def __init__(self, name="sgpars"):
        """Create the BaseSpaceGroupParameters object.

        This initializes the attributes.
        """
        RecipeContainer.__init__(self, name)
        return

    def addParameter(self, par, check=True):
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
    sg      --  The diffpy.structure.spacegroups.SpaceGroup object
                corresponding to the space group.
    sgoffset    --  Optional offset for the space group origin.
    scatterers  --  The constrained scatterer ParameterSets.
    constrainlat    --  Flag indicating whether the lattice is constrained.
    constrainadps   --  Flag indicating whether the ADPs are constrained.
    adpsymbols  --  A list of the ADP names.
    _xyzpars    --  BaseSpaceGroupParameters of free xyz Parameters that are
                    constrained to.
    xyzpars     --  Property that populates _xyzpars.
    _latpars    --  BaseSpaceGroupParameters of free lattice Parameters that
                    are constrained to.
    latpars     --  Property that populates _latpars.
    _adppars    --  BaseSpaceGroupParameters of free ADPs that are constrained
                    to.
    adppars     --  Property that populates _adppars.
    """

    def __init__(
        self,
        phase,
        sg,
        scatterers,
        sgoffset,
        constrainlat,
        constrainadps,
        adpsymbols,
        isosymbol,
    ):
        """Create the SpaceGroupParameters object.

        Arguments:
        phase   --  A BaseStructure object to be constrained.
        sg      --  The space group number or symbol (compatible with
                    diffpy.structure.spacegroups.GetSpaceGroup.
        sgoffset    --  Optional offset for sg origin.
        scatterers  --  The scatterer ParameterSets to constrain. If scatterers
                    is None, then all scatterers accessible from
                    phase.getScatterers will be constrained.
        constrainlat    --  Flag indicating whether to constrain the lattice.
        constrainadps   --  Flag indicating whether to constrain the ADPs.
        adpsymbols  --  A list of the ADP names. The names must be given in the
                    same order as
                    diffpy.structure.symmetryutilities.stdUsymbols.
        isosymbol   --  Symbol for isotropic ADP (default "Uiso"). If None,
                    isotropic ADPs will be constrained via the anisotropic
                    ADPs.
        """
        BaseSpaceGroupParameters.__init__(self)
        self._latpars = None
        self._xyzpars = None
        self._adppars = None

        self._parsets = {}
        self._manage(self._parsets)

        self.phase = phase
        self.sg = sg
        self.sgoffset = sgoffset
        self.scatterers = scatterers
        self.constrainlat = constrainlat
        self.constrainadps = constrainadps
        self.adpsymbols = adpsymbols
        self.isosymbol = isosymbol

        return

    def __iter__(self):
        """Iterate over top-level parameters."""
        if (
            self._latpars is None
            or self._xyzpars is None
            or self._adppars is None
        ):
            self._makeConstraints()
        return RecipeContainer.__iter__(self)

    latpars = property(lambda self: self._getLatPars())

    def _getLatPars(self):
        """Accessor for _latpars."""
        if self._latpars is None:
            self._constrainLattice()
        return self._latpars

    xyzpars = property(lambda self: self._getXYZPars())

    def _getXYZPars(self):
        """Accessor for _xyzpars."""
        positions = []
        for scatterer in self.scatterers:
            xyz = [scatterer.x, scatterer.y, scatterer.z]
            positions.append([p.value for p in xyz])
        if self._xyzpars is None:
            self._constrainXYZs(positions)
        return self._xyzpars

    adppars = property(lambda self: self._getADPPars())

    def _getADPPars(self):
        """Accessor for _adppars."""
        positions = []
        for scatterer in self.scatterers:
            xyz = [scatterer.x, scatterer.y, scatterer.z]
            positions.append([p.value for p in xyz])
        if self._adppars is None:
            self._constrainADPs(positions)
        return self._adppars

    def _makeConstraints(self):
        """Constrain the structure to the space group.

        This works as described by the constrainAsSpaceGroup method.
        """

        # Start by clearing the constraints
        self._clearConstraints()

        scatterers = self.scatterers

        # Prepare positions
        positions = []
        for scatterer in scatterers:
            xyz = [scatterer.x, scatterer.y, scatterer.z]
            positions.append([p.value for p in xyz])

        self._constrainLattice()
        self._constrainXYZs(positions)
        self._constrainADPs(positions)

        return

    def _clearConstraints(self):
        """Clear old constraints.

        This only clears constraints where new ones are going to be
        applied.
        """
        phase = self.phase
        scatterers = self.scatterers
        isosymbol = self.isosymbol
        adpsymbols = self.adpsymbols

        # Clear xyz
        for scatterer in scatterers:

            for par in [scatterer.x, scatterer.y, scatterer.z]:
                if scatterer.isConstrained(par):
                    scatterer.unconstrain(par)
                par.setConst(False)

        # Clear the lattice
        if self.constrainlat:

            lattice = phase.getLattice()
            latpars = [
                lattice.a,
                lattice.b,
                lattice.c,
                lattice.alpha,
                lattice.beta,
                lattice.gamma,
            ]
            for par in latpars:
                if lattice.isConstrained(par):
                    lattice.unconstrain(par)
                par.setConst(False)

        # Clear ADPs
        if self.constrainadps:
            for scatterer in scatterers:
                if isosymbol:
                    par = scatterer.get(isosymbol)
                    if par is not None:
                        if scatterer.isConstrained(par):
                            scatterer.unconstrain(par)
                        par.setConst(False)

                for pname in adpsymbols:
                    par = scatterer.get(pname)
                    if par is not None:
                        if scatterer.isConstrained(par):
                            scatterer.unconstrain(par)
                        par.setConst(False)

        return

    def _constrainLattice(self):
        """Constrain the lattice parameters."""

        if not self.constrainlat:
            return

        phase = self.phase
        sg = self.sg

        lattice = phase.getLattice()
        system = sg.crystal_system
        if not system:
            system = "Triclinic"
        system = system.title()
        # This makes the constraints
        f = _constraintMap[system]
        f(lattice)

        # Now get the unconstrained, non-constant lattice pars and store them.
        self._latpars = BaseSpaceGroupParameters("latpars")
        latpars = [
            lattice.a,
            lattice.b,
            lattice.c,
            lattice.alpha,
            lattice.beta,
            lattice.gamma,
        ]
        pars = [p for p in latpars if not p.const and not p.constrained]
        for par in pars:
            # FIXME - the original parameter will still appear as
            # constrained.
            newpar = self.__addPar(par.name, par)
            self._latpars.addParameter(newpar)

        return

    def _constrainXYZs(self, positions):
        """Constrain the positions.

        positions   --  The coordinates of the scatterers.
        """

        from diffpy.structure.symmetryutilities import SymmetryConstraints

        sg = self.sg
        sgoffset = self.sgoffset

        # We do this without ADPs here so we can skip much complication. See
        # the _constrainADPs method for details.
        g = SymmetryConstraints(sg, positions, sgoffset=sgoffset)

        scatterers = self.scatterers
        self._xyzpars = BaseSpaceGroupParameters("xyzpars")

        # Make proxies to the free xyz parameters
        xyznames = [name[:1] + "_" + name[1:] for name, val in g.pospars]
        for pname in xyznames:
            name, idx = pname.rsplit("_", 1)
            idx = int(idx)
            par = scatterers[idx].get(name)
            newpar = self.__addPar(pname, par)
            self._xyzpars.addParameter(newpar)

        # Constrain non-free xyz parameters
        fpos = g.positionFormulas(xyznames)
        for idx, tmp in enumerate(zip(scatterers, fpos)):
            scatterer, fp = tmp

            # Extract the constraint equation from the formula
            for parname, formula in fp.items():
                _makeconstraint(
                    parname, formula, scatterer, idx, self._parameters
                )

        return

    def _constrainADPs(self, positions):
        """Constrain the ADPs.

        positions   --  The coordinates of the scatterers.
        """

        from diffpy.structure.symmetryutilities import (
            SymmetryConstraints,
            stdUsymbols,
        )

        if not self.constrainadps:
            return

        sg = self.sg
        sgoffset = self.sgoffset
        scatterers = self.scatterers
        isosymbol = self.isosymbol
        adpsymbols = self.adpsymbols
        adpmap = dict(zip(stdUsymbols, adpsymbols))
        self._adppars = BaseSpaceGroupParameters("adppars")

        # Prepare ADPs. Note that not all scatterers have constrainable ADPs.
        # For example, MoleculeParSet from objcryststructure does not. We
        # discard those.
        nonadps = []
        Uijs = []
        for sidx, scatterer in enumerate(scatterers):

            pars = [scatterer.get(symb) for symb in adpsymbols]

            if None in pars:
                nonadps.append(sidx)
                continue

            Uij = numpy.zeros((3, 3), dtype=float)
            for idx, par in enumerate(pars):
                i, j = _idxtoij[idx]
                Uij[i, j] = Uij[j, i] = par.getValue()

            Uijs.append(Uij)

        # Discard any positions for the nonadps
        positions = list(positions)
        nonadps.reverse()
        [positions.pop(idx) for idx in nonadps]

        # Now we can create symmetry constraints without having to worry about
        # the nonadps
        g = SymmetryConstraints(sg, positions, Uijs, sgoffset=sgoffset)

        adpnames = [adpmap[name[:3]] + "_" + name[3:] for name, val in g.Upars]

        # Make proxies to the free adp parameters. We start by filtering out
        # the isotropic ones so we can use the isotropic parameter.
        isoidx = []
        isonames = []
        for pname in adpnames:
            name, idx = pname.rsplit("_", 1)
            idx = int(idx)
            # Check for isotropic ADPs
            scatterer = scatterers[idx]
            if isosymbol and g.Uisotropy[idx] and idx not in isoidx:
                isoidx.append(idx)
                par = scatterer.get(isosymbol)
                if par is not None:
                    parname = "%s_%i" % (isosymbol, idx)
                    newpar = self.__addPar(parname, par)
                    self._adppars.addParameter(newpar)
                    isonames.append(newpar.name)
            else:
                par = scatterer.get(name)
                if par is not None:
                    newpar = self.__addPar(pname, par)
                    self._adppars.addParameter(newpar)

        # Constrain dependent isotropics
        for idx, isoname in zip(isoidx[:], isonames):
            for j in g.coremap[idx]:
                if j == idx:
                    continue
                isoidx.append(j)
                scatterer = scatterers[j]
                scatterer.constrain(isosymbol, isoname, ns=self._parameters)

        fadp = g.UFormulas(adpnames)

        # Constrain dependent anisotropics. We use the fact that an
        # anisotropic cannot be dependent on an isotropic.
        for idx, tmp in enumerate(zip(scatterers, fadp)):
            if idx in isoidx:
                continue
            scatterer, fa = tmp
            # Extract the constraint equation from the formula
            for stdparname, formula in fa.items():
                pname = adpmap[stdparname]
                _makeconstraint(
                    pname, formula, scatterer, idx, self._parameters
                )

    def __addPar(self, parname, par):
        """Constrain a parameter via proxy with a specified name.

        par     --  Parameter to constrain
        idx     --  Index to identify scatterer from which par comes
        """
        newpar = ParameterProxy(parname, par)
        self.addParameter(newpar)
        return newpar


# End class SpaceGroupParameters

# crystal system rules
# ref: Benjamin, W. A., Introduction to crystallography,
# New York (1969), p.60


def _constrainTriclinic(lattice):
    """Make constraints for Triclinic systems."""
    return


def _constrainMonoclinic(lattice):
    """Make constraints for Monoclinic systems.

    alpha and beta are fixed to 90 unless alpha != beta and alpha ==
    gamma, in which case alpha and gamma are constrained to 90.
    """
    afactor = 1
    if lattice.angunits == "rad":
        afactor = deg2rad
    ang90 = 90.0 * afactor
    lattice.alpha.setConst(True, ang90)
    beta = lattice.beta.getValue()
    gamma = lattice.gamma.getValue()

    if ang90 != beta and ang90 == gamma:
        lattice.gamma.setConst(True, ang90)
    else:
        lattice.beta.setConst(True, ang90)
    return


def _constrainOrthorhombic(lattice):
    """Make constraints for Orthorhombic systems.

    alpha, beta and gamma are constrained to 90
    """
    afactor = 1
    if lattice.angunits == "rad":
        afactor = deg2rad
    ang90 = 90.0 * afactor
    lattice.alpha.setConst(True, ang90)
    lattice.beta.setConst(True, ang90)
    lattice.gamma.setConst(True, ang90)
    return


def _constrainTetragonal(lattice):
    """Make constraints for Tetragonal systems.

    b is constrained to a and alpha, beta and gamma are constrained to
    90.
    """
    afactor = 1
    if lattice.angunits == "rad":
        afactor = deg2rad
    ang90 = 90.0 * afactor
    lattice.alpha.setConst(True, ang90)
    lattice.beta.setConst(True, ang90)
    lattice.gamma.setConst(True, ang90)
    lattice.constrain(lattice.b, lattice.a)
    return


def _constrainTrigonal(lattice):
    """Make constraints for Trigonal systems.

    If gamma == 120, then b is constrained to a, alpha and beta are
    constrained to 90 and gamma is constrained to 120. Otherwise, b and
    c are constrained to a, beta and gamma are constrained to alpha.
    """
    afactor = 1
    if lattice.angunits == "rad":
        afactor = deg2rad
    ang90 = 90.0 * afactor
    ang120 = 120.0 * afactor
    if lattice.gamma.getValue() == ang120:
        lattice.constrain(lattice.b, lattice.a)
        lattice.alpha.setConst(True, ang90)
        lattice.beta.setConst(True, ang90)
        lattice.gamma.setConst(True, ang120)
    else:
        lattice.constrain(lattice.b, lattice.a)
        lattice.constrain(lattice.c, lattice.a)
        lattice.constrain(lattice.beta, lattice.alpha)
        lattice.constrain(lattice.gamma, lattice.alpha)
    return


def _constrainHexagonal(lattice):
    """Make constraints for Hexagonal systems.

    b is constrained to a, alpha and beta are constrained to 90 and
    gamma is constrained to 120.
    """
    afactor = 1
    if lattice.angunits == "rad":
        afactor = deg2rad
    ang90 = 90.0 * afactor
    ang120 = 120.0 * afactor
    lattice.constrain(lattice.b, lattice.a)
    lattice.alpha.setConst(True, ang90)
    lattice.beta.setConst(True, ang90)
    lattice.gamma.setConst(True, ang120)
    return


def _constrainCubic(lattice):
    """Make constraints for Cubic systems.

    b and c are constrained to a, alpha, beta and gamma are constrained
    to 90.
    """
    afactor = 1
    if lattice.angunits == "rad":
        afactor = deg2rad
    ang90 = 90.0 * afactor
    lattice.constrain(lattice.b, lattice.a)
    lattice.constrain(lattice.c, lattice.a)
    lattice.alpha.setConst(True, ang90)
    lattice.beta.setConst(True, ang90)
    lattice.gamma.setConst(True, ang90)
    return


# This is used to map the correct crystal system to the proper constraint
# function.
_constraintMap = {
    "Triclinic": _constrainTriclinic,
    "Monoclinic": _constrainMonoclinic,
    "Orthorhombic": _constrainOrthorhombic,
    "Tetragonal": _constrainTetragonal,
    "Trigonal": _constrainTrigonal,
    "Hexagonal": _constrainHexagonal,
    "Cubic": _constrainCubic,
}


def _makeconstraint(parname, formula, scatterer, idx, ns={}):
    """Constrain a parameter according to a formula.

    parname     --  Name of parameter
    formula     --  Constraint formula
    scatterer   --  scatterer containing par of parname
    idx         --  Index to identify scatterer from which par comes
    ns          --  namespace to draw extra names from (default {})

    Returns the parameter if it is free.
    """
    par = scatterer.get(parname)

    if par is None:
        return

    compname = "%s_%i" % (parname, idx)

    # Check to see if this parameter is free
    pat = r"%s *([+-] *\d+)?$" % compname
    if re.match(pat, formula):
        return par

    # Check to see if it is a constant
    fval = _getFloat(formula)
    if fval is not None:
        par.setConst()
        return

    # If we got here, then we have a constraint equation
    # Fix any division issues
    formula = formula.replace("/", "*1.0/")
    scatterer.constrain(par, formula, ns=ns)
    return


def _getFloat(formula):
    """Get a float from a formula string, or None if this is not possible."""
    try:
        return eval(formula)
    except NameError:
        return None


# Constants needed above
_idxtoij = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
deg2rad = numpy.pi / 180
rad2deg = 1.0 / deg2rad


# End of file
