#!/usr/bin/env python
########################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2008 Trustees of the Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
########################################################################
"""Simple FitRecipe class that includes a FitContribution and Profile.

"""
from diffpy.srfit.fitbase.fitrecipe import FitRecipe
from diffpy.srfit.fitbase.fitcontribution import FitContribution
from diffpy.srfit.fitbase.profile import Profile

class SimpleRecipe(FitRecipe):
    """SimpleRecipe class.

    This is a FitRecipe with a built-in Profile (the 'profile' attribute) and
    FitContribution (the 'contribution' attribute). Unique methods from each of
    these are exposed through this class to facilitate the creation of a simple
    fit recipe.

    Attributes
    profile         --  The built-in Profile object.
    contribution    --  The built-in FitContribution object.
    name            --  A name for this FitRecipe.
    fithook         --  An object to be called whenever within the residual
                        (default FitHook()) that can pass information out of
                        the system during a refinement.
    _constraints    --  A dictionary of Constraints, indexed by the constrained
                        Parameter. Constraints can be added using the
                        'constrain' method.
    _oconstraints   --  An ordered list of the constraints from this and all
                        sub-components.
    _calculators    --  A managed dictionary of Calculators.
    _contributions  --  A managed OrderedDict of FitContributions.
    _parameters     --  A managed OrderedDict of parameters (in this case the
                        parameters are varied).
    _parsets        --  A managed dictionary of ParameterSets.
    _eqfactory      --  A diffpy.srfit.equation.builder.EquationFactory
                        instance that is used to create constraints and
                        restraints from string
    _fixed          --  A set of parameters that are not actually varied.
    _restraintlist  --  A list of restraints from this and all sub-components.
    _restraints     --  A set of Restraints. Restraints can be added using the
                        'restrain' or 'confine' methods.
    _ready          --  A flag indicating if all attributes are ready for the
                        calculation.
    _tagdict        --  A dictionary of tags to variables.
    _weights        --  List of weighing factors for each FitContribution. The
                        weights are multiplied by the residual of the
                        FitContribution when determining the overall residual.

    """

    def __init__(self, name = "fit"):
        """Initialization."""
        FitRecipe.__init__(self, name)
        contribution = FitContribution("contribution")
        self.profile = Profile()

        contribution.setProfile(self.profile)
        self.addContribution(contribution)
        return

    def vary(self, pname, value = None, name = None, fixed = False, tag = None,
            tags = []):
        """Vary a Parameter from the contribution.
        
        pname   --  Name of a contribution parameter to vary

        Other arguments are as addVar.

        Returns the new variable.

        Raises ValueError if the contribution does not have a parameter named
        pname.

        """
        par = self.contribution.get(pname)
        if par is None:
            raise ValueError("No parameter named '%s'"%pname)
        return self.addVar(par, value, name, fixed, tag, tags)
    
    # Profile methods
    def loadParsedData(self, parser):
        """See Profile class"""
        return self.profile.loadParsedData(parser)

    def setObservedProfile(self, xobs, yobs, dyobs = None):
        """See Profile class"""
        return self.profile.setObservedProfile(xobs, yobs, dyobs)

    def setCalculationRange(self, xmin = None, xmax = None, dx = None):
        """See Profile class"""
        return self.profile.setCalculationRange(xmin = None, xmax = None,
                dx = None)

    def setCalculationPoints(self, x):
        """See Profile class"""
        return self.profile.setCalculationPoints(x)

    def loadtxt(self, *args, **kw):
        """See Profile class"""
        return self.profile.loadtxt(*args, **kw)

    # FitContribution
    def addProfileGenerator(self, gen, name = None):
        """See ProfileGenerator class"""
        return self.contribution.addProfileGenerator(gen, name = None)

    def setEquation(self, eqstr, ns = {}):
        """See ProfileGenerator class"""
        return self.contribution.setEquation(eqstr, ns = {})

    def setResidualEquation(self, eqstr = None):
        """See ProfileGenerator class"""
        return self.contribution.setResidualEquation(eqstr = None)

# End class SimpleRecipe

# version
__id__ = "$Id$"

#
# End of file
