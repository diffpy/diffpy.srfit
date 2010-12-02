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
from diffpy.srfit.fitbase.fitresults import FitResults
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
    results         --  The built-in FitResults object.
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

    Properties
    names           --  Variable names (read only). See getNames.
    values          --  Variable values (read only). See getValues.

    """

    def __init__(self, name = "fit", conclass = FitContribution):
        """Initialization."""
        FitRecipe.__init__(self, name)
        self.fithooks[0].verbose = 3
        contribution = conclass("contribution")
        self.profile = Profile()
        contribution.setProfile(self.profile)
        self.addContribution(contribution)
        self.results = FitResults(self, update = False)
        return

    # Profile methods
    def loadParsedData(self, parser):
        """Load parsed data from a ProfileParser.

        This sets the xobs, yobs, dyobs arrays as well as the metadata.

        """
        return self.profile.loadParsedData(parser)

    def setObservedProfile(self, xobs, yobs, dyobs = None):
        """Set the observed profile.

        Arguments
        xobs    --  Numpy array of the independent variable
        yobs    --  Numpy array of the observed signal.
        dyobs   --  Numpy array of the uncertainty in the observed signal. If
                    dyobs is None (default), it will be set to 1 at each
                    observed xobs.

        Raises ValueError if len(yobs) != len(xobs)
        Raises ValueError if dyobs != None and len(dyobs) != len(xobs)

        """
        return self.profile.setObservedProfile(xobs, yobs, dyobs)

    def setCalculationRange(self, xmin = None, xmax = None, dx = None):
        """Set the calculation range

        Arguments
        xmin    --  The minimum value of the independent variable.
                    If xmin is None (default), the minimum observed value will
                    be used. This is clipped to the minimum observed x.
        xmax    --  The maximum value of the independent variable.
                    If xmax is None (default), the maximum observed value will
                    be used. This is clipped to the maximum observed x.
        dx      --  The sample spacing in the independent variable. If dx is
                    None (default), then the spacing in the observed points
                    will be preserved.

        Note that xmin is always inclusive (unless clipped). xmax is inclusive
        if it is within the bounds of the observed data.

        raises AttributeError if there is no observed profile
        raises ValueError if xmin > xmax
        raises ValueError if dx > xmax-xmin
        raises ValueError if dx <= 0

        """
        return self.profile.setCalculationRange(xmin = None, xmax = None,
                dx = None)

    def setCalculationPoints(self, x):
        """Set the calculation points.

        Arguments
        x   --  A non-empty numpy array containing the calculation points. If
                xobs exists, the bounds of x will be limited to its bounds.

        This will create y and dy on the specified grid if xobs, yobs and
        dyobs exist.

        """
        return self.profile.setCalculationPoints(x)

    def loadtxt(self, *args, **kw):
        """Use numpy.loadtxt to load data.

        Arguments are passed to numpy.loadtxt. 
        unpack = True is enforced. 
        The first two arrays returned by numpy.loadtxt are assumed to be x and
        y.  If there is a third array, it is assumed to by dy. Any other arrays
        are ignored. These are passed to setObservedProfile.

        Raises ValueError if the call to numpy.loadtxt returns fewer than 2
        arrays.

        Returns the x, y and dy arrays loaded from the file

        """
        return self.profile.loadtxt(*args, **kw)

    # FitContribution
    def addProfileGenerator(self, gen, name = None):
        """Add a ProfileGenerator to be used by this FitContribution.

        The ProfileGenerator is given a name so that it can be used as part of
        the profile equation (see setEquation). This can be different from the
        name of the ProfileGenerator used for attribute access.
        FitContributions should not share ProfileGenerator instances. Different
        ProfileGenerators can share Parameters and ParameterSets, however.
        
        Calling addProfileGenerator sets the profile equation to call the
        calculator and if there is not a profile equation already.

        gen     --  A ProfileGenerator instance
        name    --  A name for the calculator. If name is None (default), then
                    the ProfileGenerator's name attribute will be used.

        Raises ValueError if the ProfileGenerator has no name.
        Raises ValueError if the ProfileGenerator has the same name as some
        other managed object.

        """
        return self.contribution.addProfileGenerator(gen, name = None)

    def setEquation(self, eqstr, ns = {}):
        """Set the profile equation for the FitContribution.

        This sets the equation that will be used when generating the residual.
        The equation will be usable within setResidualEquation as "eq", and it
        takes no arguments.

        eqstr   --  A string representation of the equation. Variables will be
                    extracted from this equation and be given an initial value
                    of 0.
        ns      --  A dictionary of Parameters, indexed by name, that are used
                    in the eqstr, but not registered (default {}).
        
        Raises ValueError if ns uses a name that is already used for a
        variable.

        """
        self.contribution.setEquation(eqstr, ns = {})
        # Extract variables
        for par in self.contribution:
            # Skip Profile  Parameters
            if par.name in ("x", "y", "dy"): continue
            if par.value is None:
                par.value = 0
            if par.name not in self._parameters:
                self.addVar(par)
        return

    def setResidualEquation(self, eqstr = None):
        """Set the residual equation for the FitContribution.

        eqstr   --  A string representation of the residual. If eqstr is None
                    (default), then the previous residual equation will be
                    used, or the chi2 residual will be used if that does not
                    exist.

        Two residuals are preset for convenience, "chiv" and "resv".
        chiv is defined such that dot(chiv, chiv) = chi^2.
        resv is defined such that dot(resv, resv) = Rw^2.
        You can call on these in your residual equation. Note that the quantity
        that will be optimized is the summed square of the residual equation.
        Keep that in mind when defining a new residual or using the built-in
        ones.

        Raises AttributeError if the Profile is not yet defined.
        Raises ValueError if eqstr depends on a Parameter that is not part of
        the FitContribution.

        """
        return self.contribution.setResidualEquation(eqstr = None)

    def evaluate(self):
        """Evaluate the contribution equation."""
        return self.contribution.evaluate()

    def __call__(self):
        """Evaluate the contribution equation."""
        return self.contribution.evaluate()

    def printResults(self, header = "", footer = ""):
        """Format and print the results.

        header  --  A header to add to the output (default "")
        footer  --  A footer to add to the output (default "")

        """
        self.results.printResults(header, footer, True)
        return

    def saveResults(self, filename, header = "", footer = ""):
        """Format and save the results.

        filename -  Name of the save file.
        header  --  A header to add to the output (default "")
        footer  --  A footer to add to the output (default "")

        """
        self.results.saveResults(filename, header, footer, True)

# End class SimpleRecipe

# version
__id__ = "$Id$"

#
# End of file
