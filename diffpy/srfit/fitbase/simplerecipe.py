#!/usr/bin/env python
##############################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2008 The Trustees of Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################

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

        # Adopt all the FitContribution methods
        public = [aname for aname in dir(contribution) if aname not in
                dir(self) and not aname.startswith("_")]
        for mname in public:
            method = getattr(contribution, mname)
            setattr(self, mname, method)
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


    def setCalculationRange(self, xmin=None, xmax=None, dx=None):
        """Set epsilon-inclusive calculation range.

        Adhere to the observed ``xobs`` points when ``dx`` is the same
        as in the data.  ``xmin`` and ``xmax`` are clipped at the bounds
        of the observed data.

        Parameters
        ----------

        xmin : float or "obs", optional
            The minimum value of the independent variable.  Keep the
            current minimum when not specified.  If specified as "obs"
            reset to the minimum observed value.
        xmax : float or "obs", optional
            The maximum value of the independent variable.  Keep the
            current maximum when not specified.  If specified as "obs"
            reset to the maximum observed value.
        dx : float or "obs", optional
            The sample spacing in the independent variable.  When different
            from the data, resample the ``x`` as anchored at ``xmin``.

        Note that xmin is always inclusive (unless clipped). xmax is inclusive
        if it is within the bounds of the observed data.

        Raises
        ------
        AttributeError
            If there is no observed data.
        ValueError
            When xmin > xmax or if dx <= 0.  Also if dx > xmax - xmin.
        """
        return self.profile.setCalculationRange(xmin, xmax, dx)


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

    def __call__(self):
        """Evaluate the contribution equation."""
        return self.contribution.evaluate()

    # FitResults methods

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

# End of file
