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
"""FitContribution class.

FitContributions generate a residual function for a FitRecipe. A
FitContribution associates an Equation for generating a signal,
optionally one or more ProfileGenerators or Calculators that help in
this, and a Profile that holds the observed and calculated signals.

See the examples in the documentation for how to use a FitContribution.
"""

__all__ = ["FitContribution"]

from diffpy.srfit.exceptions import SrFitError
from diffpy.srfit.fitbase.parameter import ParameterProxy
from diffpy.srfit.fitbase.parameterset import ParameterSet
from diffpy.srfit.fitbase.profile import Profile
from diffpy.srfit.fitbase.recipeorganizer import equationFromString


class FitContribution(ParameterSet):
    """FitContribution class.

    FitContributions organize an Equation that calculates the signal, and a
    Profile that holds the signal. ProfileGenerators and Calculators can be
    used as well.  Constraints and Restraints can be created as part of a
    FitContribution.

    Attributes
    name            --  A name for this FitContribution.
    profile         --  A Profile that holds the measured (and calculated)
                        signal.
    _calculators    --  A managed dictionary of Calculators, indexed by name.
    _constraints    --  A set of constrained Parameters. Constraints can be
                        added using the 'constrain' methods.
    _generators     --  A managed dictionary of ProfileGenerators.
    _parameters     --  A managed OrderedDict of parameters.
    _restraints     --  A set of Restraints. Restraints can be added using the
                        'restrain' method.
    _parsets        --  A managed dictionary of ParameterSets.
    _eqfactory      --  A diffpy.srfit.equation.builder.EquationFactory
                        instance that is used to create constraints and
                        restraints from string
    _eq             --  The FitContribution equation that will be optimized.
    _reseq          --  The residual equation.
    _xname          --  Name of the x-variable
    _yname          --  Name of the y-variable
    _dyname         --  Name of the dy-variable

    Properties
    names           --  Variable names (read only). See getNames.
    values          --  Variable values (read only). See getValues.
    """

    def __init__(self, name):
        """Initialization."""
        ParameterSet.__init__(self, name)
        self._eq = None
        self._reseq = None
        self.profile = None
        self._xname = None
        self._yname = None
        self._dyname = None

        self._generators = {}
        self._manage(self._generators)
        return

    def setProfile(self, profile, xname=None, yname=None, dyname=None):
        """Assign the Profile for this FitContribution.

        profile --  A Profile that specifies the calculation points and that
                    will store the calculated signal.
        xname   --  The name of the independent variable from the Profile. If
                    this is None (default), then the name specified by the
                    Profile for this parameter will be used.  This variable is
                    usable within string equations with the specified name.
        yname   --  The name of the observed Profile.  If this is None
                    (default), then the name specified by the Profile for this
                    parameter will be used.  This variable is usable within
                    string equations with the specified name.
        dyname  --  The name of the uncertainty in the observed Profile. If
                    this is None (default), then the name specified by the
                    Profile for this parameter will be used.  This variable is
                    usable within string equations with the specified name.
        """
        # Enforce type of the profile argument
        if not isinstance(profile, Profile):
            emsg = "Argument must be an instance of the Profile class."
            raise TypeError(emsg)

        # Set the Profile and add its parameters to this organizer.
        self.profile = profile

        if xname is None:
            xname = self.profile.xpar.name
        if yname is None:
            yname = self.profile.ypar.name
        if dyname is None:
            dyname = self.profile.dypar.name

        self._xname = xname
        self._yname = yname
        self._dyname = dyname

        xpar = ParameterProxy(xname, self.profile.xpar)
        ypar = ParameterProxy(yname, self.profile.ypar)
        dypar = ParameterProxy(dyname, self.profile.dypar)
        self.addParameter(xpar, check=False)
        self.addParameter(ypar, check=False)
        self.addParameter(dypar, check=False)

        # If we have ProfileGenerators, set their Profiles.
        for gen in self._generators.values():
            gen.setProfile(profile)

        # If we have _eq, but not _reseq, set the residual
        if self._eq is not None and self._reseq is None:
            self.setResidualEquation("chiv")

        return

    def addProfileGenerator(self, gen, name=None):
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
        if name is None:
            name = gen.name

        # Register the generator with the equation factory and add it as a
        # managed object.
        self._eqfactory.registerOperator(name, gen)
        self._addObject(gen, self._generators, True)

        # If we have a profile, set the profile of the generator.
        if self.profile is not None:
            gen.setProfile(self.profile)

        # Make this our equation if we don't have one. This will set the
        # residual equation if necessary.
        if self._eq is None:
            self.setEquation(name)

        return

    def setEquation(self, eqstr, ns={}):
        """Set the profile equation for the FitContribution.

        This sets the equation that will be used when generating the residual
        for this FitContribution.  The equation will be usable within
        setResidualEquation as "eq", and it takes no arguments.

        eqstr   --  A string representation of the equation. Any Parameter
                    registered by addParameter or setProfile, or function
                    registered by setCalculator, registerFunction or
                    registerStringFunction can be can be used in the equation
                    by name. Other names will be turned into Parameters of this
                    FitContribution.
        ns      --  A dictionary of Parameters, indexed by name, that are used
                    in the eqstr, but not registered (default {}).

        Raises ValueError if ns uses a name that is already used for a
        variable.
        """
        # Build the equation instance.
        eq = equationFromString(eqstr, self._eqfactory, buildargs=True, ns=ns)
        eq.name = "eq"

        # Register any new Parameters.
        for par in self._eqfactory.newargs:
            self._addParameter(par)

        # Register eq as an operator
        self._eqfactory.registerOperator("eq", eq)
        self._eqfactory.wipeout(self._eq)
        self._eq = eq

        # Set the residual if we need to
        if self.profile is not None and self._reseq is None:
            self.setResidualEquation("chiv")

        return

    def getEquation(self):
        """Get math expression string for the active profile equation.

        Return normalized math expression or an empty string if profile
        equation has not been set yet.
        """
        from diffpy.srfit.equation.visitors import getExpression

        rv = ""
        if self._eq is not None:
            rv = getExpression(self._eq)
        return rv

    def setResidualEquation(self, eqstr):
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

        Raises SrFitError if the Profile is not yet defined.
        Raises ValueError if eqstr depends on a Parameter that is not part of
        the FitContribution.
        """
        if self.profile is None:
            raise SrFitError("Assign the Profile first")
        if self._eq is None:
            raise SrFitError("Assign the Equation first")

        chivstr = "(eq - %s)/%s" % (self._yname, self._dyname)
        resvstr = "(eq - %s)/sum(%s**2)**0.5" % (self._yname, self._yname)

        # Get the equation string if it is not defined
        if eqstr == "chiv":
            eqstr = chivstr
        elif eqstr == "resv":
            eqstr = resvstr

        reseq = equationFromString(eqstr, self._eqfactory)
        self._eqfactory.wipeout(self._reseq)
        self._reseq = reseq

        return

    def getResidualEquation(self):
        """Get math expression string for the active residual equation.

        Return normalized math formula or an empty string if residual
        equation has not been configured yet.
        """
        from diffpy.srfit.equation.visitors import getExpression

        rv = ""
        if self._reseq is not None:
            rv = getExpression(self._reseq, eqskip="eq$")
        return rv

    def residual(self):
        """Calculate the residual for this fitcontribution.

        When this method is called, it is assumed that all parameters have been
        assigned their most current values by the FitRecipe. This will be the
        case when being called as part of a FitRecipe refinement.

        The residual is by default an array chiv:
        chiv = (eq() - self.profile.y) / self.profile.dy
        The value that is optimized is dot(chiv, chiv).

        The residual equation can be changed with the setResidualEquation
        method.
        """
        # Assign the calculated profile
        self.profile.ycalc = self._eq()
        # Note that equations only recompute when their inputs are modified, so
        # the following will not recompute the equation.
        return self._reseq()

    def evaluate(self):
        """Evaluate the contribution equation and update profile.ycalc."""
        yc = self._eq()
        if self.profile is not None:
            self.profile.ycalc = yc
        return yc

    def _validate(self):
        """Validate my state.

        This performs profile validations. This performs
        ProfileGenerator validations. This validates _eq. This validates
        _reseq and residual.

        Raises SrFitError if validation fails.
        """
        self.profile._validate()
        ParameterSet._validate(self)

        # Try to get the value of eq.
        from diffpy.srfit.equation.visitors import validate

        try:
            validate(self._eq)
        except ValueError as e:
            raise SrFitError(e)
        if self._eq is None:
            raise SrFitError("_eq is None")

        val = None
        try:
            val = self._eq()
        except TypeError as e:
            raise SrFitError("_eq cannot be evaluated: %s" % e)

        if val is None:
            raise SrFitError("_eq evaluates to None")

        # Try to get the value for residual
        try:
            validate(self._reseq)
        except ValueError as e:
            raise SrFitError(e)
        try:
            val = self.residual()
        except TypeError:
            raise SrFitError("residual cannot be evaluated")
        if val is None:
            raise SrFitError("residual evaluates to None")
        return


# End of file
