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
"""FitContribution class. 

FitContributions are generate a residual function for a FitRecipe. A FitContribution
associates an Equation for generating a signal, optionally a Calculator that
helps in this, and a Profile that holds the observed and calculated signals.  

See the examples in the documention for how to use a FitContribution.

"""

from numpy import concatenate, sqrt, inf, dot

from diffpy.srfit.equation import Equation
from diffpy.srfit.equation.builder import EquationFactory

from .recipeorganizer import RecipeOrganizer, equationFromString
from .parameter import ParameterProxy

class FitContribution(RecipeOrganizer):
    """FitContribution class.

    FitContributions organize an Equation that calculates the signal, and a
    Profile that holds the signal. A Calculator can be used as well.
    Contraints and Restraints can be created as part of a FitContribution.

    Attributes
    clicker         --  A Clicker instance for recording changes in the
                        Parameters or the residual components.
    name            --  A name for this FitContribution.
    calculator      --  A Calculator instance for generating a signal
                        (optional). If a calculator is not defined, the equation
                        to refine must be set with the setEquation method.
    profile         --  A Profile that holds the measured (and calcuated)
                        signal.
    _constraints    --  A dictionary of Constraints, indexed by the constrained
                        Parameter. Constraints can be added using the
                        'constrain' method.

    _eq             --  The FitContribution equation that will be optimized.
    _eqfactory      --  A diffpy.srfit.equation.builder.EquationFactory
                        instance that is used to create constraints and
                        restraints from strings.
    _organizers     --  A reference to the Calcualtor's _organizers attribute.
    _orgdict        --  A reference to the Calculator's _orgdict attribute.
    _parameters     --  A reference to the Calculator's _parameters attribute.
    _restraints     --  A set of Restraints. Restraints can be added using the
                        'restrain' or 'confine' methods.
    _xname          --  The name of of the independent variable from the
                        profile (default None). 
    _yname          --  The name of of the observed profile (default None). 
    _dyname         --  The name of of the uncertainty in the observed profile
                        (default None).
    """

    def __init__(self, name):
        """Initialization."""
        RecipeOrganizer.__init__(self, name)
        self._eq = None
        self._reseq = None
        self.profile = None
        self.calculator = None
        self._calcname = None
        self._xname = None
        self._yname = None
        self._dyname = None
        return
    
    # Make some methods public that were protected
    addParameter = RecipeOrganizer._addParameter
    newParameter = RecipeOrganizer._newParameter
    removeParameter = RecipeOrganizer._removeParameter

    def setProfile(self, profile, xname = None, yname = None, dyname = None):
        """Assign the profile for this fitcontribution.

        This resets the current residual (see setResidualEquation).
        
        profile --  A Profile that specifies the calculation points and that
                    will store the calculated signal.
        xname   --  The name of the independent variable from the Profile. If
                    this is None (default), then the name specified by the
                    Profile for this parameter will be used.  This variable is
                    usable within the Equation with the specified name.
        yname   --  The name of the observed profile.  If this is None
                    (default), then the name specified by the Profile for this
                    parametere will be used.  This variable is usable within
                    the Equation with the specified name.
        dyname  --  The name of the uncertainty in the observed profile. If
                    this is None (default), then the name specified by the
                    Profile for this parametere will be used.  This variable is
                    usable within the Equation with the specified name.

        """
        seteq = self.profile is None

        # Clear previously watched parameters
        if self.profile is not None:
            self.removeParameter(self._orgdict[self._xname])
            self.removeParameter(self._orgdict[self._yname])
            self.removeParameter(self._orgdict[self._dyname])

        # Set the profile and add its parameters to this organizer.

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

        xpar = ParameterProxy(self._xname, self.profile.xpar)
        ypar = ParameterProxy(self._yname, self.profile.ypar)
        dypar = ParameterProxy(self._dyname, self.profile.dypar)
        self.addParameter(xpar, check=False)
        self.addParameter(ypar, check=False)
        self.addParameter(dypar, check=False)

        # If we have a calculator, set its profile as well, and assign the
        # default residual equation if we're setting the profile for the first
        # time.
        if self.calculator is not None:
            self.calculator.setProfile(profile)

        if self._eq is not None and seteq:
            self.setResidualEquation()

        return

    def setCalculator(self, calc, name = None):
        """Set the Calculator to be used by this FitContribution.

        The Calculator is given a name so that it can be used as part of the
        profile equation (see setEquation). This can be different from the name
        of the Calculator used for attribute access. Each fitcontribution should
        have its own calculator instance. Those calculators can share
        Parameters and ParameterSets, however.
        
        Calling setCalculator sets the profile equation to call the calculator
        and resets the residual equation (see setResidualEquation).

        calc    --  A Calculator instance
        name    --  A name for the calculator. If name is None (default), then
                    the Calculator's name attribute will be used.

        """
        self.calculator = calc

        if name is None:
            name = calc.name

        # Let the RecipeOrganizer structure know of the calculator
        self._addOrganizer(calc)

        # Register the calculator with the equation factory
        self._eqfactory.registerGenerator(name, calc)

        # Set the default fitting equation, which is just the calculator
        self.setEquation(name)

        # If we have a profile already, let the calculator know about it.
        if self.profile is not None:
            calc.setProfile(self.profile)
        return

    def setEquation(self, eqstr, makepars = True, ns = {}):
        """Set the profile equation for the FitContribution.

        This sets the equation that will be used when generating the residual
        for this FitContribution.  The equation will be usable within
        setResidualEquation as "eq", and it takes no arguments.  Calling
        setEquation resets the residual equation.

        eqstr   --  A string representation of the equation. Any Parameter
                    registered by addParameter or setProfile, or function
                    registered by setCalculator, registerFunction or
                    registerStringFunction can be can be used in the equation
                    by name.
        makepars    --  A flag indicating whether missing Parameters can be
                    created by the Factory (default True). If False, then the a
                    ValueError will be raised if there are undefined arguments
                    in the eqstr. 
        ns      --  A dictionary of Parameters, indexed by name, that are used
                    in the eqstr, but not part of the FitRecipe (default {}).
        
        Raises ValueError if ns uses a name that is already used for a
        variable.
        Raises ValueError if makepars is false and eqstr depends on a Parameter
        that is not in ns or part of the FitContribution.

        """
        seteq = self._eq is None
        self._eq = self.registerStringFunction(eqstr, "eq", makepars, ns)
        self._eq.clicker.addSubject(self.clicker)

        if seteq and self.profile is not None:
            self.setResidualEquation()
        return

    def setResidualEquation(self, eqstr = None):
        """Set the residual equation for the FitContribution.

        eqstr   --  A string representation of the residual. If eqstr is None
                    (default), then the chi2 residual will be used.

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
        if self.profile is None:
            raise AttributeError("Assign the Profile first")

        # Register some convenient residuals
        chivstr = "(eq - %s)/%s" % (self._yname, self._dyname)
        chiv = equationFromString(chivstr, self._eqfactory)
        self._swapAndRegister("chiv", chiv)

        resvstr = "(eq - %s)/sum(%s**2)**0.5" % (self._yname, self._yname)
        resv = equationFromString(resvstr, self._eqfactory)
        self._swapAndRegister("resv", resv)

        # Now set the residual to one of these or create a new one
        if eqstr is None:
            self._reseq = chiv
        else:
            self._reseq = equationFromString(eqstr, self._eqfactory)

        return

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



# version
__id__ = "$Id$"

#
# End of file
