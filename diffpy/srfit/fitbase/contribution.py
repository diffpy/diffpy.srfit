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
"""Contribution class. 

Contributions organize an Equation and Calculator that calculate the signal,
and a Profile that holds the signal.
"""

from numpy import concatenate, sqrt, inf, dot

from diffpy.srfit.equation import Equation
from diffpy.srfit.equation.literals import Generator

from .parameter import Parameter
from .modelorganizer import ModelOrganizer, equationFromString


class Contribution(ModelOrganizer):
    """Contribution class.

    Contributions organize an Equation that calculates the signal, and a
    Profile that holds the signal. Contraints and Restraints can be created as
    part of a Contribution.

    Attributes
    clicker         --  A Clicker instance for recording changes in the
                        Parameters or the residual components.
    name            --  A name for this Contribution.
    profile         --  A Profile that holds the measured (and calcuated)
                        signal.
    _calcname       --  A name for the Calculator.
    _calculator     --  A Calculator instance for generating a signal.
                        Contributions can share a Calculator instance.
    _constraints    --  A dictionary of Constraints, indexed by the constrained
                        Parameter. Constraints can be added using the
                        'constrain' method.
    _eq             --  An Equation instance that generates a modified profile
                        with the Calculator.
    _eqfactory      --  A diffpy.srfit.equation.builder.EquationFactory
                        instance that is used to create constraints and
                        restraints from strings.
    _organizers     --  A reference to the Calcualtor's _organizers attribute.
    _orgdict        --  A reference to the Calculator's _orgdict attribute.
    _parameters     --  A reference to the Calculator's _parameters attribute.
    _eq             --  An Equation instance that generates the residual
                        equation.
    _restraints     --  A set of Restraints. Restraints can be added using the
                        'restrain' method.
    _xname          --  The name of of the independent variable from the
                        profile (default None). 
    _yname          --  The name of of the observed profile (default None). 
    _dyname         --  The name of of the uncertainty in the observed profile
                        (default None).
    """

    def __init__(self, name):
        """Initialization."""
        ModelOrganizer.__init__(self, name)
        self._eq = None
        self._reseq = None
        self.profile = None
        self._calculator = None
        self._calcname = None
        self._xname = None
        self._yname = None
        self._dyname = None
        return
    
    # Make some methods public that were protected
    addParameter = ModelOrganizer._addParameter
    newParameter = ModelOrganizer._newParameter
    removeParameter = ModelOrganizer._removeParameter

    def setProfile(self, profile, xname = None, yname = None, dyname = None):
        """Assign the profile for this contribution.

        This resets the current residual.
        
        profile --  A Profile that specifies the calculation points and which
                    will store the calculated signal.
        xname   --  The name of the independent variable from the Profile. If
                    this is None (default), then the name specified by the
                    Profile for this parametere will be used.  This variable is
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
        self.profile = profile

        # Clear the previous profile information
        self._eqfactory.deRegisterBuilder(self._xname)
        self._eqfactory.deRegisterBuilder(self._yname)
        self._eqfactory.deRegisterBuilder(self._dyname)

        if xname is None:
            xname = self.profile.xpar.name
        if yname is None:
            yname = self.profile.ypar.name
        if dyname is None:
            dyname = self.profile.dypar.name

        self._xname = xname
        self._eqfactory.registerArgument(xname, self.profile.xpar)
        self._yname = yname
        self._eqfactory.registerArgument(yname, self.profile.ypar)
        self._dyname = dyname
        self._eqfactory.registerArgument(dyname, self.profile.dypar)

        if self._calculator is not None:
            self.setResidualEquation()
        return

    def setCalculator(self, calc, name = None):
        """Set the Calculator to be used by this Contribution.

        The Calculator is given a name so that it can be used as part of the
        equation that is used to generate the signal. This can be different
        from the name of the Calculator for attribute purposes. This resets the
        current equation and residual.

        calc    --  A Calculator instance
        name    --  A name for the calculator. If name is None (default), then
                    the Calculator's name attribute will be used.

        """
        self._calculator = calc

        if name is None:
            name = calc.name

        self._calcname = name

        # Let the ModelOrganizer structure know of the calculator
        self._addOrganizer(calc)

        # Register the calculator with the equation factory
        self._eqfactory.registerGenerator(name, calc)

        # Create the default equation
        self._eq = self._eqfactory.makeEquation(name)

        self.setEquation(name)
        return

    def setEquation(self, eqstr, ns = {}):
        """Set the refinement equation for the Contribution.

        eqstr   --  A string representation of the equation. Any quantity
                    registered by setProfile and setCalculator can be can be
                    used in the equation by name.
        ns      --  A dictionary of Parameters, indexed by name, that are used
                    in the eqstr, but not part of the Contribution (default
                    {}).

        The equation will be usable within setResidualEquation by calling
        "eq()". Calling this function resets the residual equation.

        Raises AttributeError if the Calculator is not yet defined.
        Raises ValueError if ns uses a name that is already used for a
        Parameter.
        Raises ValueError if eqstr depends on a Parameter that is not part of
        the Contribution and that is not defined in ns.
        """
        if self._calculator is None:
            raise AttributeError("Define the calculator first")

        self._eq = equationFromString(eqstr, self._eqfactory, ns)

        # Now register this with the equation factory
        self._eqfactory.registerEquation("eq", self._eq)

        if self._calculator is not None and self.profile is not None:
            self.setResidualEquation()

        return

    def setResidualEquation(self, eqstr = None, ns = {}):
        """Set the residual equation for the Contribution.

        eqstr   --  A string representation of the residual. If eqstr is None
                    (default), then the chi2 residual will be used (see the
                    residual method.)
        ns      --  A dictionary of Parameters, indexed by name, that are used
                    in the eqstr, but not part of the Contribution (default
                    {}).

        Two residuals are preset for convenience, "chiv" and "resv".
        chiv is defined such that dot(chiv, chiv) = chi^2.
        resv is defined such that dot(resv, resv) = Rw.
        You can call on these in your residual equation. Note that the quantity
        that will be optimized is the summed square of the residual equation.
        Keep that in mind when defining a new residual or using the built-in
        ones.

        Raises AttributeError if either the Calculator or the Profile is not
        yet defined.
        Raises ValueError if ns uses a name that is already used for a
        Parameter.
        Raises ValueError if eqstr depends on a Parameter that is not part of
        the Contribution and that is not defined in ns.
        """
        if self._calculator is None:
            raise AttributeError("Define the calculator first")
        if self.profile is None:
            raise AttributeError("Define the profile first")

        # Register some convenient residuals
        chivstr = "(eq() - %s)/%s" % (self._yname, self._dyname)
        chiv = equationFromString(chivstr, self._eqfactory)
        self._eqfactory.registerEquation("chiv", chiv)

        resvstr = "(eq() - %s)/sum(%s**2)**0.5" % (self._yname, self._yname)
        resv = equationFromString(resvstr, self._eqfactory)
        self._eqfactory.registerEquation("resv", resv)

        # Now set the residual to one of these or create a new one
        if eqstr is None:
            self._reseq = chiv
        else:
            self._reseq = equationFromString(eqstr, self._eqfactory, ns)

        return

    def residual(self):
        """Calculate the residual for this contribution.

        It is assumed that all parameters have been assigned their most current
        values by the FitModel.

        The residual is by default an array chiv:
        chiv = (eq() - self.profile.y) / self.profile.dy
        The value that is optimized is dot(residual, residual).

        The residual equation can be changed with the setResidualEquation
        method.
        
        """
        # Make sure the calculator is working on my profile
        self._calculator.setProfile(self.profile)
        # Assign the calculated profile
        self.profile.ycalc = self._eq()
        # Note that equations only recompute when their inputs are modified, so
        # the following will not recompute the equation.
        return self._reseq()


# version
__id__ = "$Id$"

#
# End of file
