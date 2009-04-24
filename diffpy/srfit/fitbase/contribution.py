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
    clicker         --  A Clicker instance for recording changes in contained
                        Parameters and Contributions.
    name            --  A name for this Contribution.
    _calcname       --  A name for the Calculator.
    _calculator     --  A Calculator instance for generating a signal.
                        Contributions can share a Calculator instance.
    _constraints    --  A dictionary of Constraints, indexed by the constrained
                        Parameter. Constraints can be added using the
                        'constrain' method.
    _profile        --  A Profile that holds the measured (and calcuated)
                        signal.
    _eq             --  An Equation instance that generates a modified profile
                        with the Calculator.
    _eqfactory      --  A diffpy.srfit.equation.builder.EquationFactory
                        instance that is used to create constraints and
                        restraints from strings.
    _organizers     --  A reference to the Calcualtor's _organizers attribute.
    _orgdict        --  A reference to the Calculator's _orgdict attribute.
    _parameters     --  A reference to the Calculator's _parameters attribute.
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
        self._profile = None
        self._calculator = None
        self._calcname = None
        self._xname = None
        self._yname = None
        self._dyname = None
        return

    def setProfile(self, profile, xname = None, yname = None, dyname = None):
        """Assign the profile for this contribution.
        
        profile --  A Profile that specifies the calculation points and which
                    will store the calculated signal.
        xname   --  The name of the independent variable from the Profile
                    (default None). If this is provided, then the variable will
                    be usable within the Equation with the specified name.
        yname   --  The name of the observed profile (default None). If this is
                    provided, then the observed profile will be usable within
                    the Equation with the specified name.
        dyname  --  The name of the uncertainty in the observed profile
                    (default None). If this is provided, then the uncertainty
                    in the observed profile will be usable within the Equation
                    with the specified name.

        """
        self._profile = profile

        # Clear the previous profile information
        if self._xname is not None:
            self._eqfactory.deRegisterGenerator(self._xname)
            self._xname = None

        if self._yname is not None:
            self._eqfactory.deRegisterGenerator(self._yname)
            self._yname = None

        if self._dyname is not None:
            self._eqfactory.deRegisterGenerator(self._dyname)
            self._dyname = None

        # This creates a Generator that will produce the desired array in an
        # equation. This cannot be done with a Parameter alone because the user
        # may change the calculation points in the profile after it has been
        # added to the Contribution.
        def registerArray(name, attrname):
            g = Generator(name)
            g.literal = Parameter(name)

            def generate(clicker):
                a = getattr(self._profile, attrname)
                g.literal.setValue(a)
                return

            # Set the generate method of the generator to produce the array.
            g.generate = generate

            self._eqfactory.registerGenerator(name, g)
            return

        if xname is not None:
            self._xname = xname
            registerArray(xname, "x")

        if yname is not None:
            self._yname = yname
            registerArray(yname, "y")

        if dyname is not None:
            self._dyname = dyname
            registerArray(dyname, "dy")

        return

    def setCalculator(self, calc, name = None):
        """Set the Calculator to be used by this Contribution.

        The Calculator is given a name so that it can be used as part of the
        equation that is used to generate the signal. This can be different
        from the name of the Calculator for attribute purposes.

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

        return

    def setEquation(self, eqstr, ns = {}):
        """Set the refinement equation for the Contribution.

        eqstr   --  A string representation of the equation. The name of the
                    Calculator can be used in the equation. Any variables that
                    appear within eqstr will be added to the Contribution, and
                    will be accessible as attributes.
        ns      --  A dictionary of Parameters, indexed by name, that are used
                    in the eqstr, but not part of the Calculator (default
                    {}).

        Raises AttributeError if the Calculator is not yet defined.
        Raises ValueError if ns uses a name that is already used for a
        Parameter.
        Raises ValueError if eqstr depends on a Parameter that is not part of
        the Calculator and that is not defined in ns.
        """
        if self._calculator is None:
            raise AttributeError("Define the calculator first")

        self._eq = equationFromString(eqstr, self._eqfactory, ns)

        return

    def residual(self):
        """Calculate the residual for this contribution.

        It is assumed that all parameters have been assigned their most current
        values by the FitModel.

        The residual is by default an array chiv:
        chiv = (eq() - self._profile.y) / self._profile.dy
        
        """

        # Make sure that the calculator knows about the profile associated with
        # this contribution since multiple contributions may be using the same
        # calculator.
        self._calculator.setProfile(self._profile)
        self._profile.ycalc = self._eq()

        chiv = (self._profile.ycalc - self._profile.y) / self._profile.dy
        return chiv


# version
__id__ = "$Id$"

#
# End of file
