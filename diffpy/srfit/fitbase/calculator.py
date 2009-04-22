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
"""The Calculator class for generating a profile.

The Calculator has the interface of a diffpy.srfit.equation.literals.Generator
so that it can be used within a diffpy.srfit.equation.Equation.

"""

from numpy import array

from .parameter import Parameter

from .modelorganizer import ModelOrganizer

class Calculator(ModelOrganizer):
    """Base class for profile calculators.

    Attributes
    args            --  List needed by Generator interface.
    clicker         --  A Clicker instance for recording changes in contained
                        Parameters and ModelOrganizers.
    literal         --  This is literal created or modified by the generate
                        method (by default, a Parameter)
    name            --  A name for this organizer.
    _constraints    --  A dictionary of Constraints, indexed by the constrained
                        Parameter. Constraints can be added using the
                        'constrain' method.
    _orgdict        --  A dictionary containing the Parameters and
                        ModelOrganizers indexed by name.
    _parameters     --  A list of parameters that this ModelOrganizer knows
                        about.
    _profile        --  A Profile instance that contains the calculation range
                        and will contain the calculated profile.
    _restraints     --  A set of Restraints. Restraints can be added using the
                        'restrain' method.
    _organizers     --  A list of ModelOrganizers that this ModelOrganizer
                        knows about.
    _eqfactory      --  A diffpy.srfit.equation.builder.EquationFactory
                        instance that is used to create constraints and
                        restraints from string

    """

    def __init__(self, name):
        """Initialize the attributes."""
        ModelOrganizer.__init__(self, name)
        self._profile = None
        self.literal = Parameter(name)
        self.args = []
        return

    # Overload me!
    def eval(self):
        """Evaluate the profile.

        This method takes no arguments. The values of the calculator Parameters
        should be changed directly. The independent variable to calculate over
        is defined in the _profile attribute.

        This method calculates the profile and stores it in _profile.y.  This
        method needs to be overloaded.

        """
        self._profile.ycalc = array( self._profile.x )
        return

    ## No need to overload anything below here
        
    def setProfile(self, profile):
        """Assign the profile.

        profile --  A Profile that specifies the calculation points and which
                    will store the calculated signal.

        """
        self._profile = profile
        return

    # Generator methods. 

    def generate(self, clicker):
        """Generate the signal and store it in the literal attribute.

        This method conforms to the Generator interface. It does not need to be
        overloaded. By default it creates an Parameter to store the results.

        clicker --  A Clicker instance for decision making. The clicker is
                    sent by the Evaluator class to indicate its state. If the
                    clicker is greater than or equal to the Calculator's
                    clicker, then the profile should be re-evaluated.

        """
        if self.clicker >= clicker:
            self.eval()
            self.literal.setValue( self._profile.ycalc )
            self.clicker.click()
        return

    def identify(self, visitor):
        """Identify this to a visitor.

        This method conforms to the Generator interface. This does not need to
        be overloaded.

        """
        visitor.onGenerator(self)
        return


__id__ = "$Id$"

if __name__ == "__main__":
    # Check to see if everything imports correctly
    pass
