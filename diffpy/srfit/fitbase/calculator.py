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

Calculators encapsulate the evaluation and required Parameters and
ParameterSets of a profile calculator.  The Calculator class can be associated
with a FitContribution to help calculate a profile.  It implements the
diffpy.srfit.equation.literals.Generator interface so that it can be used
within a diffpy.srfit.equation.Equation.

To define a Calculator, one must implement the required Parameters and
ParameterSets as well as overload the __call__ method with the calculation. A
very simple example is
> class Gaussian(Calculator):
>
>    def __init__(self):
>        # Initialize and give this a name
>        Calculator.__init__(self, "g")
>        # Add amplitude, center and width parameters
>        self.newParameter("amp", 0)
>        self.newParameter("center", 0)
>        self.newParameter("width", 0)
>       
>    def __call__(self, x):
>        a = self.amp.getValue()
>        x0 = self.center.getValue()
>        w = self.width.getValue()
>        return a * exp(-0.5*((x-x0)/w)**2)

More examples can be found in the example directory of the documentation.


"""

from numpy import array, asarray

from .parameter import Parameter

from .recipeorganizer import RecipeOrganizer

class Calculator(RecipeOrganizer):
    """Base class for profile calculators.

    Attributes
    args            --  List needed by Generator interface.
    clicker         --  A Clicker instance for recording changes in contained
                        Parameters and RecipeOrganizers.
    literal         --  This is literal created or modified by the generate
                        method (by default, a Parameter)
    name            --  A name for this organizer.
    profile         --  A Profile instance that contains the calculation range
                        and will contain the calculated profile.
    meta            --  A dictionary of metadata needed by the calculator.
    _constraints    --  A dictionary of Constraints, indexed by the constrained
                        Parameter. Constraints can be added using the
                        'constrain' method.
    _orgdict        --  A dictionary containing the Parameters and
                        RecipeOrganizers indexed by name.
    _parameters     --  A list of parameters that this RecipeOrganizer knows
                        about.
    _restraints     --  A set of Restraints. Restraints can be added using the
                        'restrain' or 'confine' methods.
    _organizers     --  A list of RecipeOrganizers that this RecipeOrganizer
                        knows about.
    _eqfactory      --  A diffpy.srfit.equation.builder.EquationFactory
                        instance that is used to create constraints and
                        restraints from string

    """

    def __init__(self, name):
        """Initialize the attributes."""
        RecipeOrganizer.__init__(self, name)
        self.profile = None
        self.literal = Parameter(name)
        self.args = []
        self.meta = {}
        return

    # Make some methods public that were protected
    addParameter = RecipeOrganizer._addParameter
    newParameter = RecipeOrganizer._newParameter
    removeParameter = RecipeOrganizer._removeParameter
    addParameterSet = RecipeOrganizer._addOrganizer
    removeParameterSet = RecipeOrganizer._removeOrganizer

    # Overload me!
    def __call__(self, x):
        """Evaluate the profile.

        This method must be overloaded.

        This method only takes the independent variables to calculate over. The
        values of the calculator Parameters should be changed directly. 

        """
        return x

    ## No need to overload anything below here

    def eval(self):
        """Evaluate the profile.

        This method takes no arguments. The values of the calculator Parameters
        should be changed directly. The independent variable to calculate over
        is defined in the profile attribute.

        This method calculates the profile and stores it in profile.y.

        """
        y = self.__call__(self.profile.x)
        self.profile.ycalc = asarray(y)
        return

        
    def setProfile(self, profile):
        """Assign the profile.

        profile --  A Profile that specifies the calculation points and which
                    will store the calculated signal.

        """
        # FIXME - When data parsers are implemented, this should use the
        # metadata to automatically configure the Calculator.
        if profile is not self.profile:

            # Stop watching the parameters
            if self.profile is not None:
                self.clicker.removeSubject(self.profile.xpar.clicker)
                self.clicker.removeSubject(self.profile.ypar.clicker)
                self.clicker.removeSubject(self.profile.dypar.clicker)

            # Set the profile and watch its parameters
            self.profile = profile
            self.clicker.addSubject(self.profile.xpar.clicker)
            self.clicker.addSubject(self.profile.ypar.clicker)
            self.clicker.addSubject(self.profile.dypar.clicker)

            self.clicker.click()
        return

    # Generator methods. These are used when the Calculator is called from
    # within an equation.

    def generate(self, clicker):
        """Generate the signal and store it in the literal attribute.

        This method is part of the Generator interface. It does not need to be
        overloaded. By default it creates a Parameter to store the results.

        clicker --  A Clicker instance for decision making. The clicker is
                    sent by the Evaluator class to indicate its state. If the
                    clicker is greater than or equal to the Calculator's
                    clicker, or that of the calculation arrays, then the
                    profile should be re-evaluated.

        """
        if self.clicker >= clicker or self.profile.xpar.clicker >= clicker:
            self.eval()
            self.literal.setValue( self.profile.ycalc )
        return

    def identify(self, visitor):
        """Identify this to a visitor.

        This method is part of the Generator interface. This does not need to
        be overloaded.

        """
        visitor.onGenerator(self)
        return


__id__ = "$Id$"
