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
"""The Calculator for Parameter-aware functions.

Calculator is a functor class for producing a signal from embedded Parameters.
Calculators can store Parameters and ParameterSets, Constraints and Restraints.
Also, the __call__ function can be overloaded to accept external arguments.
This is useful when chaining together pieces of a forward calculation within a
FitContribution. A Calculator can be added to another RecipeOrganizer with the
'registerCalculator' method.

"""

from numpy import array, asarray

from .parameter import Parameter

from .parameterset import ParameterSet

class Calculator(ParameterSet):
    """Base class for calculators.

    A Calculator organizes Parameters and has a __call__ method that can
    calculate a generic signal.

    Attributes
    clicker         --  A Clicker instance for recording changes in contained
                        Parameters and RecipeOrganizers.
    name            --  A name for this organizer.
    meta            --  A dictionary of metadata needed by the calculator.
    _confclicker    --  A ConfigurationClicker for recording configuration
                        changes, esp.  additions and removal of managed objects.
    _calculators    --  A managed dictionary of Calculators, indexed by name.
    _constraints    --  A dictionary of Constraints, indexed by the constrained
                        Parameter. Constraints can be added using the
                        'constrain' method.
    _parameters     --  A managed OrderedDict of contained Parameters.
    _parsets        --  A managed dictionary of ParameterSets.
    _restraints     --  A set of Restraints. Restraints can be added using the
                        'restrain' or 'confine' methods.
    _eqfactory      --  A diffpy.srfit.equation.builder.EquationFactory
                        instance that is used create Equations from string.

    """

    def __init__(self, name):
        ParameterSet.__init__(self, name)
        self.meta = {}
        return

    # Overload me!
    def __call__(self, *args):
        """Calculate something.

        This method must be overloaded. When overloading, you should specify
        the arguments explicitly, otherwise the parameters must be specified
        when adding the Calculator to a RecipeOrganizer.

        """
        return 0

# End class Calculator

__id__ = "$Id$"
