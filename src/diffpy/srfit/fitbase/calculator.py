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
"""The Calculator for Parameter-aware functions.

Calculator is a functor class for producing a signal from embedded
Parameters. Calculators can store Parameters and ParameterSets,
Constraints and Restraints. Also, the __call__ function can be
overloaded to accept external arguments. Calculators are used to wrap
registered functions so that the function's Parameters are contained in
an object specific to the function.  A custom Calculator can be added to
another RecipeOrganizer with the 'registerCalculator' method.
"""

__all__ = ["Calculator"]

from diffpy.srfit.equation.literals.operators import Operator
from diffpy.srfit.fitbase.parameterset import ParameterSet


class Calculator(Operator, ParameterSet):
    """Base class for calculators.

    A Calculator organizes Parameters and has a __call__ method that can
    calculate a generic signal.

    Attributes
    name            --  A name for this organizer.
    meta            --  A dictionary of metadata needed by the calculator.
    _calculators    --  A managed dictionary of Calculators, indexed by name.
    _constraints    --  A set of constrained Parameters. Constraints can be
                        added using the 'constrain' methods.
    _parameters     --  A managed OrderedDict of contained Parameters.
    _parsets        --  A managed dictionary of ParameterSets.
    _restraints     --  A set of Restraints. Restraints can be added using the
                        'restrain' or 'confine' methods.
    _eqfactory      --  A diffpy.srfit.equation.builder.EquationFactory
                        instance that is used create Equations from string.

    Operator Attributes
    args    --  List of Literal arguments
    nin     --  Number of inputs (<1 means this is variable)
    nout    --  Number of outputs (1)
    operation   --  Function that performs the operation, self.__call__
    symbol  --  Same as name
    _value  --  The value of the Operator.
    value   --  Property for 'getValue'.

    Properties
    names           --  Variable names (read only). See getNames.
    values          --  Variable values (read only). See getValues.
    """

    # define abstract attributes from the Operator base.

    nin = -1
    nout = 1

    def __init__(self, name):
        ParameterSet.__init__(self, name)
        self.meta = {}

        # Initialize Operator attributes
        Operator.__init__(self, name)
        return

    @property
    def symbol(self):
        return self.name

    # Overload me!
    def __call__(self, *args):
        """Calculate something.

        This method must be overloaded. When overloading, you should
        specify the arguments explicitly, otherwise the parameters must
        be specified when adding the Calculator to a RecipeOrganizer.
        """
        return 0

    def operation(self, *args):
        self._value = self.__call__(*args)
        return self._value

    def _validate(self):
        """Validate my state.

        This performs ParameterSet validations. This does not validate
        the operation, since this could be costly. The operation should
        be validated with a containing equation.

        Raises AttributeError if validation fails.
        """
        ParameterSet._validate(self)

        return


# End class Calculator
