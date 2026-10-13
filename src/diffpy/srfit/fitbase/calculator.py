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
another RecipeOrganizer with the 'register_calculator' method.
"""

__all__ = ["Calculator"]

from diffpy.srfit.equation.literals.operators import Operator
from diffpy.srfit.fitbase.parameterset import ParameterSet


class Calculator(Operator, ParameterSet):
    """Base class for calculators.

    A Calculator organizes Parameters and has a ``__call__`` method that
    can calculate a generic signal. Calculator is also an Operator
    (diffpy.srfit.equation.literals.operators), so it can be used directly
    in an evaluation network.

    Parameters
    ----------
    name : str
        The name of the Calculator. It must be a valid attribute name.

    Attributes
    ----------
    name : str
        The name of the Calculator.
    meta : dict
        The metadata needed by the calculator.
    names : list of str
        The names of the top-level Parameters (read only). See
        ``get_names``.
    values : list
        The values of the top-level Parameters (read only). See
        ``get_values``.
    args : list
        The Literal arguments of the operator, added with
        ``add_literal``.
    nin : int
        The number of inputs. The class value -1 means that any number of
        inputs is allowed.
    nout : int
        The number of outputs, which is 1.
    symbol : str
        The symbolic representation of the operator, which is the same as
        ``name`` (read only).
    value : object
        The value of the operator (read only). See ``get_value``.
    _value : object
        The cached value of the operator, or None if it must be
        recomputed.
    _calculators : dict
        The managed dictionary of Calculators, indexed by name.
    _constraints : dict
        The dictionary of Constraints, indexed by the constrained
        Parameter. Constraints can be added using the ``add_constraint``
        method.
    _parameters : OrderedDict
        The managed OrderedDict of contained Parameters.
    _parsets : OrderedDict
        The managed OrderedDict of contained ParameterSets.
    _restraints : set
        The set of Restraints. Restraints can be added using the
        ``add_soft_bounds`` method.
    _eqfactory : diffpy.srfit.equation.builder.EquationFactory
        The factory that is used to create Equations from strings.
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
        """The symbolic representation, which is the Calculator name."""
        return self.name

    # Overload me!
    def __call__(self, *args):
        """Calculate the signal produced by this Calculator.

        This method must be overloaded. When overloading, you should
        specify the arguments explicitly, otherwise the parameters must
        be specified when adding the Calculator to a RecipeOrganizer.

        Parameters
        ----------
        *args
            The arguments needed to calculate the signal.

        Returns
        -------
        object
            The calculated signal.
        """
        return 0

    def operation(self, *args):
        """Calculate and cache the signal produced by this Calculator.

        This calls ``__call__`` and stores the result as the operator
        value.

        Parameters
        ----------
        *args
            The arguments needed to calculate the signal.

        Returns
        -------
        object
            The calculated signal.
        """
        self._value = self.__call__(*args)
        return self._value

    def _validate(self):
        """Validate my state.

        This performs ParameterSet validations. This does not validate
        the operation, since this could be costly. The operation should
        be validated with a containing equation.

        Raises
        ------
        AttributeError
            If validation fails.
        """
        ParameterSet._validate(self)

        return


# End class Calculator
