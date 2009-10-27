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
Calculators are used to wrap registered functions so that the function's
Parameters are contained in an object specific to the function.  A custom
Calculator can be added to another RecipeOrganizer with the
'registerCalculator' method.

"""

from .parameterset import ParameterSet

class Calculator(ParameterSet):
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

def makeCalculator(f, name, argnames):
    """Make a calculator out of a function.

    This creates a Calculator whose __call__ method calls f, and whose bound
    signature is that of f. The Calculator has a Parameter for each of the
    arguments of f.

    Returns the Calculator instance.

    """
    c = Calculator(name)
    for n in argnames:
        c.newParameter(n, None)

    # FIXME - This uses the same mechanism as Equation. Is this what want?
    def _callf(*args):
        for i, val in enumerate(args):
            c._parameters.values()[i].setValue(val)

        vals = [p.value for p in self._parameters]
        return f(*vals)

    c.__call__ = _callf

    return c



__id__ = "$Id$"
