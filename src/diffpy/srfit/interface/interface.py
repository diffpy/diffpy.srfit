#!/usr/bin/env python
##############################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2010 The Trustees of Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################
"""Interface enhancements for Parameter-type classes.

Most interface additions can thought of by considering the classes in
SrFit to be sets of parameters. "+=" adds a new parameter, when that
makes sense. "|=" is the 'union' of these sets, and in general is used
to combine different objects. See individual interface classes for
specifics.
"""

__all__ = [
    "ParameterInterface",
    "FitRecipeInterface",
    "RecipeOrganizerInterface",
]

import six

from diffpy.srfit.equation.literals.abcs import ArgumentABC


class ParameterInterface(object):
    """Mix-in class for enhancing the Parameter interface."""

    def __lshift__(self, v):
        """Set the value with ``<<``.

        Think of ``<<`` as injecting a value.

        Parameters
        ----------
        v : float or ArgumentABC
            The value, or Argument derivative, to assign.

        Returns
        -------
        ParameterInterface
            Return self after the value has been set.
        """
        if isinstance(v, ArgumentABC):
            self.value = v.value
        else:
            self.value = v
        return self


# End class ParameterInterface

# ----------------------------------------------------------------------------


class RecipeOrganizerInterface(object):
    """Mix-in class for enhancing the RecipeOrganizer interface."""

    def __imul__(self, args):
        """Constrain with ``*=``.

        Think of ``*`` as a push-pin. This accepts arguments for a
        single constraint.

        Parameters
        ----------
        args : tuple
            The arguments to pass to `constrain`.

        Returns
        -------
        RecipeOrganizerInterface
            Return self after applying the constraint.
        """
        _applyargs(args, self.constrain)
        return self

    def __imod__(self, args):
        """Restrain with ``%=``.

        Think of ``%`` as a loose rope. This accepts arguments for a
        single restraint.

        Parameters
        ----------
        args : tuple
            The arguments to pass to `restrain`.

        Returns
        -------
        RecipeOrganizerInterface
            Return self after applying the restraint.
        """
        _applyargs(args, self.restrain)
        return self

    def __iadd__(self, args):
        """Add a parameter with ``+=``.

        Think of ``+`` as addition of a Parameter. This accepts
        arguments for a single call to `_new_parameter` or
        `_add_parameter`.

        Parameters
        ----------
        args : tuple
            The arguments to pass to `_new_parameter` or
            `_add_parameter`.

        Returns
        -------
        RecipeOrganizerInterface
            Return self after adding the parameter.
        """

        # Want to detect _add_parameter or _new_parameter
        def f(*args):
            if isinstance(args[0], six.string_types):
                self._new_parameter(*args)
            else:
                self._add_parameter(*args)
            return

        _applyargs(args, f)
        return self


# End class RecipeOrganizerInterface

# ----------------------------------------------------------------------------


class FitRecipeInterface(object):
    """Mix-in class for enhancing the FitRecipe interface."""

    def __ior__(self, args):
        """Add a contribution with ``|=``.

        Think of ``|`` as the union of components. This accepts a
        single argument.

        Parameters
        ----------
        args : FitContribution
            The contribution to add.

        Returns
        -------
        FitRecipeInterface
            Return self after adding the contribution.
        """
        self.add_contribution(args)
        return self

    def __iadd__(self, args):
        """Add a variable with ``+=``.

        Think of ``+`` as addition of a variable. This accepts a
        single argument or an iterable of single arguments or
        argument tuples.

        Parameters
        ----------
        args
            The arguments to pass to `add_variable` or
            `create_new_variable`.

        Returns
        -------
        FitRecipeInterface
            Return self after adding the variable.
        """

        # Want to detect add_variable or create_new_variable
        def f(*args):
            if isinstance(args[0], six.string_types):
                self.create_new_variable(*args)
            else:
                self.add_variable(*args)
            return

        _applymanyargs(args, f)
        return self


# End class FitRecipeInterface

# Local helper functions -----------------------------------------------------


def _applymanyargs(args, f):
    """Apply arguments to a function.

    Args can be any of the following. arg (arg1, arg2, ...) ((arg1a,
    arg1b, ...), ...)
    """
    if not hasattr(args, "__iter__"):
        f(args)
        return

    for arg in args:
        if hasattr(arg, "__iter__"):
            f(*arg)
        else:
            f(arg)

    return


def _applyargs(args, f):
    """Apply arguments to a function.

    Args can be any of the following. arg (arg1, arg2, ...) ((arg1a,
    arg1b, ...), ...)
    """
    if not hasattr(args, "__iter__"):
        f(args)
    else:
        f(*args)
    return


# End of file
