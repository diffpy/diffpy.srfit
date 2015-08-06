#!/usr/bin/env python
########################################################################
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
########################################################################
"""Interface enhancements for Parameter-type classes.

Most interface additions can thought of by considering the classes in SrFit to
be sets of parameters. "+=" adds a new parameter, when that makes sense. "|="
is the 'union' of these sets, and in general is used to combine different
objects. See individual interface classes for specifics.
"""

from __future__ import print_function
import six
__all__ = ["ParameterFactory", "ParameterInterface", "FitRecipeInterface",
        "RecipeOrganizerInterface"]

from diffpy.srfit.equation.literals.abcs import isinstance, ArgumentABC

class ParameterInterface(object):
    """Mix-in class for enhancing the Parameter interface."""

    def __lshift__(self, v):
        """setValue with <<

        Think of '<<' as injecting a value

        v   --  value or Argument derivative
        """
        if isinstance(v, ArgumentABC):
            self.value = v.value
        else:
            self.value = v
        return self

# End class ParameterInterface

class RecipeOrganizerInterface(object):
    """Mix-in class for enhancing the RecipeOrganizer interface."""

    def __imul__(self, args):
        """constrain with *=

        Think of '*' as a push-pin.

        This accepts arguments for a single constraint.
        """
        _applyargs(args, self.constrain)
        return self

    def __imod__(self, args):
        """restrain with %=

        This of '%' as a loose rope.

        This accepts arguments for a single restraint.
        """
        _applyargs(args, self.restrain)
        return self

    def __iadd__(self, args):
        """_newParameter or _addParameter with +=

        Think of "+" as addition of a Parameter.

        This accepts arguments for a single function call.
        """
        # Want to detect _addParameter or _newParameter
        def f(*args):
            if isinstance(args[0], six.string_types):
                self._newParameter(*args)
            else:
                self._addParameter(*args)
            return

        _applyargs(args, f)
        return self

# End class RecipeOrganizerInterface

class FitContributionInterface(object):
    """Mix-in class for enhancing the FitContribution interface."""

    def __ior__(self, args):
        """setProfile, addProfileGenerator or setEquation with |=

        Think of "|" as the union of components.

        This accepts arguments for a single Profile or ProfileGenerator.
        """
        def f(*args):
            from diffpy.srfit.fitbase.profile import Profile
            from diffpy.srfit.fitbase.profilegenerator import ProfileGenerator
            if isinstance(args[0], Profile):
                self.setProfile(*args)
            elif isinstance(args[0], ProfileGenerator):
                self.addProfileGenerator(*args)
            elif isinstance(args[0], six.string_types):
                self.setEquation(*args)
            else:
                raise TypeError("Invalid argument")
            return
        _applyargs(args, f)
        return self


# End class FitContributionInterface

class FitRecipeInterface(object):
    """Mix-in class for enhancing the FitRecipe interface."""

    def __ior__(self, args):
        """addContribution with |=

        Think of "|" as the union of components.

        This accepts a single argument.
        """
        self.addContribution(args)
        return self

    def __iadd__(self, args):
        """addVar or newVar with +=

        Think of "+" as addition of a variable.

        This accepts a single argument or an iterable of single arguments or
        argument tuples.
        """
        # Want to detect addVar or newVar
        def f(*args):
            if isinstance(args[0], six.string_types):
                self.newVar(*args)
            else:
                self.addVar(*args)
            return

        _applymanyargs(args, f)
        return self

# End class FitRecipeInterface

class ParameterFactory(object):
    """This class is used to create Parameters on demand."""

    def __init__(self, parclass):
        """Initialize the factory.

        parclass    --  The class used to create new parameters (default
                        diffpy.srfit.fitbase.parameter.Parameter)
        """
        self._pclass = parclass

    def __getattr__(self, name):
        pcls = object.__getattribute__(self, '_pclass')
        par = pcls(name=name)
        setattr(self, name, par)
        return par

# End class ParameterFactory

def _applymanyargs(args, f):
    """Apply arguments to a function.

    Args can be any of the following.
    arg
    (arg1, arg2, ...)
    ((arg1a, arg1b, ...), ...)
    """
    if not hasattr(args, '__iter__'):
        f(args)
        return

    for arg in args:
        if hasattr(arg, '__iter__'):
            f(*arg)
        else:
            f(arg)

    return

def _applyargs(args, f):
    """Apply arguments to a function.

    Args can be any of the following.
    arg
    (arg1, arg2, ...)
    ((arg1a, arg1b, ...), ...)
    """
    if not hasattr(args, '__iter__'):
        f(args)
    else:
        f(*args)
    return

# End of file
