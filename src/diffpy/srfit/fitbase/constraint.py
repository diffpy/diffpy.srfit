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
"""Constraint class.

Constraints are used by a FitRecipe (and other RecipeOrganizers) to
organize constraint equations. They store a Parameter object and an
Equation object that is used to compute its value. The
Constraint.constrain method is used to create this association.
"""

__all__ = ["Constraint"]

from diffpy.srfit.exceptions import SrFitError
from diffpy.srfit.fitbase.validatable import Validatable
from diffpy.utils._deprecator import build_deprecation_message, deprecated

base = "diffpy.srfit.fitbase.constraint.Constraint"
removal_version = "4.0.0"

constrain_deprecation_msg = build_deprecation_message(
    base,
    "constrain",
    "add_constraint",
    removal_version,
)

unconstrain_deprecation_msg = build_deprecation_message(
    base,
    "unconstrain",
    "remove_constraint",
    removal_version,
)


class Constraint(Validatable):
    """Constraint class.

    Constraints are designed to be stored in only one place. (The holder of the
    constraint owns it).

    Attributes
    ----------
    par
        A Parameter that is the subject of the constraint.
    eq
        An equation whose evaluation is used to set the value of the
        constraint.
    """

    def __init__(self):
        """Initialization."""
        self.par = None
        self.eq = None
        return

    def add_constraint(self, par, eq):
        """Constrain a Parameter according to an Equation.

        The parameter will be set constant once it is constrained. This
        will keep it from being constrained multiple times.

        Parameters
        ----------
        par : Parameter
            The Parameter to constrain.
        eq : Equation
            The Equation to use to constrain the Parameter.

        Raises
        ------
        ValueError
            If par is constant or already constrained.
        """

        if par.const:
            raise ValueError("The parameter '%s' is constant" % par)

        if par.constrained:
            raise ValueError("The parameter '%s' is already constrained" % par)

        par.constrained = True

        self.par = par
        self.eq = eq
        self.update()
        return

    @deprecated(constrain_deprecation_msg)
    def constrain(self, par, eq):
        """This function has been deprecated and will be removed in
        version 4.0.0.

        Please use
        diffpy.srfit.fitbase.constraint.Constraint.add_constraint
        instead.
        """
        self.add_constraint(par, eq)
        return

    def remove_constraint(self):
        """Clear the constraint from a Parameter."""
        self.par.constrained = False
        self.par = None
        self.eq = None
        return

    @deprecated(unconstrain_deprecation_msg)
    def unconstrain(self):
        """This function has been deprecated and will be removed in
        version 4.0.0.

        Please use
        diffpy.srfit.fitbase.constraint.Constraint.remove_constraint
        instead.
        """
        self.remove_constraint()
        return

    def update(self):
        """Update the parameter according to the equation."""
        # This will be evaluated quickly thanks to the Equation class.
        val = self.eq()
        # This will only change the Parameter if val is different from the
        # currently stored value.
        self.par.set_value(val)
        return

    def _validate(self):
        """Validate my state.

        This validates that par is not None. This validates eq.

        Raises SrFitError if validation fails.
        """
        if self.par is None:
            raise SrFitError("par is None")
        if self.eq is None:
            raise SrFitError("eq is None")
        self.par._validate()
        from diffpy.srfit.equation.visitors import validate

        try:
            validate(self.eq)
        except ValueError as e:
            raise SrFitError(e)

        # Try to get the value of eq.
        try:
            val = self.eq()
            self.par.set_value(val)
        except TypeError:
            raise SrFitError("eq cannot be evaluated")
        finally:
            if val is None:
                raise SrFitError("eq evaluates to None")

        return


# End of file
