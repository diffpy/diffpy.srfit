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

"""Validatable class.

A Validatable has state that must be validated before a FitRecipe can first
calculate the residual.
"""

__all__ = ["Validatable"]


class Validatable(object):
    """Validatable class.

    A Validatable has state that must be validated by a FitRecipe.

    """

    def _validateOthers(self, iterable):
        """Method to validate configuration of Validatables in iterable.

        This is provided as a convenience for derived classes.  No need to
        overload this. Call this method from overloaded _validate with an
        iterable of other Validatables.

        """
        for obj in iterable:
            if obj is self: continue
            if isinstance(obj, Validatable):
                obj._validate()

        return

    def _validate(self):
        """Validate self and Validatables.

        Overload this in a derived class.

        Raises AttributeError if validation fails.

        """
        # Validate self here.
        # Then validate others.
        return

# End class Validatable

# End of file
