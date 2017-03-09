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

"""Configurable class.

A Configurable has state of which a FitRecipe must be aware.
"""

__all__ = ["Configurable"]


class Configurable(object):
    """Configurable class.

    A Configurable has state of which a FitRecipe must be aware.

    Attributes
    _configobjs     --  Set of Configureables in a hierarcy or instances.
                        Messasges get passed up the hierarcy to a FitReciple
                        via these objects.

    """

    def __init__(self):
        self._configobjs = set()
        return

    def _updateConfiguration(self):
        """Notify Configurables in hierarchy of configuration change."""
        for obj in self._configobjs:
            obj._updateConfiguration()
        return

    def _storeConfigurable(self, obj):
        """Store a Configurable.

        The passed obj is only stored if it is a a Configurable, otherwise this
        method quietly exits.

        """
        if isinstance(obj, Configurable):
            self._configobjs.add(obj)
        return

# End class Configurable

# End of file
