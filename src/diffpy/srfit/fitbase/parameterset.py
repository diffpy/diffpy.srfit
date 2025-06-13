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
"""ParameterSet class.

ParameterSets organize Parameters, Constraints, Restraints and other
ParameterSets. They provide attribute-access of other ParameterSets and
embedded Parameters.
"""

__all__ = ["ParameterSet"]


from collections import OrderedDict

from diffpy.srfit.fitbase.recipeorganizer import RecipeOrganizer


class ParameterSet(RecipeOrganizer):
    """Class for organizing Parameters and other ParameterSets.

    ParameterSets are hierarchical organizations of Parameters, Constraints,
    Restraints and other ParameterSets.

    Contained Parameters and other ParameterSets can be accessed by name as
    attributes in order to facilitate multi-level constraints and restraints.
    These constraints and restraints can be placed at any level and a flattened
    list of them can be retrieved with the getConstraints and getRestraints
    methods.

    Attributes
    name            --  A name for this organizer.
    _calculators    --  A managed dictionary of Calculators, indexed by name.
    _constraints    --  A set of constrained Parameters. Constraints can be
                        added using the 'constrain' methods.
    _parameters     --  A managed OrderedDict of parameters.
    _restraints     --  A set of Restraints. Restraints can be added using the
                        'restrain' methods.
    _parsets        --  A managed dictionary of ParameterSets.
    _eqfactory      --  A diffpy.srfit.equation.builder.EquationFactory
                        instance that is used to create constraints and
                        restraints from string equations.

    Properties
    names           --  Variable names (read only). See getNames.
    values          --  Variable values (read only). See getValues.
    """

    def __init__(self, name):
        """Initialize.

        name    --  The name of this ParameterSet.
        """
        RecipeOrganizer.__init__(self, name)

        self._parsets = OrderedDict()
        self._manage(self._parsets)
        return

    # Alias Parameter accessors.
    addParameter = RecipeOrganizer._addParameter
    newParameter = RecipeOrganizer._newParameter
    removeParameter = RecipeOrganizer._removeParameter

    def addParameterSet(self, parset):
        """Add a ParameterSet to the hierarchy.

        parset  --  The ParameterSet to be stored.

        Raises ValueError if the ParameterSet has no name.
        Raises ValueError if the ParameterSet has the same name as some other
        managed object.
        """
        self._addObject(parset, self._parsets, True)
        return

    def removeParameterSet(self, parset):
        """Remove a ParameterSet from the hierarchy.

        Raises ValueError if parset is not managed by this object.
        """
        self._removeObject(parset, self._parsets)
        return

    def setConst(self, const=True):
        """Set every parameter within the set to a constant.

        const   --  Flag indicating if the parameter is constant (default
                    True).
        """
        for par in self.iterPars():
            par.setConst(const)

        return


# End class ParameterSet

# End of file
