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

from diffpy.utils._deprecator import build_deprecation_message, deprecated

from diffpy.srfit.fitbase.recipeorganizer import RecipeOrganizer

base = "diffpy.srfit.fitbase.parameterset.ParameterSet"
removal_version = "4.0.0"

addparset_dep_msg = build_deprecation_message(
    base, "addParameterSet", "add_parameter_set", removal_version
)

removeParameterSet_dep_msg = build_deprecation_message(
    base, "removeParameterSet", "remove_parameter_set", removal_version
)

setConst_dep_msg = build_deprecation_message(
    base, "setConst", "set_constant", removal_version
)


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
    ----------
    name
        A name for this organizer.
    _calculators
        A managed dictionary of Calculators, indexed by name.
    _constraints
        A set of constrained Parameters. Constraints can be
        added using the 'constrain' methods.
    _parameters
        A managed OrderedDict of parameters.
    _restraints
        A set of Restraints. Restraints can be added using the
        'restrain' methods.
    _parsets
        A managed dictionary of ParameterSets.
    _eqfactory
        A diffpy.srfit.equation.builder.EquationFactory
        instance that is used to create constraints and
        restraints from string equations.

    Properties
    ----------
    names
        Variable names (read only). See get_names.
    values
        Variable values (read only). See get_values.
    """

    def __init__(self, name):
        """Initialize.

        Parameters
        ----------
        name
            The name of this ParameterSet.
        """
        RecipeOrganizer.__init__(self, name)

        self._parsets = OrderedDict()
        self._manage(self._parsets)
        return

    # Alias Parameter accessors.
    addParameter = RecipeOrganizer._add_parameter
    newParameter = RecipeOrganizer._new_parameter
    removeParameter = RecipeOrganizer._remove_parameter

    def add_parameter_set(self, parset):
        """Add a ParameterSet to the hierarchy.

        Parameters
        ----------
        parset
            The ParameterSet to be stored.


        Raises ValueError if the ParameterSet has no name.
        Raises ValueError if the ParameterSet has the same name as some other
        managed object.
        """
        self._add_object(parset, self._parsets, True)
        return

    @deprecated(addparset_dep_msg)
    def addParameterSet(self, parset):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use
        diffpy.srfit.fitbase.parameterset.ParameterSet.add_parameter_set
        instead.
        """
        self.add_parameter_set(parset)
        return

    def remove_parameter_set(self, parset):
        """Remove a ParameterSet from the hierarchy.

        Raises ValueError if parset is not managed by this object.
        """
        self._remove_object(parset, self._parsets)
        return

    @deprecated(removeParameterSet_dep_msg)
    def removeParameterSet(self, parset):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use
        diffpy.srfit.fitbase.parameterset.ParameterSet.remove_parameter_set
        instead.
        """
        self.remove_parameter_set(parset)
        return

    def set_constant(self, is_constant=True):
        """Set every parameter within the set to a constant.

        Parameters
        ----------
        is_constant : bool, optional
            The flag indicating if the parameter is constant (default
            True).
        """
        for par in self.iterate_over_parameters():
            par.set_constant(is_constant)
        return

    @deprecated(setConst_dep_msg)
    def setConst(self, const=True):
        """This function has been deprecated and will be removed in
        version 4.0.0.

        Please use
        diffpy.srfit.fitbase.parameterset.ParameterSet.set_constant
        instead.
        """
        self.set_constant(const)
        return


# End class ParameterSet

# End of file
