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
"""ParameterSet class. 

ParameterSets organize Parameters, Constraints, Restraints and other
ParameterSets.
"""

from .modelorganizer import ModelOrganizer

class ParameterSet(ModelOrganizer):
    """Class for organizing Parameters and other ParameterSets.

    ParameterSets are hierarchical organizations of Parameters, Constraints,
    Restraints and other ParameterSets. 

    Contained parameters and other ParameterSets can be accessed by name as
    attributes in order to facilitate multi-level constraints and restraints.
    These constraints and restraints can be placed at any level and a flattened
    list of them can be retrieved with the getConstraints and getRestraints
    methods.

    Attributes
    clicker         --  A Clicker instance for recording changes in contained
                        Parameters and ModelOrganizers.
    name            --  A name for this organizer.
    _constraints    --  A dictionary of Constraints, indexed by the constrained
                        Parameter. Constraints can be added using the
                        'constrain' method.
    _orgdict        --  A dictionary containing the Parameters and
                        ModelOrganizers indexed by name.
    _parameters     --  A list of parameters that this ModelOrganizer knows
                        about.
    _restraints     --  A set of Restraints. Restraints can be added using the
                        'restrain' or 'confine' methods.
    _organizers     --  A list of ParameterSets that this ParameterSet knows
                        about.
    _eqfactory      --  A diffpy.srfit.equation.builder.EquationFactory
                        instance that is used to create constraints and
                        restraints from string

    """

    def __init__(self, name):
        """Initialize.

        name    --  The name of this ParameterSet.
        """
        ModelOrganizer.__init__(self, name)
        return

    # Alias the _addParameter to addParameter
    addParameter = ModelOrganizer._addParameter
    # Alias the addOrganizer method to addParameterSet
    addParameterSet = ModelOrganizer._addOrganizer

    def setConst(self, const = True):
        """Set every parameter within the set to a constant.

        const   --  Flag indicating if the parameter is constant (default
                    True).

        """
        for par in self._iterPars():
            par.setConst(const)

        return

# End class ParameterSet

# version
__id__ = "$Id$"

#
# End of file
