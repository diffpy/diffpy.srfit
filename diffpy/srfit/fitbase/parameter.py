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
"""Parameter class. 

Parameters encapsulate an adjustable parameter of a Calculator.
"""

from diffpy.srfit.equation.literals import Argument
from diffpy.srfit.equation import Clicker
from diffpy.srfit.equation.builder import EquationFactory

from .modelorganizer import ModelOrganizer

class Parameter(Argument):
    """Parameter class.
    
    Attributes
    name    --  A name for this Parameter. Names should be unique within a
                ParameterSet.
    clicker --  A Clicker instance for recording change in the value.
    value   --  The value of the Parameter. Modified with setValue.
    """

    def __init__(self, name, value = 0):
        """Initialization.
        
        value   --  The initial value of this Parameter (default 0).
        name    --  The name of this Parameter (optional, default None).
        """
        Argument.__init__(self, value, name)
        return

    def getValue(self):
        """Get the value of this Parameter.

        This is provided as an alternative to reading the value from the
        'value' attribute. Collaborators should interact through this method
        whenever possible in order to allow for easy extension of the method
        through inheritance.
        """
        return self.value

    def __str__(self):
        if self.name:
            return "Parameter(" + self.name + ")"
        return self.__repr__()

# End class Parameter

class ParameterSet(ModelOrganizer):
    """A collection of Parameters and ParameterSets.

    This class organizes Parameters and other ParameterSets. Contained objects
    can be accessed within a ParameterSet as if they were attributes according
    to their name.

    ParameterSets inherit from dict, where entries are indexed by their name.
    It is recommended to add objects using the 'insert' method, which takes
    care of assigning a key and associating the clickers of the objects.

    Constraints and Restraints can be applied within a ParameterSet.

    Attributes
    name    --  A name for this ParameterSet. Names should be unique within a
                ParameterSet.
    clicker --  A Clicker instance for recording changes in contained objects.
    constraints     --  A dictionary of Constraints, indexed by the constrained
                        Parameter. 
    pardict --  A dictionary containing the Parameters and ParameterSets held
                herein.
    restraints      --  A set of Restraints.
    suborganizers   --  A list of ParameterSets that this ParameterSets knows
                        about (for quick access).
    _eqfactory      --  A diffpy.srfit.equation.builder.EquationFactory
                        instance that is used to create constraints and
                        restraints from strings.

    """

    def __init__(self, name):
        """Initialize.

        name    --  The name of this ParameterSet.
        """
        ModelOrganizer.__init__(self)
        self.name = name
        self.clicker = Clicker()
        self.constraints = {}
        self.restraints = set()
        self._eqfactory = EquationFactory()
        self.pardict = {}
        return

    def __getattr__(self, name):
        """Gives access to the contained objects as attributes."""
        arg = self.pardict.get(name)
        if arg is None:
            raise AttributeError(name)
        return arg

    def insert(self, par, check=True):
        """Store a Parameter or ParameterSet.

        par     --  The Parameter or ParameterSet to be stored.
        check   --  If True (default), an ValueError is raised if the parameter
                    has an invalid name, or if an object of that name has
                    already been inserted.
        """
        if check:
            message = ""
            if not par.name:
                message = "Object '%s' has no name"%par
            elif par.name in self.pardict:
                message = "Object with name '%s' already exists"%par.name

            if message:
                raise ValueError(message)

        self.pardict[par.name] = par
        self.clicker.addSubject(par.clicker)

        if isinstance(par, ParameterSet):
            self.suborganizers.add(par)
        else:
            self._eqfactory.registerArgument(par.name, par)

        return

# End class ParameterSet

# version
__id__ = "$Id$"

#
# End of file
