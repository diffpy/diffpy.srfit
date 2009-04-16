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


class ParameterSet(dict):
    """A collection of Parameters and ParameterSets.

    This class organizes Parameters and other ParameterSets. Contained objects
    can be accessed within a ParameterSet as if they were attributes according
    to their name.

    ParameterSets inherit from dict, where entries are indexed by their name.
    It is recommended to add objects using the 'insert' method, which takes
    care of assigning a key and associating the clickers of the objects.

    There is no attempt to preserve the order of objects within a ParameterSet.
    When an order is required, assume that contained objects are ordered
    alphabetically.

    Attributes
    name    --  A name for this ParameterSet. Names should be unique within a
                ParameterSet.
    clicker --  A Clicker instance for recording changes in contained objects.

    """

    def __init__(self, name):
        """Initialize.

        name    --  The name of this ParameterSet.
        """
        dict.__init__(self)
        self.name = name
        self.clicker = Clicker()
        return

    def __getattr__(self, name):
        """Gives access to the contained objects as attributes."""
        arg = self.get(name)
        if arg is None:
            raise AttributeError(name)
        return arg

    def insert(self, par, check=True):
        """Store a Parameter or ParameterSet.

        par     --  The Parameter or ParameterSet to be stored.
        check   --  If True (default), an ValueError is raised if the name
                    of par computes as False, or if an object of that name has
                    already been inserted.
        """
        if check:
            message = ""
            if not par.name:
                message = "Object '%s' has no name"%par
            elif par.name in self:
                message = "Object with name '%s' already exists"%par.name

            if message:
                raise ValueError(message)

        self[par.name] = par
        self.clicker.addSubject(par.clicker)
        return

# version
__id__ = "$Id$"

#
# End of file
