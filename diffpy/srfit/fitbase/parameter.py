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
        
        name    --  The name of this Parameter.
        value   --  The initial value of this Parameter (default 0).
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

class ParameterProxy(object):
    """A Parameter proxy for another parameter. This allows for the same
    parameter to have multiple names.
    """
    
    # this will let the proxy register as a real Parameter
    __class__ = Parameter

    def __init__(self, name, par):
        """Initialization.

        name    --  The name of this ParameterProxy.
        par     --  The Parameter this is a proxy for.
        """
        self.name = name
        self.par = par
        return

    def __getattr__(self, name):
        """Redirect accessors and attributes to the reference Parameter."""
        return getattr(self.par, name)

        

# version
__id__ = "$Id$"

#
# End of file
