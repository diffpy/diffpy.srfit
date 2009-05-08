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

    Right now this adds no significant functionality to Argument. It is
    subclassed for future extensibility.
    """

    def __init__(self, name, value = 0, const = False):
        """Initialization.
        
        name    --  The name of this Parameter.
        value   --  The initial value of this Parameter (default 0).
        const   --  A flag inticating whether the Parameter is a constant (like
                    pi).
        """
        Argument.__init__(self, value, name, const)
        return

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

# End class ParameterProxy

class ParameterWrapper(Parameter):
    """An adapter for parameter-like objects.

    This class wraps an object as a Paramter. The getValue and setValue methods
    of Parameter directly modify the appropriate attribute of the paramter-like
    object.
    """

    def __init__(self, obj, name, getter = None, setter = None, attr = None):
        """Wrap an object as a Parameter.

        obj     --  The object to be wrapped.
        name    --  The name of the object.
        getter  --  The unbound function that can be used to access the
                    attribute containing the paramter value. getter(obj) should
                    return the Paramter value.  If getter is None (default), it
                    is assumed that an attribute is accessed directly.
        setter  --  The unbound function that can be used to modify the
                    attribute containing the paramter value. setter(obj, value)
                    should set the attribute to the passed value. If getter is
                    None (default), it is assumed that an attribute is accessed
                    directly.
        attr    --  The name of the attribute that contains the value of the
                    parameter. If attr is None (default), then both getter and
                    setter must be specified.

        raises ValueError if exactly one of getter or setter is not None, or if
        getter, setter and attr ar all None.
        """
        if getter is None and setter is None and attr is None:
            raise ValueError("Specify attribute access")
        if [getter, setter].count(None) == 1:
            raise ValueError("Specify both getter and setter")


        self.obj = obj
        self.getter = getter
        self.setter = setter
        self.attr = attr

        if attr is None:
            self.getValue = self.getAccessValue
            self.setValue = self.setAccessValue
        else:
            self.getValue = self.getAttrValue
            self.setValue = self.setAttrValue

        value = self.getValue()
        Parameter.__init__(self, name, value)
        return

    def getAttrValue(self):
        """Get the value of the attribute."""
        return getattr(self.obj, self.attr)

    def setAttrValue(self, value):
        """Set the value of the attribute."""
        if value != self.getValue():
            setattr(self.obj, self.attr, value)
            self.clicker.click()
        return

    def getAccessValue(self):
        """Get the value from the accessor."""
        return self.getter(self.obj)

    def setAccessValue(self, value):
        """Set the value from the accessor."""
        if value != self.getValue():
            self.setter(self.obj, value)
            self.clicker.click()
        return
                    
# End class ParameterWrapper

# version
__id__ = "$Id$"

#
# End of file
