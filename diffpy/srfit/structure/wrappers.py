#!/usr/bin/env python
"""This module contains general wrappers for interfacing structure containers
as a hierarchy of ParameterSets and Parameters.
"""

from diffpy.srfit.fitbase.parameter import Parameter

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

        Parameter.__init__(self, name)
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
                    

__id__ = "$Id$"
