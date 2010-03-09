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
"""Parameter classes.

Parameters encapsulate an adjustable parameter within SrFit.

"""
# IDEA - Add onConstrain, onRestrain, onVary so that adaptors to Parameters
# can have more fine control over the construction of FitRecipes.
# IDEA - Add tags to parameters so they can be easily retrieved.
# IDEA - Consider scaling parameters to avoid precision issues in optimizers.

__all__ = [ "Parameter", "ParameterProxy", "ParameterAdapter"]

from numpy import inf

from diffpy.srfit.equation.literals import Argument
from diffpy.srfit.equation.literals.abcs import ArgumentABC
from diffpy.srfit.util.nameutils import validateName


class Parameter(Argument):
    """Parameter class.
    
    Attributes
    name    --  A name for this Parameter.
    const   --  A flag indicating whether this is considered a constant.
    _value  --  The value of the Parameter. Modified with 'setValue'.
    value   --  Property for 'getValue' and 'setValue'.
    constrained --  A flag indicating if the Parameter is constrained
                (default False).
    bounds  --  A 2-list defining the bounds on the Parameter. This can be
                used by some optimizers when the Parameter is varied. These
                bounds are unrelated to restraints on a Parameter.

    """

    def __init__(self, name, value = None, const = False):
        """Initialization.
        
        name    --  The name of this Parameter (must be a valid attribute
                    identifier)
        value   --  The initial value of this Parameter (default 0).
        const   --  A flag inticating whether the Parameter is a constant (like
                    pi).
        constrained --  A flag indicating if the Parameter is constrained
                    (default False).

        Raises ValueError if the name is not a valid attribute identifier
        
        """
        self.constrained = False
        self.bounds = [-inf, inf]
        validateName(name)
        Argument.__init__(self, name, value, const)
        return

    def setConst(self, const = True, value = None):
        """Toggle the Parameter as constant.

        const   --  Flag indicating if the parameter is constant (default
                    True).
        value   --  An optional value for the parameter (default None). If this
                    is not None, then the parameter will get a new value,
                    constant or otherwise.

        """
        self.const = bool(const)
        if value is not None:
            self.setValue(value)
        return

# End class Parameter

class ParameterProxy(object):
    """A Parameter proxy for another parameter. 
    
    This allows for the same parameter to have multiple names.

    Attributes
    name    --  A name for this ParameterProxy. Names should be unique within a
                RecipeOrganizer and should be valid attribute names.
    par     --  The Parameter this is a proxy for.

    """


    def __init__(self, name, par):
        """Initialization.

        name    --  The name of this ParameterProxy.
        par     --  The Parameter this is a proxy for.

        Raises ValueError if the name is not a valid attribute identifier

        """
        validateName(name)

        self.name = name
        self.par = par
        return

    def __getattr__(self, attrname):
        """Redirect accessors and attributes to the reference Parameter."""
        return getattr(self.par, attrname)

    value = property( lambda self: self.par.getValue(),
            lambda self, val: self.par.setValue(val) )

# End class ParameterProxy

# Make sure that this is registered as an Argument class
ArgumentABC.register(ParameterProxy)

class ParameterAdapter(Parameter):
    """An adapter for parameter-like objects.

    This class wraps an object as a Paramter. The getValue and setValue methods
    of Parameter directly modify the appropriate attribute of the
    parameter-like object.

    """

    def __init__(self, name, obj, getter = None, setter = None, attr = None):
        """Wrap an object as a Parameter.

        name    --  The name of this Parameter.
        obj     --  The object to be wrapped.
        getter  --  The unbound function that can be used to access the
                    attribute containing the paramter value. getter(obj) should
                    return the Parameter value.  If getter is None (default),
                    it is assumed that an attribute is accessed via attr. If
                    attr is also specified, then the Parameter value will be
                    accessed via getter(obj, attr).
        setter  --  The unbound function that can be used to modify the
                    attribute containing the paramter value. setter(obj, value)
                    should set the attribute to the passed value. If setter is
                    None (default), it is assumed that an attribute is accessed
                    via attr. If attr is also specified, then the Parameter
                    value will be set via setter(obj, attr, value).
        attr    --  The name of the attribute that contains the value of the
                    parameter. If attr is None (default), then both getter and
                    setter must be specified.

        Raises ValueError if exactly one of getter or setter is not None, or if
        getter, setter and attr are all None.

        """
        if getter is None and setter is None and attr is None:
            raise ValueError("Specify attribute access")
        if [getter, setter].count(None) == 1:
            raise ValueError("Specify both getter and setter")

        self.obj = obj
        self.getter = getter
        self.setter = setter
        self.attr = attr

        if attr is not None:
            if getter is None:
                self.getter = lambda obj: getattr(obj, self.attr)
            else:
                self.getter = lambda obj: getter(obj, self.attr)

            if setter is None:
                self.setter = lambda obj, val: setattr(obj, self.attr, val)
            else:
                self.setter = lambda obj, val: setter(obj, self.attr, val)

        value = self.getValue()
        Parameter.__init__(self, name, value)
        return

    def getValue(self):
        """Get the value of the Parameter."""
        return self.getter(self.obj)

    def setValue(self, value):
        """Set the value of the Parameter."""
        if value != self.getValue():
            self.setter(self.obj, value)
            self.notify()

        return
                    
# End class ParameterAdapter

# version
__id__ = "$Id$"

#
# End of file
