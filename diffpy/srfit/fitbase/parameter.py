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

Parameters encapsulate an adjustable parameter of a Calculator. There are three
varieties of Parameter here.
Parameter           --  Base parameter class.
ParameterWrapper    --  Wrap attributes or methods of an object as a Parameter.
ParameterProxy      --  A proxy for another Parameter, but with a different
                        name.

"""
# IDEA - Add onConstrain, onRestrain, onVary so that adaptors to Parameters
# can have more fine control over the construction of FitRecipes.
# IDEA - Add tags to parameters so they can be easily retrieved.
# IDEA - Consider scaling parameters to avoid precision issues in optimizers.

from diffpy.srfit.equation.literals import Argument
from diffpy.srfit.equation.literals.abcs import ArgumentABC
from diffpy.srfit.util.clicker import Clicker
from diffpy.srfit.util.nameutils import validateName

class Parameter(Argument):
    """Parameter class.
    
    Attributes
    name    --  A name for this Parameter. Names should be unique within a
                RecipeOrganizer and should be valid attribute names.
    clicker --  A Clicker instance for recording change in the value.
    value   --  The value of the Parameter. Modified with setValue.
    const   --  A flag indicating whether the Parameter is constant. A constant
                parameter cannot be made into a variable, nor can it be
                constrained to something else.

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
        validateName(name)

        Argument.__init__(self, value, name, const)

        self.constrained = False
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

    def __str__(self):
        if self.name:
            return "Parameter(" + self.name + ")"
        return self.__repr__()

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

# End class ParameterProxy

# Make sure that this is registered as an Argument class
ArgumentABC.register(ParameterProxy)

class ParameterWrapper(Parameter):
    """An adapter for parameter-like objects.

    This class wraps an object as a Paramter. The getValue and setValue methods
    of Parameter directly modify the appropriate attribute of the paramter-like
    object.

    """

    def __init__(self, name, obj, getter = None, setter = None, attr = None):
        """Wrap an object as a Parameter.

        name    --  The name of this Parameter.
        obj     --  The object to be wrapped.
        getter  --  The unbound function that can be used to access the
                    attribute containing the paramter value. getter(obj) should
                    return the Parameter value.  If getter is None (default),
                    it is assumed that an attribute is accessed directly. If
                    attr is also specified, then the Parameter value will be
                    accessed via getter(obj, attr).
        setter  --  The unbound function that can be used to modify the
                    attribute containing the paramter value. setter(obj, value)
                    should set the attribute to the passed value. If setter is
                    None (default), it is assumed that an attribute is accessed
                    directly. If attr is also specified, then the Parameter
                    value will be set via setter(obj, attr, value).
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
            self.clicker.click()

        return
                    
# End class ParameterWrapper

# version
__id__ = "$Id$"

#
# End of file
