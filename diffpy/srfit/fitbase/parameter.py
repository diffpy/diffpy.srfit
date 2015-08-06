#!/usr/bin/env python
########################################################################
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
########################################################################
"""Parameter classes.

Parameters encapsulate an adjustable parameter within SrFit.

"""
from __future__ import print_function
import six
import sys
# IDEA - Add onConstrain, onRestrain, onVary so that adaptors to Parameters
#        can have more fine control over the construction of FitRecipes.
# IDEA - Add tags to parameters so they can be easily retrieved.
# IDEA - Consider scaling parameters to avoid precision issues in optimizers.

__all__ = [ "Parameter", "ParameterProxy", "ParameterAdapter"]

import numpy as np

from diffpy.srfit.equation.literals import Argument
from diffpy.srfit.equation.literals.abcs import ArgumentABC
from diffpy.srfit.util.nameutils import validateName
from diffpy.srfit.util.argbinders import bind2nd
from diffpy.srfit.interface import _parameter_interface
from diffpy.srfit.fitbase.validatable import Validatable


class Parameter(_parameter_interface, Argument, Validatable):
    """Parameter class.

    Attributes
    name    --  A name for this Parameter.
    const   --  A flag indicating whether this is considered a constant.
    _value  --  The value of the Parameter. Modified with 'setValue'.
    value   --  Property for 'getValue' and 'setValue'.
    constrained --  A flag indicating if the Parameter is constrained
                (default False).
    bounds  --  A 2-list defining the bounds on the Parameter. This can be
                used by some optimizers when the Parameter is varied. See
                FitRecipe.getBounds and FitRecipe.boundsToRestraints.

    """

    def __init__(self, name, value = None, const = False):
        """Initialization.

        name    --  The name of this Parameter (must be a valid attribute
                    identifier)
        value   --  The initial value of this Parameter (default 0).
        const   --  A flag inticating whether the Parameter is a constant (like
                    pi).

        Raises ValueError if the name is not a valid attribute identifier

        """
        self.constrained = False
        self.bounds = [-np.inf, np.inf]
        validateName(name)
        Argument.__init__(self, name, value, const)
        return

    def setValue(self, val, lb = None, ub = None):
        """Set the value of the Parameter and the bounds.

        val     --  The value to assign.
        lb      --  The lower bounds for the bounds list. If this is None
                    (default), then the lower bound will not be alterered.
        ub      --  The upper bounds for the bounds list. If this is None
                    (default), then the upper bound will not be alterered.

        Returns self so that mutators can be chained.

        """
        Argument.setValue(self, val)
        if lb is not None: self.bounds[0] = lb
        if ub is not None: self.bounds[1] = ub
        return self

    def setConst(self, const = True, value = None):
        """Toggle the Parameter as constant.

        const   --  Flag indicating if the parameter is constant (default
                    True).
        value   --  An optional value for the parameter (default None). If this
                    is not None, then the parameter will get a new value,
                    constant or otherwise.

        Returns self so that mutators can be chained.

        """
        self.const = bool(const)
        if value is not None:
            self.setValue(value)
        return self

    def boundWindow(self, lr = 0, ur = None):
        """Create bounds centered on the current value of the Parameter.

        lr  --  The radius of the lower bound (default 0). The lower bound is
                computed as value - lr.
        ur  --  The radius of the upper bound. The upper bound is computed as
                value + ur. If this is None (default), then the value of the
                lower radius is used.

        Returns self so that mutators can be chained.

        """
        val = self.getValue()
        lb = val - lr
        if ur is None:
            ur = lr
        ub = val + ur
        self.bounds = [lb, ub]
        return self

    def _validate(self):
        """Validate my state.

        This validates that value is not None.

        Raises AttributeError if validation fails.

        """
        if self.value is None:
            raise AttributeError("value of '%s' is None"%self.name)
        return

# End class Parameter

class ParameterProxy(_parameter_interface, Validatable):
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
        par = object.__getattribute__(self, 'par')
        return getattr(par, attrname)

    value = property( lambda self: self.par.getValue(),
            lambda self, val: self.par.setValue(val) )

    def __str__(self):
        return "%s(%s)"%(self.__class__.__name__, self.name)

    def _validate(self):
        """Validate my state.

        This validates that value and par are not None.

        Raises AttributeError if validation fails.

        """
        if self.value is None:
            raise AttributeError("value is None")
        if self.par is None:
            raise AttributeError("par is None")
        return

# End class ParameterProxy

# I'm not sure that this is needed in python 3...
if sys.version_info[0] == 2:
    # Make sure that this is registered as an Argument class
    ArgumentABC.register(ParameterProxy)

class ParameterAdapter(Parameter):
    """An adapter for parameter-like objects.

    This class wraps an object as a Parameter. The getValue and setValue
    methods defer to the data of the wrapped object.

    """

    def __init__(self, name, obj, getter = None, setter = None, attr = None):
        """Wrap an object as a Parameter.

        name    --  The name of this Parameter.
        obj     --  The object to be wrapped.
        getter  --  The unbound function that can be used to access the
                    attribute containing the parameter value. getter(obj) should
                    return the Parameter value.  If getter is None (default),
                    it is assumed that an attribute is accessed via attr. If
                    attr is also specified, then the Parameter value will be
                    accessed via getter(obj, attr).
        setter  --  The unbound function that can be used to modify the
                    attribute containing the parameter value. setter(obj, value)
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
                self.getter = bind2nd(getattr, self.attr)
            else:
                self.getter = bind2nd(getter, self.attr)

            if setter is None:
                self.setter = bind2nd(setattr, self.attr)
            else:
                self.setter = bind2nd(setter, self.attr)

        value = self.getValue()
        Parameter.__init__(self, name, value)
        return

    def getValue(self):
        """Get the value of the Parameter."""
        return self.getter(self.obj)

    def setValue(self, value, lb = None, ub = None):
        """Set the value of the Parameter."""
        if value != self.getValue():
            self.setter(self.obj, value)
            self.notify()

        if lb is not None: self.bounds[0] = lb
        if ub is not None: self.bounds[1] = ub

        return self

# End class ParameterAdapter

# End of file
