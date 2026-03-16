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
"""Parameter classes.

Parameters encapsulate an adjustable parameter within SrFit.
"""

# IDEA - Add onConstrain, onRestrain, onVary so that adaptors to Parameters
# can have more fine control over the construction of FitRecipes.
# IDEA - Add tags to parameters so they can be easily retrieved.
# IDEA - Consider scaling parameters to avoid precision issues in optimizers.

__all__ = ["Parameter", "ParameterProxy", "ParameterAdapter"]

from functools import wraps

import numpy

from diffpy.srfit.equation.literals import Argument
from diffpy.srfit.exceptions import SrFitError
from diffpy.srfit.fitbase.validatable import Validatable
from diffpy.srfit.interface import _parameter_interface
from diffpy.srfit.util.argbinders import bind2nd
from diffpy.srfit.util.nameutils import validateName
from diffpy.utils._deprecator import build_deprecation_message, deprecated

parameter_base = "diffpy.srfit.fitbase.Parameter"
removal_version = "4.0.0"
setValue_dep_msg = build_deprecation_message(
    parameter_base, "setValue", "set_value", removal_version
)

setConst_dep_msg = build_deprecation_message(
    parameter_base, "setConst", "set_constant", removal_version
)

boundRange_dep_msg = build_deprecation_message(
    parameter_base, "boundRange", "bound_range", removal_version
)


class Parameter(_parameter_interface, Argument, Validatable):
    """Parameter class.

    Attributes
    ----------
    name
        A name for this Parameter.
    const
        A flag indicating whether this is considered a constant.
    _value
        The value of the Parameter. Modified with 'set_value'.
    value
        Property for 'getValue' and 'set_value'.
    constrained
        A flag indicating if the Parameter is constrained
        (default False).
    bounds
        A 2-list defining the bounds on the Parameter. This can be
        used by some optimizers when the Parameter is varied. See
        FitRecipe.get_bounds_pairs and FitRecipe.convert_bounds_to_restraints.
    """

    def __init__(self, name, value=None, const=False):
        """Initialization.

        Parameters
        ----------
        name
            The name of this Parameter (must be a valid attribute
            identifier)
        value
            The initial value of this Parameter (default 0).
        const
            A flag inticating whether the Parameter is a constant (like
            pi).


        Raises ValueError if the name is not a valid attribute identifier
        """
        self.constrained = False
        self.bounds = [-numpy.inf, +numpy.inf]
        validateName(name)
        Argument.__init__(self, name, value, const)
        return

    def set_value(self, val):
        """Set the value of the Parameter and the bounds.

        Parameters
        ----------
        val
            The value to assign.
        lower_bound : float
            The lower bounds for the bounds list. If this is None
            (default), then the lower bound will not be alterered.
        upper_bound : float
            The upper bounds for the bounds list. If this is None
            (default), then the upper bound will not be alterered.

        Returns
        -------
        self
            Returns self so that mutators can be chained.
        """
        Argument.set_value(self, val)
        return self

    @deprecated(setValue_dep_msg)
    def setValue(self, val):
        """This function has been deprecated and will be removed in
        version 4.0.0.

        Please use diffpy.srfit.fitbase.Parameter.set_value instead.
        """
        return self.set_value(val)

    def set_constant(self, is_constant=True, value=None):
        """Toggle the Parameter as constant.

        Parameters
        ----------
        is_constant : bool, optional
            The flag indicating if the parameter is constant (default
            True).
        value : float, optional
            The value value for the parameter to be set to (default None).
            If this is not None, then the parameter will get a new value,
            constant or otherwise.

        Returns
        -------
        self
            Returns self so that mutators can be chained.
        """
        self.const = bool(is_constant)
        if value is not None:
            self.set_value(value)
        return self

    @deprecated(setConst_dep_msg)
    def setConst(self, const=True, value=None):
        """This function has been deprecated and will be removed in
        version 4.0.0.

        Please use diffpy.srfit.fitbase.Parameter.set_constant instead.
        """
        self.set_constant(const, value)
        return self

    def bound_range(self, lower_bound=None, upper_bound=None):
        """Set lower and upper bound of the Parameter.

        Parameters
        ----------
        lower_bound : float
            The lower bound for the bounds list.
        upper_bound : float
            The upper bound for the bounds list.

        Returns
        -------
        self
            Returns self so that mutators can be chained.
        """
        if lower_bound is not None:
            self.bounds[0] = lower_bound
        if upper_bound is not None:
            self.bounds[1] = upper_bound
        return self

    @deprecated(boundRange_dep_msg)
    def boundRange(self, lower_bound=None, upper_bound=None):
        """This function has been deprecated and will be removed in
        version 4.0.0.

        Please use diffpy.srfit.fitbase.Parameter.bound_range instead.
        """
        self.bound_range(lower_bound, upper_bound)
        return self

    def boundWindow(self, lr=0, ur=None):
        """Create bounds centered on the current value of the Parameter.

        Parameters
        ----------
        lr
            The radius of the lower bound (default 0). The lower bound is
            computed as value - lr.
        ur
            The radius of the upper bound. The upper bound is computed as
            value + ur. If this is None (default), then the value of the
            lower radius is used.

        Returns
        -------
        self
            Returns self so that mutators can be chained.
        """
        val = self.getValue()
        lower_bound = val - lr
        if ur is None:
            ur = lr
        upper_bound = val + ur
        self.bounds = [lower_bound, upper_bound]
        return self

    def _validate(self):
        """Validate my state.

        This validates that value is not None.

        Raises SrFitError if validation fails.
        """
        if self.value is None:
            raise SrFitError("value of '%s' is None" % self.name)
        return


# End class Parameter


class ParameterProxy(Parameter):
    """A Parameter proxy for another parameter.

    This allows for the same parameter to have multiple names.

    Attributes
    ----------
    name
        A name for this ParameterProxy. Names should be unique within a
        RecipeOrganizer and should be valid attribute names.
    par
        The Parameter this is a proxy for.
    """

    def __init__(self, name, par):
        """Initialization.

        Parameters
        ----------
        name
            The name of this ParameterProxy.
        par
            The Parameter this is a proxy for.


        Raises ValueError if the name is not a valid attribute identifier
        """
        validateName(name)

        self.name = name
        self.par = par
        return

    # define properties to use attributes of the proxied Parameter -----------

    @property
    def constrained(self):
        """A flag indicating if the proxied Parameter is constrained."""
        return self.par.constrained

    @constrained.setter
    def constrained(self, value):
        self.par.constrained = bool(value)
        return

    @property
    def bounds(self):
        """List of lower and upper bounds of the proxied Parameter.

        This can be used by some optimizers when the Parameter is
        varied. See FitRecipe.get_bounds_pairs and
        FitRecipe.convert_bounds_to_restraints.
        """
        return self.par.bounds

    @bounds.setter
    def bounds(self, value):
        self.par.bounds = value
        return

    @property
    def _observers(self):
        return self.par._observers

    # wrap Parameter methods to use the target object ------------------------

    @wraps(Parameter.set_value)
    def set_value(self, val):
        return self.par.set_value(val)

    @wraps(Parameter.getValue)
    def getValue(self):
        return self.par.getValue()

    @wraps(Parameter.set_constant)
    def set_constant(self, const=True, value=None):
        return self.par.set_constant(const, value)

    @wraps(Parameter.bound_range)
    def bound_range(self, lower_bound=None, upper_bound=None):
        return self.par.bound_range(lower_bound, upper_bound)

    @wraps(Parameter.boundWindow)
    def boundWindow(self, lr=0, ur=None):
        return self.par.boundWindow(lr, ur)

    def _validate(self):
        """Validate my state.

        This validates that value and par are not None.

        Raises SrFitError if validation fails.
        """
        if self.par is None:
            raise SrFitError("par is None")
        self.par._validate()
        return


# End class ParameterProxy


class ParameterAdapter(Parameter):
    """An adapter for parameter-like objects.

    This class wraps an object as a Parameter. The getValue and
    set_value methods defer to the data of the wrapped object.
    """

    def __init__(self, name, obj, getter=None, setter=None, attr=None):
        """Wrap an object as a Parameter.

        Parameters
        ----------
        name
            The name of this Parameter.
        obj
            The object to be wrapped.
        getter
            The unbound function that can be used to access the
            attribute containing the parameter value. getter(obj)
            should return the Parameter value.  If getter is None
            (default), it is assumed that an attribute is accessed
            via attr. If attr is also specified, then the Parameter
            value will be accessed via getter(obj, attr).
        setter
            The unbound function that can be used to modify the
            attribute containing the parameter value.
            setter(obj, value) should set the attribute to the
            passed value. If setter is None (default), it is assumed
            that an attribute is accessed via attr. If attr is also
            specified, then the Parameter value will be set via
            setter(obj, attr, value).
        attr
            The name of the attribute that contains the value of the
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

    def set_value(self, value):
        """Set the value of the Parameter."""
        if value != self.getValue():
            self.setter(self.obj, value)
            self.notify()
        return self


# End class ParameterAdapter

# End of file
