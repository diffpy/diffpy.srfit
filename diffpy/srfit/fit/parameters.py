#!/usr/bin/env python

from diffpy.srfit.adapters.nodes import Parameter
from diffpy.srfit.util import makeName

__all__ = ["Parameter", "Par", "Pars", "Var", "Vars", "Fixed"]

# Shorthand for Parameter
Par = Parameter

# FIXME - check types!
def Pars(*args): 
    """Factory for multiple parameters.

    args should consist of (name, value) pairs. Each pair is passed to Par and
    the list of outputs is returned in a list.

    """
    plist = []
    for name, value in args:
        par = Par(name, value)
        plist.append(par)
    return plist

def Var(name, value = None):
    """Factory for variable parameters.

    name    --  Name of the parameter
    value   --  Value of the parameter

    The returned parameter is varied by default, but can be fixed.

    """
    p = Par(name, value)
    p.vary()
    return p

def Vars(*args):
    """Factory for multiple variable parameters.

    args should consist of (name, value) pairs. Each pair is passed to Var and
    the list of outputs is returned in a list.

    """
    plist = []
    for name, value in args:
        par = Var(name, value)
        plist.append(par)
    return plist

class FixedParameter(Parameter):
    """Parameter container for a reference value.

    This parameter cannot be varied.

    """

    def vary(self, val = None, dovary = True):
        """Raises AttributeError: cannot vary this parameter"""
        raise AttributeError("Parameter cannot be varied")

    def fix(self, val = None):
        """Set this as fixed during a fit. This is the default state.

        val --  New value for the parameter. If this is None , the
                parameter value will not change.

        Returns self so that mutator methods can be chained.
        
        """
        Parameter.vary(self, val, False)
        return self

    def varied(self):
        """Indicate if this Parameter is varied. Always False."""
        return False

# End class FixedParameter

# Shorthand for FixedParameter
Fixed = FixedParameter
