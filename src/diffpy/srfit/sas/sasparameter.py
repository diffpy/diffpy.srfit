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
"""SAS profile generator.

The SASGenerator class wraps a sas.models.BaseModel object as a
ProfileGenerator.
"""

__all__ = ["SASParameter"]

from diffpy.srfit.fitbase.parameter import Parameter


class SASParameter(Parameter):
    """Class adapting a sasmodel parameter to srfit Parameter.

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
    _model  --  The BaseModel to which the underlying parameter belongs.
    _parname    --  The name of the underlying BaseModel parameter.
    """

    def __init__(self, name, model, parname=None):
        """Create the Parameter.

        name    --  Name of the Parameter
        model   --  The BaseModel to which the underlying parameter belongs
        parname --  Name of parameter used by the model. If this is None
                    (default), then name is used.
        """
        self._parname = parname or name
        val = model.getParam(self._parname)
        self._model = model
        Parameter.__init__(self, name, val)
        return

    def getValue(self):
        """Get the value of the Parameter."""
        value = self._model.getParam(self._parname)
        return value

    def setValue(self, value):
        """Set the value of the Parameter."""
        if value != self.getValue():
            self._model.setParam(self._parname, value)
            self.notify()
        return self


# End of class SASParameter
