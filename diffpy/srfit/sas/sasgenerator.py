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
"""SAS profile generator.

The SASGenerator class wraps a sans.models.BaseModel object as a
ProfileGenerator.

"""

import numpy

from diffpy.srfit.fitbase import ProfileGenerator
from diffpy.srfit.fitbase.parameter import Parameter, ParameterAdapter
from diffpy.srfit.fitbase.parameterset import ParameterSet

class SASParameter(Parameter):
    """Class adapting a sansmodel parameter to srfit Parameter.
    
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

    class _SASBounds(object):

        def __init__(self, name, model):
            self.name = name
            self._model = model
            # We use -inf..inf for unbounded
            if self._model.details[name][1] is None:
                self._model.details[name][1] = -numpy.inf
            if self._model.details[name][2] is None:
                self._model.details[name][2] = numpy.inf
            return

        def __getitem__(self, i):
            return self._model.details[self.name][i+1]

        def __setitem__(self, i, val):
            self._model.details[self.name][i+1] = float(val)
            return

        def __str__(self):
            return self._model.details[self.name][1:].__str__()

    # End class _SASBounds

    def __init__(self, name, model):
        """Create the Parameter.

        name    --  Name of the Parameter
        model   --  The BaseModel to which the underlying parameter belongs

        """
        val = model.getParam(name)
        self._model = model
        Parameter.__init__(self, name, val)
        self.bounds = self.__class__._SASBounds(name, model)
        return

    def getValue(self):
        """Get the value of the Parameter."""
        return self._model.getParam(self.name)

    def setValue(self, value):
        """Set the value of the Parameter."""
        if value != self.getValue():
            self._model.setParam(self.name, value)
            self.notify()
        return

class SASGenerator(ProfileGenerator):
    """A class for calculating I(Q) from a scattering type.

    Attributes:
    _model      --  BaseModel object this adapts.

    Managed Parameters:
    These depend on the parameters of the BaseModel object held by _model. They
    are created from the 'params' attribute of the BaseModel. If a dispersion
    is set for the BaseModel, the dispersion "width" will be accessible under
    "<parname>width", where <parname> is the name a parameter adjusted by
    dispersion.

    """

    def __init__(self, name, model):
        """Initialize the generator.

        name    --  A name for the SASGenerator
        model   --  SASModel object this adapts.
        
        """
        ProfileGenerator.__init__(self, name)

        self._model = model

        # Wrap normal parameters
        for pname in model.params:
            par = SASParameter(pname, model)
            self.addParameter(par)

        # Wrap dispersion parameters
        def _dispgetter(obj, _name):
            return obj.getParam(_name)
        def _dispsetter(obj, _name, val):
            return obj.setParam(_name, val)

        for pname in model.dispersion:
            attrname = "%s.width"%pname
            pname = "%swidth"%pname

            par = ParameterAdapter(pname, self._model, _dispgetter,
                    _dispsetter, attrname)
            self.addParameter(par)

        return

    def __call__(self, q):
        """Calculate I(Q) for the BaseModel."""
        return self._model.evalDistribution(q)

# End class SASGenerator
