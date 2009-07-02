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
"""Adapter classes for translating a srfit FitModel to a PARK Fitness.

The FitnessAdapter class adapts a FitModel object (from diffpy.srfit.fitbase)
to a Fitness object (from the 'pak' branch of PARK). This will allow for
FitModels to be refined using PARK optimizers and services.

The adaptation exposes residual from the
FitModel as the residuals function in the Fitness. Variables in the FitModel
are exposes as a parameters in the Fitness using the ParkParameterProxy
class.

"""
from diffpy.srfit.fitbase.parameter import Parameter
import park

class ParkParameterProxy(park.Parameter):
    """A park.Parameter subclass that serves as a proxy to a SrFit Parameter.

    This class is designed to expose fitted variables from a FitModel to the
    PARK interface. Changes made to either of the srfit Parameter or its proxy
    will reflect in the other.

    """

    def __init__(self, par):
        """Initialize with a srfit parameter."""
        park.Parameter.__init__(self, name=par.name)
        self._par = par
        self.status = "fitted"
        return

    def _getValue(self):
        """Retrieves value of srfit parameter."""
        return self._par.getValue()

    def _setValue(self, val):
        """Sets value of srfit parameter."""
        self._par.setValue(val)
        return

    value = property(_getValue, _setValue)

# End class ParkParameterProxy

class FitnessAdapter(park.Fitness):
    """A park.Fitness subclass that serves as an adapter to a FitModel.

    """

    def __init__(self, fitmodel):
        """Initialize with a FitModel instance."""
        park.Fitness.__init__(self)
        self._fitmodel = fitmodel
        self._parmap = {}
        # Make a parameter set and add the variables from the FitModel
        self.__parameterset = park.ParameterSet()
        for par in fitmodel._parameters:
            parkpar = ParkParameterProxy(par)
            self.__parameterset.append(parkpar)
            self._parmap[par] = parkpar

        return

    def _parameterset(self):
        return self.__parameterset
    parameterset = property(_parameterset)

    def derivs(self):
        """Parameters with analytic derivatives.
        
        Does nothing as of now.
        """
        return
        # The FitModel residual requires input paramters

    def residuals(self):
        """Calculate the residuals."""

        return self._fitmodel.residual()

    def residuals_deriv(self, pars=[]):
        """No derivatives here."""
        return []

    def set(self, **kw):
        """Set parameters in the model."""
        for k,v in kw.items():
            self.parameterset[k].set(v)
        return

# End class FitnessAdapter


# version
__id__ = "$Id$"

#
# End of file
