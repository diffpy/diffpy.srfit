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

__all__ = ["SASGenerator"]

from diffpy.srfit.fitbase import ProfileGenerator
from diffpy.srfit.sas.sasparameter import SASParameter


class SASGenerator(ProfileGenerator):
    """A class for calculating I(Q) from a scattering type.

    Attributes:
    _model      --  BaseModel object this adapts.

    Managed Parameters:
    These depend on the parameters of the BaseModel object held by _model. They
    are created from the 'params' attribute of the BaseModel. If a dispersion
    is set for the BaseModel, the dispersion "width" will be accessible under
    "<parname>_width", where <parname> is the name a parameter adjusted by
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
        for parname in model.params:
            par = SASParameter(parname, model)
            self.addParameter(par)

        # Wrap dispersion parameters
        for parname in model.dispersion:
            name = parname + "_width"
            parname += ".width"
            par = SASParameter(name, model, parname)
            self.addParameter(par)

        return

    def __call__(self, q):
        """Calculate I(Q) for the BaseModel."""
        return self._model.evalDistribution(q)


# End class SASGenerator
