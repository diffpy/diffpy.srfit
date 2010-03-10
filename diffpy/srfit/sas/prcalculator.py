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
"""Nanoparticle form factor P(r) calculator.

The PrCalculator class wraps a sans.pr.invertor.Invertor object as a
Calculator. This is not wrapped as a ProfileGenerator because it will be used
to share information between SAS I(Q) to PDF G(r), but it does not use the same
profile as the PDF, which is where the calculator will be applied.

"""

__all__ = ["PrCalculator"]

import numpy

from diffpy.srfit.fitbase import Calculator

class PrCalculator(Calculator):
    """A class for calculating P(r) from data.

    This is provided so that SAS data can be used to calculate a nanoparticle
    form factor (characteristic function) for nanoparticle PDF refinements. For
    a crystal-like nanoparticle:
    Gnano(r) = f(r)Gcrystal(r),
    where f(r) is the nanoparticle form factor. 
    This is obtained from P(r) as
    P(r) = 4 pi r**2 f(r).

    Attributes:
    _invertor   --  sans.pr.invertor.Invertor object. This object is configured
                    by the user and added during initialization. This allows
                    the user to specify a background, regularization, etc.
                    using the Invertor interface.

    Managed Parameters:
    scale       --  The scale factor (default 1).
    """

    def __init__(self, name, invertor):
        """Initialize the generator.

        name        --  A name for the PrCalculator
        invertor    --  Configured invertor. This may or may not be configured
                        with data. If it is not, use the setProfile method.
        
        """
        Calculator.__init__(self, name)

        self._invertor = invertor

        self._newParameter("scale", 1)
        return

    def __call__(self, r, q, iq, diq = None):
        """Calculate P(r) from the data or calculated signal."""
        self._invertor.d_max = max(r) + 5.0
        # Assume profile doesn't include 0. It's up to the user to make this
        # happen.
        self._invertor.x = q
        self._invertor.y = iq
        if diq is None:
            diq = numpy.ones_like(q)
        self._invertor.err = diq
        c, c_cov = self._invertor.invert_optimize()
        l = lambda x: self._invertor.pr(c, x)
        pr = map(l, r)

        pr = numpy.array(pr)
        return self.scale.value * pr

# End class PrCalculator
