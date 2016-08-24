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

"""Nanoparticle form factor P(r) calculator.

The PrCalculator class wraps a sas.pr.invertor.Invertor object as a
Calculator. This is not wrapped as a ProfileGenerator because it will be used
to share information between SAS I(Q) to PDF G(r), but it does not use the same
profile as the PDF, which is where the calculator will be applied.
"""

__all__ = ["PrCalculator", "CFCalculator"]

import numpy

from diffpy.srfit.fitbase import Calculator

Invertor = None

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
    _invertor   --  sas.pr.invertor.Invertor object. This object is internal,
                    but can be configured by the user after initialization.
                    Note that the 'x', 'y' and 'err' attributes get overwritten
                    every time the invertor is used.

    Managed Parameters:
    scale       --  The scale factor (default 1).
    q           --  The q-values of the I(q) signal
    iq          --  The I(q) signal
    diq         --  The uncertainty in I(q)

    """

    def __init__(self, name):
        """Initialize the generator.

        name        --  A name for the PrCalculator

        """
        Calculator.__init__(self, name)

        # delayed import of Invertor
        global Invertor
        if Invertor is None:
            from diffpy.srfit.sas.sasimport import sasimport
            Invertor = sasimport('sas.pr.invertor').Invertor

        self._invertor = Invertor()

        self._newParameter("scale", 1)
        self._newParameter("q", None)
        self._newParameter("iq", None)
        self._newParameter("diq", None)
        return

    def __call__(self, r):
        """Calculate P(r) from the data or calculated signal."""
        q = self.q.value
        iq = self.iq.value
        diq = self.diq.value
        if diq is None:
            diq = numpy.ones_like(q)
        self._invertor.d_max = max(r) + 5.0
        # Assume profile doesn't include 0. It's up to the user to make this
        # happen.
        self._invertor.x = q
        self._invertor.y = iq
        self._invertor.err = diq
        c, c_cov = self._invertor.invert_optimize()
        l = lambda x: self._invertor.pr(c, x)
        pr = map(l, r)

        pr = numpy.array(pr)
        return self.scale.value * pr

# End class PrCalculator

class CFCalculator(PrCalculator):
    """A class for calculating the characteristic function (CF) from data.

    This calculator produces
    f(r) = P(r) / 4 pi r**2
    which is the nanoparticle form factor scaled by density.

    Attributes:
    _invertor   --  sas.pr.invertor.Invertor object. This object is internal,
                    but can be configured by the user after initialization.
                    Note that the 'x', 'y' and 'err' attributes get overwritten
                    every time the invertor is used.

    Managed Parameters:
    scale       --  The scale factor (default 1).
    q           --  The q-values of the I(q) signal
    iq          --  The I(q) signal
    diq         --  The uncertainty in I(q)

    """

    def __call__(self, r):
        """Calculate P(r) from the data or calculated signal."""
        fr = PrCalculator.__call__(self, r)
        fr /= 4 * numpy.pi * r**2
        if r[0] == 0:
            # Assume the scale makes fr properly normalized. We don't have much
            # other choice.
            fr[0] = 1
        return fr

# End class CFCalculator
