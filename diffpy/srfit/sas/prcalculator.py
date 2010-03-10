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
    _profile    --  Profile containing I(q).
    _usecalc    --  Use the calculated signal in the profile to calculate P(r)
                    (default False).

    """

    def __init__(self, name, invertor):
        """Initialize the generator.

        name        --  A name for the PrCalculator
        invertor    --  Configured invertor. This may or may not be configured
                        with data. If it is not, use the setProfile method.
        
        """
        Calculator.__init__(self, name)

        self._invertor = invertor
        self._profile = None
        self._usecalc = False
        return

    def setProfile(self, profile, usecalc = False):
        """Set the profile that is used to calculate P(r).


        profile     --  The Profile instance to refer to when calculating P(r).
                        This should contain the I(q) signal. By default, the x,
                        y and dy arrays are used to calculate I(q).
        usecalc     --  Use the calculated signal in the profile to calculate
                        P(r) (default False). This makes repetitive calculation
                        more time consuming, but allows for the simultaneous
                        characterization of I(q) and G(r).
        
        """

        if self._profile is not None:
            if self._usecalc:
                self._profile.ycpar.removeObserver(self._flush)
            else:
                self._profile.removeObserver(self._flush)

        self._usecalc = usecalc
        self._profile = profile

        if usecalc:
            profile.ycpar.addObserver(self._flush)
        else:
            profile.addObserver(self._flush)
        self._flush(self)

        return

    def __call__(self, r):
        """Calculate P(r) from the data or calculated signal."""
        print "Calculate P(r)"
        self._invertor.d_max = max(r) + 5.0
        # Assume profile doesn't include 0. It's up to the user to make this
        # happen.
        self._invertor.x = self._profile.x
        if self._usecalc:
            if self._profile.ycalc is not None:
                self._invertor.y = self._profile.ycalc
            else:
                self._invertor.y = self._profile.y
        else:
            self._invertor.y = self._profile.y
        self._invertor.err = self._profile.dy
        c, c_cov = self._invertor.invert_optimize()
        l = lambda x: self._invertor.pr(c, x)
        pr = map(l, r)

        pr = numpy.array(pr)
        return pr

# End class PrCalculator
