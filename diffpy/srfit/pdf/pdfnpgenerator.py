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
"""Prototype PDF generator for nanoparticles.

"""

__all__ = ["PDFNPGenerator"]

import numpy

from periodictable import elements

from diffpy.srfit.fitbase import ProfileGenerator, Profile
from diffpy.srfit.fitbase import FitContribution, FitRecipe
from diffpy.srfit.fitbase import FitResults
from diffpy.srfit.structure.objcryststructure import ObjCrystParSet

class PDFNPGenerator(ProfileGenerator):
    """A class for calculating the PDF for an isolated scatterer.""" 

    def __init__(self, name):
        """Initialize the generator.
        
        """
        ProfileGenerator.__init__(self, name)

        # Add any non-structural parameters here
        self.newParameter("delta2", 0)
        self.newParameter("qbroad", 0)

        # Non-Parameters that are needed by the generator.
        self.meta["qmin"] = 1
        self.meta["qmax"] = 20
        self._phase = None
        return

    def setCrystal(self, cryst):
        """Set the pyobjcryst.Crystal instance for the generator.

        This converts the Crystal to a ObjCrystParSet that is used to organize
        the fitting parameters for the generator.  The generator will have its
        own parameters, each of which will be a proxy for some part of the
        crystal. The parameters will be accessible by name under the 'phase'
        attribute of this generator.
        
        """

        # Create a custom ParameterSet designed to interface with
        # pyobjcryst.Crystal
        parset = ObjCrystParSet(cryst, "phase")
        # Put this ParameterSet in the ProfileGenerator.
        self.addParameterSet(parset)
        self._phase = parset
        return

    def __call__(self, r):
        """Calculate the PDF.

        This ProfileGenerator will be used in a fit equation that will be
        optimized to fit some data.  By the time this function is evaluated,
        the crystal has been updated by the optimizer via the ObjCrystParSet
        created in setCrystal. Thus, we need only call pdf with the internal
        structure object.

        """
        qmin = self.meta["qmin"]
        qmax = self.meta["qmax"]
        return pdf(self.phase.stru, r, self.delta2.getValue(),
                self.qbroad.getValue(), qmin, qmax)

# End class PDFNPGenerator

def pdf(cryst, r, delta2 = 0, qbroad = 0, qmin = 1, qmax = 20):
    """Calculate the PDF from a diffpy.Structure

    cryst   --  A pyobjcryst.Crystal instance. It is assumed that the structure
                is that of an isolated scatterer. Periodic boundary conditions
                are not applied.
    r       --  The r-points to calculate over.
    delta2  --  The correlation term in the Debye-Waller factor (default 0).
    qbroad  --  Resolution-based broadening term.
    qmin    --  The minimum observed q-value (default 1).
    qmax    --  The maximum observed q-value (default 20).
    
    """
    pi = numpy.pi

    # First we must calculate F(Q). To do this, we need to create the q-points
    # corresponding to the r-points. It is assumed that the r points are
    # equidistant.
    rmax = r[-1]
    dr = r[1] - r[0]

    qmaxcalc = 2*pi/dr
    # Note that this should be 2*pi/rmax, but the fft calculates a signal for
    # both positive and negative values of frequency, and we're going to throw
    # away half of the transformed values.
    #dq = pi/rmax
    rmaxcalc = rmax + 10
    dq = pi / rmaxcalc

    # We must create the q-points that correspond to the r-points
    q = numpy.arange(0, qmaxcalc, dq)

    # Calculate F(Q)
    fq = fofq(cryst, q, delta2, qbroad)

    # and enforce the qmin and qmax bounds
    fq[numpy.logical_or(q < qmin, q > qmax)] = 0

    # Now we transform this to gr
    from numpy.fft import ifft, fftfreq

    # sin FT = imaginary part of ifft
    grcalc = ifft(fq).imag
    rspan = 2 * pi / dq
    rcalc = fftfreq(len(fq)) * rspan

    # Properly normalize this. 
    grcalc *= (2 / pi) * qmaxcalc

    # Now we need to pick out the positive frequency terms. 
    #rmaxidx = int(numpy.ceil(r[-1]/dr))
    # We also need to discard values for r < rmin.
    #rminidx = int(numpy.floor(r[0]/dr))

    #gr = grcalc[rminidx:rmaxidx]
    gr = grcalc[ numpy.logical_and(rcalc >= r[0], rcalc <= r[-1]) ]

    assert(len(r) == len(gr))

    return gr


def fofq(cryst, q, delta2, qbroad):
    """Calculate F(Q) (X-ray) using the Debye Equation.

    F(Q) = 2/(N<f>^2) 
           sum(i, j > i) f_i(Q) f_j(Q) sin(rij Q)/rij exp(-0.5 ssij Q**2)
    (The exponential term is the Debye-Waller factor.)

    cryst   --  A pyobjcryst.Crystal instance. It is assumed that the structure
                is that of an isolated scatterer. Periodic boundary conditions
                are not applied.
    q       --  The q-points to calculate over.
    delta2  --  The correlation term in the Debye-Waller factor.
    qbroad  --  Resolution-based broadening term.

    The generator uses cctbx for the calculation of the f_i if it is
    available, otherwise f_i = 1.

    """
    # The functions we need
    sin = numpy.sin
    exp = numpy.exp
    pi = numpy.pi

    # The brute-force calculation is very slow. Thus we optimize a little bit.

    # The precision of distance measurements.
    deltad = 1e-8
    dmult = int(1/deltad)
    deltau = deltad**2
    umult = int(1/deltau)

    pairdict = {}
    elcount = {}

    # Loop over scattering components in the crystal.

    scl = cryst.GetScatteringComponentList()
    n = len(scl)


    for i in xrange(n):

        # Get the ScatteringPower
        scli = scl[i]
        spi = scli.mpScattPow

        # Check for dummy element
        if spi is None: continue

        # Count the element
        eli = spi.GetSymbol().title()
        m = elcount.get(eli, 0)
        elcount[eli] = m + 1

        xif = scli.X
        yif = scli.Y
        zif = scli.Z
        xi, yi, zi = cryst.FractionalToOrthonormalCoords(xif, yif, zif)

        for j in xrange(i+1,n):

            sclj = scl[j]
            # Get the ScatteringPower
            spj = sclj.mpScattPow

            # Check for dummy element
            if spj is None: continue

            elj = spj.GetSymbol().title()

            # Get the pair
            els = [eli, elj]
            els.sort()

            # Get the coordinates
            xjf = sclj.X
            yjf = sclj.Y
            zjf = sclj.Z
            xj, yj, zj = cryst.FractionalToOrthonormalCoords(xjf, yjf, zjf)

            # Get the distance to the desired precision
            d = ((xi - xj)**2 + (yi - yj)**2 + (zi - zj)**2)**0.5
            D = int(d*dmult)

            # Get the DW factor to the same precision (calculated from Biso,
            # not Uiso)
            ss = 0
            for sp in [spi, spj]:
                if spi.IsIsotropic():
                    ss += numpy.fabs(sp.Biso)
                else:
                    B = numpy.array([
                        [sp.B11, sp.B12, sp.B13],
                        [sp.B12, sp.B22, sp.B23],
                        [sp.B13, sp.B23, sp.B33]], dtype=float)
                    r = numpy.array([xi-xj, yi-yj, zi-zj], dtype=float)
                    br = numpy.dot(r, numpy.dot(B, r))
                    br /= numpy.dot(r, r)
                    ss += numpy.fabs(br)
            corfact = 1 - delta2/d**2 + (qbroad*d)**2
            if corfact > 0:
                ss *= corfact
            SS = int(ss*umult)

            # Record the multiplicity of this pair
            key = (els[0], els[1], D, SS)
            mult = pairdict.get(key, 0)
            pairdict[key] = mult + 1


    # Now we can calculate G(r) from the pair dictionary. Making the dictionary
    # first reduces the amount of calls to sin and exp we have to make.

    # First we must cache the scattering factors
    fdict = {}
    nfavg = 0
    for el, m in elcount.iteritems():
        fdict[el] = f = getXScatteringFactor(el, q)
        nfavg += m * f

    # Now we can compute F(Q) for the i != j pairs
    y = 0
    BtoU = 1.0 / (8 * pi**2)
    for key, mult in pairdict.items():
        eli = key[0]
        elj = key[1]
        fi = fdict[eli]
        fj = fdict[elj]
        D = key[2]
        SS = key[3]

        d = D * deltad

        # Add in the fitcontribution
        y += fi * fj * mult * sin(q * d) / d *\
                exp(-0.5 * (SS * deltau) * BtoU * q**2)


    # We must multiply by 2 since we only counted j > i pairs.
    y *= 2

    # Now divide by N <f>^2

    y *= n / nfavg**2

    return y

def getXScatteringFactor(el, q):
    """Get the x-ray scattering factor for an element over the q range.
    
    This uses q-independent x-ray scattering factor.

    """
    elobj = getattr(elements, el)
    f = elobj.number
    fq = numpy.ones_like(q) * f
    return fq

# End of file
