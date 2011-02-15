#!/usr/bin/env python
########################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2011 Trustees of the Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
########################################################################
from diffpy.Structure import Structure
from diffpy.srreal.pdfcalculator import PDFCalculator
import diffpy.srfit.pdf
from diffpy.srfit import *

def main(ciffile, datname):
    """Create a fitting recipe for crystalline PDF data."""

    # Load the data. The PDFParser simplifies this task. In short, there is no
    # reason to have both a PDFParser and a Profile in this example. In the
    # future, Parsers will inherit from Profile so we can save a few lines of
    # code here.
    profile = loadProfile(datname, "PDF")
    profile.setRange(xmax = 20)
    r, gr, dgr = profile
    _r = r.get()

    # Create a PDFCalculator. We can use information from the metadata loaded
    # from file to initialize it. We can't pass the whole meta dictionary since
    # it contains entries that PDFCalculator can't use. In addition, passing a
    # non-zero qmin value to PDFCalculator cretes a "bad" PDF.
    calc = PDFCalculator(qmax = profile.meta["qmax"])
    calc.setScatteringFactorTableByType(profile.meta["stype"])
    calc.rmin = _r[0]
    calc.rstep = _r[1] - _r[0]
    calc.rmax = _r[-1] + 0.5 * calc.rstep

    # Now we adapt the PDFCalculator. There is no special code written to do
    # this. The default ContainerAdapter is being used to wrap the calculator.
    # I will probably write factor functions that streamline the creation of
    # the calculator and adapter.
    g = adapt(calc)

    # Create the structure and adapt it. Again, this could be simplified with a
    # factory function.
    stru = Structure(filename = ciffile)
    stru.Uisoequiv = 0.005
    s = adapt(stru, "nickel")
    
    # Here retrieve the parameters that we want to refine. We write the
    # space group constraints manually. 
    a = s.lattice.a
    a.vary()
    dummy = Par("dummy")
    #s.lattice.b.constrain(a)
    #s.lattice.c.constrain(a)
    # FIXME - need an addConstraint method, or something similar.
    # lattice.addConstraint(s.lattice.setLatPar(a = a, b = a, c = a))

    # Retrive the independent Uisoequiv parameter. Note that mutators, such as
    # 'vary' return 'self' so that they can be chained.
    Uisoequiv = s[0].Uisoequiv.vary()
    for atom in s[1:]:
        atom.Uisoequiv.constrain(Uisoequiv)

    # Now vary parameters from the calculator.
    g.scale.vary(1)
    g.qdamp.vary(0.01)
    g.delta2.vary(5)

    def talker(val):
        print "************ called"
        return val
    atalker = adapt(talker, "talker")

    # Create the fit equation.
    out = gcalc = g(s)
    out = atalker(gcalc)
    print s._cache
    print s._cache._neighbors, gcalc._cache
    print gcalc._cache, g._cache
    print out._cache

    rcalc, gcalc = out
    rcalc.rename("rcalc")
    gcalc.rename("gcalc")

    fiteq = interp(r, rcalc, gcalc)
    # Create the residual equation. Note that 'chi' creates a vector residual
    # that can be dotted into itself to generate 'chi^2'.
    reseq = chi(fiteq, gr, dgr)
    # Create the residual object. This step and the previous can be performed
    # at once with the 'reschi2' function.
    res = residual(reseq)

    print res.names, res.values

    # Optimize. 
    from scipy.optimize import leastsq, fmin
    # XXX Why aren't the atoms valid after the first couple of calls?
    leastsq(res.vec, res.values)
    #fmin(res, res.values)

    # Get the results
    results = FitResults(res, showfixed=True)
    results.show()
    return

    # and plot 
    from pylab import plot, show
    # Note that 'value' is a property for 'get' and 'set'.
    plot(r.get(), gr.get(), 'bo')
    plot(r.get(), fiteq.get(), 'r-')
    show()

if __name__ == "__main__":

    # Make the data and the recipe
    ciffile = "data/ni.cif"
    data = "data/ni-q27r100-neutron.gr"
    main(ciffile, data)

# End of file
