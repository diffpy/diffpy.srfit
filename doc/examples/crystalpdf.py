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
from diffpy.srfit.pdf import PDFParser
from diffpy.srfit import *

def main(ciffile, datname):
    """Create a fitting recipe for crystalline PDF data."""

    # Load the data. The PDFParser simplifies this task. In short, there is no
    # reason to have both a PDFParser and a Profile in this example. In the
    # future, Parsers will inherit from Profile so we can save a few lines of
    # code here.
    parser = PDFParser()
    parser.parseFile(datname)
    profile = Profile()
    profile.load(parser)
    profile.setRange(xmax = 20)
    r, gr, dgr = profile
    r = r.get()

    # Create a PDFCalculator. We can use information from the metadata loaded
    # from file to initialize it. We can't pass the whole meta dictionary since
    # it contains entries that PDFCalculator can't use. In addition, passing a
    # non-zero qmin value to PDFCalculator cretes a "bad" PDF.
    calc = PDFCalculator(qmax = profile.meta["qmax"])
    calc.setScatteringFactorTableByType(profile.meta["stype"])
    calc.rmin = r[0]
    calc.rstep = r[1] - r[0]
    calc.rmax = r[-1] + 0.5 * calc.rstep
    r = adapt(r, "r")

    # Now we adapt the PDFCalculator. There is no special code written to do
    # this. The default ContainerAdapter is being used to wrap the calculator.
    # I will probably write factor functions that streamline the creation of
    # the calculator and adapter.
    g = adapt(calc)

    # Create the structure and adapt it. Again, this could be simplified with a
    # factory function.
    stru = Structure(filename = ciffile)
    for a in stru:
        a.Uisoequiv = 0.005
    s = adapt(stru, "nickel")
    
    # Here retrieve the parameters that we want to refine. We write the
    # space group constraints manually. 
    a = s.lattice.a
    a.vary()
    s.lattice.b.constrain(a)
    s.lattice.c.constrain(a)

    # Retrive the independent Uisoequiv parameter. Note that mutators, such as
    # 'vary' return 'self' so that they can be chained.
    Uisoequiv = s[0].Uisoequiv.vary(0.005)
    for atom in s[1:]:
        atom.Uisoequiv.constrain(Uisoequiv)

    # Now vary parameters from the calculator.
    g.scale.vary(1)
    g.qdamp.vary(0.01)
    g.delta2.vary(5)

    # Create the fit equation.
    out = g(s)
    # FIXME - calculation gets screwed up somehow. Changes in the parameters
    # are not seen at any level deeper than interp.
    assert(out in s._viewers)
    rcalc, gcalc = out
    assert(rcalc in out._viewers)
    assert(out in rcalc._viewers)
    assert(gcalc in out._viewers)
    assert(out in gcalc._viewers)
    rcalc.rename("rcalc")
    gcalc.rename("gcalc")
    fiteq = interp(r, rcalc, gcalc)
    assert(fiteq in r._viewers)
    assert(fiteq in rcalc._viewers)
    assert(fiteq in gcalc._viewers)
    # Create the residual equation. Note that 'chi' creates a vector residual
    # that can be dotted into itself to generate 'chi^2'.
    reseq = chi(fiteq, gr, dgr)
    # Create the residual object. This step and the previous can be performed
    # at once with the 'reschi2' function.
    res = residual(reseq)

    print res.names, res.values

    # Optimize. 
    from scipy.optimize import leastsq
    leastsq(res.vec, res.values)

    # Get the results
    results = FitResults(res)
    results.show()

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
