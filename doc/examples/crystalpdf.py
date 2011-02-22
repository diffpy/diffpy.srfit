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
"""Example of a PDF refinement using diffpy.Structure and SrReal.

This is example of fitting the fcc nickel structure to measured PDF data. The
purpose of this example is to demonstrate and describe the classes in
configuration options involved with setting up a fit in this way. The main
benefit of using SrFit for PDF refinement is the flexibility of modifying the
PDF profile function for specific needs, adding restraints to a fit and the
ability to simultaneously refine a structure to PDF data and data from other
sources. This example demonstrates only the basic configuration.

"""

from numpy import pi

from diffpy.Structure import Structure
from diffpy.srreal.pdfcalculator import PDFCalculator
# This import statement registers pdf-related adapters
import diffpy.srfit.pdf
from diffpy.srfit import *

def main(ciffile, datname):
    """Fit crystalline PDF data using SrFit."""

    # Load the data. We use 'loadProfile' as in earlier examples, but now we
    # want to pull out metadata specific to PDF data. The 
    # 'import diffpy.srfit.pdf' statement above registered a parser designed
    # just for this. We tell 'loadProfile' to use this parser by specifying the
    # 'PDF' format.
    profile = loadProfile(datname, "PDF")

    # Create a PDFCalculator. We can use information from the metadata loaded
    # from file to initialize it. We can't pass the whole meta dictionary since
    # it contains entries that PDFCalculator can't use. In addition, passing a
    # non-zero qmin value to PDFCalculator creates a "bad" PDF.
    calc = PDFCalculator(qmax = profile.meta["qmax"], 
            rmin=profile.x[0], rmax=20)
    calc.rstep = pi / calc.qmax
    calc.setScatteringFactorTableByType(profile.meta["stype"])

    # Use the same grid as the profile
    profile.setPoints(calc.rgrid)
    r, gr, dgr = profile.pars

    # Now we adapt the PDFCalculator so it can be used as a symbolic
    # calculation object. We must pass our PDFCalculator instance, and
    # optionally give it a name.
    g = adapt(calc, "g")
    # While we're at it, define which parameters we want to vary and give them
    # initial values. Since we're using the adapted calculator for attribute
    # access, these objects are adapted as parameters.
    g.scale.vary(1)
    g.qdamp.vary(0.01)
    g.delta2.vary(5)

    # Now we create the structure object we want to refine and then adapt it.
    stru = Structure(filename = ciffile)
    stru.Uisoequiv = 0.005
    s = adapt(stru, "nickel")
    
    # Now we must use the adapters to identify the parameters that we want to
    # refine and how these are inter-related. There is a
    # 'constrainAsSpaceGroup' function that simplifies this - we're doing it
    # the hard way for the sake of demonstration.
    #
    # First, the lattice parameters.
    a = s.lattice.a.vary()
    s.lattice.b.constrain(a)
    s.lattice.c.constrain(a)
    # Next, the ADPs.
    Uisoequiv = s[0].Uisoequiv.vary()
    for atom in s[1:]:
        atom.Uisoequiv.constrain(Uisoequiv)

    # Now we symbolically create the fit equation. The PDFCalculator returns
    # the calculation range and calculated PDF when called. We capture these
    # symbolically.
    rcalc, gcalc = g(s)

    # Because our data is on the grid defined by 'r', and the fit is on the
    # grid defined by 'rcalc', we have to interpolate in order to compare them.
    # The numpy function interp has been adapted for this purpose.
    fiteq = gcalc#interp_(r, rcalc, gcalc)

    # Create the residual equation. Note that 'chi' creates a vector residual
    # that can be dotted into itself to generate 'chi^2'.
    reseq = chi(fiteq, gr, dgr)
    # Create a restraint on Uisoequiv. Here we pretend that we have prior
    # knowledge that the value of Uisoequiv is 0.0055 with uncertainty 0.0005.
    # We create a restraint that effectively adds this input as data point to
    # the fit.
    rest1 = restrain(Uisoequiv, 0.0055, sig = 0.0005)
    # Create the residual object. The restraint adds to the cost of the
    # residual - we pass it as an argument.
    res = residual(reseq, rest1)

    # Optimize. 
    #from scipy.optimize import leastsq
    #leastsq(res.vec, res.values)
    from scipy.optimize import fmin
    fmin(res, res.values)

    # Get the results
    results = FitResults(res)
    results.show()

    # and plot 
    from pylab import plot, show
    plot(r.value, gr.value, 'bo')
    plot(r.value, fiteq.value, 'r-')
    show()

if __name__ == "__main__":

    # Make the data and the recipe
    ciffile = "data/ni.cif"
    data = "data/ni-q27r100-neutron.gr"
    #data = "data/ni-q27r60-xray.gr"
    main(ciffile, data)

# End of file
