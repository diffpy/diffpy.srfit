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

from diffpy.srfit.adapters.adaptersmod import ContainerAdapter
from diffpy.srfit.adapters.adaptersmod import registry, adapt

import numpy

__all__ = []

def _newcall(calc, r, stru):
    """Overloaded call for PDF calculators to set up calculator and check for
    errors.

    calc    --  The PDFCalculator or DebyePDFCalculator
    r       --  Indexible of points over which to calculate the PDF
    stru    --  The structure from which to calculate the PDF

    """
    calc.rstep = r[1] - r[0]
    calc.rmin = r[0]
    calc.rmax = r[-1] + 0.5 * calc.rstep

    rcalc, y = calc(stru)

    if numpy.isnan(y).any():
        y = numpy.zeros_like(r)
    else:
        y = numpy.interp(r, rcalc, y)
    return y

class PDFCalculatorAdapter(ContainerAdapter):
    """Adapter for PDF calculators from SrReal.

    This adapts as a normal container, but overloads __call__ so that it
    accepts the computation array and the structure used for the calculation.
    This makes it possible to interpolate the output from the PDF calculator on
    the requested grid. A 'parallel' method is also provided to easily set up a
    parallel calculation.

    """

    def __call__(self, x, stru):
        calladapter = adapt(_newcall, self.name)
        return calladapter(self, x, stru)

# End class PDFCalculatorAdapter

# Register the adapter
# FIXME - use abc to make parallel calculators look like these.
from diffpy.srreal.pdfcalculator import PDFCalculator, DebyePDFCalculator
registry[PDFCalculator] = PDFCalculatorAdapter
registry[DebyePDFCalculator] = PDFCalculatorAdapter
