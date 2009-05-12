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
"""The FitHook class for inspecting the progress of a FitModel refinement.

"""

import numpy

class FitHook(object):
    """Base class for inspecting the progress of a FitModel refinement.

    Can serve as a fithook for the FitModel class (see FitModel.setFitHook
    method.) The methods in this class are called during the preparation of the
    FitModel for refinement, and during the residual call. See the class
    methods for a description of their purpose.

    Attributes

    count   --  The number of times the residual has been called (default 0).
    verbose --  An integer telling how verbose to be (default 1).
                0   --  print nothing
                1   --  print the count during the precall
                >=2 --  print the residual during the postcall

    """

    def __init__(self):
        """Initialize the attributes."""
        self.count = 0
        self.verbose = 1
        return

    def reset(self):
        """Reset the hook data.

        This is called whenever FitModel._prepare is called.
        """
        self.count = 0
        return

    def precall(self, model):
        """This is called within FitModel.residual, before the calculation.

        model   --  The FitModel instance
        
        """
        self.count += 1
        if self.verbose > 0:
            print self.count
        return

    def postcall(self, model, chiv):
        """This is called within FitModel.residual, after the calculation.

        model   --  The FitModel instance
        chiv    --  The residual vector
        
        """
        if self.verbose < 2:
            return

        # Get the number of restraints
        numres = len(model._restraintlist)
        chi2tot = numpy.dot(chiv, chiv)
        chi2 = 0
        res = 0
        if numres > 0:
            chi = chiv[:-numres]
            resv = chiv[-numres:]
            chi2 = numpy.dot(chi, chi)
            res = numpy.dot(resv, resv)

        print "Residual:", chi2tot
        if res:
            print "Contributions:", chi2
            print "Restraints:", res

__id__ = "$Id$"
