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
"""The FitHook class for inspecting the progress of a FitRecipe refinement.

FitHooks are called by a FitRecipe during various times of the residual is
evaluation. The default FitHook simply counts the number of times the residual
is called, and reports that number every time the residual is calculated.
Depending on the verbosity, it will also report the residual and the current
variable values.

Custom FitHooks can be added to a FitRecipe with the FitRecipe.setFitHook method.

"""

import numpy

class FitHook(object):
    """Base class for inspecting the progress of a FitRecipe refinement.

    Can serve as a fithook for the FitRecipe class (see FitRecipe.setFitHook
    method.) The methods in this class are called during the preparation of the
    FitRecipe for refinement, and during the residual call. See the class
    methods for a description of their purpose.

    Attributes

    count   --  The number of times the residual has been called (default 0).
    verbose --  An integer telling how verbose to be (default 1).
                0   --  print nothing
                1   --  print the count during the precall
                2   --  print the residual during the postcall
                >=3 --  print the variables during the postcall

    """

    def __init__(self):
        """Initialize the attributes."""
        self.count = 0
        self.verbose = 1
        return

    def reset(self):
        """Reset the hook data.

        This is called whenever FitRecipe._prepare is called.
        """
        self.count = 0
        return

    def precall(self, recipe):
        """This is called within FitRecipe.residual, before the calculation.

        recipe   --  The FitRecipe instance
        
        """
        self.count += 1
        if self.verbose > 0:
            print self.count
        return

    def postcall(self, recipe, chiv):
        """This is called within FitRecipe.residual, after the calculation.

        recipe   --  The FitRecipe instance
        chiv    --  The residual vector
        
        """
        if self.verbose < 3:
            return

        # Get the number of restraints
        numres = len(recipe._restraintlist)
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
            print "FitContributions:", chi2
            print "Restraints:", res

        if self.verbose >= 3:

            print "Variables"

            vnames = recipe.getNames()
            vals = recipe.getValues()

            for name, val in zip(vnames, vals):
                print "  %s = %f" % (name, val)

__id__ = "$Id$"
