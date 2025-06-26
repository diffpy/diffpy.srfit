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
"""The FitHook class for inspecting the progress of a FitRecipe refinement.

FitHooks are called by a FitRecipe during various times of the residual
is evaluation. The default FitHook simply counts the number of times the
residual is called, and reports that number every time the residual is
calculated. Depending on the verbosity, it will also report the residual
and the current variable values.

Custom FitHooks can be added to a FitRecipe with the
FitRecipe.setFitHook method.
"""

from __future__ import print_function

__all__ = ["FitHook"]

import numpy

from diffpy.srfit.util import sortKeyForNumericString


class FitHook(object):
    """Base class for inspecting the progress of a FitRecipe refinement.

    Can serve as a fithook for the FitRecipe class (see
    FitRecipe.pushFitHook method.) The methods in this class are called
    during the preparation of the FitRecipe for refinement, and during
    the residual call. See the class methods for a description of their
    purpose.
    """

    def reset(self, recipe):
        """Reset the hook data.

        This is called whenever FitRecipe._prepare is called, which is
        whenever a configurational change to the fit hierarchy takes
        place, such as adding a new ParameterSet, constraint or
        restraint.
        """
        return

    def precall(self, recipe):
        """This is called within FitRecipe.residual, before the calculation.

        recipe  --  The FitRecipe instance
        """
        return

    def postcall(self, recipe, chiv):
        """This is called within FitRecipe.residual, after the calculation.

        recipe  --  The FitRecipe instance
        chiv    --  The residual vector
        """
        return


# End class FitHook


class PrintFitHook(FitHook):
    """Base class for inspecting the progress of a FitRecipe refinement.

    This FitHook prints out a running count of the number of times the residual
    has been called, or other information, based on the verbosity.

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

    def reset(self, recipe):
        """Reset the hook data.

        This is called whenever FitRecipe._prepare is called, which is
        whenever a configurational change to the fit hierarchy takes
        place, such as adding a new ParameterSet, constraint or
        restraint.
        """
        self.count = 0
        return

    def precall(self, recipe):
        """This is called within FitRecipe.residual, before the calculation.

        recipe  --  The FitRecipe instance
        """
        self.count += 1
        if self.verbose > 0:
            print(self.count)
        return

    def postcall(self, recipe, chiv):
        """This is called within FitRecipe.residual, after the calculation.

        recipe  --  The FitRecipe instance
        chiv    --  The residual vector
        """
        if self.verbose < 2:
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

        print("Residual:", chi2tot)
        if numres:
            print("FitContributions:", chi2)
            print("Restraints:", res)

        if self.verbose >= 3:
            print("Variables")
            vnames = recipe.getNames()
            vals = recipe.getValues()
            # byname = _byname()
            items = sorted(zip(vnames, vals), key=_byname)
            for name, val in items:
                print("  %s = %f" % (name, val))
        return


def _byname(nv):
    return sortKeyForNumericString(nv[0])


# End class PrintFitHook


# TODO - Display the chi^2 on the plot during refinement.
class PlotFitHook(FitHook):
    """This FitHook has live plotting of whatever is being refined."""

    def reset(self, recipe):
        """Set up the plot."""
        FitHook.reset(self, recipe)

        self._plots = []

        import pylab

        pylab.clf()
        pylab.ion()

        nc = len(recipe._contributions)
        if nc > 1:
            ncols = 2
            nrows = (nc + 1) / 2

        for idx, c in enumerate(recipe._contributions.values()):

            name = c.name
            xname = c._xname
            yname = c._yname
            p = c.profile

            # Create a subplot
            if nc > 1:
                pylab.subplot(nrows, ncols, idx + 1)
            pdata = pylab.plot(p.x, p.y, "bo")[0]
            pfit = pylab.plot(p.x, p.y, "r-")[0]
            self._plots.append((pdata, pfit))
            pylab.xlabel(xname)
            pylab.ylabel(yname)
            pylab.title(name)

        # Set up some event handling, so things behave nicely.
        # def redraw(event):
        #    canvas = event.canvas
        #    canvas.draw()
        #    return
        # pylab.connect('resize_event', redraw)

        return

    def postcall(self, recipe, chiv):
        """This is called within FitRecipe.residual, after the calculation.

        Find data and plot it.

        recipe  --  The FitRecipe instance
        chiv    --  The residual vector
        """
        FitHook.postcall(self, recipe, chiv)
        import pylab

        for c, plottup in zip(recipe._contributions.values(), self._plots):

            p = c.profile
            pdata = plottup[0]
            pfit = plottup[1]
            pdata.set_data(p.x, p.y)
            pfit.set_data(p.x, p.ycalc)

        pylab.draw()
        return
