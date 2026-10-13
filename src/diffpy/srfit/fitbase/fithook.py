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
"""The FitHook class for inspecting the progress of a FitRecipe
refinement.

FitHooks are called by a FitRecipe at various stages of the residual
evaluation. The default FitHook of a FitRecipe, PrintFitHook, counts the
number of times the residual is called, and reports that number every
time the residual is calculated. Depending on the verbosity, it will
also report the residual and the current variable values.

Custom FitHooks can be added to a FitRecipe with the
FitRecipe.push_fit_hook method.
"""

from __future__ import print_function

__all__ = ["FitHook"]

import numpy

from diffpy.srfit.util import sort_key_for_numeric_string


class FitHook(object):
    """Base class for inspecting the progress of a FitRecipe refinement.

    A FitHook can be added to a FitRecipe with the
    ``FitRecipe.push_fit_hook`` method. The methods in this class are
    called during the preparation of the FitRecipe for refinement, and
    during the residual call. The base implementations do nothing; derived
    classes overload them as needed.
    """

    def reset(self, recipe):
        """Reset the hook data.

        This is called by ``FitRecipe._prepare`` whenever a
        configurational change to the fit hierarchy takes place, such as
        adding a new ParameterSet, constraint or restraint.

        Parameters
        ----------
        recipe : FitRecipe
            The FitRecipe that is being prepared.
        """
        return

    def precall(self, recipe):
        """Inspect the recipe before the residual is calculated.

        This is called within ``FitRecipe.residual``, before the
        variable values are applied and the residual is calculated.

        Parameters
        ----------
        recipe : FitRecipe
            The FitRecipe whose residual is being calculated.
        """
        return

    def postcall(self, recipe, chiv):
        """Inspect the recipe after the residual is calculated.

        This is called within ``FitRecipe.residual``, after the residual
        is calculated.

        Parameters
        ----------
        recipe : FitRecipe
            The FitRecipe whose residual was calculated.
        chiv : ndarray
            The residual vector, including the restraint penalties.
        """
        return


# End class FitHook


class PrintFitHook(FitHook):
    """Print the progress of a FitRecipe refinement.

    This FitHook prints out a running count of the number of times the residual
    has been called, or other information, based on the verbosity.

    Attributes
    ----------
    count : int
        The number of times the residual has been called (default 0).
    verbose : int
        The verbosity level (default 1). Each level also prints what the
        lower levels print.

        0
            Print nothing.
        1
            Print the count during the precall.
        2
            Print the residual during the postcall.
        >=3
            Print the variables during the postcall.
    """

    def __init__(self):
        self.count = 0
        self.verbose = 1
        return

    def reset(self, recipe):
        """Reset the residual call count to 0.

        This is called by ``FitRecipe._prepare`` whenever a
        configurational change to the fit hierarchy takes place, such as
        adding a new ParameterSet, constraint or restraint.

        Parameters
        ----------
        recipe : FitRecipe
            The FitRecipe that is being prepared.
        """
        self.count = 0
        return

    def precall(self, recipe):
        """Increment the residual call count and print it.

        The count is printed when ``verbose`` is at least 1.

        Parameters
        ----------
        recipe : FitRecipe
            The FitRecipe whose residual is being calculated.
        """
        self.count += 1
        if self.verbose > 0:
            print(self.count)
        return

    def postcall(self, recipe, chiv):
        """Print the residual and, optionally, the variable values.

        When ``verbose`` is at least 2, print the total residual and, if
        the recipe has restraints, its FitContribution and restraint
        parts. When ``verbose`` is at least 3, also print the name and
        value of each variable.

        Parameters
        ----------
        recipe : FitRecipe
            The FitRecipe whose residual was calculated.
        chiv : ndarray
            The residual vector, including the restraint penalties.
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
            vnames = recipe.get_names()
            vals = recipe.get_values()
            # byname = _byname()
            items = sorted(zip(vnames, vals), key=_byname)
            for name, val in items:
                print("  %s = %f" % (name, val))
        return


def _byname(nv):
    return sort_key_for_numeric_string(nv[0])


# End class PrintFitHook


# TODO - Display the chi^2 on the plot during refinement.
class PlotFitHook(FitHook):
    """Live-plot the progress of a FitRecipe refinement.

    Each FitContribution of the recipe gets a subplot showing the
    observed profile and the calculated profile, which is updated after
    every residual call. This requires matplotlib.

    Attributes
    ----------
    _plots : list of tuple
        The (observed, calculated) line pair for each FitContribution.
        It is created by ``reset``.
    """

    def reset(self, recipe):
        """Set up the plot.

        Clear the current figure, turn on interactive mode and draw one
        subplot per FitContribution of the recipe.

        Parameters
        ----------
        recipe : FitRecipe
            The FitRecipe that is being prepared.
        """
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
        """Update the plots with the latest profiles.

        This is called within ``FitRecipe.residual``, after the residual
        is calculated.

        Parameters
        ----------
        recipe : FitRecipe
            The FitRecipe whose residual was calculated.
        chiv : ndarray
            The residual vector, including the restraint penalties.
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
