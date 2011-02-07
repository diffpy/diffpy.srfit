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
"""The FitHook class for inspecting the progress of a refinement.

FitHooks are called by a residual during various times during evaluation.  The
default FitHook simply counts the number of times the residual is called, and
reports that number every time the residual is calculated.  Depending on the
verbosity, it will also report the residual and the current variable values.

Custom FitHooks can be added to a residual with the addFitHook method. These
will be called in the order they are added to the residual.

"""
__all__ = ["FitHook"]

import numpy

class FitHook(object):
    """Base class for inspecting the progress of a refinement.

    Can serve as a fithook for residuals.  The methods in this class are called
    during the preparation of the refinement, and during the calculation of the
    residual. See the class methods for a description of their purpose.

    """

    def reset(self, res):
        """Reset the hook data.

        res --  A Residual instance

        This is called whenever the residual's _setup method is called, which
        is whenever a configurational change to the fit takes place.

        """
        return

    def precall(self, res):
        """This is called before the residual calculation.

        res --  A Residual instance
        
        """
        return

    def postcall(self, res, val):
        """This is called after the residual calculation.

        res --  A Residual instance
        val --  The value of the vector residual.
        
        """
        return

# End class FitHook

class PrintFitHook(FitHook):
    """Base class for inspecting the progress of a refinement.

    This FitHook prints out a running count of the number of times the residual
    has been called, or other information, based on the verbosity.

    Attributes

    count   --  The number of times the residual has been called (default 0).
    verbose --  An integer telling how verbose to be (default 3).
                0   --  print nothing
                1   --  print the count during the precall
                2   --  print the residual during the postcall
                >=3 --  print the variables during the postcall

    """

    def __init__(self):
        """Initialize the attributes."""
        self.count = 0
        self.verbose = 3
        return

    def reset(self, res):
        """Reset the hook data.

        res --  A Residual instance

        This is called whenever the residual's _setup method is called, which
        is whenever a configurational change to the fit takes place.

        """
        self.count = 0
        return

    def precall(self, res):
        """This is called before the residual calculation.

        res --  A Residual instance
        
        """
        self.count += 1
        if self.verbose > 0:
            print self.count
        return

    def postcall(self, res, val):
        """This is called after the residual calculation.

        res --  A Residual instance
        val --  The value of the vector residual.
        
        """
        if self.verbose < 2:
            return

        # Get the number of restraints
        print "Residual:", numpy.dot(val, val)

        if self.verbose >= 3:

            print "Variables"

            vnames = res.names
            vals = res.values
            items = zip(vnames, vals)
            for name, val in items:
                print "  %s = %s" % (name, val)

# End class PrintFitHook

# FIXME - Display the residual on the plot!
class PlotFitHook(FitHook):
    """This FitHook has live plotting of something being refined."""

    def __init__(self, x, ycalc, yobs = None):
        """Initialize the plot fit hook.

        x       --  The independent variable to plot over (Node).
        ycalc   --  An equation to evaluate to obtain the calculated profile
                    (Node).
        yobs    --  The observed profile (Node, default None).

        """
        self.x = x
        self.ycalc = ycalc
        self.yobs = yobs
        return

    def reset(self, res):
        """Set up the plot."""
        FitHook.reset(self, res)

        self._plots = []

        import pylab

        pylab.clf()
        pylab.ion()

        if self.yobs is not None:
            pdata = pylab.plot(self.x.get(), self.yobs.get(), 'bo')[0]
            self._plots.append(pdata)
        pfit = pylab.plot(self.x.get(), self.ycalc.get(), 'r-')[0]
        self._plots.append(pfit)
        pylab.xlabel(self.x.name)
        pylab.ylabel(self.ycalc.name)

        return

    def postcall(self, res, val):
        """This is called after the residual calculation.

        res --  A Residual instance
        val --  The value of the vector residual.
        
        """
        FitHook.postcall(self, res, val)
        import pylab

        if self.yobs is not None:
            pdata, pfit = self._plots
            pdata.set_data(self.x.get(), self.yobs.get())
        else:
            pfit = self._plots[0]
        pfit.set_data(self.x.get(), self.ycalc.get())

        pylab.draw()
        return


__id__ = "$Id$"
