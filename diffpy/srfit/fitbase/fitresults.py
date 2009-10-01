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
"""The FitResults and ContributionResults classes for storing results of a fit.

The FitResults class is used to display the current state of a FitRecipe. It
stores the state, and uses it to calculate useful statistics, which can be
displayed on screen or saved to file.

"""

import numpy

# FIXME - really small numbers don't print well

class FitResults(object):
    """Class for processing, presenting and storing results of a fit. 

    Attributes
    recipe       --  The recipe containing the results.
    cov         --  The covariance matrix from the recipe.
    conresults  --  A dictionary of ContributionResults for each
                    FitContribution, indexed by the FitContribution name.
    varnames    --  Names of the variables in the recipe.
    varvals     --  Values of the variables in the recipe.
    varunc      --  Uncertainties in the variable values.
    showcon     --  Show constraint values in the output (default False).
    connames    --  Names of the constrained parameters.
    convals     --  Values of the constrained parameters.
    conunc      --  Uncertainties in the constraint values.
    residual    --  The scalar residual of the recipe.
    penalty     --  The penalty to residual from the restraints.
    chi2        --  The chi2 of the recipe.
    rchi2       --  The reduced chi2 of the recipe.
    rw          --  The Rw of the recipe.
    messages    --  A list of messages about the results.
    _dcon       --  The derivatives of the constraint equations with respect to
                    the variables. This is used internally.

    Each of these attributes, except the recipe, are created or updated when the
    update method is called.

    """

    def __init__(self, recipe, update = True, showcon = False):
        """Initialize the attributes.

        recipe   --  The recipe containing the results
        update  --  Flag indicating whether to do an immediate update (default
                    True).
        showcon --  Show constraint values in the output (default False).
        
        """
        self.recipe = recipe
        self.conresults = {}
        self.varnames = []
        self.varvals = []
        self.varunc = []
        self.connames = []
        self.convals = []
        self.conunc = []
        self.cov = None
        self.residual = 0
        self.penalty = 0
        self.chi2 = 0
        self.rchi2 = 0
        self.rw = 0
        self._dcon = []
        self.messages = []

        self.showcon = showcon

        if update:
            self.update()
        return

    def update(self):
        """Update the results according to the current state of the recipe."""
        ## Note that the order of these operations are chosen to reduce
        ## computation time.

        recipe = self.recipe

        if not recipe._contributions:
            return

        # Make sure everything is ready for calculation
        recipe._prepare()

        # Store the variable names and values
        self.varnames = recipe.getNames()
        self.varvals = recipe.getValues()

        # Store the constraint information
        self.connames = [con.par.name for con in recipe._oconstraints]
        self.convals = [con.par.getValue() for con in recipe._oconstraints]

        # Calculate the covariance
        self._calculateCovariance()

        # Get the variable uncertainties
        self.varunc = [self.cov[i,i]**0.5 for i in range(len(self.varnames))]

        # Get the constraint uncertainties
        self._calculateConstraintUncertainties()

        # Store the fitting arrays and metrics for each FitContribution.
        self.conresults = {}
        for con, weight in zip(recipe._contributions.values(), recipe._weights):
            self.conresults[con.name] = ContributionResults(con, weight, self)

        # Calculate the metrics
        res = recipe.residual()
        self.residual = numpy.dot(res, res)
        self._calculateMetrics()

        # Calcualte the restraints penalty
        w = self.chi2 / len(res)
        self.penalty = sum([res.penalty(w) for res in recipe._restraintlist])

        return

    def _calculateCovariance(self):
        """Calculate the covariance matrix. This is called by update.

        This code borrowed from PARK. It finds the pseudo-inverse of the
        Jacobian using the singular value decomposition.

        """
        try:
            J = self._calculateJacobian()
            u,s,vh = numpy.linalg.svd(J,0)
            self.cov = numpy.dot(vh.T.conj()/s**2,vh)
        except numpy.linalg.LinAlgError:
            self.messages.append("Cannot compute covariance matrix.")
            l = len(self.varvals)
            self.cov = numpy.zeros((l, l), dtype=float)
        return

    def _calculateJacobian(self, step=1e-8):
        """Calculate the Jacobian for the fitting.

        Borrowed from PARK.
        Returns the derivative wrt the fit variables at point p.

        This also calculates the derivatives of the constrained parameters
        while we're at it.

        Numeric derivatives are calculated based on step, where step is the
        portion of variable value. E.g. step = dv/v.

        """
        recipe = self.recipe

        # Make sure the input vector is an array
        pvals = numpy.asarray(self.varvals)
        # Compute the numeric derivative using the center point formula.
        delta = step * pvals

        # Center point formula: 
        #     df/dv = lim_{h->0} ( f(v+h)-f(v-h) ) / ( 2h )
        #

        r = []
        # The list of constraint derivatives with respect to variables 
        # The forward difference would be faster, but perhaps not as accurate.
        conr = []
        for k,v in enumerate(pvals):
            h = delta[k]
            pvals[k] = v + h
            rk = self.recipe.residual(pvals)

            # The constraints derivatives
            cond = []
            for con in recipe._oconstraints:
                con.update()
                cond.append(con.par.getValue())

            pvals[k] = v - h
            rk -= self.recipe.residual(pvals)

            for i, con in enumerate(recipe._oconstraints):
                con.update()
                cond[i] -= con.par.getValue()
                cond[i] /= 2*h

            conr.append(cond)

            pvals[k] = v
            r.append(rk/(2*h))
            
        # Reset the constrained parameters to their original values
        for con in recipe._oconstraints:
            con.update()

        self._dcon = numpy.vstack(conr).T

        # return the jacobian
        return numpy.vstack(r).T

    def _calculateMetrics(self):
        """Calculate chi2, rchi2 and Rw for the recipe."""
        chi2 = 0
        rw = 0
        numpoints = 0
        for con in self.conresults.values():
            chi2 += con.weight * con.chi2
            rw += con.weight * con.rw
            numpoints += len(con.x)

        numpoints += len(self.recipe._restraintlist)

        rchi2 = chi2 / (numpoints - len(self.varnames))

        self.chi2 = chi2
        self.rchi2 = rchi2
        self.rw = rw
        return

    def _calculateConstraintUncertainties(self):
        """Calculate the uncertainty on the constrained parameters."""
        vu = self.varunc

        # sig^2(c) = sum_i sum_j sig(v_i) sig(v_j) (dc/dv_i)(dc/dv_j)
        # sig^2(c) = sum_i sum_j [sig(v_i)(dc/dv_i)][sig(v_j)(dc/dv_j)]
        # sig^2(c) = sum_i sum_j u_i u_j
        self.conunc = []
        for dci in self._dcon:

            # Create sig(v_i) (dc/dv_i) array.
            u = dci * vu
            # The outer product is all possible pairings of u_i and u_j
            # uu_ij = u_i u_j
            uu = numpy.outer(u, u)
            # Sum these pairings to get sig^2(c)
            sig2c = sum(uu.flatten())

            self.conunc.append(sig2c**0.5)
        return

    def formatResults(self, header = "", footer = "", update = False):
        """Format the results and return them in a string.

        This function is called by printResults and saveResults. Overloading
        the formatting here will change all three functions.

        header  --  A header to add to the output (default "")
        footer  --  A footer to add to the output (default "")
        update  --  Flag indicating whether to call update() (default False).

        Returns a string containing the formatted results.
        
        """
        if update:
            self.update()

        lines = []
        corrmin = 0.25

        # Check to see if the uncertainty values are reliable.
        certain = True
        for con in self.conresults.values():
            if numpy.array_equal(con.dy, numpy.ones_like(con.dy)):
                certain = False
                break

        # User-defined header
        if header:
            lines.append(header)

        if not certain:
            l = "Some quantities invalid due to missing profile uncertainty"
            self.messages.append(l)

        lines.extend(self.messages)

        ## Overall results
        l = "Overall"
        if not certain:
            l += " (Chi2 and Reduced Chi2 invalid)"
        lines.append(l)
        lines.append("-"*79)
        formatstr = "%-14s %-12.8f"
        lines.append(formatstr%("Residual",self.residual))
        lines.append(formatstr%("Contributions", self.residual - self.penalty))
        lines.append(formatstr%("Restraints", self.penalty))
        lines.append(formatstr%("Chi2",self.chi2))
        lines.append(formatstr%("Reduced Chi2",self.rchi2))
        lines.append(formatstr%("Rw",self.rw))

        ## Per-FitContribution results
        if len(self.conresults) > 1:
            keys = self.conresults.keys()
            numericStringSort(keys)

            lines.append("")
            l = "Contributions"
            if not certain:
                l += " (Chi2 and Reduced Chi2 invalid)"
            lines.append(l)
            lines.append("-"*79)
            formatstr = "%-10s %-42.8f"
            for name in keys:
                res = self.conresults[name]
                lines.append("")
                namestr = name + " (%f)"%res.weight
                lines.append(namestr)
                lines.append("-"*len(namestr))
                lines.append(formatstr%("Residual",res.residual))
                lines.append(formatstr%("Chi2",res.chi2))
                lines.append(formatstr%("Rw",res.rw))

        ## The variables
        lines.append("")

        l = "Variables"
        if not certain:
            m = "Uncertainties invalid"
            l += " (%s)"%m
        lines.append(l)
        lines.append("-"*79)

        varnames = self.varnames
        varvals = self.varvals
        varunc = self.varunc
        varlines = []

        w = max(map(len, varnames))
        w = str(w+1)
        # Format the lines
        formatstr = "%-"+w+"s %- 15f +/- %-15f"
        for name, val, unc in zip(varnames, varvals, varunc):
            varlines.append(formatstr%(name, val, unc))

        varlines.sort()
        lines.extend(varlines)

        ## The constraints
        if self.connames and self.showcon:
            lines.append("")
            l = "Constrained Parameters"
            if not certain:
                l += " (Uncertainties invalid)"
            lines.append(l)
            lines.append("-"*79)

            w = 0
            keys = []
            vals = {}
            for con in self.conresults.values():
                for i, loc in enumerate(con.conlocs):
                    names = [obj.name for obj in loc]
                    name = ".".join(names)
                    w = max(w, len(name))
                    val = con.convals[i]
                    unc = con.conunc[i]
                    keys.append(name)
                    vals[name] = (val, unc)
                    
            numericStringSort(keys)
            w = str(w+1)
            formatstr = "%-"+w+"s %- 15f +/- %-15f"
            for name in keys:
                val, unc = vals[name]
                lines.append(formatstr%(name, val, unc))

        ## Variable correlations
        lines.append("")
        corint = int(corrmin*100)
        l = "Variable Correlations greater than %i%%"%corint
        if not certain:
            l += " (Correlations invalid)"
        lines.append(l)
        lines.append("-"*79)
        tup = []
        cornames = []
        n = len(varnames)
        for i in xrange(n):
            for j in xrange(i+1, n):
                name = "corr(%s, %s)"%(varnames[i], varnames[j])
                val = (self.cov[i,j]/(self.cov[i,i] * self.cov[j,j])**0.5)
                if val > corrmin:
                    cornames.append(name)
                    tup.append((val, name))

        tup.sort()
        tup.reverse()

        if cornames:
            w = max(map(len, cornames))
            w = str(w + 1)
            formatstr = "%-"+w+"s %- 8.4f"
            for val, name in tup:
                lines.append(formatstr%(name, val))
        else:
            lines.append("No correlations greater than %i%%"%corint)


        # User-defined footer
        if footer:
            lines.append(footer)

        out = "\n".join(lines)
        return out

    def printResults(self, header = "", footer = "", update = False):
        """Format and print the results.

        header  --  A header to add to the output (default "")
        footer  --  A footer to add to the output (default "")
        update  --  Flag indicating whether to call update() (default False).

        """
        print self.formatResults(header, footer, update)
        return

    def saveResults(self, filename, header = "", footer = "", update = False):
        """Format and save the results.

        filename -  Name of the save file.
        header  --  A header to add to the output (default "")
        footer  --  A footer to add to the output (default "")
        update  --  Flag indicating whether to call update() (default False).

        """
        # Save the time and user
        from time import ctime
        from getpass import getuser
        myheader = "Results written: " + ctime() + "\n"
        myheader += "produced by " + getuser() + "\n"
        header = myheader + header

        res = self.formatResults(header, footer, update)
        f = file(filename, 'w')
        f.write(res)
        f.close()
        return

# End class FitResults

class ContributionResults(object):
    """Class for processing, storing FitContribution results.

    This does not store the FitContribution.
    
    Attributes
    y       --  The FitContribution's profile over the calculation range
                (default None).
    dy      --  The uncertainty in the FitContribution's profile over the
                calculation range (default None).
    x       --  A numpy array of the calculated independent variable for the
                FitContribution (default None).
    ycalc   --  A numpy array of the calculated signal for the FitContribution
                (default None).
    residual    --  The scalar residual of the FitContribution.
    chi2        --  The chi2 of the FitContribution.
    rw          --  The Rw of the FitContribution.
    weight      --  The weight of the FitContribution in the recipe.
    conlocs     --  The location of the constrained parameters in the
                    FitContribution (see the
                    RecipeContainer._locateManagedObject method).  
    convals     --  Values of the constrained parameters.
    conunc      --  Uncertainties in the constraint values.

    """

    def __init__(self, con, weight, fitres):
        """Initialize the attributes.

        con     --  The FitContribution
        weight  --  The weight of the FitContribution in the recipe
        fitres  --  The FitResults instance to contain this ContributionResults
        
        """
        self.x = None
        self.y = None
        self.dy = None
        self.ycalc = None
        self.residual = 0
        self.chi2 = 0
        self.rw = 0
        self.weight = 0
        self.conlocs = []
        self.convals = []
        self.conunc = []
        self._init(con, weight, fitres)
        return

    def _init(self, con, weight, fitres):
        """Initialize the attributes, for real."""
        ## Note that the order of these operations is chosen to reduce
        ## computation time.

        if con.profile is None:
            return

        recipe = fitres.recipe

        # Store the weight
        self.weight = weight

        # First the residual
        res = con.residual()
        self.residual = numpy.dot(res, res)

        # The arrays
        self.x = numpy.array(con.profile.x)
        self.y = numpy.array(con.profile.y)
        self.dy = numpy.array(con.profile.dy)
        self.ycalc = numpy.array(con.profile.ycalc)
        
        # The other metrics
        self._calculateMetrics()

        # Find the parameters
        for i, constraint in enumerate(recipe._oconstraints):
            par = constraint.par
            loc = con._locateManagedObject(par)
            if loc:
                self.conlocs.append(loc)
                self.convals.append(fitres.convals[i])
                self.conunc.append(fitres.conunc[i])

        return

    def _calculateMetrics(self):
        """Calculte chi2 and Rw of the recipe."""
        # We take absolute values in case the signal is complex
        num = numpy.abs(self.y - self.ycalc)
        y = numpy.abs(self.y)
        chiv = num/self.dy
        self.chi2 = numpy.dot(chiv, chiv)
        self.rw = (numpy.dot(num, num) / numpy.dot(y, y))**0.5
        return


def numericStringSort(lst):
    """Sort list of strings inplace according to general numeric value.

    Each string gets split to string and integer segments to create keys
    for comparison.  Signs, decimal points and exponents are ignored.
    
    lst  -- sorted list of strings
    
    No return value to highlight inplace sorting.

    """
    import re
    rx = re.compile(r'(\d+)')
    keys = [ rx.split(s) for s in lst ]
    for k in keys:  k[1::2] = [ int(i) for i in k[1::2] ]
    newlst = zip(keys, lst)
    newlst.sort()
    lst[:] = [kv[1] for kv in newlst]
    return

__id__ = "$Id$"
