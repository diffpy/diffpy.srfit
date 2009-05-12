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

"""

import numpy

# FIXME - add restraints??

class FitResults(object):
    """Class for processing, presenting and storing results of a fit. 

    Attributes
    model       --  The model containing the results.
    cov         --  The covariance matrix from the model.
    conresults  --  A dictionary of ContributionResults for each contribution,
                    indexed by the contribution name.
    varnames    --  Names of the variables in the model.
    varvals     --  Values of the variables in the model.
    varunc      --  Uncertainties in the variable values.
    connames    --  Names of the constrained parameters.
    convals     --  Values of the constrained parameters.
    conunc      --  Uncertainties in the constraint values.
    residual    --  The scalar residual of the model.
    chi2        --  The chi2 of the model.
    rchi2       --  The reduced chi2 of the model.
    rw          --  The Rw of the model.

    Each of these attributes, except the model, are created or updated when the
    update method is called.

    """

    def __init__(self, model, update = True):
        """Initialize the attributes.

        model   --  The model containing the results
        update  --  Flag indicating whether to do an immediate update (default
                    True).
        
        """
        self.model = model
        self.conresults = {}
        self.varnames = []
        self.varvals = []
        self.varunc = []
        self.connames = []
        self.convals = []
        self.conunc = []
        self.cov = None
        self.residual = 0
        self.chi2 = 0
        self.rchi2 = 0
        self.rw = 0
        self._conjac = None

        if update:
            self.update()
        return

    def update(self):
        """Update the results according to the current state of the model."""
        ## Note that the order of these operations are chosen to reduce
        ## computation time.

        model = self.model

        if not model._organizers:
            return

        # Store the variable names and values
        self.varnames = model.getNames()
        self.varvals = model.getValues()

        # Store the constraint information
        self.connames = [con.par.name for con in model._constraintlist]
        self.convals = [con.par.getValue() for con in model._constraintlist]

        # Calculate the covariance
        self._calculateCovariance()

        # Get the variable uncertainties
        self.varunc = [self.cov[i,i]**0.5 for i in range(len(self.varnames))]

        # Get the constraint uncertainties
        self._calculateConstraintUncertainties()

        # Store the fitting arrays and metrics for each contribution.
        self.conresults = {}
        for con, weight in zip(model._organizers, model._weights):
            self.conresults[con.name] = ContributionResults(con, weight, self)

        # Calculate the metrics
        self.residual = model.scalarResidual()
        self._calculateMetrics()

        return

    def _calculateCovariance(self):
        """Calculate the covariance matrix. This is called by update.

        This code borrowed from PARK. It finds the pseudo-inverse of the
        Jacobian using the singular value decomposition.
        """
        J = self._calculateJacobian()
        u,s,vh = numpy.linalg.svd(J,0)
        self.cov = numpy.dot(vh.T.conj()/s**2,vh)
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
        model = self.model

        # Make sure the input vector is an array
        pvals = numpy.asarray(self.varvals)
        # Compute the numeric derivative using the three point formula.
        delta = step * pvals

        # Center point formula: 
        #     df/dv = lim_{h->0} ( f(v+h)-f(v-h) ) / ( 2h )
        r = []
        # For constraints
        conr = []
        for k,v in enumerate(pvals):
            h = delta[k]
            pvals[k] = v + h
            rk = self.model.residual(pvals)

            cond = []
            for con in model._constraintlist:
                con.update()
                cond.append(con.par.getValue())

            pvals[k] = v - h
            rk -= self.model.residual(pvals)

            for i, con in enumerate(model._constraintlist):
                con.update()
                cond[i] -= con.par.getValue()
                cond[i] /= 2*h
                cond[i] = numpy.abs(cond[i])

            conr.append(cond)

            pvals[k] = v
            r.append(rk/(2*h))
        # Record the constraint derivative matrix
        for con in model._constraintlist:
            con.update()
        self._conjac = numpy.vstack(conr).T
        # return the jacobian
        return numpy.vstack(r).T

    def _calculateMetrics(self):
        """Calculate chi2, rchi2 and Rw for the model."""
        chi2 = 0
        rw = 0
        numpoints = 0
        for con in self.conresults.values():
            chi2 += con.weight * con.chi2
            rw += con.weight * con.rw
            numpoints += len(con.x)

        rchi2 = chi2 / (numpoints - len(self.varnames))

        self.chi2 = chi2
        self.rchi2 = rchi2
        self.rw = rw
        return

    def _calculateConstraintUncertainties(self):
        """Calculate the uncertainty on the constrained parameters."""

        # dp = sum_i |dp/dvi| * unc(vi)
        v = numpy.array(self.varunc)
        dp = abs(numpy.dot(self._conjac, v))
        self.conunc = dp.tolist()
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
            lines.append("")
            lines.append(l)
            lines.append("")

        ## Overall results
        l = "Overall"
        if not certain:
            l += " (Chi2 and Reduced Chi2 invalid)"
        lines.append(l)
        lines.append("-"*79)
        formatstr = "%-14s %-12.8f"
        lines.append(formatstr%("Residual",self.residual))
        lines.append(formatstr%("Chi2",self.chi2))
        lines.append(formatstr%("Reduced Chi2",self.rchi2))
        lines.append(formatstr%("Rw",self.rw))

        ## Per-contribution results
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
        d = {}
        for i, name in enumerate(varnames):
            d[name] = (varvals[i], varunc[i])
        numericStringSort(varnames)
        
        w = max(map(len, varnames))
        w = str(w+1)
        # Format the lines
        formatstr = "%-"+w+"s %- 15f +/- %-15f"
        for name in varnames:
            val, unc = d[name]
            lines.append(formatstr%(name, val, unc))

        ## The constraints
        if self.connames:
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
                    name = ".".join(loc)
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
    """Class for processing, storing contribution results.

    This does not store the contribution.
    
    Attributes
    y       --  The contribution's profile over the calculation range (default
                None).
    dy      --  The uncertainty in the contribution's profile over the
                calculation range (default None).
    x       --  A numpy array of the calculated independent variable for the
                contribution (default None).
    ycalc   --  A numpy array of the calculated signal for the contribution
                (default None).
    residual    --  The scalar residual of the contribution.
    chi2        --  The chi2 of the contribution.
    rw          --  The Rw of the contribution.
    weight      --  The weight of the contribution in the model.
    conlocs     --  The location of the constrained parameters in the
                    contribution.
    convals     --  Values of the constrained parameters.
    conunc      --  Uncertainties in the constraint values.

    """

    def __init__(self, con, weight, fitres):
        """Initialize the attributes.

        con     --  The contribution
        weight  --  The weight of the contribution in the model
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
        ## Note that the order of these operations are chosen to reduce
        ## computation time.

        if con.profile is None:
            return

        model = fitres.model

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
        for i, constraint in enumerate(model._constraintlist):
            par = constraint.par
            loc = con._findParameter(par, [])
            if loc:
                self.conlocs.append(loc)
                self.convals.append(fitres.convals[i])
                self.conunc.append(fitres.conunc[i])

        return

    def _calculateMetrics(self):
        """Calculte chi2 and Rw of the model."""
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
