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
"""The FitResults class for storing results of a fit.

The FitResults class is used to display the current state of a fit. It stores
the state, and uses it to calculate useful statistics, which can be displayed
on screen or saved to file.

"""
__all__ = ["FitResults"]

import numpy

from diffpy.srfit.util.inpututils import inputToString
from diffpy.srfit.util import getParameters, absName, isConstrained

class FitResults(object):
    """Class for processing, presenting and storing results of a fit. 

    Attributes
    res         --  The Residual containing the results.
    cov         --  The covariance matrix from the residual.
    derivstep   --  The fractional step size for calculating numeric
                    derivatives. Default 1e-8.
    varnames    --  Names of the variables in the residual.
    varvals     --  Values of the variables in the residual.
    varunc      --  Uncertainties in the variable values.
    showfixed   --  Show fixed variables (default False).
    fixednames  --  Names of the fixed variables of the residual.
    fixedvals   --  Values of the fixed variables of the residual.
    showcon     --  Show constraint values in the output (default True).
    constrained --  List of constrained parameters.
    connames    --  Names of the constrained parameters.
    convals     --  Values of the constrained parameters.
    conunc      --  Uncertainties in the constraint values.
    resval      --  The scalar residual value.
    rresval     --  Reduced residual: resval / numpoints
    numpoints   --  Estimated number of independent pieces of information. This
                    is the lenght of the residual vector less the number of
                    variables.
    resterms    --  List of names and values from the residual terms.
    messages    --  A list of messages about the results.
    precision   --  The precision of numeric output (default 8).
    _dcon       --  The derivatives of the constraint equations with respect to
                    the variables. This is used internally.

    Each of these attributes, except the residual, are created or updated when
    the update method is called.

    """

    def __init__(self, res, showcon = True, showfixed = False, update = True):
        """Initialize the attributes.

        res     --  The Residual containing the results
        showcon --  Show constraint values in the output (bool).
        showfixed --  Show fixed variables in the output (bool).
        update  --  Flag indicating whether to do an immediate update (bool).
        
        """
        self.res = res
        self.derivstep = 1e-8
        self.varnames = []
        self.varvals = []
        self.varunc = []
        self.fixednames = []
        self.fixedvals = []
        self.connames = []
        self.convals = []
        self.conunc = []
        self.cov = None
        self.resval = 0
        self.rresval = 0
        self.numpoints = 0
        self.resterms = []
        self.precision = 8
        self._dcon = []
        self.messages = []
        self._certain = True

        self.showfixed = bool(showfixed)
        self.showcon = bool(showcon)

        if update:
            self.update()
        return

    def update(self):
        """Update results according to the current state of the residual."""
        ## Note that the order of these operations are chosen to reduce
        ## computation time.

        res = self.res

        if not res.terms:
            return

        # Store the variable names and values
        self.varnames = [absName(v) for v in res.variables]
        self.varvals = res.values
        # Get the fixed parameter names and values
        fixedpars = set()
        parsets = (getParameters(term) for term in res.terms)
        map(fixedpars.update, parsets)
        fixedpars.difference_update(res.variables)
        # Pull out the constrained ones
        constrained = filter(isConstrained, fixedpars)
        fixedpars.difference_update(constrained)
        fixedpars = fixedpars
        self.fixednames = [absName(p) for p in fixedpars]
        self.fixedvals = [p.value for p in fixedpars]
        # Get the names and values of the constraints
        self.constrained = constrained
        self.connames = [absName(p) for p in constrained]
        self.convals = [p.value for p in constrained]


        if self.varnames:
            # Calculate the covariance
            self._calculateCovariance()

            # Get the variable uncertainties
            self.varunc = [self.cov[i,i]**0.5 for i in \
                    range(len(self.varnames))]

            # Get the constraint uncertainties
            self._calculateConstraintUncertainties()

        # Store the fitting arrays and metrics for each FitContribution.

        # Calculate the metrics
        self._calculateMetrics()

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

    def _calculateJacobian(self):
        """Calculate the Jacobian for the fitting.

        This also calculates the derivatives of the constrained parameters.

        Adapted from PARK, by the DANSE reflectometry team lead by Paul
        Kienzle.

        Returns the derivative wrt the fit variables at point p.

        Numeric derivatives are calculated based on step, where step is the
        portion of variable value. E.g. step = dv/v.

        """
        res = self.res
        step = self.derivstep

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
            rk = res.vec(pvals)

            # The constraints
            cond = []
            for par in self.constrained:
                cond.append(par.value)

            
            # Set the next parameter values
            pvals[k] = v - h
            rk -= res.vec(pvals)

            # FIXME - constraints are used for vectors as well!
            for i, par in enumerate(self.constrained):
                val = par.value
                if numpy.isscalar(val):
                    cond[i] -= val
                    cond[i] /= 2*h
                else:
                    cond[i] = 0.0

            conr.append(cond)

            pvals[k] = v
            r.append(rk/(2*h))
            
        self._dcon = numpy.vstack(conr).T

        # return the jacobian
        jac = numpy.vstack(r).T
        return jac

    def _calculateMetrics(self):
        """Calculate chi2, rchi2 and Rw for the recipe."""
        self.numpoints = len(self.res.vec())
        self.numpoints -= len(self.res.variables)
        self.resval = self.res()
        self.rresval = self.resval / self.numpoints

        self.resterms = []
        for term in self.res.terms:
            val = term.value
            # Make this a scalar, if it is not already
            if not numpy.isscalar(val):
                val = numpy.dot(val, val)
            self.resterms.append((term.name, val))

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

    def format(self, header = "", footer = "", update = False):
        """Format the results and return them in a string.

        This function is called by show and save. Overloading
        the formatting here will change all three functions.

        header  --  A header to add to the output.
        footer  --  A footer to add to the output.
        update  --  Flag indicating whether to call update() (bool).

        Returns a string containing the formatted results.
        
        """
        if update:
            self.update()

        lines = []
        corrmin = 0.25
        p = self.precision
        pe = "%-" + "%i.%ie" % (p+6, p)

        # User-defined header
        if header:
            lines.append(str(header))

        lines.extend(self.messages)

        ## Overall results
        l = "Overall"
        lines.append(l)
        dashedline = 79 * '-'
        lines.append(dashedline)
        formatstr = "%-18s %-12.8f"
        lines.append(formatstr%("Residual",self.resval))
        lines.append(formatstr%("Reduced residual",self.rresval))

        # See if the uncertainties are valid. Right now we do it based on the
        # names of the residual terms.
        invalid = False

        ## Per-term results
        if self.resterms:
            lines.append("")
            l = "Terms of the residual"
            lines.append(l)
            lines.append(dashedline)
            formatstr = "%-18s %-12.8f"
            for resname, resval in self.resterms:
                if resname != "chi^2": invalid = True
                lines.append(formatstr%(resname,resval))

        ## The variables
        if self.varnames:
            lines.append("")
            l = "Variables"
            if invalid:
                l += "\n(Uncertainties assume chi^2 residual, may be invalid.)"
            lines.append(l)
            lines.append(dashedline)

            varnames = self.varnames
            varvals = self.varvals
            varunc = self.varunc
            varlines = []

            # Format the lines
            w = max(map(len, varnames))
            w = str(w+1)
            formatstr = "%-"+w+"s " + pe + " +/- " + pe
            for name, val, unc in zip(varnames, varvals, varunc):
                varlines.append(formatstr%(name, val, unc))

            varlines.sort()
            lines.extend(varlines)

        ## The constraints
        if self.connames and self.showcon:
            lines.append("")
            l = "Constrained Parameters"
            if invalid:
                l += "\n(Uncertainties assume chi^2 residual, may be invalid.)"
            lines.append(l)
            lines.append(dashedline)
            conlines = []

            # Format the lines
            w = max(map(len, self.connames))
            w = str(w+1)
            formatstr = "%-"+w+"s " + pe + " +/- " + pe
            for name, val, unc in zip(self.connames, self.convals, self.conunc):
                conlines.append(formatstr%(name, val, unc))

            conlines.sort()
            lines.extend(conlines)

        # Fixed parameters
        if self.showfixed and self.fixednames:
            varlines = []
            lines.append("")
            lines.append("Fixed Parameters")
            lines.append(dashedline)
            # Format the lines
            w = max(map(len, self.fixednames))
            w = str(w+1)
            formatstr = "%-"+w+"s " + pe
            for name, val in zip(self.fixednames, self.fixedvals):
                # Skip over vector values for now. These are almost certainly
                # the data.
                if not numpy.isscalar(val): continue
                varlines.append(formatstr%(name, val))
            varlines.sort()
            lines.extend(varlines)


        ## Variable correlations
        lines.append("")
        corint = int(corrmin*100)
        l = "Variable Correlations greater than %i%%"%corint
        if invalid:
            l += "\n(Correlations assume chi^2 residual, may be invalid.)"
        lines.append(l)
        lines.append(dashedline)
        tup = []
        cornames = []
        n = len(self.varnames)
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
            lines.append(str(footer))

        out = "\n".join(lines) + '\n'
        return out

    def show(self, header = "", footer = "", update = False):
        """Format and print the results.

        header  --  A header to add to the output.
        footer  --  A footer to add to the output.
        update  --  Flag indicating whether to call update() (bool).

        """
        print self.format(header, footer, update).rstrip()
        return

    def __str__(self):
        return self.format()

    def save(self, filename, header = "", footer = "", update = False):
        """Format and save the results.

        filename -  Name of the save file.
        header  --  A header to add to the output.
        footer  --  A footer to add to the output.
        update  --  Flag indicating whether to call update() (bool).

        """
        # Save the time and user
        from time import ctime
        from getpass import getuser
        myheader = "Results written: " + ctime() + "\n"
        myheader += "produced by " + getuser() + "\n"
        header = myheader + header

        res = self.format(header, footer, update)
        f = file(filename, 'w')
        f.write(res)
        f.close()
        return

# End class FitResults

__id__ = "$Id$"
