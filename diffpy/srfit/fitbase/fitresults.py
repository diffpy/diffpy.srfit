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
"""The FitResults and ContributionResult classes for storing results of a fit.

"""

import numpy

class FitResults(object):
    """Class for processing, presenting and storing results of a fit. 

    
    Attributes
    model       --  The model containing the results.
    cov         --  The covariance matrix from the model.
    conresults  --  A dictionary of ContributionResult for each contribution,
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

        # Store the fitting arrays and metrics for each contribution.
        self.conresults = {}
        for con, weight in zip(model._organizers, model._weights):
            self.conresults[con.name] = ContributionResult(con, weight)

        # Store the variable names and values
        self.varnames = model.getNames()
        self.varvals = model.getValues()

        # Calculate the metrics
        self.residual = model.scalarResidual()
        self._calculateMetrics()

        # Store the constraint names and values
        self.connames = [con.par.name for con in model._constraintlist]
        self.convals = [con.par.getValue() for con in model._constraintlist]

        # Calculate the covariance
        self._calculateCovariance()

        # Get the variable uncertainties
        self.varunc = [self.cov[i,i]**0.5 for i in range(len(self.varnames))]

        # Get the constraint uncertainties
        self._calculateConstraintUncertainties()
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

        # dp = sum_i dp/dvi * unc(vi)
        v = numpy.array(self.varunc)
        dp = abs(numpy.dot(self._conjac, v))
        self.conunc = dp.tolist()
        return

    def formatResult(self):
        """"""
        return

    def printResult(self):
        """"""
        return

# End class FitResults

class ContributionResult(object):
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

    """

    def __init__(self, con, weight):
        """Initialize the attributes."""
        self.x = None
        self.y = None
        self.dy = None
        self.ycalc = None
        self.residual = 0
        self.chi2 = 0
        self.rw = 0
        self.weight = 0
        self._init(con, weight)
        return

    def _init(self, con, weight):
        """Initialize the attributes, for real."""
        ## Note that the order of these operations are chosen to reduce
        ## computation time.

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


__id__ = "$Id$"
