#!/usr/bin/env python

from itertools import imap

import numpy

from diffpy.srfit.adapters import adapt
from diffpy.srfit.util import getVaried, makeArray
from diffpy.srfit.util import isContained
from diffpy.srfit.util.cachemanager import CacheManager
from diffpy.srfit.fit import functions

__all__ = ["chi", "Rw", "reschi2", "resRw2", "residual"]

# FIXME - make chi and Rw classes so they can be identified in FitResults. This
# will allow us to calculate the uncertainties in FitResults even if the user
# is not using chi^2 to optimize.
def chi(fiteq, y, dy = None):
    """Make vector chi function out of a fit equation and profile.

    fiteq   --  The equation or value to be fit
    y       --  y values or parameter containing this (observed profile).
    dy      --  dy values or parameter containing this (uncertainty in observed
                profile). If this is None (default), it will not be used.

    Returns the vector chi function:
    chiv = (fiteq - y) / dy

    This will be used in a Residual to calculate chi2
    chi2 = dot(chiv, chiv)

    """
    # Make sure we're dealing with adapted objects
    fiteq = adapt(fiteq)
    y = adapt(y)
    if dy is not None:
        dy = adapt(dy)

    chiv = (fiteq - y)
    if dy is not None:
        chiv /= dy
    chiv.name = "chi^2"
    return chiv

def Rw(fiteq, y, w = None):
    """Make vector Rw function out of equation and profile.

    fiteq   --  The equation or value to be fit
    y       --  y values or parameter containing this (observed profile).
    dy      --  The weighting factor or parameter containing this.  If this is
                None (default), it will not be used.

    This returns the vector Rw function:
    rwv = (fiteq - y) * sqrt(w) / sqrt( dot(w * y, y) )

    This will be used in a residual to calculate Rw^2
    Rw^2 = dot(rwv, rwv)

    """
    # Make sure we're dealing with adapted objects
    fiteq = adapt(fiteq)
    y = adapt(y)
    if w is not None:
        w = adapt(w)

    rwv = fiteq - y
    if w is None:
        rwv /= functions.sqrt_( functions.dot_(y, y) )
    else:
        rwv *= functions.sqrt_( functions.abs_(w) )
        rwv /= functions.sqrt_( functions.dot_(w * y, y) )
    rwv.name = "Rw^2"
    return rwv

def reschi2(fiteq, y, dy = None):
    """Make the standard chi^2 residual.
    
    fiteq   --  The equation to be fit
    y       --  Parameter containing y values (observed profile).
    dy      --  Parameter containing dy values (uncertainty in observed
                profile).
    
    """
    reseq = chi(fiteq, y, dy)
    return residual(reseq)

def resRw2(fiteq, y, w = None):
    """Make the Rw^2 residual.
    
    fiteq   --  The equation to be fit
    y       --  Parameter containing y values (observed profile).
    w       --  The weighting factor. If this is None (default), it is not
                used.
    
    """
    reseq = Rw(fiteq, y, w)
    return residual(reseq)

class Residual(object):
    """Class for making residual equation.

    The 'vec' method returns the evaluated vector residual for use with vector
    optimizers, such as leastsq from scipy.optimize. The residual is a
    functor that computes the scalar residual for use with scalar optimizers,
    such as fmin from scipy.optimize.

    Attributes
    terms       --  Terms of the residual, such equations returned by 'chi',
                    'Rw' or 'restrain'.
    variables   --  Variable parameters extracted from the terms. This is
                    updated whenever the variables status of parameters
                    changes, and only when needed.

    """

    def __init__(self, *terms):
        """Initialize the residual.

        terms   --  Terms in the residual. The output of these equations will
                    be concatenated to form a vector residual, resv. This can
                    be accessed with the 'vec' method. The scalar residual,
                    dot(resv, resv), is accessible through '__call__'.
        
        """
        self._cache = CacheManager()
        self._cache.addNode(self)
        self.terms = terms
        self._networkTerms()
        # Get the fitted variables and the equation
        self._variables = None
        self.fithooks = []
        from diffpy.srfit.fit.fithook import PrintFitHook
        self.fithooks.append(PrintFitHook())
        return

    variables = property( lambda self: self._getVariables() )
    def _getVariables(self):
        """Get variables.

        If the cached variables are None, then this will scan the residual
        terms for variables and cache them. The 'reset' method of the fit hooks
        is called when the variables are recached.
        
        """
        if self._variables is not None:
            return self._variables

        # Get variables. Note that these are sorted according to name.
        varset = set()
        varsets = imap(getVaried, self.terms)
        map(varset.update, varsets)
        self._variables = sorted(varset, key = lambda v: v.name)

        # Set up the fit hooks
        for hook in self.fithooks:
            hook.reset(self)

        return self._variables


    names = property( lambda self: self._getNames() )
    def _getNames(self):
        """Get names of variables."""
        return [v.name for v in self.variables]

    values = property( lambda self: self._getValues() )
    def _getValues(self):
        """Get values of variables."""
        return [v.get() for v in self.variables]

    def _networkTerms(self):
        """Add us to the network of each term."""
        self.terms = map(adapt, self.terms)
        for term in self.terms:
            term._cache.addNode(self)
        return

    def append(self, term):
        """Add a term to the residual.

        term    --  An equation that can be computed to obtain part of the
                    vector residual.

        """
        adterm = adapt(term)
        self.terms.append(adterm)
        adterm._cache.addNode(self)
        self._variables = None
        return

    def remove(self, term):
        """Remove a term from the residual.

        term    --  An equation that can be computed to obtain part of the
                    vector residual.

        """
        self.terms.remove(term)
        term._cache.removeNode(self)
        self._variables = None
        return

    def _onVary(self):
        """Reset variables if cache manager says they have changed."""
        self._variables = None
        return

    def _setvals(self, p):
        """Update the values of all parameters."""
        for var, val in zip(self.variables, p):
            var.set(val)
        return

    def _getresv(self):
        """Get the residual vector from the residual terms."""
        vals = (term.get() for term in self.terms)
        avals = map(makeArray, vals)
        resv = numpy.concatenate(avals)
        return resv

    def vec(self, p = None):
        """Vector residual."""
        if p is not None: 
            self._setvals(p)
        for hook in self.fithooks:
            hook.precall(self)
        retval = self._getresv()
        for hook in self.fithooks:
            hook.postcall(self, retval)
        return retval

    def __call__(self, p = None):
        """Scalar residual."""
        resv = self.vec(p)
        retval = numpy.dot(resv, resv)
        return retval

# end class Residual

def residual(*terms):
    """Make a residual out of an equation.

    terms   --  Terms in the residual. The output of these equations will
                be concatenated to form a vector residual, resv.
    
    The equation should be in vector form, e.g. the output returned by chi.

    A residual has various other methods that aid in evaluation. See the
    examples for use.

    """
    res = Residual(*terms)
    return res

