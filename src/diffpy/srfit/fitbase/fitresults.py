#!/usr/bin/env python
##############################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2008 The Trustees of Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow, Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################
"""The FitResults and ContributionResults classes for storing results of
a fit.

The FitResults class is used to display the current state of a
FitRecipe. It stores the state, and uses it to calculate useful
statistics, which can be displayed on screen or saved to file.
"""

from __future__ import print_function

__all__ = ["FitResults", "ContributionResults", "initializeRecipe"]

import re
from collections import OrderedDict

import numpy
from diffpy.utils._deprecator import build_deprecation_message, deprecated

from diffpy.srfit.util import _DASHEDLINE
from diffpy.srfit.util import sortKeyForNumericString as numstr
from diffpy.srfit.util.inpututils import inputToString

fitresults_base = "diffpy.srfit.fitbase.FitResults"
removal_version = "4.0.0"

formatResults_dep_msg = build_deprecation_message(
    fitresults_base,
    "formatResults",
    "get_results_string",
    removal_version,
)

printResults_dep_msg = build_deprecation_message(
    fitresults_base,
    "printResults",
    "print_results",
    removal_version,
)

saveResults_dep_msg = build_deprecation_message(
    fitresults_base,
    "saveResults",
    "save_results",
    removal_version,
)

resultsDictionary_dep_msg = build_deprecation_message(
    "diffpy.srfit.fitbase",
    "resultsDictionary",
    "get_results_dictionary",
    removal_version,
    new_base="diffpy.srfit.fitbase.FitResults",
)

initializeRecipe_dep_msg = build_deprecation_message(
    "diffpy.srfit.fitbase",
    "initializeRecipe",
    "initialize_recipe_with_results",
    removal_version,
    new_base="diffpy.srfit.fitbase.FitRecipe",
)


class FitResults(object):
    """Class for processing, presenting and storing results of a fit.

    Attributes
    ----------
    recipe : FitRecipe
        The recipe from which the results were generated.

    cov : numpy.ndarray or None
        The covariance matrix of the refined variables. None if unavailable.

    conresults : collections.OrderedDict[str, ContributionResults]
        The ordered mapping of FitContribution name → ContributionResults.

    derivstep : float
        The fractional step size used for numerical derivatives (default 1e-8).

    varnames : list[str]
        The names of refined variables in the recipe.

    varvals : numpy.ndarray
        The optimized values of the refined variables.

    varunc : numpy.ndarray or None
        The estimated standard uncertainties of the variables. None if invalid.

    showfixed : bool
        Show the fixed variables in the formatted output
        (default True).

    fixednames : list[str]
        The names of variables held fixed during refinement.

    fixedvals : numpy.ndarray
        The values of the fixed variables.

    showcon : bool
        show the constrained parameters in the formatted output
        (default False).

    connames : list[str]
        The names of constrained parameters.

    convals : numpy.ndarray
        The values of constrained parameters.

    conunc : numpy.ndarray or None
        The uncertainties of constrained parameters. None if unavailable.

    residual : float
        The scalar residual value of the recipe.

    penalty : float
        The penalty contribution to the residual from restraints.

    chi2 : float
        The chi-squared value of the fit.

    cumchi2 : numpy.ndarray
        The cumulative chi-squared as a function of data index.

    rchi2 : float
        The reduced chi-squared of the fit.

    rw : float
        The weighted R-factor of the fit.

    cumrw : numpy.ndarray
        The cumulative weighted R-factor as a function of data index.

    messages : list[str]
        The informational or warning messages associated with the results.

    precision : int
        The number of digits used when formatting numeric output (default 8).

    _dcon : numpy.ndarray
        The jacobian of constraint equations with respect to variables.
        Used internally for uncertainty propagation.

    Each of these attributes, except the recipe, are created or updated when
    the update method is called.
    """

    def __init__(self, recipe, update=True, showfixed=True, showcon=False):
        """Initialize the attributes.

        Parameters
        ----------
        recipe : FitRecipe
            The recipe containing the results.
        update : bool
            The flag indicating whether to do an immediate update
            (default True).
        showfixed : bool
            Show fixed variables in the output (default True).
        showcon : bool
            Show constraint values in the output (default False).
        """
        self.recipe = recipe
        self.conresults = OrderedDict()
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
        self.residual = 0
        self.penalty = 0
        self.chi2 = 0
        self.rchi2 = 0
        self.rw = 0
        self.precision = 8
        self._dcon = []
        self.messages = []

        self.showfixed = bool(showfixed)
        self.showcon = bool(showcon)

        if update:
            self.update()
        return

    def update(self):
        """Update the results according to the current state of the
        recipe."""
        # Note that the order of these operations are chosen to reduce
        # computation time.

        recipe = self.recipe

        if not recipe._contributions:
            return

        # Make sure everything is ready for calculation
        recipe._prepare()

        # Store the variable names and values
        self.varnames = recipe.get_names()
        self.varvals = recipe.get_values()
        fixedpars = recipe._tagmanager.union(recipe._fixedtag)
        fixedpars = [p for p in fixedpars if not p.constrained]
        self.fixednames = [p.name for p in fixedpars]
        self.fixedvals = [p.value for p in fixedpars]

        # Store the constraint information
        self.connames = [con.par.name for con in recipe._oconstraints]
        self.convals = [con.par.getValue() for con in recipe._oconstraints]

        if self.varnames:
            # Calculate the covariance
            self._calculate_covariance()

            # Get the variable uncertainties
            self.varunc = [
                self.cov[i, i] ** 0.5 for i in range(len(self.varnames))
            ]

            # Get the constraint uncertainties
            self._calculate_constraint_uncertainties()

        # Store the fitting arrays and metrics for each FitContribution.
        self.conresults = OrderedDict()
        for con, weight in zip(
            recipe._contributions.values(), recipe._weights
        ):
            self.conresults[con.name] = ContributionResults(con, weight, self)

        # Calculate the metrics
        res = recipe.residual()
        self.residual = numpy.dot(res, res)
        self._calculate_metrics()

        # Calculate the restraints penalty
        w = self.chi2 / len(res)
        self.penalty = sum([r.penalty(w) for r in recipe._restraintlist])

        return

    def _calculate_covariance(self):
        """Calculate the covariance matrix. This is called by update.

        This code borrowed from PARK. It finds the pseudo-inverse of the
        Jacobian using the singular value decomposition.
        """
        try:
            J = self._calculate_jacobian()
            u, s, vh = numpy.linalg.svd(J, 0)
            self.cov = numpy.dot(vh.T.conj() / s**2, vh)
        except numpy.linalg.LinAlgError:
            self.messages.append("Cannot compute covariance matrix.")
            lvarvals = len(self.varvals)
            self.cov = numpy.zeros((lvarvals, lvarvals), dtype=float)
        return

    def _calculate_jacobian(self):
        """Calculate the Jacobian for the fitting.

        Adapted from PARK. Returns the derivative wrt the fit variables
        at point p.

        This also calculates the derivatives of the constrained
        parameters while we're at it.

        Numeric derivatives are calculated based on step, where step is
        the portion of variable value. E.g. step = dv/v.
        """
        recipe = self.recipe
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
        for k, v in enumerate(pvals):
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

            # FIXME - constraints are used for vectors as well!
            for i, con in enumerate(recipe._oconstraints):
                con.update()
                val = con.par.getValue()
                if numpy.isscalar(val):
                    cond[i] -= con.par.getValue()
                    cond[i] /= 2 * h
                else:
                    cond[i] = 0.0

            conr.append(cond)

            pvals[k] = v
            r.append(rk / (2 * h))

        # Reset the constrained parameters to their original values
        for con in recipe._oconstraints:
            con.update()

        self._dcon = numpy.vstack(conr).T

        # return the jacobian
        jac = numpy.vstack(r).T
        return jac

    def _calculate_metrics(self):
        """Calculate chi2, cumchi2, rchi2, rw and cumrw for the
        recipe."""
        cumchi2 = numpy.array([], dtype=float)
        # total weighed denominator for the ratio in the Rw formula
        yw2tot = 0.0
        numpoints = 0
        for con in self.conresults.values():
            cc2w = con.weight * con.cumchi2
            c2last = cumchi2[-1:].sum()
            cumchi2 = numpy.concatenate([cumchi2, c2last + cc2w])
            yw2tot += con.weight * (con.chi2 / con.rw**2)
            numpoints += len(con.x)

        chi2 = cumchi2[-1:].sum()
        cumrw = numpy.sqrt(cumchi2 / yw2tot)
        rw = cumrw[-1:].sum()

        numpoints += len(self.recipe._restraintlist)
        rchi2 = chi2 / (numpoints - len(self.varnames))

        self.chi2 = chi2
        self.rchi2 = rchi2
        self.rw = rw
        self.cumchi2 = cumchi2
        self.cumrw = cumrw
        return

    def _calculate_constraint_uncertainties(self):
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

    def get_results_string(self, header="", footer="", update=False):
        """Format the results and return them in a string.

        This function is called by print_results and save_results. Overloading
        the formatting here will change all three functions.

        Parameters
        ----------
        header : str
            The header to add to the output (default "")
        footer : str
            The footer to add to the output (default "")
        update : bool
            The flag indicating whether to call update() (default False).

        Returns
        -------
        out : str
            a string containing the formatted results.
        """
        if update:
            self.update()

        lines = []
        corrmin = 0.25
        p = self.precision
        pe = "%-" + "%i.%ie" % (p + 6, p)
        pet = "%" + ".%ie" % (p,)
        # Check to see if the uncertainty values are reliable.
        certain = True
        for con in self.conresults.values():
            if (con.dy == 1).all():
                certain = False
                break

        # User-defined header
        if header:
            lines.append(header)

        if not certain:
            err_msg = (
                "Some quantities invalid due to missing profile uncertainty"
            )
            if err_msg not in self.messages:
                self.messages.append(err_msg)

        lines.extend(self.messages)

        # Overall results
        err_msg = "Overall"
        if not certain:
            err_msg += " (Chi2 and Reduced Chi2 invalid)"
        lines.append(err_msg)
        lines.append(_DASHEDLINE)
        formatstr = "%-14s %.8f"
        lines.append(formatstr % ("Residual", self.residual))
        lines.append(
            formatstr % ("Contributions", self.residual - self.penalty)
        )
        lines.append(formatstr % ("Restraints", self.penalty))
        lines.append(formatstr % ("Chi2", self.chi2))
        lines.append(formatstr % ("Reduced Chi2", self.rchi2))
        lines.append(formatstr % ("Rw", self.rw))

        # Per-FitContribution results
        if len(self.conresults) > 1:
            keys = list(self.conresults.keys())
            keys.sort(key=numstr)

            lines.append("")
            err_msg = "Contributions"
            if not certain:
                err_msg += " (Chi2 and Reduced Chi2 invalid)"
            lines.append(err_msg)
            lines.append(_DASHEDLINE)
            formatstr = "%-10s %-42.8f"
            for name in keys:
                res = self.conresults[name]
                lines.append("")
                namestr = name + " (%f)" % res.weight
                lines.append(namestr)
                lines.append("-" * len(namestr))
                lines.append(formatstr % ("Residual", res.residual))
                lines.append(formatstr % ("Chi2", res.chi2))
                lines.append(formatstr % ("Rw", res.rw))

        # The variables
        if self.varnames:
            lines.append("")
            err_msg = "Variables"
            if not certain:
                err_msg2 = "Uncertainties invalid"
                err_msg += " (%s)" % err_msg2
            lines.append(err_msg)
            lines.append(_DASHEDLINE)

            varnames = self.varnames
            varvals = self.varvals
            varunc = self.varunc
            varlines = []

            w = max(map(len, varnames))
            w = str(w + 1)
            # Format the lines
            formatstr = "%-" + w + "s " + pe + " +/- " + pet
            for name, val, unc in zip(varnames, varvals, varunc):
                varlines.append(formatstr % (name, val, unc))

            varlines.sort()
            lines.extend(varlines)

        # Fixed variables
        if self.showfixed and self.fixednames:
            varlines = []
            lines.append("")
            lines.append("Fixed Variables")
            lines.append(_DASHEDLINE)
            w = max(map(len, self.fixednames))
            w = str(w + 1)
            formatstr = "%-" + w + "s " + pet
            for name, val in zip(self.fixednames, self.fixedvals):
                varlines.append(formatstr % (name, val))
            varlines.sort()
            lines.extend(varlines)

        # The constraints
        if self.connames and self.showcon:
            lines.append("")
            err_msg = "Constrained Parameters"
            if not certain:
                err_msg += " (Uncertainties invalid)"
            lines.append(err_msg)
            lines.append(_DASHEDLINE)

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

            keys.sort(key=numstr)
            w = str(w + 1)
            formatstr = "%-" + w + "s %- 15f +/- %-15f"
            for name in keys:
                val, unc = vals[name]
                lines.append(formatstr % (name, val, unc))

        # Variable correlations
        lines.append("")
        corint = int(corrmin * 100)
        err_msg = "Variable Correlations greater than %i%%" % corint
        if not certain:
            err_msg += " (Correlations invalid)"
        lines.append(err_msg)
        lines.append(_DASHEDLINE)
        tup = []
        cornames = []
        n = len(self.varnames)
        for i in range(n):
            for j in range(i + 1, n):
                name = "corr(%s, %s)" % (varnames[i], varnames[j])
                val = self.cov[i, j] / (self.cov[i, i] * self.cov[j, j]) ** 0.5
                if abs(val) > corrmin:
                    cornames.append(name)
                    tup.append((val, name))

        tup.sort(key=lambda vn: abs(vn[0]))
        tup.reverse()

        if cornames:
            w = max(map(len, cornames))
            w = str(w + 1)
            formatstr = "%-" + w + "s  %.4f"
            for val, name in tup:
                lines.append(formatstr % (name, val))
        else:
            lines.append("No correlations greater than %i%%" % corint)

        # User-defined footer
        if footer:
            lines.append(footer)

        out = "\n".join(lines) + "\n"
        return out

    @deprecated(formatResults_dep_msg)
    def formatResults(self, header="", footer="", update=False):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use diffpy.srfit.fitbase.FitResults.get_results_string
        instead.
        """
        return self.get_results_string(header, footer, update)

    def print_results(self, header="", footer="", update=False):
        """Format and print the results.

        Parameters
        ----------
        header
            The header to add to the output (default "")
        footer
            The footer to add to the output (default "")
        update
            The flag indicating whether to call update() (default False).
        """
        print(self.get_results_string(header, footer, update).rstrip())
        return

    @deprecated(printResults_dep_msg)
    def printResults(self, header="", footer="", update=False):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use diffpy.srfit.fitbase.FitResults.print_results
        instead.
        """
        self.print_results(header, footer, update)
        return

    def __str__(self):
        return self.get_results_string()

    def save_results(self, filename, header="", footer="", update=False):
        """Format and save the results.

        Parameters
        ----------------------------------
        filename
            The name of the save file.
        header
            The header to add to the output (default "")
        footer
            The footer to add to the output (default "")
        update
            The flag indicating whether to call update() (default False).
        """
        # Save the time and user
        from getpass import getuser
        from time import ctime

        myheader = "Results written: " + ctime() + "\n"
        myheader += "produced by " + getuser() + "\n"
        header = myheader + header

        res = self.get_results_string(header, footer, update)
        f = open(filename, "w")
        f.write(res)
        f.close()
        return

    @deprecated(saveResults_dep_msg)
    def saveResults(self, filename, header="", footer="", update=False):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use diffpy.srfit.fitbase.FitResults.save_results
        instead.
        """
        self.save_results(filename, header, footer, update)
        return

    def get_results_dictionary(self):
        """Get a dictionary of results, with variable names and values,
        and overall metrics.

        Returns
        -------
        results_dict : dict
            The dictionary containing the variable names and values,
            and overall metrics, from the FitResults.
        """
        parameter_names = self.varnames
        parameter_values = self.varvals
        results_dict = dict(zip(parameter_names, parameter_values))
        results_dict.update(
            {
                "Residual": self.residual,
                "Contributions": self.residual - self.penalty,
                "Restraints": self.penalty,
                "Chi2": self.chi2,
                "Reduced Chi2": self.rchi2,
                "Rw": self.rw,
            }
        )
        return results_dict


# End class FitResults


class ContributionResults(object):
    """Class for processing, storing FitContribution results.

    This does not store the FitContribution.

    Attributes
    ----------
    y : numpy.ndarray or None
        The FitContribution's profile over the calculation range
        (default None).
    dy : numpy.ndarray or None
        The uncertainty in the FitContribution's profile over the
        calculation range (default None).
    x : numpy.ndarray or None
        The numpy array of the calculated independent variable for the
        FitContribution (default None).
    ycalc : numpy.ndarray or None
        The numpy array of the calculated signal for the FitContribution
        (default None).
    residual : float
        The scalar residual of the FitContribution.
    chi2 : float
        The chi2 of the FitContribution.
    cumchi2 : numpy.ndarray
        The cumulative chi2 of the FitContribution.
    rw : float
        The Rw of the FitContribution.
    cumrw : numpy.ndarray
        The cumulative Rw of the FitContribution.
    weight : float
        The weight of the FitContribution in the recipe.
    conlocs : list
        The location of the constrained parameters in the
        FitContribution (see the
        RecipeContainer._locate_managed_object method).
    convals : list
        The values of the constrained parameters.
    conunc : list
        The uncertainties in the constraint values.
    """

    def __init__(self, con, weight, fitres):
        """Initialize the attributes.

        Parameters
        ----------
        con
            The FitContribution
        weight
            The weight of the FitContribution in the recipe
        fitres
            The FitResults instance to contain this ContributionResults
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
        # Note that the order of these operations is chosen to reduce
        # computation time.

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
        self._calculate_metrics()

        # Find the parameters
        for i, constraint in enumerate(recipe._oconstraints):
            par = constraint.par
            loc = con._locate_managed_object(par)
            if loc:
                self.conlocs.append(loc)
                self.convals.append(fitres.convals[i])
                self.conunc.append(fitres.conunc[i])

        return

    # FIXME: factor rw, chi2, cumrw, cumchi2 to separate functions.
    def _calculate_metrics(self):
        """Calculate chi2 and Rw of the recipe."""
        # We take absolute values in case the signal is complex
        num = numpy.abs(self.y - self.ycalc)
        y = numpy.abs(self.y)
        chiv = num / self.dy
        self.cumchi2 = numpy.cumsum(chiv**2)
        # avoid index error for empty array
        self.chi2 = self.cumchi2[-1:].sum()
        yw = y / self.dy
        yw2tot = numpy.dot(yw, yw)
        if yw2tot == 0.0:
            yw2tot = 1.0
        self.cumrw = numpy.sqrt(self.cumchi2 / yw2tot)
        # avoid index error for empty array
        self.rw = self.cumrw[-1:].sum()
        return


# End class ContributionResults


@deprecated(resultsDictionary_dep_msg)
def resultsDictionary(results):
    """**This function has been deprecated and will be** **removed in version
    4.0.0.**

    **Please use**
    **diffpy.srfit.fitbase.FitResults.get_results_dictionary instead.**

    Get dictionary of results from file.

    This reads the results from file and stores them in a dictionary to be
    returned to the caller. The dictionary may contain non-result entries.

    Parameters
    ----------
    results
        An open file-like object, name of a file that contains
        results from FitResults or a string containing fit results.
    """
    resstr = inputToString(results)

    rx = {
        "f": r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?",
        "n": r"[a-zA-Z_]\w*",
    }
    pat = r"(%(n)s)\s+(%(f)s)" % rx

    matches = re.findall(pat, resstr)
    # Prefer the first match
    matches.reverse()
    mpairs = dict(matches)
    return mpairs


@deprecated(initializeRecipe_dep_msg)
def initializeRecipe(recipe, results):
    """**This function has been deprecated and will be** **removed in
    version 4.0.0.**

    **Please use**
    **diffpy.srfit.fitbase.FitRecipe.initialize_recipe_with_results**
    **instead.**

    Initialize the variables of a recipe from a results file.

    This reads the results from file and initializes any variables (fixed or
    free) in the recipe to the results values. Note that the recipe has to be
    configured, with variables. This does not reconstruct a FitRecipe.

    Parameters
    ----------
    recipe
        A configured recipe with variables
    results
        An open file-like object, name of a file that contains
        results from FitResults or a string containing fit results.
    """

    mpairs = resultsDictionary(results)
    if not mpairs:
        raise AttributeError("Cannot find results")

    # Get variable names
    names = recipe._parameters.keys()
    for vname in names:
        value = mpairs.get(vname)
        if value is not None:
            var = recipe.get(vname)
            var.value = float(value)

    return
