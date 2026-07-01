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
"""FitRecipe class.

FitRecipes organize FitContributions, variables, Restraints and
Constraints to create a recipe of the system you wish to optimize. From
the client's perspective, the FitRecipe is a residual calculator. The
residual method does the work of updating variable values, which get
propagated to the Parameters of the underlying FitContributions via the
variables and Constraints.  This class needs no special knowledge of the
type of FitContribution or data being used. Thus, it is suitable for
combining residual equations from various types of refinements into a
single residual.

Variables added to a FitRecipe can be tagged with string identifiers.
Variables can be later retrieved or manipulated by tag. The tag name
"__fixed" is reserved.

See the examples in the documentation for how to create an optimization
problem using FitRecipe.
"""

__all__ = ["FitRecipe"]

from collections import OrderedDict
from pathlib import Path

import matplotlib.pyplot as plt
from bg_mpl_stylesheets.styles import all_styles
from numpy import array, concatenate, dot, sqrt

import diffpy.srfit.util.inpututils as utils
from diffpy.srfit.fitbase.fithook import PrintFitHook
from diffpy.srfit.fitbase.parameter import ParameterProxy
from diffpy.srfit.fitbase.recipeorganizer import RecipeOrganizer
from diffpy.srfit.interface import _fitrecipe_interface
from diffpy.srfit.util.tagmanager import TagManager
from diffpy.utils._deprecator import build_deprecation_message, deprecated

plt.style.use(all_styles["bg-style"])

base = "diffpy.srfit.fitbase.FitRecipe"
removal_version = "4.0.0"
addcontrib_dep_msg = build_deprecation_message(
    base, "addContribution", "add_contribution", removal_version
)

pushfithook_dep_msg = build_deprecation_message(
    base, "pushFitHook", "push_fit_hook", removal_version
)

popfithook_dep_msg = build_deprecation_message(
    base, "popFitHook", "pop_fit_hook", removal_version
)

getfithooks_dep_msg = build_deprecation_message(
    base, "getFitHooks", "get_fit_hooks", removal_version
)

clearfithooks_dep_msg = build_deprecation_message(
    base, "clearFitHooks", "clear_fit_hooks", removal_version
)

setweight_dep_msg = build_deprecation_message(
    base, "setWeight", "set_weight", removal_version
)

addparset_dep_msg = build_deprecation_message(
    base, "addParameterSet", "add_parameter_set", removal_version
)

removeParameterSet_dep_msg = build_deprecation_message(
    base, "removeParameterSet", "remove_parameter_set", removal_version
)

scalarResidual_dep_msg = build_deprecation_message(
    base, "scalarResidual", "scalar_residual", removal_version
)

addVar_dep_msg = build_deprecation_message(
    base, "addVar", "add_variable", removal_version
)

delVar_dep_msg = build_deprecation_message(
    base, "delVar", "delete_variable", removal_version
)

newVar_dep_msg = build_deprecation_message(
    base, "newVar", "create_new_variable", removal_version
)

isFree_dep_msg = build_deprecation_message(
    base, "isFree", "is_free", removal_version
)

getValues_dep_msg = build_deprecation_message(
    base, "getValues", "get_values", removal_version
)

getNames_dep_msg = build_deprecation_message(
    base, "getNames", "get_names", removal_version
)

getBounds_dep_msg = build_deprecation_message(
    base, "getBounds", "get_bounds_pairs", removal_version
)

getBounds2_dep_msg = build_deprecation_message(
    base, "getBounds2", "get_bounds_array", removal_version
)

boundsToRestraints_dep_msg = build_deprecation_message(
    base, "boundsToRestraints", "convert_bounds_to_restraints", removal_version
)

constrain_dep_msg = build_deprecation_message(
    base, "constrain", "add_constraint", removal_version
)

unconstrain_dep_msg = build_deprecation_message(
    base, "unconstrain", "remove_constraint", removal_version
)


class FitRecipe(_fitrecipe_interface, RecipeOrganizer):
    """FitRecipe class.

    Attributes
    ----------
    name : str
        A name for this FitRecipe.
    fithooks : list
        The list of FitHook instances that can pass information out
        of the system during a refinement. By default, this is
        populated by a PrintFitHook instance.
    _constraints : dict
        The dictionary of Constraints, indexed by the constrained
        Parameter. Constraints can be added using the
        'constrain' method.
    _oconstraints : list
        The ordered list of the constraints from this and all
        sub-components.
    _calculators : dict
        The managed dictionary of Calculators.
    _contributions : OrderedDict
        The managed OrderedDict of FitContributions.
    _parameters : OrderedDict
        The managed OrderedDict of parameters (in this case the
        parameters are varied).
    _parsets : dict
        The managed dictionary of ParameterSets.
    _eqfactory : diffpy.srfit.equation.builder.EquationFactory
        The diffpy.srfit.equation.builder.EquationFactory
        instance that is used to create constraints and
        restraints from strings.
    _restraintlist : list
        The list of restraints from this and all sub-components.
    _restraints : set
        The set of Restraints. Restraints can be added using the
        'restrain' or 'confine' methods.
    _ready : bool
        The flag indicating if all attributes are ready for the
        calculation.
    _tagmanager : TagManager
        The TagManager instance for managing tags on Parameters.
    _weights : list
        The list of weighing factors for each FitContribution. The
        weights are multiplied by the residual of the
        FitContribution when determining the overall residual.
    _fixedtag : str
        "__fixed", used for tagging variables as fixed. Don't
        use this tag unless you want issues.

    Properties
    ----------
    names : list
        The variable names (read only). See get_names.
    values : numpy.ndarray
        The variable values (read only). See get_values.
    fixednames : list
        The names of the fixed refinable variables (read only).
    fixedvalues : numpy.ndarray
        The values of the fixed refinable variables (read only).
    bounds : list of tuple
        The bounds on parameters (read only). See get_bounds_pairs.
    bounds2 : tuple of numpy.ndarray
        The bounds on parameters (read only). See get_bounds_array.
    """

    fixednames = property(
        lambda self: [
            v.name
            for v in self._parameters.values()
            if not (self.is_free(v) or self.is_constrained(v))
        ],
        doc="names of the fixed refinable variables",
    )
    fixedvalues = property(
        lambda self: array(
            [
                v.value
                for v in self._parameters.values()
                if not (self.is_free(v) or self.is_constrained(v))
            ]
        ),
        doc="values of the fixed refinable variables",
    )
    bounds = property(lambda self: self.get_bounds_pairs())
    bounds2 = property(lambda self: self.get_bounds_array())

    def __init__(self, name="fit"):
        """Initialization."""
        RecipeOrganizer.__init__(self, name)
        self.fithooks = []
        self.pushFitHook(PrintFitHook())
        self._restraintlist = []
        self._oconstraints = []
        self._ready = False
        self._fixedtag = "__fixed"

        self._weights = []
        self._tagmanager = TagManager()

        self._parsets = {}
        self._manage(self._parsets)

        self._contributions = OrderedDict()
        self._manage(self._contributions)

        self.plot_options = {
            "show_observed": True,
            "show_fit": True,
            "show_diff": True,
            "offset_scale": 1.0,
            "xmin": None,
            "xmax": None,
            "figsize": (8, 6),
            "data_style": "o",
            "fit_style": "-",
            "diff_style": "-",
            "data_color": None,
            "fit_color": None,
            "diff_color": None,
            "data_label": "Observed",
            "fit_label": "Calculated",
            "diff_label": "Difference",
            "xlabel": None,
            "ylabel": None,
            "title": None,
            "legend": True,
            "legend_loc": "best",
            "grid": False,
            "markersize": None,
            "linewidth": None,
            "alpha": 1.0,
            "show": True,
        }
        return

    def push_fit_hook(self, fithook, index=None):
        """Add a FitHook to be called within the residual method.

        The hook is an object for reporting updates, or more fundamentally,
        passing information out of the system during a refinement. See the
        diffpy.srfit.fitbase.fithook.FitHook class for the required interface.
        Added FitHooks will be called sequentially during refinement.

        Parameters
        ----------
        fithook : diffpy.srfit.fitbase.fithook.FitHook
            The FitHook instance to add to the sequence.
        index : int or None, optional
            The index for inserting fithook into the list of fit hooks. If
            this is None (default), the fithook is added to the end.
        """
        if index is None:
            index = len(self.fithooks)
        self.fithooks.insert(index, fithook)
        # Make sure the added FitHook gets its reset method called.
        self._update_configuration()
        return

    @deprecated(pushfithook_dep_msg)
    def pushFitHook(self, fithook, index=None):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use diffpy.srfit.fitbase.FitRecipe.push_fit_hook instead.
        """
        self.push_fit_hook(fithook, index)
        return

    def pop_fit_hook(self, fithook=None, index=-1):
        """Remove a FitHook by index or reference.

        Parameters
        ----------
        fithook : diffpy.srfit.fitbase.fithook.FitHook or None, optional
            The FitHook instance to remove from the sequence. If this is
            None (default), default to index.
        index : int, optional
            The index of FitHook instance to remove (default -1).

        Raises
        ------
        ValueError
            If fithook is not None, but is not present in the sequence.
        IndexError
            If the sequence is empty or index is out of range.
        """
        if fithook is not None:
            self.fithooks.remove(fithook)
            return
        self.fithook.remove(index)
        return

    @deprecated(popfithook_dep_msg)
    def popFitHook(self, fithook=None, index=-1):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use diffpy.srfit.fitbase.FitRecipe.pop_fit_hook instead.
        """
        self.pop_fit_hook(fithook, index)
        return

    def get_fit_hooks(self):
        """Get the sequence of FitHook instances."""
        return self.fithooks[:]

    @deprecated(getfithooks_dep_msg)
    def getFitHooks(self):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use diffpy.srfit.fitbase.FitRecipe.get_fit_hooks instead."""
        return self.get_fit_hooks()

    def clear_fit_hooks(self):
        """Clear the FitHook sequence."""
        del self.fithooks[:]
        return

    @deprecated(clearfithooks_dep_msg)
    def clearFitHooks(self):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use diffpy.srfit.fitbase.FitRecipe.clear_fit_hooks instead."""
        self.clear_fit_hooks()
        return

    def add_contribution(self, con, weight=1.0):
        """Add a FitContribution to the FitRecipe.

        Parameters
        ----------
        con : FitContribution
            The FitContribution to be stored.
        weight : float, optional
            The weight of the FitContribution. Default is 1.0.

        Raises
        ------
        ValueError
            If the FitContribution has no name or if the FitContribution has
            the same name as some other managed object.
        """
        self._add_object(con, self._contributions, True)
        self._weights.append(weight)
        return

    @deprecated(addcontrib_dep_msg)
    def addContribution(self, con, weight=1.0):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use diffpy.srfit.fitbase.FitRecipe.add_contribution
        instead.
        """
        self.add_contribution(con, weight)
        return

    def set_weight(self, con, weight):
        """Set the weight of a FitContribution.

        Parameters
        ----------
        con : FitContribution
            The FitContribution object whose weight is to be set.
        weight : float
            The weight value to assign to the specified FitContribution.

        Returns
        -------
        None
        """
        idx = list(self._contributions.values()).index(con)
        self._weights[idx] = weight
        return

    @deprecated(setweight_dep_msg)
    def setWeight(self, con, weight):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use diffpy.srfit.fitbase.FitRecipe.set_weight instead."""
        self.set_weight(con, weight)
        return

    def add_parameter_set(self, parset):
        """Add a ParameterSet to the hierarchy.

        Parameters
        ----------
        parset : ParameterSet
            The ParameterSet to be stored.

        Raises
        ------
        ValueError
            If the ParameterSet has no name or if the ParameterSet has the same
            name as some other managed object.
        """
        self._add_object(parset, self._parsets, True)
        return

    @deprecated(addparset_dep_msg)
    def addParameterSet(self, parset):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use diffpy.srfit.fitbase.FitRecipe.add_parameter_set instead.
        """
        self.add_parameter_set(parset)
        return

    def remove_parameter_set(self, parset):
        """Remove a ParameterSet from the hierarchy.

        This method removes the specified ParameterSet object from the internal
        hierarchy of managed ParameterSets. If the provided ParameterSet is not
        currently managed by this object, a ValueError will be raised.

        Parameters:
        -----------
        parset : ParameterSet
            The ParameterSet instance to be removed from the hierarchy.

        Raises:
        -------
        ValueError
            If the provided ParameterSet is not managed by this object.
        """
        self._remove_object(parset, self._parsets)
        return

    @deprecated(removeParameterSet_dep_msg)
    def removeParameterSet(self, parset):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use diffpy.srfit.fitbase.FitRecipe.remove_parameter_set instead.
        """
        self.remove_parameter_set(parset)
        return

    def residual(self, p=[]):
        """Calculate the vector residual to be optimized.

        The residual is by default the weighted concatenation of each
        FitContribution's residual, plus the value of each restraint. The array
        returned, denoted chiv, is such that
        dot(chiv, chiv) = chi^2 + restraints.

        Parameters
        ----------
        p : list or numpy.ndarray
            The list of current variable values, provided in the same order
            as the '_parameters' list. If p is an empty iterable (default),
            then it is assumed that the parameters have already been
            updated in some other way, and the explicit update within this
            function is skipped.

        Return
        ------
        chiv : numpy.ndarray
            The array of residuals to be optimized. The array is such that
            dot(chiv, chiv) = chi^2 + restraints.
        """

        # Prepare, if necessary
        self._prepare()

        for fithook in self.fithooks:
            fithook.precall(self)

        # Update the variable parameters.
        self._apply_values(p)

        # Update the constraints. These are ordered such that the list only
        # needs to be cycled once.
        for con in self._oconstraints:
            con.update()

        # Calculate the bare chiv
        chiv = concatenate(
            [
                wi * ci.residual().flatten()
                for wi, ci in zip(self._weights, self._contributions.values())
            ]
        )

        # Calculate the point-average chi^2
        w = dot(chiv, chiv) / len(chiv)
        # Now we must append the restraints
        penalties = [sqrt(res.penalty(w)) for res in self._restraintlist]
        chiv = concatenate([chiv, penalties])

        for fithook in self.fithooks:
            fithook.postcall(self, chiv)

        return chiv

    def scalar_residual(self, p=[]):
        """Calculate the scalar residual to be optimized.

        Parameters
        ----------
        p : list or numpy.ndarray
            The list of current variable values, provided in the same order
            as the '_parameters' list. If p is an empty iterable (default),
            then it is assumed that the parameters have already been
            updated in some other way, and the explicit update within this
            function is skipped.

        The residual is by default the weighted concatenation of each
        FitContribution's residual, plus the value of each restraint. The array
        returned, denoted chiv, is such that
        dot(chiv, chiv) = chi^2 + restraints.
        """
        chiv = self.residual(p)
        return dot(chiv, chiv)

    @deprecated(scalarResidual_dep_msg)
    def scalarResidual(self, p=[]):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use diffpy.srfit.fitbase.FitRecipe.scalar_residual
        instead.
        """
        return self.scalar_residual(p)

    def __call__(self, p=[]):
        """Same as scalar_residual method."""
        return self.scalar_residual(p)

    def _prepare(self):
        """Prepare for the residual calculation, if necessary.

        This will prepare the data attributes to be used in the residual
        calculation.

        This updates the local restraints with those of the
        contributions.

        Raises
        ------
        AttributeError
            If there are variables without a value.
        """
        # Only prepare if the configuration has changed within the recipe
        # hierarchy.
        if self._ready:
            return

        # Inform the fit hooks that we're updating things
        for fithook in self.fithooks:
            fithook.reset(self)

        # Check Profiles
        self.__verify_profiles()

        # Check parameters
        self.__verify_parameters()

        # Update constraints and restraints.
        self.__collect_constraints_and_restraints()

        # We do this here so that the calculations that take place during the
        # validation use the most current values of the parameters. In most
        # cases, this will save us from recalculating them later.
        for con in self._oconstraints:
            con.update()

        # Validate!
        self._validate()

        self._ready = True

        return

    def __verify_profiles(self):
        """Verify that each FitContribution has a Profile."""
        # Check for profile values
        for con in self._contributions.values():
            if con.profile is None:
                m = "FitContribution '%s' does not have a Profile" % con.name
                raise AttributeError(m)
            if (
                con.profile.x is None
                or con.profile.y is None
                or con.profile.dy is None
            ):

                m = "Profile for '%s' is missing data" % con.name
                raise AttributeError(m)
        return

    def __verify_parameters(self):
        """Verify that all Parameters have values."""
        # Get all parameters with a value of None
        badpars = []
        for par in self.iterate_over_parameters():
            try:
                par.getValue()
            except ValueError:
                badpars.append(par)

        # Get the bad names
        badnames = []
        for par in badpars:
            objlist = self._locate_managed_object(par)
            names = [obj.name for obj in objlist]
            badnames.append(".".join(names))

        # Construct an error message, if necessary
        m = ""
        if len(badnames) == 1:
            m = "%s is not defined or needs an initial value" % badnames[0]
        elif len(badnames) > 0:
            s1 = ",".join(badnames[:-1])
            s2 = badnames[-1]
            m = "%s and %s are not defined or need initial values" % (s1, s2)

        if m:
            raise AttributeError(m)

        return

    def __collect_constraints_and_restraints(self):
        """Collect the Constraints and Restraints from subobjects."""
        from functools import cmp_to_key
        from itertools import chain

        rset = set(self._restraints)
        cdict = {}

        for org in chain(self._contributions.values(), self._parsets.values()):
            rset.update(org._get_restraints())
            cdict.update(org._get_constraints())
        cdict.update(self._constraints)

        # The order of the restraint list does not matter
        self._restraintlist = list(rset)

        # Reorder the constraints. Constraints are ordered such that a given
        # constraint is placed before its dependencies.
        self._oconstraints = list(cdict.values())

        # Create a depth-1 map of the constraint dependencies
        depmap = {}
        for con in self._oconstraints:
            depmap[con] = set()
            # Now check the constraint's equation for constrained arguments
            for arg in con.eq.args:
                if arg in cdict:
                    depmap[con].add(cdict[arg])

        # Turn the dependency map into multi-level map.
        def _extendDeps(con):
            deps = set(depmap[con])
            for dep in depmap[con]:
                deps.update(_extendDeps(dep))

            return deps

        for con in depmap:
            depmap[con] = _extendDeps(con)

        # Now sort the constraints based on the dependency map.
        def cmp(x, y):
            # x == y if neither of them have dependencies
            if not depmap[x] and not depmap[y]:
                return 0
            # x > y if y is a dependency of x
            # x > y if y has no dependencies
            if y in depmap[x] or not depmap[y]:
                return 1
            # x < y if x is a dependency of y
            # x < y if x has no dependencies
            if x in depmap[y] or not depmap[x]:
                return -1
            # If there are dependencies, but there is no relationship, the
            # constraints are equivalent
            return 0

        self._oconstraints.sort(key=cmp_to_key(cmp))

        return

    # Variable manipulation

    def add_variable(
        self, par, value=None, name=None, fixed=False, tag=None, tags=[]
    ):
        """Add a variable to be refined.

        Parameters
        ----------
        par : diffpy.srfit.fitbase.Parameter
            The Parameter that will be varied during a fit.
        value : float or None, optional
            The initial value for the variable. If this is None
            (default), then the current value of par will be used.
        name : str or None, optional
            The name for this variable. If name is None (default), then
            the name of the parameter will be used.
        fixed : bool, optional
            Fix the variable so that it does not vary (default False).
        tag : str or None, optional
            The tag for the variable. This can be used to retrieve, fix
            or free variables by tag (default None). Note that a
            variable is automatically tagged with its name and "all".
        tags : list of str, optional
            The list of tags (default []). Both tag and tags can be
            applied.

        Returns
        -------
        ParameterProxy
            ParameterProxy (variable) for the passed Parameter.

        Raises
        ------
        ValueError
            If the name of the variable is already taken by
            another managed object.
        ValueError
            If par is constant.
        ValueError
            If par is constrained.
        """
        name = name or par.name
        if par.const:
            raise ValueError("The parameter '%s' is constant" % par)
        if par.constrained:
            raise ValueError("The parameter '%s' is constrained" % par)
        var = ParameterProxy(name, par)
        if value is not None:
            var.set_value(value)
        self._add_parameter(var)
        if fixed:
            self.fix(var)
        # Tag with passed tags and by name
        self._tagmanager.tag(var, var.name)
        self._tagmanager.tag(var, "all")
        self._tagmanager.tag(var, *tags)
        if tag is not None:
            self._tagmanager.tag(var, tag)
        return var

    @deprecated(addVar_dep_msg)
    def addVar(
        self, par, value=None, name=None, fixed=False, tag=None, tags=[]
    ):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use diffpy.srfit.fitbase.FitRecipe.add_variable instead.
        """
        return self.add_variable(par, value, name, fixed, tag, tags)

    def delete_variable(self, var):
        """Remove a variable.

        Note that constraints and restraints involving the variable are not
        modified.

        Parameters
        ----------
        var : ParameterProxy
            A variable of the FitRecipe.

        Raises
        ------
        ValueError
            If var is not part of the FitRecipe.
        """
        self._remove_parameter(var)
        self._tagmanager.untag(var)
        return

    @deprecated(delVar_dep_msg)
    def delVar(self, var):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use diffpy.srfit.fitbase.FitRecipe.delete_variable instead.
        """
        self.delete_variable(var)
        return

    def __delattr__(self, name):
        if name in self._parameters:
            self.delete_variable(self._parameters[name])
            return
        super(FitRecipe, self).__delattr__(name)
        return

    def create_new_variable(
        self, name, value=None, fixed=False, tag=None, tags=[]
    ):
        """Create a new variable of the fit.

        This method lets new variables be created that are not tied to a
        Parameter.  Orphan variables may cause a fit to fail, depending on the
        optimization routine, and therefore should only be created to be used
        in constraint or restraint equations.

        Parameters
        ----------
        name : str
            The name of the variable. The variable will be able to be
            used by this name in restraint and constraint equations.
        value : float or None, optional
            The initial value for the variable. If this is None
            (default), then the variable will be given the value of the
            first non-None-valued Parameter constrained to it. If this
            fails, an error will be thrown when 'residual' is called.
        fixed : bool, optional
            Fix the variable so that it does not vary (default False).
            The variable will still be managed by the FitRecipe.
        tag : str or None, optional
            The tag for the variable. This can be used to fix and free
            variables by tag (default None). Note that a variable is
            automatically tagged with its name and "all".
        tags : list of str, optional
            The list of tags (default []). Both tag and tags can be
            applied.

        Returns
        -------
        Parameter
            The new variable (Parameter instance).
        """
        # This will fix the Parameter
        var = self._new_parameter(name, value)
        # We may explicitly free it
        if not fixed:
            self.free(var)
        # Tag with passed tags
        self._tagmanager.tag(var, *tags)
        if tag is not None:
            self._tagmanager.tag(var, tag)

        return var

    @deprecated(newVar_dep_msg)
    def newVar(self, name, value=None, fixed=False, tag=None, tags=[]):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use diffpy.srfit.fitbase.FitRecipe.create_new_variable instead.
        """
        return self.create_new_variable(name, value, fixed, tag, tags)

    def _new_parameter(self, name, value, check=True):
        """Overloaded to tag variables.

        See RecipeOrganizer._new_parameter
        """
        par = RecipeOrganizer._new_parameter(self, name, value, check)
        # tag this
        self._tagmanager.tag(par, par.name)
        self._tagmanager.tag(par, "all")
        self.fix(par.name)
        return par

    def __get_var_and_check(self, var):
        """Get the actual variable from var.

        Parameters
        ----------
        var
            A variable of the FitRecipe, or the name of a variable.

        Returns the variable or None if the variable cannot be found in the
        _parameters list.
        """
        if isinstance(var, str):
            var = self._parameters.get(var)

        if var not in self._parameters.values():
            raise ValueError("Passed variable is not part of the FitRecipe")

        return var

    def __get_vars_from_args(self, *args, **kw):
        """Get a list of variables from passed arguments.

        This method accepts string or variable arguments. An argument of
        "all" selects all variables. Keyword arguments must be parameter
        names, followed by a value to assign to the fixed variable. This
        method is used by the fix and free methods.

        Raises ValueError if an unknown variable, name or tag is passed,
        or if a tag is passed in a keyword.
        """
        # Process args. Each variable is tagged with its name, so this is easy.
        strargs = set([arg for arg in args if isinstance(arg, str)])
        varargs = set(args) - strargs
        # Check that the tags are valid
        alltags = set(self._tagmanager.alltags())
        badtags = strargs - alltags
        if badtags:
            names = ",".join(badtags)
            raise ValueError("Variables or tags cannot be found (%s)" % names)

        # Check that variables are valid
        allvars = set(self._parameters.values())
        badvars = varargs - allvars
        if badvars:
            names = ",".join(v.name for v in badvars)
            raise ValueError("Variables cannot be found (%s)" % names)

        # Make sure that we only have parameters in kw
        kwnames = set(kw.keys())
        allnames = set(self._parameters.keys())
        badkw = kwnames - allnames
        if badkw:
            names = ",".join(badkw)
            raise ValueError("Tags cannot be passed as keywords (%s)" % names)

        # Now get all the objects referred to in the arguments.
        varargs |= self._tagmanager.union(*strargs)
        varargs |= self._tagmanager.union(*kw.keys())
        return varargs

    def fix(self, *args, **kw):
        """Fix one or more parameters by reference, name, or tag.

        This method marks specified parameters as fixed, meaning they will not
        be refined during the fitting process. By default, all parameters are
        free (not fixed). Parameters can be specified using their references,
        names, or tags. Additionally, keyword arguments can be used to assign
        specific values to the fixed parameters.

        Parameters
        ----------
            *args : str or Parameter
                The positional arguments specifying the parameters to fix.
                These can be parameter objects, their names as strings, or
                tags. The special string "all" can be used to select all
                parameters.
            **kw : dict
                The keyword arguments where the keys are parameter names and
                the values are the values to assign to the corresponding
                fixed parameters.

        Raises
        ------
            ValueError:
                If an unknown parameter, name, or tag is passed, or if a
                tag is passed as a keyword argument.

        Example
        -------

        ::

            # Fix a parameter by reference
            recipe.fix(param1)

            # Fix a parameter by name
            recipe.fix("param2")

            # Fix all parameters
            recipe.fix("all")

            # Fix parameters by tag
            recipe.fix(tag="group1")

            # Fix a parameter and assign it a value
            recipe.fix(param3=10.0)
        """
        # Check the inputs and get the variables from them
        varargs = self.__get_vars_from_args(*args, **kw)

        # Fix all of these
        for var in varargs:
            self._tagmanager.tag(var, self._fixedtag)

        # Set the kw values
        for name, val in kw.items():
            self.get(name).value = val

        return

    def free(self, *args, **kw):
        """Free one or more parameters by reference, name, or tag.

        This method marks specified parameters as free, allowing them to be
        refined during the fitting process. By default, variables
        are free unless they are constrained. Constrained variables
        cannot be freed.

        Parameters
        ----------
        *args : str or Parameter
            The positional arguments specifying the parameters to free.
            These can be:
            - Parameter objects
            - Names of parameters (as strings)
            - Tags associated with parameters (as strings)
            - The string "all" to select all parameters.
        **kw : dict
            The keyword arguments specifying parameter names as keys and
            their values to assign after freeing. This is useful
            for setting the value of a parameter while marking it as free.

        Raises
        ------
        ValueError
            If an unknown parameter, name, or tag is passed, or if a
            tag is passed as a keyword argument.

        Notes
        -----
        - Parameters that are already free will remain free.
        - Tags associated with fixed parameters will be removed when they
          are freed.
        - If keyword arguments are provided, the corresponding parameter values
          will be updated after freeing.

        Returns
        -------
        None
        """
        # Check the inputs and get the variables from them
        varargs = self.__get_vars_from_args(*args, **kw)

        # Free all of these
        for var in varargs:
            if not var.constrained:
                self._tagmanager.untag(var, self._fixedtag)

        # Set the kw values
        for name, val in kw.items():
            self.get(name).value = val

        return

    def is_free(self, var):
        """Determine if a variable is free (not fixed) in the fit
        recipe.

        This method checks whether the specified variable does not have the
        fixed tag associated with it, indicating that it is free to vary
        during the fitting process.

        Parameters
        ----------
        var : object
            The variable to check. This is typically an instance of a parameter
            or variable object used in the fit recipe.

        Returns
        -------
        bool
            True if the variable is free (not fixed), False otherwise.
        """
        return not self._tagmanager.hasTags(var, self._fixedtag)

    @deprecated(isFree_dep_msg)
    def isFree(self, var):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use diffpy.srfit.fitbase.FitRecipe.is_free instead.
        """
        return self.is_free(var)

    def remove_constraint(self, *pars):
        """Unconstrain a Parameter.

        This removes any constraints on a Parameter. If the Parameter is also a
        variable of the recipe, it will be freed as well.

        Parameters
        ----------
        *pars : str or Parameter
            The names of Parameters or Parameter objects to unconstrain.

        Raises
        ------
        ValueError
            If the Parameter is not constrained.
        """
        update = False
        for par in pars:
            if isinstance(par, str):
                name = par
                par = self.get(name)

            if par is None:
                raise ValueError("The parameter cannot be found")

            if par in self._constraints:
                self._constraints[par].remove_constraint()
                del self._constraints[par]
                update = True

            if par in self._parameters.values():
                self._tagmanager.untag(par, self._fixedtag)

        if update:
            # Our configuration changed
            self._update_configuration()

        return

    @deprecated(unconstrain_dep_msg)
    def unconstrain(self, *pars):
        """This function has been deprecated and will be removed in
        version 4.0.0.

        Please use diffpy.srfit.fitbase.FitRecipe.remove_constraint
        instead.
        """
        self.remove_constraint(*pars)
        return

    def add_constraint(self, par, con, ns={}):
        """Constrain a parameter to an equation.

        Note that only one constraint can exist on a Parameter at a time.

        This is overloaded to set the value of con if it represents a variable
        and its current value is None. A constrained variable will be set as
        fixed.

        Parameters
        ----------
        par : Parameter
            The Parameter to constrain.
        con : str or Parameter
            The string representation of the constraint equation or a
            Parameter to constrain to. A constraint equation must
            consist of numpy operators and "known" Parameters.
            Parameters are known if they are in the ns argument, or if
            they are managed by this object.
        ns : dict, optional
            The dictionary of Parameters, indexed by name, that are used
            in the eqstr, but not part of this object (default {}).

        Raises
        ------
        ValueError
            If ns uses a name that is already used for a variable.
        ValueError
            If eqstr depends on a Parameter that is not part of the FitRecipe
            and that is not defined in ns.
        ValueError
            If par is marked as constant.
        """
        if isinstance(par, str):
            name = par
            par = self.get(name)
            if par is None:
                par = ns.get(name)
            if par is None:
                raise ValueError("The parameter '%s' cannot be found" % name)

        if con in self._parameters.keys():
            con = self._parameters[con]

        if par.const:
            raise ValueError("The parameter '%s' is constant" % par)

        # This will pass the value of a constrained parameter to the initial
        # value of a parameter constraint.
        if con in self._parameters.values():
            val = con.getValue()
            if val is None:
                val = par.getValue()
                con.set_value(val)

        if par in self._parameters.values():
            self.fix(par)

        RecipeOrganizer.add_constraint(self, par, con, ns)
        return

    @deprecated(constrain_dep_msg)
    def constrain(self, par, con, ns={}):
        """This function has been deprecated and will be removed in
        version 4.0.0.

        Please use diffpy.srfit.fitbase.FitRecipe.add_constraint
        instead.
        """
        self.add_constraint(par, con, ns)
        return

    def get_values(self):
        """Retrieve the current values of all free variables in the fit
        recipe.

        This method collects the values of all parameters that are marked as
        free (i.e., adjustable during the fitting process) and returns them
        as a NumPy array.

        Returns
        -------

        values_array : numpy.ndarray
            The array containing the current values of all free
            variables in the fit recipe.
        """
        values_array = array(
            [v.value for v in self._parameters.values() if self.is_free(v)]
        )
        return values_array

    @deprecated(getValues_dep_msg)
    def getValues(self):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use diffpy.srfit.fitbase.FitRecipe.get_values instead."""
        return self.get_values()

    def get_names(self):
        """Retrieve the names of all free variables in the fit recipe.

        This method iterates through the parameters in the fit recipe and
        returns a list of names for those variables that are marked as free.

        Returns
        -------
        parameter_names :list of str
            The list containing the names of free variables.
        """
        parameter_names = [
            v.name for v in self._parameters.values() if self.is_free(v)
        ]
        return parameter_names

    @deprecated(getNames_dep_msg)
    def getNames(self):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use diffpy.srfit.fitbase.FitRecipe.get_names instead."""
        return self.get_names()

    def get_bounds_pairs(self):
        """Get the bounds on variables in a list.

        Returns
        -------
        bounds_pair_list : list of tuple of float
            The list of ``(lower, upper)`` bounds on the variables, in the same
            order as ``get_names`` and ``get_values``.
        """
        return [v.bounds for v in self._parameters.values() if self.is_free(v)]

    @deprecated(getBounds_dep_msg)
    def getBounds(self):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use diffpy.srfit.fitbase.FitRecipe.get_bounds_pairs
        instead.
        """
        return self.get_bounds_pairs()

    def get_bounds_array(self):
        """Get the bounds on variables in two numpy arrays.

        Returns
        -------
        lower_bounds : numpy.ndarray
            The numpy array of lower bounds on the variables, in the same order
            as ``get_names`` and ``get_values``.
        upper_bounds : numpy.ndarray
            The numpy array of upper bounds on the variables, in the same order
            as ``get_names`` and ``get_values``.
        """
        bounds = self.get_bounds_pairs()
        lower_bounds = array([b[0] for b in bounds])
        upper_bounds = array([b[1] for b in bounds])
        return lower_bounds, upper_bounds

    @deprecated(getBounds2_dep_msg)
    def getBounds2(self):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use diffpy.srfit.fitbase.FitRecipe.get_bounds_array instead.
        """
        return self.get_bounds_array()

    def initialize_recipe_with_recipe(self, recipe_object):
        """Initialize a FitRecipe with another FitRecipe.

        This is used to initialize a FitRecipe with the contribution(s),
        parameters, constraints and restraints of another FitRecipe.
        If a duplicate contribution, parameter, constraint, or restraint
        is added to the FitRecipe you are initializing, the value from the
        added object will be used.

        Parameters
        ----------
        recipe_object : FitRecipe
            The FitRecipe to initialize with.

        Raises
        ------
        ValueError
            If the object passed is not a FitRecipe.
        """
        if not isinstance(recipe_object, FitRecipe):
            raise ValueError(
                "The input recipe_object must be a FitRecipe, "
                f"but got {type(recipe_object)}."
            )

        for contrib_object in recipe_object._contributions.values():
            if contrib_object not in self._contributions.values():
                self.add_contribution(contrib_object)

        for param_name, param_object in recipe_object._parameters.items():
            if param_name not in self._parameters:
                self._parameters.update({param_name: param_object})

        for (
            parameter_object,
            constraint_object,
        ) in recipe_object._constraints.items():
            if parameter_object not in self._constraints:
                self._constraints.update({parameter_object: constraint_object})

        for restraint in recipe_object._restraints:
            if restraint not in self._restraints:
                self._restraints.add(restraint)

    def _pretty_print_results_dict(self, params_dict):
        """Pretty print a dictionary of parameter names and values."""
        sorted_params = sorted(params_dict.items())
        width = max(len(name) for name, _ in sorted_params)
        for name, value in sorted_params:
            if isinstance(value, float):
                value_str = f"{value:.6g}"
            else:
                value_str = str(value)
            print(f"  {name:<{width}} = {value_str}")

    def _set_parameters_from_dict(self, params_dict):
        """Set the parameters of the FitRecipe from a dictionary of
        parameter names and values."""
        for param_name, param_value in params_dict.items():
            if param_name in self._parameters:
                self._parameters[param_name].set_value(param_value)
            else:
                print(
                    f"Warning: Parameter '{param_name}' from results "
                    "not found in FitRecipe and will be ignored."
                )

    def initialize_recipe_with_results(self, results, verbose=True):
        """Initialize a FitRecipe with a FitResults object or a results
        file.

        Note that at least one FitContribution must already exist in
        the FitRecipe.

        Parameters
        ----------
        results : FitResults, pathlib.Path, or str
            The FitResults object or path to results file to initialize with.
        verbose : bool, optional
            If True, print warnings for any parameters in the results that are
            not in the FitRecipe. Default is True.

        Raises
        ------
        ValueError
            If the input results is not a FitResults object or a path to a
            results file.
        """
        if hasattr(results, "get_results_dictionary"):
            params_dict = results.get_results_dictionary()
            metrics_in_dict = [
                "Residual",
                "Contributions",
                "Restraints",
                "Chi2",
                "Reduced Chi2",
                "Rw",
            ]
            for metric in metrics_in_dict:
                params_dict.pop(metric, None)
        elif isinstance(results, (str, Path)):
            params_dict = utils.get_dict_from_results_file(results)
        else:
            raise ValueError(
                "The input results must be a FitResults object or a path to a "
                f"results file, but got {type(results)}."
            )
        self._set_parameters_from_dict(params_dict)
        if verbose:
            print()
            print("Parameters found in Results:")
            print("=" * 30)
            self._pretty_print_results_dict(params_dict)
            print()
            print("Parameters set in FitRecipe:")
            print("=" * 30)
            set_parameters_dict = {
                param.name: param.getValue()
                for param in self._parameters.values()
            }
            self._pretty_print_results_dict(set_parameters_dict)

    def set_plot_defaults(self, **kwargs):
        """Set default plotting options for all future plots.

        Any keyword argument accepted by plot_recipe() can be set here.

        Parameters
        ----------
        show_observed : bool, optional
            The observed data is plotted if True. Default is True.
        show_fit : bool, optional
            The fit to the data is plotted if True. Default is True.
        show_diff : bool, optional
            The difference curve (observed - calculated) is plotted if True.
            Default is True.
        offset_scale : float, optional
            The scaling factor for the difference curve offset. The difference
            curve is offset below the data by
            (min_y - 0.1*range) * offset_scale. Default is 1.0.
        xmin : float or None, optional
            The minimum x value to plot. If None, uses the minimum x value
            of the data. Default is None.
        xmax : float or None, optional
            The maximum x value to plot. If None, uses the maximum x value
            of the data. Default is None.
        figsize : tuple, optional
            The figure size as (width, height). Default is (8, 6).
        data_style : str, optional
            The matplotlib line/marker style for data points. Default is "o".
        fit_style : str, optional
            The matplotlib line/marker style for the calculated fit.
            Default is "-".
        diff_style : str, optional
            The matplotlib line/marker style for the difference curve.
            Default is "-".
        data_color : str or None, optional
            The color for data plot. If None, uses default matplotlib colors.
        fit_color : str or None, optional
            The color for the fit plot. If None, uses default matplotlib
            colors.
        diff_color : str or None, optional
            The color for the difference plot. If None, uses default
            matplotlib colors.
        data_label : str, optional
            The legend label for observed data. Default is "Observed".
        fit_label : str, optional
            The legend label for the calculated fit. Default is "Calculated".
        diff_label : str, optional
            The legend label for the difference curve. Default is "Difference".
        xlabel : str, optional
            The label for the x-axis.
        ylabel : str, optional
            The label for the y-axis.
        title : str or None, optional
            The plot title. Default is no title.
        legend : bool, optional
            The legend is shown if True. Default is True.
        legend_loc : str, optional
            The legend location. Default is "best".
        grid : bool, optional
            The grid is shown if True. Default is False.
        markersize : float, optional
            The size of data point markers.
        linewidth : float, optional
            The width of fit and difference lines.
        alpha : float, optional
            The transparency of all plot elements (0=transparent, 1=opaque).
            Default is 1.0.
        show : bool, optional
            The plot is displayed using `plt.show()` if True. Default is True.
        ax : matplotlib.axes.Axes or None, optional
            The axes object to plot on. If None, creates a new figure.
            Default is None.
        return_fig : bool, optional
            The figure and axes objects are returned if True. Default is False.

        Examples
        --------
        >>> recipe.set_plot_defaults(
                xlabel='r (Å)',
                ylabel='G(r) (Å⁻²)',
                data_color='black',
                fit_color='red'
            )
        """
        for key in kwargs:
            if key not in self.plot_options:
                print(
                    f"Warning: '{key}' is not a valid "
                    "plot_recipe option and will be ignored."
                )
        self.plot_options.update(kwargs)

    def _set_axes_labels_from_metadata(self, meta, plot_params):
        """Set axes labels based on filename suffix in profile metadata
        if not already set."""
        if isinstance(meta, dict):
            filename = meta.get("filename")
            if filename:
                suffix = filename.rsplit(".", 1)[-1].lower()
                if "gr" in suffix:
                    if plot_params.get("xlabel") is None:
                        plot_params["xlabel"] = r"r ($\mathrm{\AA}$)"
                    if plot_params.get("ylabel") is None:
                        plot_params["ylabel"] = r"G ($\mathrm{\AA}^{-2}$)"
        return

    def plot_recipe(self, ax=None, return_fig=False, **kwargs):
        """Plot the observed, fit, and difference curves for each
        contribution of the fit recipe.

        If the recipe has multiple contributions, a separate
        plot is created for each contribution.

        Parameters
        ----------
        ax : matplotlib.axes.Axes or None, optional
            The axes object to plot on. If None, creates a new figure.
            Default is None.
        return_fig : bool, optional
            The figure and axes objects are returned if True. Default is False.
        **kwargs : dict
            Any plotting option can be passed to override the defaults in
            `FitRecipe().plot_options`. See the
            `FitRecipe().set_plot_defaults()` method for available
            keyword arguments.

        Returns
        -------
        fig, axes : tuple of (mpl.figure.Figure, list of mpl.axes.Axes)
            The figure object and a list of axes objects (one per contribution)
            are returned if return_fig=True.

        Examples
        --------
        Plot with default settings:

        >>> recipe.plot_recipe()

        Override defaults for one plot:

        >>> recipe.plot_recipe(show_diff=False, title='My Custom Title')

        Set defaults once, use everywhere:

        >>> recipe.set_plot_defaults(xlabel='r (Å)', ylabel='G(r)')
        >>> recipe.plot_recipe()  # Uses xlabel and ylabel
        >>> recipe.plot_recipe()  # Still uses them

        Override a default for one plot:

        >>> recipe.set_plot_defaults(figsize=(10, 7))
        >>> recipe.plot_recipe()  # Uses (10, 7)
        >>> recipe.plot_recipe(figsize=(12, 8))  # Temporarily uses (12, 8)
        >>> recipe.plot_recipe()  # Back to (10, 7)

        Notes
        -----
        The default values are taken from recipe.plot_options. You can modify
        these defaults in three ways:

        1. Using set_plot_defaults():
        recipe.set_plot_defaults(xlabel='r (Å)')

        2. Direct attribute access:
        recipe.plot_options['xlabel'] = 'r (Å)'

        3. Using update():
        recipe.plot_options.update({'xlabel': 'r (Å)', 'ylabel': 'G(r)'})
        """
        plot_params = self.plot_options.copy()
        plot_params.update(kwargs)

        if not any(
            [
                plot_params["show_observed"],
                plot_params["show_fit"],
                plot_params["show_diff"],
            ]
        ):
            raise ValueError(
                "At least one of show_observed, show_fit, "
                "or show_diff must be True"
            )

        if not self._contributions:
            raise ValueError(
                "No contributions found in recipe. "
                "Add contributions before plotting."
            )
        figures = []
        axes_list = []
        for name, contrib in self._contributions.items():
            profile = contrib.profile
            x = profile.x
            yobs = profile.y
            ycalc = profile.ycalc
            if ycalc is None:
                if plot_params["show_fit"] or plot_params["show_diff"]:
                    print(
                        f"Contribution '{name}' has no calculated values "
                        "(ycalc is None). "
                        "Only observed data will be plotted."
                    )
                plot_params["show_fit"] = False
                plot_params["show_diff"] = False
            else:
                diff = yobs - ycalc
                y_min = min(yobs.min(), ycalc.min())
                y_max = max(yobs.max(), ycalc.max())
                y_range = y_max - y_min
                base_offset = y_min - 0.1 * y_range
                offset = base_offset * plot_params["offset_scale"]
            if ax is None:
                fig = plt.figure(figsize=plot_params["figsize"])
                current_ax = fig.add_subplot(111)
            else:
                current_ax = ax
                fig = current_ax.figure
            if plot_params["show_observed"]:
                current_ax.plot(
                    x,
                    yobs,
                    plot_params["data_style"],
                    label=plot_params["data_label"],
                    color=plot_params["data_color"],
                    markersize=plot_params["markersize"],
                    alpha=plot_params["alpha"],
                )
            if plot_params["show_fit"]:
                current_ax.plot(
                    x,
                    ycalc,
                    plot_params["fit_style"],
                    label=plot_params["fit_label"],
                    color=plot_params["fit_color"],
                    linewidth=plot_params["linewidth"],
                    alpha=plot_params["alpha"],
                )
            if plot_params["show_diff"]:
                current_ax.plot(
                    x,
                    diff + offset,
                    plot_params["diff_style"],
                    label=plot_params["diff_label"],
                    color=plot_params["diff_color"],
                    linewidth=plot_params["linewidth"],
                    alpha=plot_params["alpha"],
                )
                current_ax.axhline(
                    offset,
                    color="black",
                )
            meta = getattr(profile, "meta", None)
            if meta:
                self._set_axes_labels_from_metadata(meta, plot_params)
            if plot_params["xlabel"] is not None:
                current_ax.set_xlabel(plot_params["xlabel"])
            if plot_params["ylabel"] is not None:
                current_ax.set_ylabel(plot_params["ylabel"])
            if plot_params["title"] is not None:
                current_ax.set_title(plot_params["title"])
            if plot_params["legend"]:
                current_ax.legend(loc=plot_params["legend_loc"], frameon=True)
            if plot_params["grid"]:
                current_ax.grid(True)
            if (
                plot_params["xmin"] is not None
                or plot_params["xmax"] is not None
            ):
                current_ax.set_xlim(
                    left=plot_params["xmin"], right=plot_params["xmax"]
                )
            fig.tight_layout()
            figures.append(fig)
            axes_list.append(current_ax)
            if plot_params["show"] and ax is None:
                plt.show()
        if return_fig:
            if len(figures) == 1:
                return figures[0], axes_list[0]
            else:
                return figures, axes_list

    def convert_bounds_to_restraints(self, sig=1, scaled=False):
        """Turn all bounded parameters into restraints.

        The bounds become limits on the restraint.

        Parameters
        ----------
        sig : float or iterable of float, optional
            The number of standard deviations associated with each bound.
            Smaller values produce stronger restraints. If a scalar is given,
            the same value is applied to all parameters. If an iterable is
            provided, it must match the number of parameters. Default is 1.

        scaled : bool, optional
            If True, scale each restraint by the magnitude of the corresponding
            parameter, consistent with the behavior of :meth:`restrain`.
            Default is False.
        """
        pars = self._parameters.values()
        if not hasattr(sig, "__iter__"):
            sig = [sig] * len(pars)
        for par, x in zip(pars, sig):
            self.add_soft_bounds(
                par, par.bounds[0], par.bounds[1], sig=x, scaled=scaled
            )
        return

    @deprecated(boundsToRestraints_dep_msg)
    def boundsToRestraints(self, sig=1, scaled=False):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use diffpy.srfit.fitbase.FitRecipe.convert_bounds_to_restraints
        instead.
        """
        self.convert_bounds_to_restraints(sig, scaled)
        return

    def _apply_values(self, p):
        """Apply variable values to the variables."""
        if len(p) == 0:
            return
        vargen = (v for v in self._parameters.values() if self.is_free(v))
        for var, pval in zip(vargen, p):
            var.set_value(pval)
        return

    def _update_configuration(self):
        """Notify RecipeContainers in hierarchy of configuration
        change."""
        self._ready = False
        return


# End of file
