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

import matplotlib.pyplot as plt
from bg_mpl_stylesheets.styles import all_styles
from numpy import array, concatenate, dot, sqrt

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


class FitRecipe(_fitrecipe_interface, RecipeOrganizer):
    """FitRecipe class.

    Attributes
    ----------
    name
        A name for this FitRecipe.
    fithooks
        List of FitHook instances that can pass information out
        of the system during a refinement. By default, the is
        populated by a PrintFitHook instance.
    _constraints
        A dictionary of Constraints, indexed by the constrained
        Parameter. Constraints can be added using the
        'constrain' method.
    _oconstraints
        An ordered list of the constraints from this and all
        sub-components.
    _calculators
        A managed dictionary of Calculators.
    _contributions
        A managed OrderedDict of FitContributions.
    _parameters
        A managed OrderedDict of parameters (in this case the
        parameters are varied).
    _parsets
        A managed dictionary of ParameterSets.
    _eqfactory
        A diffpy.srfit.equation.builder.EquationFactory
        instance that is used to create constraints and
        restraints from string
    _restraintlist
        A list of restraints from this and all sub-components.
    _restraints
        A set of Restraints. Restraints can be added using the
        'restrain' or 'confine' methods.
    _ready
        A flag indicating if all attributes are ready for the
        calculation.
    _tagmanager
        A TagManager instance for managing tags on Parameters.
    _weights
        List of weighing factors for each FitContribution. The
        weights are multiplied by the residual of the
        FitContribution when determining the overall residual.
    _fixedtag
        "__fixed", used for tagging variables as fixed. Don't
        use this tag unless you want issues.

    Properties
    ----------
    names
        Variable names (read only). See get_names.
    values
        Variable values (read only). See get_values.
    fixednames
        Names of the fixed refinable variables (read only).
    fixedvalues
        Values of the fixed refinable variables (read only).
    bounds
        Bounds on parameters (read only). See get_bounds_pairs.
    bounds2
        Bounds on parameters (read only). See get_bounds_array.
    """

    fixednames = property(
        lambda self: [
            v.name
            for v in self._parameters.values()
            if not (self.is_free(v) or self.isConstrained(v))
        ],
        doc="names of the fixed refinable variables",
    )
    fixedvalues = property(
        lambda self: array(
            [
                v.value
                for v in self._parameters.values()
                if not (self.is_free(v) or self.isConstrained(v))
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

        Attributes
        ----------
        fithook
            FitHook instance to add to the sequence
        index
            Index for inserting fithook into the list of fit hooks.  If
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

        Attributes
        ----------
        fithook
            FitHook instance to remove from the sequence. If this is
            None (default), default to index.
        index
            Index of FitHook instance to remove (default -1).


        Raises ValueError if fithook is not None, but is not present in the
        sequence.
        Raises IndexError if the sequence is empty or index is out of range.
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

        Attributes
        ----------
        con
            The FitContribution to be stored.


        Raises ValueError if the FitContribution has no name
        Raises ValueError if the FitContribution has the same name as some
        other managed object.
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
        """Set the weight of a FitContribution."""
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

        Attributes
        ----------
        parset
            The ParameterSet to be stored.


        Raises ValueError if the ParameterSet has no name.
        Raises ValueError if the ParameterSet has the same name as some other
        managed object.
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

        Raises ValueError if parset is not managed by this object.
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

        Parameters
        ----------
        p
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
        p
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

        Raises AttributeError if there are variables without a value.
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
        for par in self.iterPars():
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

        Attributes
        ----------
        par
            A Parameter that will be varied during a fit.
        value
            An initial value for the variable. If this is None
            (default), then the current value of par will be used.
        name
            A name for this variable. If name is None (default), then
            the name of the parameter will be used.
        fixed
            Fix the variable so that it does not vary (default False).
        tag
            A tag for the variable. This can be used to retrieve, fix
            or free variables by tag (default None). Note that a
            variable is automatically tagged with its name and "all".
        tags
            A list of tags (default []). Both tag and tags can be
            applied.


        Returns
        -------
        vars
            ParameterProxy (variable) for the passed Parameter.


        Raises ValueError if the name of the variable is already taken by
        another managed object.
        Raises ValueError if par is constant.
        Raises ValueError if par is constrained.
        """
        name = name or par.name
        if par.const:
            raise ValueError("The parameter '%s' is constant" % par)
        if par.constrained:
            raise ValueError("The parameter '%s' is constrained" % par)
        var = ParameterProxy(name, par)
        if value is not None:
            var.setValue(value)
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

        Attributes
        ----------
        var
            A variable of the FitRecipe.


        Raises ValueError if var is not part of the FitRecipe.
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

        Attributes
        ----------
        name
            The name of the variable. The variable will be able to be
            used by this name in restraint and constraint equations.
        value
            An initial value for the variable. If this is None
            (default), then the variable will be given the value of the
            first non-None-valued Parameter constrained to it. If this
            fails, an error will be thrown when 'residual' is called.
        fixed
            Fix the variable so that it does not vary (default False).
            The variable will still be managed by the FitRecipe.
        tag
            A tag for the variable. This can be used to fix and free
            variables by tag (default None). Note that a variable is
            automatically tagged with its name and "all".
        tags
            A list of tags (default []). Both tag and tags can be
            applied.


        Returns the new variable (Parameter instance).
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

        Attributes
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
        """Fix a parameter by reference, name or tag.

        A fixed variable is not refined.  Variables are free by default.

        This method accepts string or variable arguments. An argument of
        "all" selects all variables. Keyword arguments must be parameter
        names, followed by a value to assign to the fixed variable.

        Raises ValueError if an unknown Parameter, name or tag is
        passed, or if a tag is passed in a keyword.
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
        """Free a parameter by reference, name or tag.

        A free variable is refined.  Variables are free by default.
        Constrained variables are not free.

        This method accepts string or variable arguments. An argument of
        "all" selects all variables. Keyword arguments must be parameter
        names, followed by a value to assign to the fixed variable.

        Raises ValueError if an unknown Parameter, name or tag is
        passed, or if a tag is passed in a keyword.
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
        """Check if a variable is fixed."""
        return not self._tagmanager.hasTags(var, self._fixedtag)

    @deprecated(isFree_dep_msg)
    def isFree(self, var):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use diffpy.srfit.fitbase.FitRecipe.is_free instead.
        """
        return self.is_free(var)

    def unconstrain(self, *pars):
        """Unconstrain a Parameter.

        This removes any constraints on a Parameter. If the Parameter is also a
        variable of the recipe, it will be freed as well.

        Attributes
        ----------
        *pars
            The names of Parameters or Parameters to unconstrain.


        Raises ValueError if the Parameter is not constrained.
        """
        update = False
        for par in pars:
            if isinstance(par, str):
                name = par
                par = self.get(name)

            if par is None:
                raise ValueError("The parameter cannot be found")

            if par in self._constraints:
                self._constraints[par].unconstrain()
                del self._constraints[par]
                update = True

            if par in self._parameters.values():
                self._tagmanager.untag(par, self._fixedtag)

        if update:
            # Our configuration changed
            self._update_configuration()

        return

    def constrain(self, par, con, ns={}):
        """Constrain a parameter to an equation.

        Note that only one constraint can exist on a Parameter at a time.

        This is overloaded to set the value of con if it represents a variable
        and its current value is None. A constrained variable will be set as
        fixed.

        Attributes
        ----------
        par
            The Parameter to constrain.
        con
            A string representation of the constraint equation or a
            Parameter to constrain to.  A constraint equation must
            consist of numpy operators and "known" Parameters.
            Parameters are known if they are in the ns argument, or if
            they are managed by this object.
        ns
            A dictionary of Parameters, indexed by name, that are used
            in the eqstr, but not part of this object (default {}).


        Raises ValueError if ns uses a name that is already used for a
        variable.
        Raises ValueError if eqstr depends on a Parameter that is not part of
        the FitRecipe and that is not defined in ns.
        Raises ValueError if par is marked as constant.
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
                con.setValue(val)

        if par in self._parameters.values():
            self.fix(par)

        RecipeOrganizer.constrain(self, par, con, ns)
        return

    def get_values(self):
        """Get the current values of the variables in a list."""
        return array(
            [v.value for v in self._parameters.values() if self.is_free(v)]
        )

    @deprecated(getValues_dep_msg)
    def getValues(self):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use diffpy.srfit.fitbase.FitRecipe.get_values instead."""
        return self.get_values()

    def get_names(self):
        """Get the names of the variables in a list."""
        return [v.name for v in self._parameters.values() if self.is_free(v)]

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
            A list of ``(lower, upper)`` bounds on the variables, in the same
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
            A numpy array of lower bounds on the variables, in the same order
            as ``get_names`` and ``get_values``.
        upper_bounds : numpy.ndarray
            A numpy array of upper bounds on the variables, in the same order
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
        """Set axes labels based on filename suffix in profile metadata if not
        already set."""
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
        """Plot the observed, fit, and difference curves for each contribution
        of the fit recipe.

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

        Attributes
        ----------
        sig
            The uncertainty on the bounds (scalar or iterable,
            default 1).
        scaled
            Scale the restraints, see restrain.
        """
        pars = self._parameters.values()
        if not hasattr(sig, "__iter__"):
            sig = [sig] * len(pars)
        for par, x in zip(pars, sig):
            self.restrain(
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
            var.setValue(pval)
        return

    def _update_configuration(self):
        """Notify RecipeContainers in hierarchy of configuration change."""
        self._ready = False
        return


# End of file
