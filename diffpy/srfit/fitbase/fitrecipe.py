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
"""FitRecipe class. 

FitRecipes organize FitContributions, variables, Restraints and Constraints to
create a recipe of the system you wish to optimize. From the client's
perspective, the FitRecipe is a residual calculator. The residual method does
the work of updating variable values, which get propagated to the Parameters of
the underlying FitContributions via the varibles and Constraints.  This class
needs no special knowledge of the type of FitContribution or data being used.
Thus, it is suitable for combining residual equations from various types of
refinements into a single residual.

See the examples in the documentation for how to create an optimization problem
using FitRecipe.

"""
__all__ = ["FitRecipe"]

from numpy import array, concatenate, sqrt, dot

from diffpy.srfit.interface import _fitrecipe_interface
from diffpy.srfit.util.ordereddict import OrderedDict
from .parameter import ParameterProxy
from .recipeorganizer import RecipeOrganizer
from .fithook import FitHook

class FitRecipe(_fitrecipe_interface, RecipeOrganizer):
    """FitRecipe class.

    Attributes
    name            --  A name for this FitRecipe.
    fithook         --  An object to be called whenever within the residual
                        (default FitHook()) that can pass information out of
                        the system during a refinement.
    _constraints    --  A dictionary of Constraints, indexed by the constrained
                        Parameter. Constraints can be added using the
                        'constrain' method.
    _oconstraints   --  An ordered list of the constraints from this and all
                        sub-components.
    _calculators    --  A managed dictionary of Calculators.
    _contributions  --  A managed OrderedDict of FitContributions.
    _parameters     --  A managed OrderedDict of parameters (in this case the
                        parameters are varied).
    _parsets        --  A managed dictionary of ParameterSets.
    _eqfactory      --  A diffpy.srfit.equation.builder.EquationFactory
                        instance that is used to create constraints and
                        restraints from string
    _fixed          --  A set of parameters that are not actually varied.
    _restraintlist  --  A list of restraints from this and all sub-components.
    _restraints     --  A set of Restraints. Restraints can be added using the
                        'restrain' or 'confine' methods.
    _ready          --  A flag indicating if all attributes are ready for the
                        calculation.
    _tagdict        --  A dictionary of tags to variables.
    _weights        --  List of weighing factors for each FitContribution. The
                        weights are multiplied by the residual of the
                        FitContribution when determining the overall residual.

    """

    def __init__(self, name = "fit"):
        """Initialization."""
        RecipeOrganizer.__init__(self, name)
        self.fithook = FitHook()
        self._restraintlist = []
        self._oconstraints = []
        self._ready = False

        self._weights = []
        self._tagdict = {}

        self._fixed = set()

        self._parsets = {}
        self._manage(self._parsets)

        self._contributions = OrderedDict()
        self._manage(self._contributions)

        return

    def setFitHook(self, fithook):
        """Set a hook to be called within the residual method.

        The hook is an object for reporting updates, or more fundamentally,
        passing information out of the system during a refinement. It must have
        'precall' and 'postcall' methods, which are called at the start and at
        the end of the residual calculation. The precall method must accept a
        single argument, which is this FitRecipe object. The postcall method
        must accept the recipe and the chiv, vector residual.  It must also
        have a reset method that takes no arguments, which is called whenver
        the FitRecipe is prepared for a refinement.

        See the FitHook class for the interface.

        """
        self.fithook = fithook
        self._updateConfiguration()
        return

    def addContribution(self, con, weight = 1.0):
        """Add a FitContribution to the FitRecipe.

        con     --  The FitContribution to be stored.

        Raises ValueError if the FitContribution has no name
        Raises ValueError if the FitContribution has the same name as some
        other managed object.
        
        """
        self._addObject(con, self._contributions, True)
        self._weights.append(weight)
        return

    def setWeight(self, con, weight):
        """Set the weight of a FitContribution."""
        idx = self._contributions.values().index(con)
        self._weights[idx] = weight
        return

    def addParameterSet(self, parset):
        """Add a ParameterSet to the hierarchy.

        parset  --  The ParameterSet to be stored.

        Raises ValueError if the ParameterSet has no name.  
        Raises ValueError if the ParameterSet has the same name as some other
        managed object.

        """
        self._addObject(parset, self._parsets, True)
        return

    def removeParameterSet(self, parset):
        """Remove a ParameterSet from the hierarchy.

        Raises ValueError if parset is not managed by this object.

        """
        self._removeObject(parset, self._parsets)
        return

    def residual(self, p = []):
        """Calculate the residual to be optimized.

        Arguments
        p   --  The list of current variable values, provided in the same order
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

        if self.fithook:
            self.fithook.precall(self)

        # Update the variable parameters.
        self.__applyValues(p)

        # Update the constraints. These are ordered such that the list only
        # needs to be cycled once.
        for con in self._oconstraints:
            con.update()

        # Calculate the bare chiv
        chiv = concatenate([ 
            sqrt(self._weights[i])*\
                    self._contributions.values()[i].residual().flatten() \
                    for i in range(len(self._contributions))])

        # Calculate the point-average chi^2
        w = dot(chiv, chiv)/len(chiv)
        # Now we must append the restraints
        penalties = [ sqrt(res.penalty(w)) for res in self._restraintlist ]
        chiv = concatenate( [ chiv, penalties ] )

        if self.fithook:
            self.fithook.postcall(self, chiv)

        return chiv

    def scalarResidual(self, p = []):
        """A scalar version of the residual.

        See the residual method. This returns dot(chiv, chiv).

        """
        chiv = self.residual(p)
        return dot(chiv, chiv)

    def _prepare(self):
        """Prepare for the residual calculation, if necessary.

        This will prepare the data attributes to be used in the residual
        calculation.

        This updates the local restraints with those of the contributions.

        Raises AttributeError if there are variables without a value.

        """

        # Only prepare if the configuration has changed within the recipe
        # hierarchy.
        if self._ready:
            return

        # Validate!
        self._validate()

        # Inform the fit hook that we're updating things
        if self.fithook:
            self.fithook.reset(self)

        # Check Profiles
        self.__verifyProfiles()

        # Check parameters
        self.__verifyParameters()

        # Update constraints and restraints. 
        self.__collectConstraintsAndRestraints()

        self._ready = True

        return

    def __verifyProfiles(self):
        """Verify that each FitContribution has a Profile."""
        # Check for profile values
        for con in self._contributions.values():
            if con.profile is None:
                m = "FitContribution '%s' does not have a Profile"%con.name
                raise AttributeError(m)
            if con.profile.x is None or\
                con.profile.y is None or\
                con.profile.dy is None:

                    m = "Profile for '%s' is missing data"%con.name
                    raise AttributeError(m)

        return

    def __verifyParameters(self):
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
            objlist = self._locateManagedObject(par)
            names = [obj.name for obj in objlist]
            badnames.append( ".".join(names) )

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

    def __collectConstraintsAndRestraints(self):
        """Collect the Constraints and Restraints from subobjects."""
        rset = set(self._restraints)
        cdict = {}

        for org in self._contributions.values() + self._parsets.values():
            rset.update( org._getRestraints() )
            cdict.update( org._getConstraints() )
        cdict.update(self._constraints)

        # The order of the restraint list does not matter
        self._restraintlist = list(rset)

        # Reorder the constraints. Constraints are ordered such that a given
        # constraint is placed before its dependencies.
        self._oconstraints = cdict.values()

        # Create a depth-1 map of the constraint dependencies
        depmap = {}
        for con in self._oconstraints:
            depmap[con] = set()
            # Now check the constraint's equation for constrained arguments
            for arg in con.eq.args:
                if arg in cdict:
                    depmap[con].add( cdict[arg] )

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

        self._oconstraints.sort(cmp)

        return

    # Variable manipulation

    def addVar(self, par, value = None, name = None, fixed = False, tag = None,
            tags = []):
        """Add a variable to be refined.

        par     --  A Parameter that will be varied during a fit.
        value   --  An initial value for the variable. If this is None
                    (default), then the current value of par will be used.
        name    --  A name for this variable. If name is None (default), then
                    the name of the parameter will be used.
        fixed   --  Fix the variable so that it does not vary (default False).
        tag     --  A tag for the variable. This can be used to fix and free
                    variables by tag (default None).
        tags    --  A list of tags (default []). Both tag and tags can be
                    applied.

        Returns the ParameterProxy (variable) for the passed Parameter.

        Raises ValueError if the name of the variable is already taken by
        another managed object.
        Raises ValueError if par is constant.
        Raises ValueError if par is constrained.

        """
        name = name or par.name

        if par.const:
            raise ValueError("The parameter '%s' is constant"%par)

        if par.constrained:
            raise ValueError("The parameter '%s' is constrained"%par)

        var = ParameterProxy(name, par)
        if value is not None:
            var.setValue(value)

        self._addParameter(var)

        if fixed:
            self.fixVar(var)
          
        # Deal with tags
        self.__applyTags(var, tag, tags)

        return var

    def __applyTags(self, var, tag, tags):
        """Apply tags to a variable.

        tag     --  A tag for the variable. This can be used to fix and free
                    variables by tag (default None).
        tags    --  A list of tags (default []). Both tag and tags can be
                    applied.

        """
        if tag:
            self.__tagVar(var, tag)

        for t in tags:
            self.__tagVar(var, tag)

        return

    def __tagVar(self, var, tag):
        """Private function to tag a variable."""
        vset = self._tagdict.get(tag, set())
        vset.add(var)
        self._tagdict[tag] = vset
        return

    def delVar(self, var):
        """Remove a variable.

        Note that constraints and restraints involving the variable are not
        modified.

        var     --  A variable of the FitRecipe.

        Raises ValueError if var is not part of the FitRecipe.

        """

        self._removeParameter(var)
        self._fixed.discard(var)

        # Remove tags
        for vset in self._tagdict.items():
            vset.discard(var)
        
        return

    def __delattr__(self, name):
        if name in self._parameters:
            self.delVar( self._parameters[name] )
            return
        super(FitRecipe, self).__delattr__(name)
        return


    def newVar(self, name, value = None, fixed = False, tag = None, tags = []):
        """Create a new variable of the fit.

        This method lets new variables be created that are not tied to a
        Parameter.  Orphan variables may cause a fit to fail, depending on the
        optimization routine, and therefore should only be created to be used
        in contraint or restraint equations.

        name    --  The name of the variable. The variable will be able to be
                    used by this name in restraint and constraint equations.
        value   --  An initial value for the variable. If this is None
                    (default), then the variable will be given the value of the
                    first non-None-valued Parameter constrained to it. If this
                    fails, an error will be thrown when 'residual' is called.
        fixed   --  Fix the variable so that it does not vary (default False).
                    The variable will still be managed by the FitRecipe.
        tag     --  A tag for the variable. This can be used to fix and free
                    variables by tag (default None).
        tags    --  A list of tags (default []). Both tag and tags can be
                    applied.

        Returns the new variable (Parameter instance).

        """
        var = self._newParameter(name, value)

        if fixed:
            self.fixVar(var)

        # Deal with tags
        self.__applyTags(var, tag, tags)

        return var

    def fixVar(self, var, value = None):
        """Fix a variable so that it doesn't change.

        The variable will still be managed by the FitRecipe.

        var     --  A variable of the FitRecipe, or the name of a variable.
        value   --  A new value for the variable. If this is None
                    (default), then the value will not be changed.

        Raises ValueError if var is not part of the FitRecipe.

        """
        if isinstance(var, str):
            var = self._parameters.get(var)

        if var not in self._parameters.values():
            raise ValueError("Passed variable is not part of the FitRecipe")

        self._fixed.add(var)

        if value is not None:
            var.setValue(value)

        return

    def freeVar(self, var, value = None):
        """Free a variable so that it is refined.

        Variables are free by default.

        var     --  A variable of the FitRecipe, or the name of a variable.
        value   --  A new value for the variable. If this is None
                    (default), then the value will not be changed.

        Raises ValueError if var is not part of the FitRecipe.

        """
        if isinstance(var, str):
            var = self._parameters.get(var)

        if var not in self._parameters.values():
            raise ValueError("Passed variable is not part of the FitRecipe")

        self._fixed.discard(var)

        if value is not None:
            var.setValue(value)

        return

    def fixAll(self, *tags):
        """Fix all variables.

        Extra arguments are assumed to be tags. If present, only variables with
        the given tag will be fixed.

        Raises ValueError when passed tags do not rever to any variables.

        """
        for tag in tags:
            if tag not in self._tagdict:
                raise ValueError("Tag '%s' not found" % tag)

        tagset = self.__getTagSet(tags)

        if tagset:
            self._fixed.update(tagset)
        else:
            self._fixed.update(self._parameters.values())

        return

    def freeAll(self, *tags):
        """Free all variables.

        Extra arguments are assumed to be tags. If present, only variables with
        the given tag will be freed.

        Raises ValueError when passed tags do not rever to any variables.

        """
        for tag in tags:
            if tag not in self._tagdict:
                raise ValueError("Tag '%s' not found" % tag)

        tagset = self.__getTagSet(tags)

        if tagset:
            self._fixed.difference_update(tagset)
        else:
            self._fixed.clear()

        return

    def unconstrain(self, par, free = True):
        """Unconstrain a Parameter.

        This removes any constraints on a Parameter. 

        par     --  The name of a Parameter or a Parameter to unconstrain.
        free    --  Flag indicating whether to free the Parameter after
                    removing the constraint (bool, default True)
        
        Raises ValueError if the Parameter is not constrained.

        """
        if isinstance(par, str):
            par = self.get(par)
        RecipeOrganizer.unconstrain(self, par)
        if free:
            self._fixed.discard(par)
        return

    def constrain(self, par, con, ns = {}):
        """Constrain a parameter to an equation.

        Note that only one constraint can exist on a Parameter at a time.

        This is overloaded to set the value of con if it represents a variable
        and its current value is None. A constrained variable will be set as
        fixed.

        par     --  The Parameter to constrain.
        con     --  A string representation of the constraint equation or a
                    Parameter to constrain to.  A constraint equation must
                    consist of numpy operators and "known" Parameters.
                    Parameters are known if they are in the ns argument, or if
                    they are managed by this object.
        ns      --  A dictionary of Parameters, indexed by name, that are used
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
                raise ValueError("The parameter '%s' cannot be found"%name)

        if con in self._parameters.keys():
            con = self._parameters[con]

        if par.const:
            raise ValueError("The parameter '%s' is constant"%par)

        # This will pass the value of a constrained parameter to the initial
        # value of a parameter constraint.
        if con in self._parameters.values(): 
            val = con.getValue()
            if val is None:
                val = par.getValue()
                con.setValue(val)

        if par in self._parameters.values():
            self.fixVar(par)

        RecipeOrganizer.constrain(self, par, con, ns)
        return

    def __getTagSet(self, tags):
        """Get all variables with the given tags."""
        tagset = set()
        for tag in tags:
            tagset.update( self._tagdict.get(tag, []) )
        return tagset

    def getValues(self):
        """Get the current values of the variables in a list."""
        return array([v.getValue() for v in self._parameters.values() if v
            not in self._fixed])

    def getNames(self):
        """Get the names of the variables in a list."""
        return [par.name for par in self._parameters.values() if par not in
                self._fixed]

    def getBounds(self):
        """Get the bounds on variables in a list.

        Returns a list of (lb, ub) pairs, where lb is the lower bound and ub is
        the upper bound.

        """
        return [v.bounds for v in self._parameters.values() if v not in
                self._fixed]

    def getBounds2(self):
        """Get the bounds on variables in two lists.

        Returns lower- and upper-bound lists of variable bounds.

        """
        bounds = self.getBounds()
        lb = array([b[0] for b in bounds])
        ub = array([b[1] for b in bounds])
        return lb, ub

    def boundsToRestraints(self, prefactor = 1, scaled = False):
        """Turn all bounded parameters into restraints.

        The bounds become limits on the restraint.

        prefactor   --  prefactor for each parameter (scalar or iterable)
        scaled      --  Scale the restraints, see restrain.

        """
        pars = self._parameters.values()
        if not hasattr(prefactor, "__iter__"):
            prefactor = [prefactor] * len(pars)
        for par, x in zip(pars, prefactor):
            self.restrain(par, par.bounds[0], par.bounds[1], prefactor = x,
                    scaled = scaled)
        return

    def __applyValues(self, p):
        """Apply variable values to the variables."""
        if len(p) == 0: return
        vars = [v for v in self._parameters.values() if v not in self._fixed]
        for var, pval in zip(vars, p):
            var.setValue(pval)
        return

    def _updateConfiguration(self):
        """Notify RecipeContainers in hierarchy of configuration change."""
        self._ready = False
        return

# version
__id__ = "$Id$"

#
# End of file
