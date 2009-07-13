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

FitRecipes organize FitContributions, Parameters, Restraints and Constraints to
create a recipe of the system you wish to optimize. From the client's
perspective, the FitRecipe is a residual calculator. The residual method does
the work of updating variable values, which get propagated to the Parameters of
the underlying FitContributions via the varibles, Restraints and Constraints. As a
result, this class can be used without subclassing.

See the examples in the documentation for how to create an optimization problem
using FitRecipe.

"""

from numpy import concatenate, sqrt, inf, dot
from functools import update_wrapper

from .parameter import Parameter, ParameterProxy
from .recipeorganizer import RecipeOrganizer
from .fithook import FitHook

class FitRecipe(RecipeOrganizer):
    """FitRecipe class.

    Attributes
    clicker         --  A Clicker instance for recording changes in contained
                        Parameters and FitContributions.
    name            --  A name for this FitRecipe.
    fithook         --  An object to be called whenever within the residual
                        (default FitHook()).
    _constraintlist --  An ordered list of the constraints from this and all
                        sub-components.
    _constraints    --  A dictionary of Constraints, indexed by the constrained
                        Parameter. Constraints can be added using the
                        'constrain' method.
    _doprepare      --  A flag indicating that the the attributes need to be
                        prepared for calculation.
    _eqfactory      --  A diffpy.srfit.equation.builder.EquationFactory
                        instance that is used to create constraints and
                        restraints from strings.
    _fixed          --  A list of Parameters that are fixed, but still managed
                        by the FitRecipe.
    _organizers     --  A list of FitContributions to the recipe. Modified by the
                        addContribution method.
    _orgdict        --  A dictionary containing the Parameters and
                        FitContributions indexed by name.
    _parameters     --  A list of variable Parameters.
    _restraintlist  --  A list of restraints from this and all sub-components.
    _restraints     --  A set of Restraints. Restraints can be added using the
                        'restrain' or 'confine' methods.
    _tagdict        --  A dictionary of tags to variables.
    _weights        --  The weighing factor for each fitcontribution. This value
                        is multiplied by the residual of the fitcontribution when
                        determining the overall residual.

    """

    def __init__(self, name = "fit"):
        """Initialization."""
        RecipeOrganizer.__init__(self, name)
        self.fithook = FitHook()
        self._constraintlist = []
        self._restraintlist = []
        self._fixed = []
        self._weights = []
        self._doprepare = True
        self._tagdict = {}
        return

    def setFitHook(self, fithook):
        """Set a hook to be called within the residual method.

        The hook is an object for reportind updates, or doing whatever else. It
        must have 'precall' and 'postcall' methods, which are called at the
        start and at the end of the residual calculation. The precall method
        must accept a single argument, which is this FitRecipe object. The
        postcall method must accept the recipe and the chiv, vector residual.
        It must also have a reset method that takes no arguments, which is
        called whenver the FitRecipe is prepared for a refinement.

        See the FitHook class for the interface.

        """
        self.fithook = fithook
        self._doprepare = True
        return

    def addContribution(self, con, weight = 1.0):
        """Add a fitcontribution to the FitRecipe."""
        self._addOrganizer(con, check=True)
        self._weights.append(weight)
        self._doprepare = True
        return

    def setWeight(self, con, weight):
        """Set the weight of a fitcontribution."""
        idx = self._organizers.index(con)
        self._weights[idx] = weight
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
        fitcontribution's residual, plus the value of each restraint. The array
        returned, denoted chiv, is such that 
        dot(chiv, chiv) = chi^2 + restraints.

        """

        if self._doprepare:
            self._prepare()

        if self.fithook:
            self.fithook.precall(self)

        # Update the variable parameters.
        for i, val in enumerate(p):
            self._parameters[i].setValue(val)

        # Update the constraints. These are ordered such that the list only
        # needs to be cycled once.
        for con in self._constraintlist:
            con.update()

        # Calculate the bare chiv
        chiv = concatenate([ 
            sqrt(self._weights[i])*self._organizers[i].residual() \
                    for i in range(len(self._organizers))])

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
        """Prepare for the calculation.

        This will prepare the data attributes to be used in the residual
        calculation.

        This updates the local constraints and restraints with those of the
        contributions.

        Constraints can have inter-dependencies that require that they are
        updated in a specific order. This will set the proper order.

        Raises AttributeError if there are variables without a value.
        Raises AttributeError if there are multiple constraints on the same
        parameter defined in different places within the recipe hierarchy.

        """
        # Inform the fit hook that we're updating things
        if self.fithook:
            self.fithook.reset()

        # Check for variable values
        varvals = self.getValues()
        names = self.getNames()
        m = ""
        badvars = []
        for idx, val in enumerate(self.getValues()):
            if val is not None: continue
            badvars.append( "'%s'"%names[idx] )

        if len(badvars) == 1:
            m = "Variable %s needs an initial value" % badvars[0]
        elif len(badvars) > 0:
            s1 = ",".join(badvars[:-1])
            s2 = badvars[-1]
            m = "Variables %s and %s need initial values" % (s1, s2)

        if m:
            raise AttributeError(m)

        # Update constraints and restraints. 
        rset = set(self._restraints)
        cdict = {}
        # We let constraints closer to the FitRecipe override all others.
        # Constraints on the same parameter in different organizers cannot be
        # resolved without some guesswork, so throw an error instead.
        for con in self._organizers:
            rset.update( con._getRestraints() )
            constraints = con._getConstraints()
            pars = set(cdict.keys()).intersection( constraints.keys() )
            if pars:
                m = "There are multiple internal constraints on '%s'"%\
                        pars.pop()
                raise AttributeError(m)
            cdict.update(constraints)
        # Update with local constraints last
        cdict.update(self._constraints)

        # Reorder the constraints. Constraints are ordered such that a given
        # constraint is placed before its dependencies.

        self._restraintlist = list(rset)
        self._constraintlist = cdict.values()

        # Create a depth-1 map of the constraint dependencies
        depmap = {}
        for con in self._constraintlist:
            depmap[con] = set()
            # Now check the constraint's equation for constrained arguments
            for arg in con.eq.args:
                if arg in self._constraintlist:
                    depmap[con].add( arg )

        # Turn the dependency map into multi-level map.
        def _extendDeps(con):
            deps = set(depmap[con])
            for dep in depmap[con]:
                deps.update(_extendDeps(dep))

            return deps

        for con in depmap:
            depmap[con] = _extendDeps(con)

        # Now sort the constraints based on the dependency map
        def cmp(x, y):
            # If y is a dependency of x, then y must come first
            if y in depmap[x]:
                return 1
            # If x is a dependency of y, then x must come first
            if x in depmap[y]:
                return -1
            # Otherwise, these are equivalant
            return 0

        self._constraintlist.sort(cmp)

        self._doprepare = False
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

        Returns the variable (ParameterProxy instance).

        Raises ValueError if the name of the variable is already taken by
        another variable or a fitcontribution.
        Raises ValueError if par is constant.

        """
        name = name or par.name

        if par.const:
            raise ValueError("The parameter '%s' is constant"%par)

        var = ParameterProxy(name, par)
        if value is not None:
            var.setValue(value)

        self._addParameter(var)

        if fixed:
            self.fixVar(var)
          
        # Deal with tags
        if tag:
            self.__tagVar(tag, var)

        for t in tags:
            self.__tagVar(t, var)

        self._doprepare = True

        return var

    def __tagVar(self, tag, var):
        """Private function to tag a variable."""
        vset = self._tagdict.get(tag, set())
        vset.add(var)
        self._tagdict[tag] = vset
        return

    def delVar(self, var):
        """Remove a variable.

        var     --  A variable of the FitRecipe.

        Raises ValueError if var is not part of the FitRecipe.

        """
        if var in self._parameters:
            self._parameters.remove(var)
        elif var in self._fixed:
            self._fixed.remove(var)
        else:
            raise ValueError("'%s' is not part of the FitRecipe"%var)

        # De-register the Parameter with the equation factory
        self._eqfactory.deRegisterBuilder(var.name)

        # Remove this from the organizer dictionary
        del self._orgdict[var.name]
        self.clicker.removeSubject(var.clicker)

        # Remove tags
        for vset in self._tagdict.items():
            vset.discard(var)
        
        return

    def newVar(self, name, value, fixed = False, tag = None, tags = []):
        """Create a new variable of the fit.

        This method lets new variables be created that are not tied to a
        Parameter.  Orphan Parameters may cause a fit to fail, depending on the
        optimization routine, and therefore should only be created to be used
        in contraint or restraint equations.

        name    --  The name of the variable. The variable will be able to be
                    used by this name in restraint and constraint equations.
        value   --  An initial value for the variable. 
        fixed   --  Fix the variable so that it does not vary (default False).
        tag     --  A tag for the variable. This can be used to fix and free
                    variables by tag (default None).
        tags    --  A list of tags (default []). Both tag and tags can be
                    appolied.

        Returns the new variable (Parameter instance).

        """
        var = self._newParameter(name, value)

        if fixed:
            self.fixVar(var)

        # Deal with tags
        if tag:
            self.__tagVar(tag, var)

        for t in tags:
            self.__tagVar(t, var)

        self._doprepare = True
        return var

    def fixVar(self, var, value = None):
        """Fix a variable so that it doesn't change.

        var     --  A variable of the FitRecipe.
        value   --  A new value for the variable. If this is None
                    (default), then the value will not be changed.

        Raises ValueError if var is not part of the FitRecipe.
        
        """
        if var in self._parameters:
            self._parameters.remove(var)
            self._fixed.append(var)
        elif var in self._fixed:
            pass
        else:
            raise ValueError("'%s' is not part of the FitRecipe"%var)

        if value is not None:
            var.setValue(value)

        return

    def freeVar(self, var, value = None):
        """Free a variable so that it is refined.

        Variables are free by default.

        var     --  A variable of the FitRecipe.
        value   --  A new value for the variable. If this is None
                    (default), then the value will not be changed.

        This will disturb the order of the variables.

        Raises ValueError if var is not part of the FitRecipe.
        
        """
        if var in self._parameters:
            pass
        elif var in self._fixed:
            self._fixed.remove(var)
            self._parameters.append(var)
        else:
            raise ValueError("'%s' is not part of the FitRecipe"%var)

        if value is not None:
            var.setValue(value)

        return

    def fixAll(self, *args):
        """Fix all variables.

        Extra arguments are assumed to be tags. If present, only variables with
        the given tag will be fixed.

        """

        fixset = set()
        for tag in args:
            vset = self._tagdict.get(tag)
            if vset:
                fixset.update(vset)

        if not fixset:
            fixset = self._parameters[:]

        for var in fixset:
            self.fixVar(var)

        return

    def freeAll(self, *args):
        """Free all variables.

        Extra arguments are assumed to be tags. If present, only variables with
        the given tag will be fixed.

        This will disturb the order of the variables.
        
        """

        freeset = set()
        for tag in args:
            vset = self._tagdict.get(tag)
            if vset:
                freeset.update(vset)

        if not freeset:
            freeset = self._fixed[:]

        for var in freeset:
            self.freeVar(var)

        return

    def getValues(self):
        """Get the current values of the variables in a list."""
        return [par.getValue() for par in self._parameters]

    def getNames(self):
        """Get the names of the variables in a list."""
        return [par.name for par in self._parameters]


    # Overloaded

    def constrain(self, par, con, ns = {}):
        """Constrain a parameter to an equation.

        Note that only one constraint can exist on a Parameter at a time. The
        most recent constraint override all other user-defined constraints.
        Built-in constraints override all other constraints.

        par     --  The Parameter to constrain. It does not need to be a
                    variable.
        con     --  A string representation of the constraint equation or a
                    parameter to constrain to.  A constraint equation must
                    consist of numpy operators and "known" Parameters.
                    Parameters are known if they are in the ns argument, or if
                    they have been added to this FitRecipe with the 'add' or
                    'new' methods.
        ns      --  A dictionary of Parameters, indexed by name, that are used
                    in the eqstr, but not part of the FitRecipe (default {}).

        Raises ValueError if ns uses a name that is already used for a
        variable.
        Raises ValueError if eqstr depends on a Parameter that is not part of
        the FitRecipe and that is not defined in ns.

        """
        RecipeOrganizer.constrain(self, par, con, ns)
        self._doprepare = True
        return

    def unconstrain(self, par):
        """Unconstrain a Parameter.

        par     --  The Parameter to unconstrain.

        This removes any constraints on a parameter, including built-in
        constraints.

        """
        RecipeOrganizer.unconstrain(self, par)
        self._doprepare = True
        return

    def restrain(self, res, lb = -inf, ub = inf, prefactor = 1, power = 2,  
            scaled = False, ns = {}):
        """Restrain an expression to specified bounds

        res     --  An equation string or Parameter to restrain.
        lb      --  The lower bound on the restraint evaluation (default -inf).
        ub      --  The lower bound on the restraint evaluation (default inf).
        prefactor   --  A multiplicative prefactor for the restraint 
                        (default 1).
        power   --  The power of the penalty (default 2).
        scaled  --  A flag indicating if the restraint is scaled (multiplied)
                    by the unrestrained point-average chi^2 (chi^2/numpoints)
                    (default False).
        ns      --  A dictionary of Parameters, indexed by name, that are used
                    in the eqstr, but not part of the FitRecipe 
                    (default {}).

        The penalty is calculated as 
        prefactor * max(0, lb - val, val - ub) ** power
        and val is the value of the calculated equation. This is multipled by
        the average chi^2 if scaled is True.

        Raises ValueError if ns uses a name that is already used for a
        Parameter.
        Raises ValueError if eqstr depends on a Parameter that is not part of
        the FitRecipe and that is not defined in ns.

        Returns the Restraint selfect for use with the 'unrestrain' method.

        """
        res = RecipeOrganizer.restrain(self, res, lb, ub, prefactor, power,
                scaled, ns)
        self._doprepare = True
        return res

    def confine(self, res, lb = -inf, ub = inf, ns = {}):
        """Confine an expression to hard bounds.

        res     --  An equation string or Parameter to restrain.
        lb      --  The lower bound on the restraint evaluation (default -inf).
        ub      --  The lower bound on the restraint evaluation (default inf).
        ns      --  A dictionary of Parameters, indexed by name, that are used
                    in the eqstr, but not part of the FitRecipe 
                    (default {}).

        The penalty is infinite if the value of the calculated equation is
        outside the bounds.

        Raises ValueError if ns uses a name that is already used for a
        Parameter.
        Raises ValueError if eqstr depends on a Parameter that is not part of
        the FitRecipe and that is not defined in ns.

        Returns the BoundsRestraint object for use with the 'unrestrain' method.

        """
        res = RecipeOrganizer.confine(self, res, lb, ub, ns)
        self._doprepare = True
        return res

    def unrestrain(self, res):
        """Remove a restraint from the FitRecipe.
        
        res     --  A Restraint returned from the 'restrain' method.

        """
        RecipeOrganizer.unrestrain(self, res)
        self._doprepare = True
        return

# version
__id__ = "$Id$"

#
# End of file
