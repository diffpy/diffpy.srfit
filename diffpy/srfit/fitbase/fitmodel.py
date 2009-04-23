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
"""FitModel class. 

FitModels organize Contributions, Restraints and Constraints. The main purpose
of a FitModel is to serve as an error calculator for an exernal fitting
framework.
"""

from numpy import concatenate, sqrt, inf, dot

from .parameter import Parameter
from .modelorganizer import ModelOrganizer

class FitModel(ModelOrganizer):
    """FitModel class.

    Attributes
    clicker         --  A Clicker instance for recording changes in contained
                        Parameters and Contributions.
    name            --  A name for this FitModel.
    _aliasmap       --  A map from Parameters to their aliases.
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
                        by the FitModel.
    _organizers     --  A list of Contributions to the model. Modified by the
                        addContribution method.
    _orgdict        --  A dictionary containing the Parameters and
                        Contributions indexed by name.
    _parameters     --  A list of variable Parameters.
    _restraintlist  --  A list of restraints form this and all sub-components.
    _restraints     --  A set of Restraints. Restraints can be added using the
                        'restrain' method.
    _weights        --  The weighing factor for each contribution. This value
                        is multiplied by the residual of the contribution when
                        determining the overall residual.

    """

    def __init__(self, name = "fit"):
        """Initialization."""
        ModelOrganizer.__init__(self, name)
        self._constraintlist = []
        self._restraintlist = []
        self._fixed = []
        self._aliasmap = {}
        self._weights = []
        self._doprepare = True
        return

    def addContribution(self, con, weight = 1.0):
        """Add a contribution to the FitModel."""
        self._addOrganizer(con, check=True)
        self._weights.append(weight)
        self._doprepare = True
        return

    def setWeight(con, weight):
        """Set the weight of a contribution."""
        idx = self._organizers.index(con)
        self._weights[idx] = weight
        return

    def residual(self, p):
        """Calculate the residual to be optimized.

        Arguments
        p   --  The list of current variable values, provided in the same order
                as the '_parameters' list. If p is an empty iterable, then it
                is assumed that the parameters have already been updated in
                some other way, and the update is skipped.

        The residual is by default the weighted concatenation of each 
        contribution's residual, plus the value of each restraint. The array
        returned, denoted chiv, is such that 
        dot(chiv, chiv) = chi^2 + restraints.
        """

        if self._doprepare:
            self._prepare()

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

        return chiv

    def _prepare(self):
        """Prepare for the calculation.

        This will prepare the data attributes to be used in the residual
        calculation.

        It updates the constraints and restraints with those of its
        contributions.

        Constraints can have inter-dependencies that require that they are
        updated in a specific order. This will set the proper order.
        """
        # Update constraints and restraints.
        rset = set(self._restraints)
        cdict = dict(self._constraints)
        for con in self._organizers:
            rset.update( con._getRestraints() )
            constraints = con._getConstraints()
            pars = set(cdict.keys()).intersection( constraints.keys() )
            if pars:
                message = "There are multiple constraints on '%s'"%pars.pop()
                raise AttributeError(message)
            cdict.update(constraints)

        # Reorder the constraints. Constraints are ordered such that a given
        # constraint is placed before its dependencies.

        self._restraintlist = list(rset)
        self._constraintlist = cdict.values()

        # Create a depth-1 map of the constraint dependencies
        depmap = {}
        for con in self._constraintlist:
            depmap[con] = set()
            # Now check the constraint's equation for constrained arguments
            for arg in con.eq.args.values():
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

    def addVar(self, par, value = None, alias = None, fixed = False):
        """Add a variable to be refined.

        par     --  A Parameter that will be varied during a fit.
        value   --  An initial value for the variable. If this is None
                    (default), then the current value of par will be used.
        alias   --  A name for this variable. The variable will be able to be
                    used by this name in restraint and constraint equations.
        fixed   --  Fix the variable so that it does not vary (default False).

        """
        self.clicker.addSubject(par.clicker)

        if value is not None:
            par.setValue(value)

        if fixed:
            self._fixed.append(par)
        else:
            self._addParameter(par)

        # Record the aliases
        self._aliasmap[par] = []
        self._aliasmap[par].append(par.name)
                
        # Register the parameter aliases with the eqfactory
        if alias is not None:
            self._aliasmap[par].append(alias)
            self._eqfactory.registerArgument(alias, par)

        self._doprepare = True

        return

    def delVar(self, par):
        """Remove a variable.

        par     --  A Parameter that has already been added to the FitModel.

        Raises ValueError if par is not part of the FitModel.

        """
        if par in self._parameters:
            self._parameters.remove(par)
        elif par in self._fixed:
            self._fixed.remove(par)
        else:
            raise ValueError("'%s' is not part of the FitModel"%par)

        # De-register the Parameter with the equation factory
        self._eqfactory.deRegisterBuilder(par.name)
        for name in self._aliasmap[par]:
            self._eqfactory.deRegisterBuilder(name)

        # Remove this from the alias map and the organizer dictionary
        del self._aliasmap[par]
        del self._orgdict[par.name]
        self.clicker.removeSubject(par.clicker)
        
        return

    def newVar(self, name, value, fixed = False):
        """Create a new Parameter that will be varied in the fit.

        This method lets new Parameters be created that can be used in
        constraint and restraint equations. A Parameter created by this method
        that is not involved with a restraint or constraint is called an
        orphan. Orphan Parameters may cause a fit to fail, depending on the
        optimization routine.

        name    --  The name of the variable. The variable will be able to be
                    used by this name in restraint and constraint equations.
        value   --  An initial value for the variable. 
        fixed   --  Fix the variable so that it does not vary (default False).

        Returns the new Parameter.

        """
        par = Parameter(name, value=value)
        self.addVar(par, fixed=fixed)
        return par

    def fixVar(self, par, value = None):
        """Fix a variable so that it doesn't change.

        par     --  A Parameter that has already been added to the FitModel.
        value   --  A new value for the variable. If this is None
                    (default), then the value will not be changed.

        Raises ValueError if par is not part of the FitModel.
        
        """
        if par in self._parameters:
            self._parameters.remove(par)
            self._fixed.append(par)
        elif par in self._fixed:
            pass
        else:
            raise ValueError("'%s' is not part of the FitModel"%par)

        if value is not None:
            par.setValue(value)

        return

    def freeVar(self, par, value = None):
        """Free a variable so that it is refined.

        Note that variables are free by default.

        par     --  A Parameter that has already been added to the FitModel.
        value   --  A new value for the variable. If this is None
                    (default), then the value will not be changed.

        Raises ValueError if par is not part of the FitModel.
        
        """
        if par in self._parameters:
            pass
        elif par in self._fixed:
            self._fixed.remove(par)
            self._parameters.append(par)
        else:
            raise ValueError("'%s' is not part of the FitModel"%par)

        if value is not None:
            par.setValue(value)

        return

    # Overloaded
    def constrain(self, par, eqstr, ns = {}):
        """Constrain a parameter to an equation.

        Note that only one constraint can exist on a Parameter at a time. The
        most recent constraint override all other user-defined constraints.
        Built-in constraints override all other constraints.

        par     --  The Parameter to constrain. It does not need to be a
                    variable.
        eqstr   --  A string representation of the constraint equation. The
                    constraint equation must consist of numpy operators and
                    "known" Parameters. Parameters are known if they are in the
                    ns argument, or if they have been added to this FitModel
                    with the 'add' or 'new' methods.
        ns      --  A dictionary of Parameters, indexed by name, that are used
                    in the eqstr, but not part of the FitModel (default {}).

        Raises ValueError if ns uses a name that is already used for a
        variable.
        Raises ValueError if eqstr depends on a Parameter that is not part of
        the FitModel and that is not defined in ns.

        """
        ModelOrganizer.constrain(self, par, eqstr, ns)
        self._doprepare = True
        return

    def unconstrain(self, par):
        """Unconstrain a Parameter.

        par     --  The Parameter to unconstrain.

        This removes any constraints on a parameter, including built-in
        constraints.
        """
        ModelOrganizer.unconstrain(self, par)
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
                    in the eqstr, but not part of the FitModel 
                    (default {}).

        The penalty is calculated as 
        prefactor * max(0, lb - val, val - ub) ** power
        and val is the value of the calculated equation. This is multipled by
        the average chi^2 if scaled is True.

        Raises ValueError if ns uses a name that is already used for a
        Parameter.
        Raises ValueError if eqstr depends on a Parameter that is not part of
        the FitModel and that is not defined in ns.

        Returns the Restraint selfect for use with the 'unrestrain' method.

        """
        res = ModelOrganizer.restrain(self, res, lb, ub, prefactor, power,
                scaled, ns)
        self._doprepare = True
        return res

    def unrestrain(self, res):
        """Remove a restraint from the FitModel.
        
        res     --  A Restraint returned from the 'restrain' method.
        """
        ModelOrganizer.unrestrain(self, res)
        self._doprepare = True
        return

# version
__id__ = "$Id$"

#
# End of file
