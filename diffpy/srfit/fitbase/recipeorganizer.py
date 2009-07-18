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
"""The RecipeOrganizer class and equationFromString method.

RecipeOrganizer is the base class for various classes within the FitRecipe
hierarchy. equationFromString creates an Equation instance from a string. It
checks for specific conditions on the string, as defined in the method.

"""

from numpy import inf

from .constraint import Constraint
from .restraint import Restraint, BoundsRestraint
from .parameter import Parameter
from .utils import validateName

from diffpy.srfit.equation.builder import EquationFactory
from diffpy.srfit.equation import Clicker
from diffpy.srfit.equation import Equation

class RecipeOrganizer(object):
    """A base class for organizing pieces of a FitRecipe.

    RecipeOrganizers are hierarchical organizations of Parameters, Constraints,
    Restraints and other RecipeOrganizers. This class is used throughout the
    hierarcy of a FitRecipe and provides attributes and members that help
    organize these objects at any level of the hierarcy.

    Contained parameters and other RecipeOrganizers can be accessed by name as
    attributes in order to facilitate multi-level constraints and restraints.
    These constraints and restraints can be placed at any level and a flattened
    list of them can be retrieved with the _getConstraints and _getRestraints
    methods. Parameters and other organizers can be found within the hierarchy
    with the _locateChild method.

    Attributes
    clicker         --  A Clicker instance for recording changes in contained
                        Parameters and RecipeOrganizers.
    name            --  A name for this organizer. Names should be unique
                        within a RecipeOrganizer and should be valid attribute
                        names.
    _constraints    --  A dictionary of Constraints, indexed by the constrained
                        Parameter. Constraints can be added using the
                        'constrain' method.
    _organizers     --  A list of RecipeOrganizers that this RecipeOrganizer
                        knows about.
    _orgdict        --  A dictionary containing the Parameters and
                        RecipeOrganizers indexed by name.
    _parameters     --  A list of parameters that this RecipeOrganizer knows
                        about.
    _restraints     --  A set of Restraints. Restraints can be added using the
                        'restrain' or 'confine' methods.
    _eqfactory      --  A diffpy.srfit.equation.builder.EquationFactory
                        instance that is used to create constraints and
                        restraints from string
    _equations      --  Equation instances built with the _eqfactory
                        (dictionary).

    Raises ValueError if the name is not a valid attribute identifier

    """

    def __init__(self, name):
        validateName(name)
        self.name = name
        self.clicker = Clicker()
        self._organizers = []
        self._parameters = []
        self._constraints = {}
        self._restraints = set()
        self._orgdict = {}
        self._eqfactory = EquationFactory()
        self._equations = {}
        return

    def __getattr__(self, name):
        """Gives access to the contained objects as attributes."""
        arg = self._orgdict.get(name)
        if arg is None:
            raise AttributeError(name)
        return arg

    def registerCalculator(self, f, argnames = None, makepars = True):
        """Register a Calculator so it can be used within equation strings.

        A Calculator is an elaborate function that can organize parameters.
        This creates a function with this class that can be used within string
        equations.  The resulting equation does not require the arguments to be
        passed in the equation string, as this will be handled automatically.

        f           --  The Calculator to register.
        argnames    --  The names of the arguments to f (list or None). 
                        If this is None, then the argument names will be
                        extracted from the Calculator.
        makepars    --  Flag indicating whether to make parameters from the
                        argnames if they don't already appear in this object
                        (default True). Parameters created in this way are
                        added using the _newParameter method.  If makepars is
                        False , then it is necessary that the parameters are
                        already part of this object in order to make the
                        function.

        Raises TypeError if argnames cannot be automatically extracted.
        Raises AttributeError if makepars is False and the parameters are not
        part of this object.

        """
        self._addOrganizer(f)
        self.registerFunction(f, f.name, argnames, makepars)
        return

    def registerFunction(self, f, name = None, argnames = None, 
            makepars = True):
        """Register a function so it can be used within equation strings.

        This creates a function with this class that can be used within string
        equations.  The resulting equation does not require the arguments to be
        passed in the equation string, as this will be handled automatically.

        f           --  The callable to register. If this is an Equation
                        instance, then all that needs to be provied is a name.
        name        --  The name of the function to be used in equations. If
                        this is None (default), the method will try to
                        determine the name of the function automatically.
        argnames    --  The names of the arguments to f (list or None). 
                        If this is None, then the argument names will be
                        extracted from the function.
        makepars    --  Flag indicating whether to make parameters from the
                        argnames if they don't already appear in this object
                        (default True). Parameters created in this way are
                        added using the _newParameter method.  If makepars is
                        False , then it is necessary that the parameters are
                        already part of this object in order to make the
                        function.

        Note that name and argnames can be extracted from regular python
        functions (of type 'function'), bound class methods and callable
        classes.

        Raises TypeError if name or argnames cannot be automatically
        extracted.
        Raises TypeError if an automatically extracted name is '<lambda>'.
        Raises AttributeError if makepars is False and the parameters are not
        part of this object.
        Raises ValueError if f is an Equation object and name is None.

        Returns the callable Equation object.

        """

        if isinstance(f, Equation):
            if name is None:
                m = "Equation must be given a name"
                raise ValueError(m)

            self._swapAndRegister(name, f)
            return f

        # Extract the name and argnames if necessary
        if name is None or argnames is None:

            import inspect

            func_code = None

            # This will let us offset the argument list to eliminate 'self'
            offset = 0

            # check regular functions
            if inspect.isfunction(f):
                func_code = f.func_code
            # check class method
            elif inspect.ismethod(f):
                    func_code = f.im_func.func_code
                    offset = 1
            # check functor
            elif hasattr(f, "__call__") and hasattr(f.__call__, 'im_func'):
                    func_code = f.__call__.im_func.func_code
                    offset = 1
            else:
                m = "Cannot extract name or argnames"
                raise TypeError(m)

            # Extract the name
            if name is None:
                name = func_code.co_name
                if name == '<lambda>':
                    m = "You must supply a name name for a lambda function"
                    raise TypeError(m)

            # Extract the arguments
            if argnames is None:
                argnames = list(func_code.co_varnames)
                argnames = argnames[offset:func_code.co_argcount]

        if makepars:
            for pname in argnames:
                if pname not in self._eqfactory.builders:
                    self._newParameter(pname, 0)

        # In order to make the function callable without plugging in the
        # arguments, we will build an Equation instance and register it. We'll
        # build it in another EquationFactory so we don't clutter up our
        # namespace.
        factory = EquationFactory()

        for pname in argnames:
            parbuilder = self._eqfactory.builders.get(pname)
            if parbuilder is None:
                m = "Function requires unspecified parameters (%s)."%pname
                raise AttributeError(m)

            factory.registerBuilder(pname, parbuilder)

        factory.registerFunction(name, f, len(argnames))

        argstr = ",".join(argnames)
        eq = equationFromString("%s(%s)"%(name,argstr), factory)

        self._swapAndRegister(name, eq)

        return eq

    def registerStringFunction(self, fstr, name, makepars = True, ns = {}):
        """Register a string function.

        This creates a function with this class that can be used within string
        equations.  The resulting equation does not require the arguments to be
        passed in the function string, as this will be handled automatically.

        fstr        --  A string equation to register.
        name        --  The name of the function to be used in equations. 
        makepars    --  Flag indicating whether to make parameters from the
                        arguments of fstr if they don't already appear in this
                        object (default True). Parameters created in this way
                        are added to the fitcontribution using the _newParameter
                        method.  If makepars is False , then it is necessary
                        that the parameters are already part of this object in
                        order to make the function.
        ns          --  A dictionary of Parameters, indexed by name, that are
                        used in fstr, but not part of the FitRecipe (default
                        {}).

        Raises ValueError if ns uses a name that is already used for a
        variable.
        Raises AttributeError if makepars is False and the parameters are not
        part of this object.

        Returns the callable Equation object.

        """

        # Build the equation instance.
        eq = equationFromString(fstr, self._eqfactory, buildargs = makepars)

        # Register any new Parameters. Note that these are ParameterReferences,
        # so we must create an actual Parameter to register.
        for par in self._eqfactory.newargs:
            self._addParameter(par)

        # Register the equation by name and do any necessary swapping.
        self._swapAndRegister(name, eq)

        return eq

    def evaluateEquation(self, eqstr, ns = {}):
        """Evaluate a string equation.

        eqstr   --  A string equation to evaluate. The equation is evaluated at
                    the current value of the registered parameters. The
                    equation can be specified as described in the setEquation
                    method.
        ns      --  A dictionary of Parameters, indexed by name, that are
                    used in fstr, but not part of the FitRecipe (default {}).

        Raises ValueError if ns uses a name that is already used for a
        variable.

        """
        eq = equationFromString(eqstr, self._eqfactory, ns = {})
        return eq()

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
        Raises ValueError if par is marked as constant.

        """
        if par.const:
            raise ValueError("The parameter '%s' is constant"%par)

        if isinstance(con, str):
            eq = equationFromString(con, self._eqfactory, ns)
        else:
            eq = Equation(root = con)


        # Make and store the constraint
        con = Constraint()
        con.constrain(par, eq)
        self._constraints[par] = con

        # Store this equation for future swapping
        self._equations[repr(con)] = eq
        return

    def unconstrain(self, par):
        """Unconstrain a Parameter.

        par     --  The Parameter to unconstrain.

        This removes any constraints on a parameter.

        """
        if par in self._constraints:
            con = self._constraints[par]
            del self._equations[repr(con)]
            del self._constraints[par]
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
                    in the eqstr, but not part of the RecipeOrganizer 
                    (default {}).

        The penalty is calculated as 
        prefactor * max(0, lb - val, val - ub) ** power
        and val is the value of the calculated equation. This is multipled by
        the average chi^2 if scaled is True.

        Raises ValueError if ns uses a name that is already used for a
        Parameter.
        Raises ValueError if eqstr depends on a Parameter that is not part of
        the RecipeOrganizer and that is not defined in ns.

        Returns the Restraint object for use with the 'unrestrain' method.

        """

        if isinstance(res, str):
            eq = equationFromString(res, self._eqfactory, ns)
        else:
            eq = Equation(root = res)

        # Make and store the constraint
        res = Restraint()
        res.restrain(eq, lb, ub, prefactor, power, scaled)
        self._restraints.add(res)

        self._equations[repr(res)] = eq

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
        the RecipeOrganizer and that is not defined in ns.
        Raises ValueError if lb == ub.

        Returns the BoundsRestraint object for use with the 'unrestrain' method.

        """

        if lb == ub:
            m = "Bounds must be different"
            raise ValueError(m)

        if isinstance(res, str):
            eq = equationFromString(res, self._eqfactory, ns)
        else:
            eq = Equation(root = res)
            val = res.getValue()
            # Reset the value of the parameter if it is out of bounds.
            # Optimizers might choke on getting an infinte residual right out
            # of the gate.
            if val < lb or val > ub:
                if lb is not -inf:
                    res.setValue(lb)
                elif ub is not inf:
                    res.setValue(ub)
                else:
                    res.setValue(0)

        # Make and store the constraint
        res = BoundsRestraint()
        res.confine(eq, lb, ub)
        self._restraints.add(res)

        self._equations[repr(res)] = eq
        return res

    def unrestrain(self, res):
        """Remove a Restraint or BoundsRestraint from the RecipeOrganizer.
        
        res     --  A Restraint returned from the 'restrain' method.

        """
        if res in self._restraints:
            del self._equations[repr(res)]
            self._restraints.remove(res)

        return

    def _addParameter(self, par, check=True):
        """Store a Parameter.

        Parameters added in this way are registered with the _eqfactory.

        par     --  The Parameter to be stored.
        check   --  If True (default), a ValueError is raised a Parameter or
                    RecipeOrganizer of the specified name has already been
                    inserted.

        Raises ValueError if the Parameter has no name.

        """
        
        message = ""
        if not par.name:
            message = "Parameter has no name"%par

        if check:
            if par.name in self._orgdict:
                message = "Object with name '%s' already exists"%par.name

        if message:
            raise ValueError(message)

        # Swap parameters
        oldpar = self._orgdict.get(par.name)
        if oldpar is not None:
            self._swapEquationObject(oldpar, par)
            self._removeParameter(oldpar)
        
        self._orgdict[par.name] = par
        self._parameters.append(par)
        self._eqfactory.registerArgument(par.name, par)
        self.clicker.addSubject(par.clicker)

        return

    def _newParameter(self, name, value, check=True):
        """Add a new parameter that can be used in the equation.

        Returns the parameter.

        """
        p = Parameter(name, value)
        self._addParameter(p, check)
        return p

    def _removeParameter(self, par):
        """Remove a parameter.
        
        raises ValueError if par is not part of the RecipeOrganizer.

        """
        if par not in self._parameters:
            m = "'%s' is not part of the %s" % (par, self.__class__.__name__)
            raise ValueError(m)

        self._parameters.remove(par)
        del self._orgdict[par.name]
        self._eqfactory.deRegisterBuilder(par.name)
        self.clicker.removeSubject(par.clicker)
        return

    def _addOrganizer(self, org, check=True):
        """Store a RecipeOrganizer.

        org     --  The RecipeOrganizer to be stored.
        check   --  If True (default), an ValueError is raised if the
                    RecipeOrganizer has an invalid name, or if a Parameter or
                    RecipeOrganizer of that name has already been inserted.

        """
        if check:
            message = ""
            if not org.name:
                message = "RecipeOrganizer has no name"%org
            elif org.name in self._orgdict:
                message = "Object with name '%s' already exists"%org.name

            if message:
                raise ValueError(message)

        self._orgdict[org.name] = org
        self._organizers.append(org)
        self.clicker.addSubject(org.clicker)
        return

    def _removeOrganizer(self, org):
        """Remove an organizer.
        
        raises ValueError if organizer is not part of the RecipeOrganizer.

        """
        if org not in self._organizers:
            m = "'%s' is not part of the %s" % (org, self.__class__.__name__)
            raise ValueError(m)

        self._organizers.remove(org)
        del self._orgdict[org.name]
        self.clicker.removeSubject(org.clicker)

        return

    def _getConstraints(self):
        """Get the Constraints for this and embedded ParameterSets."""
        constraints = {}
        for org in self._organizers:
            constraints.update( org._getConstraints() )
        # Update with local constraints last. These override the others.
        constraints.update(self._constraints)

        return constraints

    def _getRestraints(self):
        """Get the Restraints for this and embedded ParameterSets."""
        restraints = set(self._restraints)
        for org in self._organizers:
            restraints.update( org._getRestraints() )

        return restraints

    def _locateChild(self, obj, loc = None):
        """Find the location of a parameter or organizer within the organizer

        obj     --  The Parameter or RecipeOrganizer to locate
        loc     --  A list containing the path to the object. The 
                    The name of this RecipeOrganizer gets appended to the list,
                    which gets passed on util the parameter is located. If the
                    parameter is not located herein, the name of this
                    RecipeOrganizer is not appended.  This defaults to None, in
                    which case a new list is created and passed along.

        Returns a list of objects. Each entry in the list is the object
        containing the next object in the list. The last object is obj, if it
        can be found, otherwise, the list is empty.

        """
        if loc is None:
            loc = []
        loc.append(self)

        loclen = len(loc)

        if obj in self._orgdict.itervalues():
            loc.append(obj)
            return loc

        for org in self._organizers:
            org._locateChild(obj, loc)
            if len(loc) > loclen:
                break;

        if len(loc) == loclen:
            loc.pop()

        return loc

    def _iterPars(self):
        """Iterate over parameters."""
        for par in self._parameters:
            yield par

        for org in self._organizers:
            org.iterParameters()

        return

    def _swapEquationObject(self, oldobj, newobj):
        """Swap all instances of an Equation object for another.

        This method is used to curate Equations belonging to this object. If
        one overloads an object that can appear in an Equation, it is the
        responsibility of this method to replace all instances of the original
        object with the new object.

        Note that this method is likely to be slow if many Equations are
        registered.
        
        """
        for eq in self._equations.values():
            eq.swap(oldobj, newobj)
        return

    def _swapAndRegister(self, name, newobj):
        """Swap and register an equation.

        This registers an equation, and performs any swapping necessary with
        the new equation.
        
        """
        self._eqfactory.registerEquation(name, newobj)
        oldobj = self._equations.get(name)
        if oldobj is not None:
            self._swapEquationObject(oldobj, newobj)
        self._equations[name] = newobj
        return

# End RecipeOrganizer

def equationFromString(eqstr, factory, ns = {}, buildargs = False,
        argclass = Parameter, argkw = {}):
    """Make an equation from a string.

    eqstr   --  A string representation of the equation. The equation must
                consist of numpy operators and "known" Parameters. Parameters
                are known if they are in ns, or already defined in the factory.
    factory --  An EquationFactory instance.
    ns      --  A dictionary of Parameters indexed by name that are used
                in the eqstr but not already defined in the factory 
                (default {}).
    buildargs   --  A flag indicating whether missing Parameters can be created
                by the Factory (default False). If False, then the a ValueError
                will be raised if there are undefined arguments in the eqstr. 
    argclass    --  Class to use when creating new Arguments (default
                Parameter). The class constructor must accept the 'name' key
                word.
    argkw   --  Key word dictionary to pass to the argclass constructor
                (default {}).

    Raises ValueError if ns uses a name that is already defined in the factory.
    Raises ValueError if the equation has undefined parameters.

    """

    defined = set(factory.builders.keys())

    # Check if ns overloads any parameters.
    if defined.intersection(ns.keys()):
        raise ValueError("ns contains defined names")

    # Register the ns parameters in the equation factory
    for name, arg in ns.items():
        factory.registerArgument(name, arg)

    eq = factory.makeEquation(eqstr, buildargs, argclass, argkw)

    # Clean the ns parameters
    for name in ns:
        factory.deRegisterBuilder(name)

    return eq

__id__ = "$Id$"
