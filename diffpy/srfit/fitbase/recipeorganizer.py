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
"""Base classes and tools for constructing a FitRecipe.

RecipeContainer is the base class for organizing Parameters, and other
RecipeContainers.  RecipeOrganizer is an extended RecipeContainer that
incorporates equation building, constraints and Restraints.  equationFromString
creates an Equation instance from a string.

"""
__all__ = [ "RecipeContainer", "RecipeOrganizer", "equationFromString"]

from numpy import inf
from itertools import chain, ifilter
import re

from .constraint import Constraint
from .restraint import Restraint
from .parameter import Parameter

from diffpy.srfit.util.observable import Observable
from diffpy.srfit.equation import Equation
from diffpy.srfit.equation.builder import EquationFactory
from diffpy.srfit.util.nameutils import validateName
from diffpy.srfit.util.ordereddict import OrderedDict
from diffpy.srfit.interface import _recipeorganizer_interface

class RecipeContainer(Observable):
    """Base class for organizing pieces of a FitRecipe.

    RecipeContainers are hierarchical organizations of Parameters and other
    RecipeContainers. This class provides attribute-access to these contained
    objects.  Parameters and other RecipeContainers can be found within the
    hierarchy with the _locateManagedObject method.

    A RecipeContainer can manage dictionaries for that store various objects.
    These dictionaries can be added to the RecipeContainer using the _manage
    method. RecipeContainer methods that add, remove or retrieve objects will
    work with any managed dictionary. This makes it easy to add new types of
    objects to be contained by a RecipeContainer. By default, the
    RecipeContainer is configured to manage an OrderedDict of Parameter
    objects.

    RecipeContainer is an Observable, and observes its managed objects and
    Parameters. This allows hierarchical calculation elements, such as
    ProfileGenerator, to detect changes in Parameters and Restraints on which
    it may depend.

    Attributes
    name            --  A name for this RecipeContainer. Names should be unique
                        within a RecipeContainer and should be valid attribute
                        names.
    _parameters     --  A managed OrderedDict of contained Parameters.
    __managed       --  A list of managed dictionaries. This is used for
                        attribute access, addition and removal.
    _configobjs     --  A set of configurable objects that must know of
                        configuration changes within this object.

    """

    def __init__(self, name):
        Observable.__init__(self)
        validateName(name)
        self.name = name
        self._parameters = OrderedDict()
        self._configobjs = set()

        self.__managed = []
        self._manage(self._parameters)

        return

    def _manage(self, d):
        """Manage a dictionary of objects.

        This adds the dictionary to the __managed list. Dictionaries in
        __managed are used for attribute access, addition, and removal.

        """
        self.__managed.append(d)
        return

    def _iterManaged(self):
        """Get iterator over managed objects."""
        return chain(*(d.values() for d in self.__managed))

    def iterPars(self, name = ".", recurse = True):
        """Iterate over Parameters.
        
        name    --  Select parameters with this name (regular expression,
                    default "."). 
        recurse --  Recurse into managed objects (default True)

        """
        for par in self._parameters.itervalues():
            if re.match(name, par.name):
                yield par

        if not recurse:
            return

        # Iterate over objects within the managed dictionaries.
        managed = self.__managed[:]
        managed.remove(self._parameters)
        for m in managed:
            for obj in m.values():
                if hasattr(obj, "iterPars"):
                    for par in obj.iterPars(name = name):
                        yield par

        return

    def __iter__(self):
        """Iterate over top-level parameters."""
        return self._parameters.itervalues()

    def __getattr__(self, name):
        """Gives access to the contained objects as attributes."""
        arg = self.get(name)
        if arg is None:
            raise AttributeError(name)
        return arg
     
    def __delattr__(self, name):
        if name in self._parameters:
            self._removeParameter( self._parameters[name] )
            return

        m = self.get(name)
        if m is not None:
            raise AttributeError("Cannot delete '%s'"%name)

        super(RecipeContainer, self).__delattr__(name)
        return

    def get(self, name, default = None):
        """Get a managed object."""
        for d in self.__managed:
            arg = d.get(name)
            if arg is not None:
                return arg

        return default

    def _addObject(self, obj, d, check = True):
        """Add an object to a managed dictionary.

        obj     --  The object to be stored.
        d       --  The managed dictionary to store the object in.
        check   --  If True (default), a ValueError is raised an object of the
                    given name already exists.

        Raises ValueError if the object has no name.
        Raises ValueError if the object has the same name as some other managed
        object.

        """

        # Check name
        if not obj.name:
            message = "%s has no name" % par.__class__.__name__
            raise ValueError(message)

        # Check for extant object in d with same name
        oldobj = d.get(obj.name)
        if check and oldobj is not None:
            message = "%s with name '%s' already exists"%\
                    (obj.__class__.__name__, obj.name)
            raise ValueError(message)

        # Check for object with same name in other dictionary.
        if oldobj is None and self.get(obj.name) is not None:
            message = "Non-%s with name '%s' already exists"%\
                    (obj.__class__.__name__, obj.name)
            raise ValueError(message)

        # Detach the old object, if there is one
        if oldobj is not None:
            oldobj.removeObserver(self._flush)

        # Add the object
        d[obj.name] = obj

        # Observe the object
        obj.addObserver(self._flush)

        # Store this as a configurable object
        if hasattr(obj, "_updateConfiguration"):
            self._configobjs.add(obj)

        return

    def _removeObject(self, obj, d):
        """Remove an object from a managed dictionary.
        
        Raises ValueError if obj is not part of the dictionary.

        """
        if obj not in d.values():
            m = "'%s' is not part of the %s" % (obj, self.__class__.__name__)
            raise ValueError(m)

        del d[obj.name]
        obj.removeObserver(self._flush)

        return

    def _locateManagedObject(self, obj):
        """Find the location a managed object within the hierarchy.

        obj     --  The object to find.

        Returns a list of objects. Each entry in the list contains the next
        entry. The last object is obj, if it can be found, otherwise, the list
        is empty.

        """
        loc = [self]

        # This handles the case that an object is asked to locate itself.
        if obj is self:
            return loc

        for m in self._iterManaged():

            # Check locally for the object
            if m is obj:
                loc.append(obj)
                return loc

            # Check within managed objects
            if hasattr(m, "_locateManagedObject"):

                subloc = m._locateManagedObject(obj)
                if subloc:
                    return loc + subloc

        return []

    def _flush(self, other):
        """Invalidate cached state.

        This will force any observer to invalidate its state. By default this
        does nothing.

        """
        self.notify()
        return

    def _updateConfiguration(self):
        """Notify RecipeContainers in hierarchy of configuration change."""
        for obj in self._configobjs:
            obj._updateConfiguration()
        return


# End class RecipeContainer

class RecipeOrganizer(_recipeorganizer_interface, RecipeContainer):
    """Extended base class for organizing pieces of a FitRecipe.

    This class extends RecipeContainer by organizing constraints and
    Restraints, as well as Equations that can be used in Constraint and
    Restraint equations.  These constraints and Restraints can be placed at any
    level and a flattened list of them can be retrieved with the
    _getConstraints and _getRestraints methods. 

    Attributes
    name            --  A name for this organizer. Names should be unique
                        within a RecipeOrganizer and should be valid attribute
                        names.
    _calculators    --  A managed dictionary of Calculators, indexed by name.
    _parameters     --  A managed OrderedDict of contained Parameters.
    _constraints    --  A dictionary of Constraints, indexed by the constrained
                        Parameter. Constraints can be added using the
                        'constrain' method.
    _restraints     --  A set of Restraints. Restraints can be added using the
                        'restrain' or 'confine' methods.
    _eqfactory      --  A diffpy.srfit.equation.builder.EquationFactory
                        instance that is used create Equations from string.

    Raises ValueError if the name is not a valid attribute identifier

    """

    def __init__(self, name):
        RecipeContainer.__init__(self, name)
        self._restraints = set()
        self._constraints = {}
        self._eqfactory = EquationFactory()

        self._calculators = {}
        self._manage(self._calculators)
        return

    # Parameter management

    def _newParameter(self, name, value, check=True):
        """Add a new Parameter to the container.

        This creates a new Parameter and adds it to the container using the
        _addParameter method.

        Returns the Parameter.

        """
        p = Parameter(name, value)
        self._addParameter(p, check)
        return p

    def _addParameter(self, par, check=True):
        """Store a Parameter.

        Parameters added in this way are registered with the _eqfactory.

        par     --  The Parameter to be stored.
        check   --  If True (default), a ValueError is raised a Parameter of
                    the specified name has already been inserted.

        Raises ValueError if the Parameter has no name.
        Raises ValueError if the Parameter has the same name as a contained
        RecipeContainer.

        """

        # Store the Parameter
        RecipeContainer._addObject(self, par, self._parameters, check)

        # Register the Parameter
        self._eqfactory.registerArgument(par.name, par)
        return

    def _removeParameter(self, par):
        """Remove a parameter.

        This de-registers the Parameter with the _eqfactory. The Parameter will
        remain part of built equations.

        Note that constraints and restraints involving the Parameter are not
        modified.
        
        Raises ValueError if par is not part of the RecipeOrganizer.

        """
        self._removeObject(par, self._parameters)
        self._eqfactory.deRegisterBuilder(par.name)
        return

    def registerCalculator(self, f, argnames = None):
        """Register a Calculator so it can be used within equation strings.

        A Calculator is an elaborate function that can organize Parameters.
        This creates a function with this class that can be used within string
        equations. The resulting equation can be used in a string with
        arguments like a function or without, in which case the values of the
        Parameters created from argnames will be be used to compute the value.

        f           --  The Calculator to register.
        argnames    --  The names of the arguments to f (list or None). 
                        If this is None, then the argument names will be
                        extracted from the function.

        """
        self._eqfactory.registerOperator(f.name, f)
        self._addObject(f, self._calculators)
        # Register arguments of the calculator
        if argnames is None:
            func_code = f.__call__.im_func.func_code
            argnames = list(func_code.co_varnames)
            argnames = argnames[1:func_code.co_argcount]

        for pname in argnames:
            if pname not in self._eqfactory.builders:
                par = self._newParameter(pname, 0)
            else:
                par = self.get(pname)
            f.addLiteral(par)

        # Now return an equation object
        eq = self._eqfactory.makeEquation(f.name)
        return eq

    def registerFunction(self, f, name = None, argnames = None):
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

        Note that name and argnames can be extracted from regular python
        functions (of type 'function'), bound class methods and callable
        classes.

        Raises TypeError if name or argnames cannot be automatically
        extracted.
        Raises TypeError if an automatically extracted name is '<lambda>'.
        Raises ValueError if f is an Equation object and name is None.

        Returns the callable Equation object.

        """

        # If the function is an equation, we treat it specially. This is
        # required so that the objects observed by the root get observed if the
        # Equation is used within another equation. It is assumed that a plain
        # function is not observable.
        if isinstance(f, Equation):
            if name is None:
                m = "Equation must be given a name"
                raise ValueError(m)
            self._eqfactory.registerOperator(name, f)
            return f

        #### Introspection code
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
                raise ValueError(m)

            # Extract the name
            if name is None:
                name = func_code.co_name
                if name == '<lambda>':
                    m = "You must supply a name name for a lambda function"
                    raise ValueError(m)

            # Extract the arguments
            if argnames is None:
                argnames = list(func_code.co_varnames)
                argnames = argnames[offset:func_code.co_argcount]

        #### End introspection code

        # Make missing Parameters
        for pname in argnames:
            if pname not in self._eqfactory.builders:
                self._newParameter(pname, 0)

        # Initialize and register
        from .calculator import Calculator
        if isinstance(f, Calculator):
            for pname in argnames:
                par = self.get(pname)
                f.addLiteral(par)
            self._eqfactory.registerOperator(name, f)
        else:
            self._eqfactory.registerFunction(name, f, argnames)

        # Now we can create the Equation and return it to the user.
        eq = self._eqfactory.makeEquation(name)

        return eq

    def registerStringFunction(self, fstr, name, ns = {}):
        """Register a string function.

        This creates a function with this class that can be used within string
        equations.  The resulting equation does not require the arguments to be
        passed in the function string, as this will be handled automatically.

        fstr        --  A string equation to register.
        name        --  The name of the function to be used in equations. 
        ns          --  A dictionary of Parameters, indexed by name, that are
                        used in fstr, but not part of the FitRecipe (default
                        {}).

        Raises ValueError if ns uses a name that is already used for another
        managed object.
        Raises ValueError if the function name is the name of another managed
        object.

        Returns the callable Equation object.

        """

        # Build the equation instance.
        eq = equationFromString(fstr, self._eqfactory, ns = ns, buildargs =
                True)
        eq.name = name

        # Register any new Parameters.
        for par in self._eqfactory.newargs:
            self._addParameter(par)

        # Register the equation as a callable function. 
        argnames = eq.argdict.keys()
        return self.registerFunction(eq, name, argnames)

    def evaluateEquation(self, eqstr, ns = {}):
        """Evaluate a string equation.

        eqstr   --  A string equation to evaluate. The equation is evaluated at
                    the current value of the registered Parameters.
        ns      --  A dictionary of Parameters, indexed by name, that are
                    used in fstr, but not part of the FitRecipe (default {}).

        Raises ValueError if ns uses a name that is already used for a
        variable.

        """
        eq = equationFromString(eqstr, self._eqfactory, ns)
        return eq()

    def constrain(self, par, con, ns = {}):
        """Constrain a parameter to an equation.

        Note that only one constraint can exist on a Parameter at a time.

        par     --  The name of a Parameter or a Parameter to constrain.
        con     --  A string representation of the constraint equation or a
                    Parameter to constrain to.  A constraint equation must
                    consist of numpy operators and "known" Parameters.
                    Parameters are known if they are in the ns argument, or if
                    they are managed by this object.
        ns      --  A dictionary of Parameters, indexed by name, that are used
                    in the eqstr, but not part of this object (default {}).

        Raises ValueError if ns uses a name that is already used for a
        variable.
        Raises ValueError if par is a string but not part of this object or in
        ns.
        Raises ValueError if eqstr depends on a Parameter that is not part of
        this object and that is not defined in ns.
        Raises ValueError if par is marked as constant.

        """
        if isinstance(par, str):
            name = par
            par = self.get(name)
            if par is None:
                par = ns.get(name)

        if par is None:
            raise ValueError("The parameter cannot be found")

        if par.const:
            raise ValueError("The parameter '%s' is constant"%par)

        if isinstance(con, str):
            eq = equationFromString(con, self._eqfactory, ns)
        else:
            eq = Equation(root = con)

        eq.name = "_constraint_%s"%par.name

        # Make and store the constraint
        con = Constraint()
        con.constrain(par, eq)
        self._constraints[par] = con

        # Our configuration changed
        self._updateConfiguration()

        return

    def unconstrain(self, par):
        """Unconstrain a Parameter.

        par     --  The Parameter to unconstrain.

        This removes any constraints on a Parameter. This does nothing if the
        Parameter is not constrained.

        """
        if par in self._constraints:
            self._constraints[par].unconstrain()
            del self._constraints[par]

            # Our configuration changed
            self._updateConfiguration()

        return

    def getConstrainedPars(self, recurse = False):
        """Get a list of constrained managed Parameters in this object.

        recurse --  Recurse into managed objects and retrive their constrained
                    Parameters as well (default False).

        """
        const = self._getConstraints(recurse)
        return const.keys()

    def clearConstraints(self):
        """Clear all constraints."""
        for par in self._parameters:
            self.unconstrain(par)
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

        # Make and store the restraint
        res = Restraint()
        res.restrain(eq, lb, ub, prefactor, power, scaled)
        self._restraints.add(res)

        # Our configuration changed. Notify observers.
        self._updateConfiguration()

        return res

    def unrestrain(self, res):
        """Remove a Restraint from the RecipeOrganizer.
        
        res     --  A Restraint returned from the 'restrain' method.

        """
        if res in self._restraints:
            self._restraints.remove(res)

            # Our configuration changed
            self._updateConfiguration()

        return

    def clearRestraints(self):
        """Clear all restraints."""
        for res in self._restraints:
            self.unrestrain(res)
        return

    def _getConstraints(self, recurse = True):
        """Get the constrained Parameters for this and managed sub-objects."""
        constraints = {}
        if recurse:
            f = lambda m : hasattr(m, "_getConstraints")
            for m in ifilter(f, self._iterManaged()):
                constraints.update( m._getConstraints(recurse) )
        
        constraints.update( self._constraints)

        return constraints

    def _getRestraints(self, recurse = True):
        """Get the Restraints for this and embedded ParameterSets.
        
        This returns a set of Restraint objects.

        """
        restraints = set(self._restraints)
        if recurse:
            f = lambda m : hasattr(m, "_getRestraints")
            for m in ifilter(f, self._iterManaged()):
                restraints.update( m._getRestraints(recurse) )

        return restraints

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
