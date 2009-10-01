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
incorporates Equations, Constraints and Restraints.  equationFromString creates
an Equation instance from a string. It checks for specific conditions on the
string, as defined in the method. 

"""

from numpy import inf
from itertools import chain, ifilter
import re

from .constraint import Constraint
from .restraint import Restraint
from .parameter import Parameter

from diffpy.srfit.util.clicker import Clicker
from diffpy.srfit.equation import Equation
from diffpy.srfit.equation.builder import EquationFactory
from diffpy.srfit.util.nameutils import validateName
from diffpy.srfit.util.ordereddict import OrderedDict

class RecipeContainer(object):
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

    Attributes
    clicker         --  A Clicker instance for recording changes in contained
                        Parameters and RecipeOrganizers.
    name            --  A name for this RecipeContainer. Names should be unique
                        within a RecipeContainer and should be valid attribute
                        names.
    _confclicker    --  A Clicker for recording configuration changes, esp.
                        additions and removal of managed objects.
    _parameters     --  A managed OrderedDict of contained Parameters.
    __managed       --  A list of managed dictionaries. This is used for
                        attribute access, addition and removal.

    """

    def __init__(self, name):
        validateName(name)
        self.name = name
        self.clicker = Clicker()
        self._confclicker = Clicker()
        self._parameters = OrderedDict()

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
                else:
                    break

        return

    def __getattr__(self, name):
        """Gives access to the contained objects as attributes."""
        arg = self.get(name)
        if arg is None:
            raise AttributeError(name)
        return arg

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
            self.clicker.removeSubject(oldobj.clicker)
            if hasattr(oldobj, "_confclicker"):
                self._confclicker.removeSubject(oldobj._confclicker)

        # Add the object
        d[obj.name] = obj

        # Attach the object to the Clicker
        self.clicker.addSubject(obj.clicker)
        if hasattr(obj, "_confclicker"):
            self._confclicker.addSubject(obj._confclicker)

        # Our configuration changed
        self._confclicker.click()

        return

    def _removeObject(self, obj, d):
        """Remove an object from a managed dictionary.
        
        Raises ValueError if obj is not part of the dictionary.

        """
        if obj not in d.values():
            m = "'%s' is not part of the %s" % (obj, self.__class__.__name__)
            raise ValueError(m)

        del d[obj.name]
        self.clicker.removeSubject(obj.clicker)
        if hasattr(obj, "_confclicker"):
            self._confclicker.removeSubject(obj._confclicker)

        # Our configuration changed
        self._confclicker.click()

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

# End class RecipeContainer

class RecipeOrganizer(RecipeContainer):
    """Extended base class for organizing pieces of a FitRecipe.

    This class extends RecipeContainer by organizing Constraints and
    Restraints, as well as Equations that can be used in Constraint and
    Restraint equations.  These Constraints and Restraints can be placed at any
    level and a flattened list of them can be retrieved with the
    _getConstraints and _getRestraints methods. 

    Attributes
    clicker         --  A Clicker instance for recording changes in contained
                        Parameters and RecipeOrganizers.
    name            --  A name for this organizer. Names should be unique
                        within a RecipeOrganizer and should be valid attribute
                        names.
    _confclicker    --  A Clicker for recording configuration
                        changes, esp.  additions and removal of managed
                        objects.
    _calculators    --  A managed dictionary of Calculators, indexed by name.
    _constraints    --  A dictionary of Constraints, indexed by the constrained
                        Parameter. Constraints can be added using the
                        'constrain' method.
    _parameters     --  A managed OrderedDict of contained Parameters.
    _restraints     --  A set of Restraints. Restraints can be added using the
                        'restrain' or 'confine' methods.
    _eqfactory      --  A diffpy.srfit.equation.builder.EquationFactory
                        instance that is used create Equations from string.

    Raises ValueError if the name is not a valid attribute identifier

    """

    def __init__(self, name):
        RecipeContainer.__init__(self, name)
        self._constraints = {}
        self._restraints = set()
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
        
        Raises ValueError if par is not part of the RecipeContainer.

        """
        self._removeObject(par, self._parameters)
        self._eqfactory.deRegisterBuilder(par.name)
        return

    def registerCalculator(self, f, argnames = None, makepars = True):
        """Register a Calculator so it can be used within equation strings.

        A Calculator is an elaborate function that can organize Parameters.
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
        Raises ValueError if the calculator's name is already in use by a
        managed object or an Equation object.

        """
        self._addObject(f, self._calculators, True)
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
        Raises ValueError if the function name is already in use by an Equation
        object.

        Returns the callable Equation object.

        """

        if isinstance(f, Equation):
            if name is None:
                m = "Equation must be given a name"
                raise ValueError(m)

            self._registerEquation(name, f)
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

        # Verify the name 
        if name in self._eqfactory.builders:
            m = "Equation object named '%s' already exists"%name
            raise ValueError(m)

        # Make missing Parameters, if requested
        if makepars:
            for pname in argnames:
                if pname not in self._eqfactory.builders:
                    self._newParameter(pname, 0)

        # In order to make the function callable without plugging in the
        # arguments, we will build an Equation instance and register it. We'll
        # build it in another EquationFactory so we don't clutter up our
        # namespace.
        tempfactory = EquationFactory()

        for pname in argnames:
            parbuilder = self._eqfactory.builders.get(pname)
            if parbuilder is None:
                m = "Function requires unspecified parameters (%s)."%pname
                raise AttributeError(m)

            tempfactory.registerBuilder(pname, parbuilder)

        tempfactory.registerFunction(name, f, len(argnames))

        argstr = ",".join(argnames)
        eq = equationFromString("%s(%s)"%(name,argstr), tempfactory)

        eq.name = name

        self._registerEquation(name, eq)

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

        Raises ValueError if ns uses a name that is already used for another
        managed object.
        Raises AttributeError if makepars is False and the Parameters are not
        part of this object.
        Raises ValueError if the function name is the name of another managed
        object.

        Returns the callable Equation object.

        """

        # Build the equation instance.
        eq = equationFromString(fstr, self._eqfactory, ns = ns, buildargs =
                makepars)

        eq.name = name

        # Register the equation
        self._registerEquation(name, eq)

        # FIXME - this will have to be changed when proper swapping is
        # implemented.
        # Register any new Parameters.
        for par in self._eqfactory.newargs:
            self._addParameter(par)

        return eq

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
                raise ValueError("The parameter '%s' cannot be found"%name)

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

        # Our configuration changed
        self._confclicker.click()

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
            self._confclicker.click()

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

        # Make and store the constraint
        res = Restraint()
        res.restrain(eq, lb, ub, prefactor, power, scaled)
        self._restraints.add(res)

        # Our configuration changed
        self._confclicker.click()

        return res

    def unrestrain(self, res):
        """Remove a Restraint from the RecipeOrganizer.
        
        res     --  A Restraint returned from the 'restrain' method.

        """
        if res in self._restraints:
            self._restraints.remove(res)

            # Our configuration changed
            self._confclicker.click()

        return

    def clearRestraints(self):
        """Clear all constraints."""
        for res in self._restraints:
            self.unrestrain(res)
        return

    def _getConstraints(self, recurse = True):
        """Get the Constraints for this and managed sub-objects.

        This returns a {Parameter : Constraint} dictionary.

        """

        constraints = {}

        if recurse:
            f = lambda m : hasattr(m, "_getConstraints")
            for m in ifilter(f, self._iterManaged()):
                constraints.update( m._getConstraints(recurse) )

        # Update with local constraints last. These override the others.
        constraints.update(self._constraints)

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

    def _registerEquation(self, name, newobj, check = True):
        """Register an Equation with the _eqfactory

        If check is True this raises ValueError if _eqfactory already has a
        builder of the given name.
        
        """
        if check and name in self._eqfactory.builders:
            m = "Equation object named '%s' already exists"%name
            raise ValueError(m)

        self._eqfactory.registerEquation(name, newobj)
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
