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
"""Base classes and tools for constructing a FitRecipe.

RecipeContainer is the base class for organizing Parameters, and other
RecipeContainers.  RecipeOrganizer is an extended RecipeContainer that
incorporates equation building, constraints and Restraints.
equationFromString creates an Equation instance from a string.
"""

__all__ = ["RecipeContainer", "RecipeOrganizer", "equationFromString"]

import re
from collections import OrderedDict
from functools import partial
from itertools import chain, groupby

import six
from numpy import inf

from diffpy.srfit.equation import Equation
from diffpy.srfit.equation.builder import EquationFactory
from diffpy.srfit.fitbase.configurable import Configurable
from diffpy.srfit.fitbase.constraint import Constraint
from diffpy.srfit.fitbase.parameter import Parameter
from diffpy.srfit.fitbase.restraint import Restraint
from diffpy.srfit.fitbase.validatable import Validatable
from diffpy.srfit.interface import _recipeorganizer_interface
from diffpy.srfit.util import _DASHEDLINE
from diffpy.srfit.util import sortKeyForNumericString as numstr
from diffpy.srfit.util.nameutils import validateName
from diffpy.srfit.util.observable import Observable


class RecipeContainer(Observable, Configurable, Validatable):
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

    Properties
    names           --  Variable names (read only). See getNames.
    values          --  Variable values (read only). See getValues.
    """

    names = property(lambda self: self.getNames())
    values = property(lambda self: self.getValues())

    def __init__(self, name):
        Observable.__init__(self)
        Configurable.__init__(self)
        validateName(name)
        self.name = name
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

    def iterPars(self, pattern="", recurse=True):
        """Iterate over the Parameters contained in this object.

        Parameters
        ----------
        pattern : str
            Iterate over parameters with names matching this regular
            expression (all parameters by default).
        recurse : bool
            Recurse into managed objects when True (default).
        """
        regexp = re.compile(pattern)
        for par in list(self._parameters.values()):
            if regexp.search(par.name):
                yield par
        if not recurse:
            return
        # Iterate over objects within the managed dictionaries.
        managed = self.__managed[:]
        managed.remove(self._parameters)
        for m in managed:
            for obj in m.values():
                if hasattr(obj, "iterPars"):
                    for par in obj.iterPars(pattern=pattern):
                        yield par
        return

    def __iter__(self):
        """Iterate over top-level parameters."""
        return iter(self._parameters.values())

    def __len__(self):
        """Get number of top-level parameters."""
        return len(self._parameters)

    def __getitem__(self, idx):
        """Get top-level parameters by index."""
        # need to wrap this in a list for python 3 compatibility.
        return list(self._parameters.values())[idx]

    def __getattr__(self, name):
        """Gives access to the contained objects as attributes."""
        arg = self.get(name)
        if arg is None:
            raise AttributeError(name)
        return arg

    # Ensure there is no __dir__ override in the base class.
    assert (
        getattr(Observable, "__dir__", None)
        is getattr(Configurable, "__dir__", None)
        is getattr(Validatable, "__dir__", None)
        is getattr(object, "__dir__", None)
    )

    def __dir__(self):
        "Return sorted list of attributes for this object."
        rv = set(dir(type(self)))
        rv.update(self.__dict__)
        # self.get fetches looks up for items in all managed dictionaries.
        # Add keys from each dictionary in self.__managed.
        rv.update(*self.__managed)
        rv = sorted(rv)
        return rv

    # Needed by __setattr__
    _parameters = OrderedDict()
    __managed = []

    def __setattr__(self, name, value):
        """Parameter access and object checking."""
        if name in self._parameters:
            par = self._parameters[name]
            if isinstance(value, Parameter):
                par.value = value.value
            else:
                par.value = value
            return

        m = self.get(name)
        if m is not None:
            raise AttributeError("Cannot set '%s'" % name)

        super(RecipeContainer, self).__setattr__(name, value)
        return

    def __delattr__(self, name):
        """Delete parameters with del.

        This does not allow deletion of non-parameters, as this may
        require configuration changes that are not yet handled in a
        general way.
        """
        if name in self._parameters:
            self._removeParameter(self._parameters[name])
            return

        m = self.get(name)
        if m is not None:
            raise AttributeError("Cannot delete '%s'" % name)

        super(RecipeContainer, self).__delattr__(name)
        return

    def get(self, name, default=None):
        """Get a managed object."""
        for d in self.__managed:
            arg = d.get(name)
            if arg is not None:
                return arg

        return default

    def getNames(self):
        """Get the names of managed parameters."""
        return [p.name for p in self._parameters.values()]

    def getValues(self):
        """Get the values of managed parameters."""
        return [p.value for p in self._parameters.values()]

    def _addObject(self, obj, d, check=True):
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
            message = "%s has no name" % obj.__class__.__name__
            raise ValueError(message)

        # Check for extant object in d with same name
        oldobj = d.get(obj.name)
        if check and oldobj is not None:
            message = "%s with name '%s' already exists" % (
                obj.__class__.__name__,
                obj.name,
            )
            raise ValueError(message)

        # Check for object with same name in other dictionary.
        if oldobj is None and self.get(obj.name) is not None:
            message = "Non-%s with name '%s' already exists" % (
                obj.__class__.__name__,
                obj.name,
            )
            raise ValueError(message)

        # Detach the old object, if there is one
        if oldobj is not None:
            oldobj.removeObserver(self._flush)

        # Add the object
        d[obj.name] = obj

        # Observe the object
        obj.addObserver(self._flush)

        # Store this as a configurable object
        self._storeConfigurable(obj)
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

        Returns a list of objects. The first member of the list is this object,
        and each subsequent member is a sub-object of the previous one.  The
        last entry in the list is obj. If obj cannot be found, the list is
        empty.
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

        This will force any observer to invalidate its state. By default
        this does nothing.
        """
        self.notify(other)
        return

    def _validate(self):
        """Validate my state.

        This validates that contained Parameters and managed objects are
        valid.

        Raises AttributeError if validation fails.
        """
        iterable = chain(self.__iter__(), self._iterManaged())
        self._validateOthers(iterable)
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
                        'restrain' method.
    _eqfactory      --  A diffpy.srfit.equation.builder.EquationFactory
                        instance that is used create Equations from string.

    Properties
    names           --  Variable names (read only). See getNames.
    values          --  Variable values (read only). See getValues.

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

        This creates a new Parameter and adds it to the container using
        the _addParameter method.

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

        This de-registers the Parameter with the _eqfactory. The
        Parameter will remain part of built equations.

        Note that constraints and restraints involving the Parameter are
        not modified.

        Raises ValueError if par is not part of the RecipeOrganizer.
        """
        self._removeObject(par, self._parameters)
        self._eqfactory.deRegisterBuilder(par.name)
        return

    def registerCalculator(self, f, argnames=None):
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
            fncode = f.__call__.__func__.__code__
            argnames = list(fncode.co_varnames)
            argnames = argnames[1 : fncode.co_argcount]

        for pname in argnames:
            if pname not in self._eqfactory.builders:
                par = self._newParameter(pname, 0)
            else:
                par = self.get(pname)
            f.addLiteral(par)

        # Now return an equation object
        eq = self._eqfactory.makeEquation(f.name)
        return eq

    def registerFunction(self, f, name=None, argnames=None):
        """Register a function so it can be used within equation strings.

        This creates a function with this class that can be used within string
        equations.  The resulting equation does not require the arguments to be
        passed in the equation string, as this will be handled automatically.

        f           --  The callable to register. If this is an Equation
                        instance, then all that needs to be provided is a name.
        name        --  The name of the function to be used in equations. If
                        this is None (default), the method will try to
                        determine the name of the function automatically.
        argnames    --  The names of the arguments to f (list or None).
                        If this is None (default), then the argument names will
                        be extracted from the function.

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

        # Introspection code
        if name is None or argnames is None:

            import inspect

            fncode = None

            # This will let us offset the argument list to eliminate 'self'
            offset = 0

            # check regular functions
            if inspect.isfunction(f):
                fncode = f.__code__
            # check class method
            elif inspect.ismethod(f):
                fncode = f.__func__.__code__
                offset = 1
            # check functor
            elif hasattr(f, "__call__") and hasattr(f.__call__, "__func__"):
                fncode = f.__call__.__func__.__code__
                offset = 1
            else:
                m = "Cannot extract name or argnames"
                raise ValueError(m)

            # Extract the name
            if name is None:
                name = fncode.co_name
                if name == "<lambda>":
                    m = "You must supply a name name for a lambda function"
                    raise ValueError(m)

            # Extract the arguments
            if argnames is None:
                argnames = list(fncode.co_varnames)
                argnames = argnames[offset : fncode.co_argcount]

        # End introspection code

        # Make missing Parameters
        for pname in argnames:
            if pname not in self._eqfactory.builders:
                self._newParameter(pname, 0)

        # Initialize and register
        from diffpy.srfit.fitbase.calculator import Calculator

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

    def registerStringFunction(self, fstr, name, ns={}):
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
        eq = equationFromString(fstr, self._eqfactory, ns=ns, buildargs=True)
        eq.name = name

        # Register any new Parameters.
        for par in self._eqfactory.newargs:
            self._addParameter(par)

        # Register the equation as a callable function.
        argnames = eq.argdict.keys()
        return self.registerFunction(eq, name, argnames)

    def evaluateEquation(self, eqstr, ns={}):
        """Evaluate a string equation.

        eqstr   --  A string equation to evaluate. The equation is evaluated at
                    the current value of the registered Parameters.
        ns      --  A dictionary of Parameters, indexed by name, that are
                    used in fstr, but not part of the FitRecipe (default {}).

        Raises ValueError if ns uses a name that is already used for a
        variable.
        """
        eq = equationFromString(eqstr, self._eqfactory, ns)
        try:
            rv = eq()
        finally:
            self._eqfactory.wipeout(eq)
        return rv

    def constrain(self, par, con, ns={}):
        """Constrain a parameter to an equation.

        Note that only one constraint can exist on a Parameter at a time.

        par     --  The name of a Parameter or a Parameter to constrain.
        con     --  A string representation of the constraint equation or a
                    Parameter to constrain to.  A constraint equation must
                    consist of numpy operators and "known" Parameters.
                    Parameters are known if they are in the ns argument, or if
                    they are managed by this object.
        ns      --  A dictionary of Parameters, indexed by name, that are used
                    in the parameter, but not part of this object (default {}).

        Raises ValueError if ns uses a name that is already used for a
        variable.
        Raises ValueError if par is a string but not part of this object or in
        ns.
        Raises ValueError if par is marked as constant.
        """
        if isinstance(par, six.string_types):
            name = par
            par = self.get(name)
            if par is None:
                par = ns.get(name)

        if par is None:
            raise ValueError("The parameter cannot be found")

        if par.const:
            raise ValueError("The parameter '%s' is constant" % par)

        if isinstance(con, six.string_types):
            eqstr = con
            eq = equationFromString(con, self._eqfactory, ns)
        else:
            eq = Equation(root=con)
            eqstr = con.name

        eq.name = "_constraint_%s" % par.name

        # Make and store the constraint
        con = Constraint()
        con.constrain(par, eq)
        # Store the equation string so it can be shown later.
        con.eqstr = eqstr
        self._constraints[par] = con

        # Our configuration changed
        self._updateConfiguration()

        return

    def isConstrained(self, par):
        """Determine if a Parameter is constrained in this object.

        par     --  The name of a Parameter or a Parameter to check.
        """
        if isinstance(par, six.string_types):
            name = par
            par = self.get(name)

        return par in self._constraints

    def unconstrain(self, *pars):
        """Unconstrain a Parameter.

        This removes any constraints on a Parameter.

        *pars   --  The names of Parameters or Parameters to unconstrain.


        Raises ValueError if the Parameter is not constrained.
        """
        update = False
        for par in pars:
            if isinstance(par, six.string_types):
                name = par
                par = self.get(name)

            if par is None:
                raise ValueError("The parameter cannot be found")

            if par in self._constraints:
                self._constraints[par].unconstrain()
                del self._constraints[par]
                update = True

        if update:
            # Our configuration changed
            self._updateConfiguration()

        else:

            raise ValueError("The parameter is not constrained")

        return

    def getConstrainedPars(self, recurse=False):
        """Get a list of constrained managed Parameters in this object.

        recurse --  Recurse into managed objects and retrieve their constrained
                    Parameters as well (default False).
        """
        const = self._getConstraints(recurse)
        return const.keys()

    def clearConstraints(self, recurse=False):
        """Clear all constraints managed by this organizer.

        recurse --  Recurse into managed objects and clear all constraints
                    found there as well.

        This removes constraints that are held in this organizer, no matter
        where the constrained parameters are from.
        """
        if self._constraints:
            self.unconstrain(*self._constraints)

        if recurse:
            for m in filter(_has_clear_constraints, self._iterManaged()):
                m.clearConstraints(recurse)
        return

    def restrain(self, res, lb=-inf, ub=inf, sig=1, scaled=False, ns={}):
        """Restrain an expression to specified bounds.

        res     --  An equation string or Parameter to restrain.
        lb      --  The lower bound on the restraint evaluation (default -inf).
        ub      --  The lower bound on the restraint evaluation (default inf).
        sig     --  The uncertainty on the bounds (default 1).
        scaled  --  A flag indicating if the restraint is scaled (multiplied)
                    by the unrestrained point-average chi^2 (chi^2/numpoints)
                    (default False).
        ns      --  A dictionary of Parameters, indexed by name, that are used
                    in the equation string, but not part of the RecipeOrganizer
                    (default {}).

        The penalty is calculated as
        (max(0, lb - val, val - ub)/sig)**2
        and val is the value of the calculated equation.  This is multiplied by
        the average chi^2 if scaled is True.

        Raises ValueError if ns uses a name that is already used for a
        Parameter.
        Raises ValueError if res depends on a Parameter that is not part of
        the RecipeOrganizer and that is not defined in ns.

        Returns the Restraint object for use with the 'unrestrain' method.
        """

        if isinstance(res, six.string_types):
            eqstr = res
            eq = equationFromString(res, self._eqfactory, ns)
        else:
            eq = Equation(root=res)
            eqstr = res.name

        # Make and store the restraint
        res = Restraint(eq, lb, ub, sig, scaled)
        res.eqstr = eqstr
        self.addRestraint(res)
        return res

    def addRestraint(self, res):
        """Add a Restraint instance to the RecipeOrganizer.

        res     --  A Restraint instance.
        """
        self._restraints.add(res)
        # Our configuration changed. Notify observers.
        self._updateConfiguration()
        return

    def unrestrain(self, *ress):
        """Remove a Restraint from the RecipeOrganizer.

        *ress   --  Restraints returned from the 'restrain' method or added
                    with the 'addRestraint' method.
        """
        update = False
        restuple = tuple(self._restraints)
        for res in ress:
            if res in restuple:
                self._restraints.remove(res)
                update = True

        if update:
            # Our configuration changed
            self._updateConfiguration()

        return

    def clearRestraints(self, recurse=False):
        """Clear all restraints.

        recurse --  Recurse into managed objects and clear all restraints
                    found there as well.
        """
        self.unrestrain(*self._restraints)
        if recurse:
            for msg in filter(_has_clear_restraints, self._iterManaged()):
                msg.clearRestraints(recurse)
        return

    def _getConstraints(self, recurse=True):
        """Get the constrained Parameters for this and managed sub-objects."""
        constraints = {}
        if recurse:
            for m in filter(_has_get_constraints, self._iterManaged()):
                constraints.update(m._getConstraints(recurse))

        constraints.update(self._constraints)

        return constraints

    def _getRestraints(self, recurse=True):
        """Get the Restraints for this and embedded ParameterSets.

        This returns a set of Restraint objects.
        """
        restraints = set(self._restraints)
        if recurse:
            for m in filter(_has_get_restraints, self._iterManaged()):
                restraints.update(m._getRestraints(recurse))

        return restraints

    def _validate(self):
        """Validate my state.

        This performs RecipeContainer validations. This validates
        contained Restraints and Constraints.

        Raises AttributeError if validation fails.
        """
        RecipeContainer._validate(self)
        iterable = chain(self._restraints, self._constraints.values())
        self._validateOthers(iterable)
        return

    # For printing the configured recipe to screen

    def _formatManaged(self, prefix=""):
        """Format hierarchy of managed parameters for showing.

        Parameters
        ----------
        prefix : str
            The leading string to be prefixed to each parameter name.

        Returns
        -------
        list
            List of formatted lines, one per each Parameter.
        """
        lines = []
        formatstr = "{:<W}{}"
        # Format own parameters.
        if self._parameters:
            w0 = max(len(n) for n in self._parameters)
            w1 = ((w0 + len(prefix) + 1) // 4 + 1) * 4
            fmt = formatstr.replace("W", str(w1))
            lines.extend(
                fmt.format(prefix + n, p.value)
                for n, p in self._parameters.items()
            )
        # Recurse into managed objects.
        for obj in self._iterManaged():
            if hasattr(obj, "_formatManaged"):
                oprefix = prefix + obj.name + "."
                tlines = obj._formatManaged(prefix=oprefix)
                lines.extend([""] if lines and tlines else [])
                lines.extend(tlines)
        return lines

    def _formatConstraints(self):
        """Format constraints for showing.

        This collects constraints on all levels of the hierarchy and displays
        them with respect to this level.

        Returns
        -------
        list
            List of formatted lines displaying the defined constraints.
            Return empty list when no constraints were defined.
        """
        cdict = self._getConstraints()
        # Find each constraint and format the equation
        clines = []
        for par, con in cdict.items():
            loc = self._locateManagedObject(par)
            if loc:
                locstr = ".".join(o.name for o in loc[1:])
                clines.append("%s <-- %s" % (locstr, con.eqstr))
            else:
                clines.append("%s <-- %s" % (par.name, con.eqstr))
        clines.sort(key=numstr)
        return clines

    def _formatRestraints(self):
        """Format restraints for showing.

        This collects restraints on all levels of the hierarchy and displays
        them with respect to this level.

        Returns
        -------
        list
            List of formatted lines displaying the defined restraints.
            Return empty list when no restraints were defined.
        """
        rset = self._getRestraints()
        rlines = []
        for res in rset:
            line = "%s: lb = %f, ub = %f, sig = %f, scaled = %s" % (
                res.eqstr,
                res.lb,
                res.ub,
                res.sig,
                res.scaled,
            )
            rlines.append(line)
        rlines.sort(key=numstr)
        return rlines

    def show(self, pattern="", textwidth=78):
        """Show the configuration hierarchy on the screen.

        This will print out a summary of all contained objects.

        Parameters
        ----------
        pattern : str, optional
            Limit output to only those parameters that match this regular
            expression (match all by default).
        textwidth : int, optional
            Trim formatted lines at this text width to avoid folding at
            the screen width.  Do not trim when negative or 0.
        """
        regexp = re.compile(pattern)
        _pmatch_with_re = partial(_pmatch, regexp=regexp)
        # Show sub objects and their parameters
        lines = []
        tlines = self._formatManaged()
        if tlines:
            lines.extend(["Parameters", _DASHEDLINE])
            linesok = filter(_pmatch_with_re, tlines)
            lastnotblank = False
            # squeeze repeated blank lines
            for lastnotblank, g in groupby(linesok, bool):
                lines.extend(g if lastnotblank else [""])
            # remove trailing blank line
            if not lastnotblank:
                lines.pop(-1)

        # FIXME - parameter names in equations not particularly informative
        # Show constraints
        cmatch = regexp.search
        tlines = self._formatConstraints()
        if tlines:
            if lines:
                lines.append("")
            lines.extend(["Constraints", _DASHEDLINE])
            lines.extend(filter(cmatch, tlines))

        # FIXME - parameter names in equations not particularly informative
        # Show restraints
        tlines = self._formatRestraints()
        if tlines:
            if lines:
                lines.append("")
            lines.extend(["Restraints", _DASHEDLINE])
            lines.extend(filter(_pmatch_with_re, tlines))

        # Determine effective text width tw.
        tw = textwidth if (textwidth is not None and textwidth > 0) else None
        # Avoid outputting "\n" when there is no output.
        if lines:
            print("\n".join(s[:tw] for s in lines))
        return


# End RecipeOrganizer


def equationFromString(
    eqstr, factory, ns={}, buildargs=False, argclass=Parameter, argkw={}
):
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


def _has_clear_constraints(msg):
    return hasattr(msg, "clearConstraints")


def _has_clear_restraints(msg):
    return hasattr(msg, "clearRestraints")


def _has_get_restraints(msg):
    return hasattr(msg, "_getRestraints")


def _has_get_constraints(msg):
    return hasattr(msg, "_getConstraints")


def _pmatch(inp_str, regexp):
    parts = inp_str.split(None, 1)
    return len(parts) < 2 or regexp.search(parts[0])
