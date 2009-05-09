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
"""Contribution class. 

Contributions organize an Equation and Calculator that calculate the signal,
and a Profile that holds the signal.
"""

from numpy import concatenate, sqrt, inf, dot

from diffpy.srfit.equation import Equation
from diffpy.srfit.equation.literals import Generator
from diffpy.srfit.equation.builder import EquationFactory

from .parameter import Parameter
from .modelorganizer import ModelOrganizer, equationFromString


class Contribution(ModelOrganizer):
    """Contribution class.

    Contributions organize an Equation that calculates the signal, and a
    Profile that holds the signal. Contraints and Restraints can be created as
    part of a Contribution.

    Attributes
    clicker         --  A Clicker instance for recording changes in the
                        Parameters or the residual components.
    name            --  A name for this Contribution.
    calculator      --  A Calculator instance for generating a signal
                        (optional). Contributions can share a Calculator
                        instance. If a calculator is not defined, the equation
                        to refine must be set with the setEquation method.
    profile         --  A Profile that holds the measured (and calcuated)
                        signal.
    _constraints    --  A dictionary of Constraints, indexed by the constrained
                        Parameter. Constraints can be added using the
                        'constrain' method.

    _eq             --  The Contribution equation that will be optimized.
    _eqfactory      --  A diffpy.srfit.equation.builder.EquationFactory
                        instance that is used to create constraints and
                        restraints from strings.
    _organizers     --  A reference to the Calcualtor's _organizers attribute.
    _orgdict        --  A reference to the Calculator's _orgdict attribute.
    _parameters     --  A reference to the Calculator's _parameters attribute.
    _restraints     --  A set of Restraints. Restraints can be added using the
                        'restrain' or 'confine' methods.
    _xname          --  The name of of the independent variable from the
                        profile (default None). 
    _yname          --  The name of of the observed profile (default None). 
    _dyname         --  The name of of the uncertainty in the observed profile
                        (default None).
    """

    def __init__(self, name):
        """Initialization."""
        ModelOrganizer.__init__(self, name)
        self._eq = None
        self._reseq = None
        self.profile = None
        self.calculator = None
        self._calcname = None
        self._xname = None
        self._yname = None
        self._dyname = None
        return
    
    # Make some methods public that were protected
    addParameter = ModelOrganizer._addParameter
    newParameter = ModelOrganizer._newParameter
    removeParameter = ModelOrganizer._removeParameter

    def setProfile(self, profile, xname = None, yname = None, dyname = None):
        """Assign the profile for this contribution.

        This resets the current residual.
        
        profile --  A Profile that specifies the calculation points and which
                    will store the calculated signal.
        xname   --  The name of the independent variable from the Profile. If
                    this is None (default), then the name specified by the
                    Profile for this parametere will be used.  This variable is
                    usable within the Equation with the specified name.
        yname   --  The name of the observed profile.  If this is None
                    (default), then the name specified by the Profile for this
                    parametere will be used.  This variable is usable within
                    the Equation with the specified name.
        dyname  --  The name of the uncertainty in the observed profile. If
                    this is None (default), then the name specified by the
                    Profile for this parametere will be used.  This variable is
                    usable within the Equation with the specified name.
        

        """
        self.profile = profile

        # Clear the previous profile information
        self._eqfactory.deRegisterBuilder(self._xname)
        self._eqfactory.deRegisterBuilder(self._yname)
        self._eqfactory.deRegisterBuilder(self._dyname)

        if xname is None:
            xname = self.profile.xpar.name
        if yname is None:
            yname = self.profile.ypar.name
        if dyname is None:
            dyname = self.profile.dypar.name

        self._xname = xname
        self._eqfactory.registerArgument(xname, self.profile.xpar)
        self._yname = yname
        self._eqfactory.registerArgument(yname, self.profile.ypar)
        self._dyname = dyname
        self._eqfactory.registerArgument(dyname, self.profile.dypar)

        if self.calculator is not None:
            self.setResidualEquation()
        return

    def setCalculator(self, calc, name = None):
        """Set the Calculator to be used by this Contribution.

        The Calculator is given a name so that it can be used as part of the
        equation that is used to generate the signal. This can be different
        from the name of the Calculator for attribute purposes. 
        
        Calling setCalculator sets the equation to call the calculator and
        resets the residual.

        calc    --  A Calculator instance
        name    --  A name for the calculator. If name is None (default), then
                    the Calculator's name attribute will be used.

        """
        self.calculator = calc

        if name is None:
            name = calc.name

        # Let the ModelOrganizer structure know of the calculator
        self._addOrganizer(calc)

        # Register the calculator with the equation factory
        self._eqfactory.registerGenerator(name, calc)

        self.setEquation(name)
        return

    def setEquation(self, eqstr, makepars = True):
        """Set the refinement equation for the Contribution.

        eqstr   --  A string representation of the equation. Any Parameter
                    registered by addParameter or setProfile, or function
                    registered by setCalculator or registerFunction  can be can
                    be used in the equation by name.
        makepars    --  A flag indicating whether missing Parameters can be
                    created by the Factory (default True). If False, then the a
                    ValueError will be raised if there are undefined arguments
                    in the eqstr. 

        The equation will be usable within setResidualEquation by calling
        "eq()" or "eq". 
        
        Calling setEquation resets the residual equation.

        Raises ValueError if makepars is false and eqstr depends on a Parameter
        that is not part of the Contribution.

        """
        self._eq = self.registerStringFunction(eqstr, "eq", makepars)

        if self.profile is not None:
            self.setResidualEquation()
        return

    def setResidualEquation(self, eqstr = None):
        """Set the residual equation for the Contribution.

        eqstr   --  A string representation of the residual. If eqstr is None
                    (default), then the chi2 residual will be used (see the
                    residual method.)

        Two residuals are preset for convenience, "chiv" and "resv".
        chiv is defined such that dot(chiv, chiv) = chi^2.
        resv is defined such that dot(resv, resv) = Rw.
        You can call on these in your residual equation. Note that the quantity
        that will be optimized is the summed square of the residual equation.
        Keep that in mind when defining a new residual or using the built-in
        ones.

        Raises AttributeError if the Profile is not yet defined.
        Raises ValueError if ns uses a name that is already used for a
        Parameter.
        Raises ValueError if eqstr depends on a Parameter that is not part of
        the Contribution and that is not defined in ns.
        """
        if self.profile is None:
            raise AttributeError("Define the profile first")

        # Register some convenient residuals
        chivstr = "(eq() - %s)/%s" % (self._yname, self._dyname)
        chiv = equationFromString(chivstr, self._eqfactory)
        self._eqfactory.registerEquation("chiv", chiv)

        resvstr = "(eq() - %s)/sum(%s**2)**0.5" % (self._yname, self._yname)
        resv = equationFromString(resvstr, self._eqfactory)
        self._eqfactory.registerEquation("resv", resv)

        # Now set the residual to one of these or create a new one
        if eqstr is None:
            self._reseq = chiv
        else:
            self._reseq = equationFromString(eqstr, self._eqfactory)

        return

    def residual(self):
        """Calculate the residual for this contribution.

        When this method is called, it is assumed that all parameters have been
        assigned their most current values by the FitModel. This will be the
        case when being called as part of a FitModel refinement.

        The residual is by default an array chiv:
        chiv = (eq() - self.profile.y) / self.profile.dy
        The value that is optimized is dot(residual, residual).

        The residual equation can be changed with the setResidualEquation
        method.
        
        """
        # If we have a calculator, make sure it is working on my profile
        if self.calculator is not None:
            self.calculator.setProfile(self.profile)
        # Assign the calculated profile
        self.profile.ycalc = self._eq()
        # Note that equations only recompute when their inputs are modified, so
        # the following will not recompute the equation.
        return self._reseq()

    def registerFunction(self, f, name = None, argnames = None, 
            makepars = True):
        """Register a function so it can be used in the Contribution equation.

        This creates a function useable within setEquation and
        setResidualEquation. The function does not require the arguments to be
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
                        argnames if they don't already appear in the
                        Contribution (default True). Parameters created in this
                        way are added to the contribution using the
                        newParameter method.  If makepars is False , then it is
                        necessary that the parameters are already part of this
                        object in order to make the function.

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
            self._eqfactory.registerEquation(name, f)
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
            # ehck functor
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
                    self.newParameter(pname, 0)

        # In order to make the function callable without plugging in the
        # arguments, we will build an Equation instance and register it. We'll
        # build it in another EquationFactory so we don't clutter up our
        # namespace.
        factory = EquationFactory()

        for pname in argnames:
            par = self._eqfactory.builders.get(pname)
            if par is None:
                m = "Function requires unspecified parameters (%s)."%pname
                raise AttributeError(m)

            factory.registerBuilder(pname, par)

        factory.registerFunction(name, f, len(argnames))

        argstr = ",".join(argnames)
        eq = equationFromString("%s(%s)"%(name,argstr), factory)

        self._eqfactory.registerEquation(name, eq)

        return eq

    def registerStringFunction(self, fstr, name, makepars = True):
        """Register a string function.

        This creates a function useable within setEquation and
        setResidualEquation. The function does not require the arguments to be
        passed in the equation string, as this will be handled automatically.

        fstr        --  A string equation to register.
        name        --  The name of the function to be used in equations. 
        makepars    --  Flag indicating whether to make parameters from the
                        arguments of fstr if they don't already appear in the
                        Contribution (default True). Parameters created in this
                        way are added to the contribution using the
                        newParameter method.  If makepars is False , then it is
                        necessary that the parameters are already part of this
                        object in order to make the function.

        Raises AttributeError if makepars is False and the parameters are not
        part of this object.

        Returns the callable Equation object.

        """

        # Build the equation instance.
        eq = equationFromString(fstr, self._eqfactory, buildargs =
                makepars)

        # Register any new parameters
        for par in self._eqfactory.newargs:
            self._addParameter(par)

        # Register the equation by name
        self._eqfactory.registerEquation(name, eq)

        return eq

    def evaluateEquation(self, eqstr):
        """Evaluate a string equation.

        eqstr   --  A string equation to evaluate. The equation is evaluated at
                    the current value of the registered parameters. The
                    equation can be specified as described in the setEquation
                    method.

        """
        eq = equationFromString(eqstr, self._eqfactory)
        return eq()



# version
__id__ = "$Id$"

#
# End of file
