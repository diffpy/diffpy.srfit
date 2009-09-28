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

"""The base fit classes for diffpy.srfit.

This package contains modules and subpackages that are used to define a fit
problem in SrFit. Unaltered, these classes will help set up a fit problem that
can be optimized within a fitting framework. They provide the basic framework
for defining a forward calculator, and from that defining a fit problem with
data, constraints and restraints. The classes involved in this collaboration
can be tied to a fitting framework at various levels through inheritance. One
can create a fitting problem using the FitContribution, Profile and FitRecipe
classes.

Modules:
    calculator      --  Home to Calculator class (see below).
    constraint      --  Home to Constraint class. Constraints are used by
                        FitRecipe to associate Parameters and variables of a
                        fit.
    fitcontribution    --  Home to FitContribution class (see below).
    fithook         --  Home to FitHook class. The FitHook is used by the
                        FitRecipe to display information or otherwise operate on
                        the recipe during a running fit.
    fitrecipe       --  Home to FitRecipe class (see below).
    fitresults      --  Home to FitResults class (see below).
    parameter       --  Home to Parameter classes. Parameters encapsulate
                        the name, value and other properties of an equation
                        parameter.
    parameterset    --  Home to ParameterSet class. ParameterSets organize
                        Parameters and other ParameterSets within a FitRecipe.
    profile         --  Home to Profile class (see below).
    profilegenerator    --  Home to ProfileGenerator class (see below.)
    recipeorganizer --  Home to RecipeOrganizer class. This is the base class
                        for other classes that organize Parameters,
                        Constraints, Restraints and other RecipeOrganizers.
    restraint       --  Home to the Restraint class. Restraints are used by the
                        FitRecipe to coerce variables within a specified range
                        of values during a fit.
    utils           --  Various utilities used within the diffpy.srfit.fitbase
                        package.

Classes:

    Calculator      --  A class that can be used with a FitContribution (or any
                        RecipyOrganizer) to help define a calculation.
    FitContribution --  FitContributions organize Calculators and Profiles and
                        use them to calculate a residual function for a
                        FitRecipe.
    FitRecipe       --  FitRecipes are used to define a system (or recipe) to
                        be optimized (or fit).  FitRecipes combine one or more
                        FitContributions and Restraints to create a residual
                        function for the system. Parameters from the various
                        FitContributions can be turned into or constrained to
                        fitting variables.
    FitResults      --  FitResults interrogate a FitRecipe and stores all
                        relavent fitting results, such as the goodness of fit,
                        variables and constraint values. It estimates
                        uncertainties based on these values and can display
                        this information on screen or save it to a file.
    Profile         --  Class for holding observed and calculated profiles.
                        Profiles are used by ProfileGenerators to provide the
                        set of evaluation points, and by FitContributions to
                        calculate a residual.
    ProfileGenerator   --  Base class for profile generators. This class is
                        used by a FitContribution to produce a residual
                        function to be compared with data.  It can be used on
                        its own for the organziation of calculation parameters
                        and evaluation of a profile function.

Various code and design taken from Paul Kienzle's PARK package.
http://www.reflectometry.org/danse/park.html

"""

# package version
from diffpy.srfit.version import __version__

from .calculator import Calculator
from .fitcontribution import FitContribution
from .fitrecipe import FitRecipe
from .fitresults import FitResults
from .profile import Profile
from .profilegenerator import ProfileGenerator

__id__ = "$Id$"
# End of file
