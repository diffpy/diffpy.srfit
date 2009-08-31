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
"""Example of fitting a gaussian to simulated data.

This is an example of building a FitRecipe in order to fit experimental data.

The makeRecipe function shows how to build a FitRecipe, which describes what we
want to fit, and how to fit it. The scipyOptimize and parkOptimize functions
show two different ways of refining the recipe, using scipy and PARK,
respectively.

"""

import numpy

from diffpy.srfit.fitbase import FitContribution, FitRecipe, Profile, FitResults

####### Example Code

def makeRecipe():
    """Make a FitRecipe for fitting a Gaussian curve to data.

    The instructions for what we want to refine, and how to refine it will be
    defined within a FitRecipe instance. The job of a FitRecipe is to collect
    and associate all the data, the fitting equations, fitting variables,
    constraints and restraints. We will demonstrate each of these within the
    code. 
    
    Once we define the FitRecipe, we can send it an optimizer to be optimized.
    See the scipyOptimize and parkOptimize functions.
    
    """

        
    ## The Profile
    # Create a Profile to hold the experimental and calculated signal.
    profile = Profile()

    # Load data and add it to the profile. It is our responsibility to get our
    # data into the profile. Here we read the data from file.
    x, y = numpy.loadtxt("data/gaussian.dat", unpack=1)
    profile.setObservedProfile(x, y)

    ## The FitContribution
    # The FitContribution associates the Profile with a fitting equation. The
    # FitContribution also stores the parameters of the fitting equation. We
    # give our FitContribution then name "g1". We will be able to access the
    # FitContribution by that name within the FitRecipe.
    contribution = FitContribution("g1")
    # Tell the FitContribution about the Profile. The FitContribution will give
    # us access to the data held within the Profile. Here, we can tell it what
    # name we want to use for the independent variable. We tell it to use the
    # name "x".
    contribution.setProfile(profile, xname="x")

    # Now we need to create a fitting equation. We do that by writing out the
    # equation as a string. The FitContribution will turn this into a callable
    # function internally. In the process, it extracts all the parameters from
    # the equation (A, x, x0, sigma) and turns them into Parameter objects
    # internally. These objects can be accessed as attributes of the
    # contribution.  Since we told the contribution that our independent
    # variable is named "x", this value will be substituted into the fitting
    # equation whenever it is called.
    contribution.setEquation("A * exp(-0.5*(x-x0)**2/sigma**2)")

    # To demonstrate how these parameters are used, we will give "A" an initial
    # value. Note that Parameters are not numbers, but are containers for
    # numbers. To change the value of a parameter, use its 'setValue' method.
    # To get its value, use the 'getValue' method. Parameters also have a
    # 'name' attribute.
    contribution.A.setValue(1.0)

    ## The FitRecipe
    # The FitRecipe lets us define what we want to fit. It is where we can
    # create variables, constraints and restraints.
    recipe = FitRecipe()

    # Here we tell the FitRecipe to use our FitContribution. When the FitRecipe
    # calculates its residual function, it will call on the FitContribution to
    # do part of the work.
    recipe.addContribution(contribution)

    # Specify which Parameters we want to vary in the fit.  This will add
    # Variables to the FitRecipe that directly modify the Parameters of the
    # FitContribution.
    # 
    # Here we create a Variable for the 'A' Parameter from our fit equation.
    # The resulting Variable will be named 'A' as well, but it will be accessed
    # via the FitRecipe.
    recipe.addVar(contribution.A)
    # Here we create the Variable for 'x0' and give it an initial value of 5.
    recipe.addVar(contribution.x0, 5)
    # Here we create a Variable named 'sig', which is tied to the 'sigma'
    # Parameter of our FitContribution. We give it an initial value through the
    # FitRecipe instance.
    recipe.addVar(contribution.sigma, name = "sig")
    recipe.sig.setValue(1)

    return recipe

def scipyOptimize(recipe):
    """Optimize the recipe created above using scipy.

    The FitRecipe we created in makeRecipe has a 'residual' method that we can be
    minimized using a scipy optimizer. The details are described in the source.

    """

    # We're going to use the least-squares (Levenberg-Marquardt) optimizer from
    # scipy. We simply have to give it the function to minimize
    # (recipe.residual) and the starting values of the Variables
    # (recipe.getValues()).
    from scipy.optimize.minpack import leastsq
    print "Fit using scipy's LM optimizer"
    leastsq(recipe.residual, recipe.getValues())

    return

def parkOptimize(recipe):
    """Optimize the recipe created above using PARK.
    
    This requires the 'pak' branch of PARK to be installed on your system.

    """
    from diffpy.srfit.park import FitnessAdapter

    # We have to turn the recipe into something that PARK can use. In PARK, a
    # Fitness object is the equivalent of a SrFit FitContribution. However, we
    # want to use the varibles, constraints and restraints, defined in our
    # FitRecipe, so we will turn it into a Fitness object. To do this, we have
    # written a special FitnessAdapter class in the diffpy.srfit.park package.
    f = FitnessAdapter(recipe)

    # Now we can fit this using park.
    from park.fitting.fit import fit
    print "Fit using the default PARK optimizer"
    result = fit([f])

    return


def plotResults(recipe):
    """Plot the results contained within a refined FitRecipe."""

    # We can access the data and fit profile through the Profile we created
    # above. We get to it through our FitContribution, which we named "g1".
    #
    # The independent variable. This is always under the "x" attribute.
    x = recipe.g1.profile.x
    # The observed profile that we loaded earlier, the "y" attribute.
    y = recipe.g1.profile.y
    # The calculated profile, the "ycalc" attribute.
    ycalc = recipe.g1.profile.ycalc

    # This stuff is specific to pylab from the matplotlib distribution.
    import pylab
    pylab.plot(x, y, 'b.', label = "observed Gaussian")
    pylab.plot(x, ycalc, 'g-', label = "calculated Gaussian")
    pylab.legend(loc = (0.0,0.8))
    pylab.xlabel("x")
    pylab.ylabel("y")

    pylab.show()
    return


if __name__ == "__main__":

    # Create the recipe
    recipe = makeRecipe()

    # Refine using the optimizer of your choice
    scipyOptimize(recipe)
    #parkOptimize(recipe)

    # Get the results in a FitResults object. The FitResults object stores the
    # current state of the recipe, and uses it to calculate useful statistics
    # about the fit.  If you later modify the recipe, the FitResults object
    # will hold the recipe values from when it was created. You can tell it to
    # update its values by calling its 'update' method.
    res = FitResults(recipe)

    # Print the results
    res.printResults()

    # Plot the results
    plotResults(recipe)


# End of file
