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

This is an example of building a FitModel in order to fit experimental data.

The makeModel function shows how to build a FitModel, which describes what we
want to fit, and how to fit it. The scipyOptimize and parkOptimize functions
show two different ways of refining the model, using scipy and PARK,
respectively.

"""

import numpy

from diffpy.srfit.fitbase import Contribution, FitModel, Profile, FitResults

####### Example Code

def makeModel():
    """Make a FitModel for fitting a Gaussian curve to data.

    Our model will be defined within a FitModel instance. The job of a FitModel
    is to collect and associate the fitting equations, data, fitting
    variables, constraints and restraints, which collectively define a fit.
    
    Once we define the FitModel, we can send it an optimizer to be optimized.
    See the scipyOptimize and parkOptimize functions.
    
    """

        
    ## The Profile
    # Create a Profile to hold the experimental and calculated signal.
    profile = Profile()

    # Load data and add it to the profile. It is our responsibility to get our
    # data into the profile. Here we read the data from file.
    x, y = numpy.loadtxt("data/gaussian.dat", unpack=1)
    profile.setObservedProfile(x, y)

    ## The Contribution
    # The Contribution associates the Profile with a fitting equation. The
    # Contribution also stores the parameters of the fitting equation. We give
    # our Contribution then name "g1". We will be able to access the
    # Contribution by that name within the FitModel.
    contribution = Contribution("g1")
    # Tell the Contribution about the Profile. The Contribution will give us
    # access to the data held within the Profile. Here, we can tell it what
    # name we want to use for the independent variable. We tell it to use the
    # name "x".
    contribution.setProfile(profile, xname="x")

    # Now we need to create a fitting equation. We do that by writing out the
    # equation as a string. The contribution will turn this into a callable
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

    ## The FitModel
    # The FitModel lets us define what we want to fit. It is where we can
    # create variables, constraints and restraints.
    model = FitModel()

    # Here we tell the FitModel to use our Contribution. When the FitModel
    # calculates its residual function, it will call on the Contribution to do
    # part of the work.
    model.addContribution(contribution)

    # Specify which Parameters we want to vary in the fit.  This will add
    # Parameters to the FitModel that directly modify the Parameters of the
    # Contribution. We refer to the FitModel's Parameters as variables.
    # 
    # Here we create a variable for the 'A' Parameter from our fit equation.
    # The resulting variable will be named 'A' as well, but it will be accessed
    # via the FitModel.
    model.addVar(contribution.A)
    # Here we create the variable for 'x0' and give it an initial value of 5.
    model.addVar(contribution.x0, 5)
    # Here we create a variable named 'sig', which is tied to the 'sigma'
    # Parameter of our Contribution. We give it an initial value through the
    # model.
    model.addVar(contribution.sigma, name = "sig")
    model.sig.setValue(1)

    return model

def scipyOptimize(model):
    """Optimize the model created above using scipy.

    The FitModel we created in makeModel has a 'residual' method that we can be
    minimized using a scipy optimizer. The details are described in the source.

    """

    # We're going to use the least-squares (Levenberg-Marquardt) optimizer from
    # scipy. We simply have to give it the function to minimize
    # (model.residual) and the starting values of the variables
    # (model.getValues()).
    from scipy.optimize.minpack import leastsq
    print "Fit using scipy's LM optimizer"
    leastsq(model.residual, model.getValues())

    return

def parkOptimize(model):
    """Optimize the model created above using PARK.
    
    This requires the 'pak' branch of PARK to be installed on your system.
    """
    from diffpy.srfit.park import FitnessAdapter

    # We have to turn the model into something that PARK can use. In PARK, a
    # Fitness object is the equivalent of a SrFit Contribution. However, we
    # want to use the varibles, constraints and restraints, defined in our
    # FitModel, so we will turn it into a Fitness object. To do this, we have
    # written a special FitnessAdapter class in the diffpy.srfit.park package.
    f = FitnessAdapter(model)

    # Now we can fit this using park.
    from park.fitting.fit import fit
    print "Fit using the default PARK optimizer"
    result = fit([f])

    return


def plotResults(model):
    """Plot the results contained within a refined FitModel."""

    # Plot this.

    # We can access the data and fit profile through the Profile we created
    # above. We get to it through our Contribution, which we named "g1".
    #
    # The independent variable. This is always named "x"; the name we gave it
    # when we added the Profile to the Contribution was for the purpose of
    # writing the fit equation.
    x = model.g1.profile.x
    # The observed profile that we loaded earlier.
    y = model.g1.profile.y
    # The calculated profile.
    ycalc = model.g1.profile.ycalc

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

    # Create the model
    model = makeModel()

    # Refine using the optimizer of your choice
    scipyOptimize(model)
    #parkOptimize(model)

    # Get the results in a FitResults object. The FitResults object stores the
    # current state of the model, and uses it to calculate useful statistics
    # about the fit.  If you later modify the model, the FitResults object will
    # hold the model values from when it was created. You can tell it to update
    # its values by calling its 'update' method.
    res = FitResults(model)

    # Print the results
    res.printResults()

    # Plot the results
    plotResults(model)


# End of file
