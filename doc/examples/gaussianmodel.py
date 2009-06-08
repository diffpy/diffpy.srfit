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
"""Example of fitting A gaussian to experimental Debye-Waller factors.

This is an example of building a FitModel in order to fit experimental data.

The makeModel function shows how to build a FitModel that will fit our model to
the data. The scipyOptimize and parkOptimize functions show two different ways
of refining the model, using scipy and park, respectively.

Once you understand this, move on to the debyemodelII.py example.
"""

import numpy

from diffpy.srfit.fitbase import Contribution, FitModel, Profile, FitResults
from diffpy.srfit.park import FitnessAdapter

####### Example Code

def makeModel():
    """Make the model for our problem.

    Our model will be defined within a FitModel instance. The job of a FitModel
    is to collect and associate all the data, the fitting equations, fitting
    variables, constraints and restrations. We will demonstrate each of these
    within the code. 

    Data is held within a Profile object. The Profile is simply a container
    that holds the data, and the theoretical profile once it has been
    calculated.

    Data is associated with a fitting equation within a Contribution. The
    Contribution defines the equation and parameters that will be adjusted to
    fit the data. The fitting equation can be defined within a function or
    optionally within the Calculator class. We won't need the Calculator class
    in this example since the signature of the fitting equationis so simple.
    The contribution also defines the residual function to optimize for the
    data/equation pair. This can be modified, but we won't do that here.

    Once we define the FitModel, we can send it an optimizer to be optimized.
    See the scipyOptimize and parkOptimize functions.
    
    """

        
    ## The Profile
    # Create a Profile to hold the experimental and calculated signal.
    profile = Profile()

    # Load data and add it to the profile. It is our responsibility to get our
    # data into the profile.
    x, y = numpy.loadtxt("data/gaussian.dat", unpack=1)
    profile.setObservedProfile(x, y)

    ## The Contribution
    # The Contribution associates the profile with the Debye Calculator. 
    contribution = Contribution("g1")
    # Tell the contribution about the Profile. We will need to use the
    # independent variable (the temperature) from the data to calculate the
    # theoretical signal, so give it an informative name ('T') that we can use
    # later.
    contribution.setProfile(profile, xname="x")

    # We have a function registered to the contribution, but we have yet to
    # define the fitting equation.
    contribution.setEquation("A * exp(-0.5*(x-x0)**2/sigma**2)")

    ## The FitModel
    # The FitModel lets us define what we want to fit. It is where we can
    # create variables, constraints and restraints. If we had multiple profiles
    # to fit simultaneously, the contribution from each could be added to the
    # model.
    model = FitModel()
    model.addContribution(contribution)

    # Specify which parameters we want to refine. We can give them initial
    # values in the process. This tells the model to vary the offset and to
    # give it an initial value of 0.
    model.addVar(contribution.A, 1)
    model.addVar(contribution.x0, 5)
    model.addVar(contribution.sigma, 1)

    return model

def scipyOptimize(model):
    """Optimize the model created above using scipy.

    The FitModel we created in makeModel has a 'residual' method that we can be
    minimized using a scipy optimizer. The details are described in the source.

    """

    # We're going to use the least-squares (Levenberg-Marquardt) optimizer from
    # scipy. We simply have to give it the function to minimize
    # (model.residual) and the starting values of the variables
    # (model.getValues()). We defined offset and tvar as variables, so that is
    # what will be adjusted by the optimizer.
    from scipy.optimize.minpack import leastsq
    print "Fit using scipy's LM optimizer"
    leastsq(model.residual, model.getValues())

    return

def parkOptimize(model):
    """Optimize the model created above using PARK."""

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
    # Note that since the contribution was given the name "pb", it is
    # accessible from the model with this name. This is a useful way to
    # organize multiple contributions to a fit.
    x = model.g1.profile.x
    y = model.g1.profile.y
    ycalc = model.g1.profile.ycalc

    import pylab
    pylab.plot(x, y, x, ycalc)

    pylab.show()
    return


if __name__ == "__main__":

    #model = makeModel()
    #scipyOptimize(model)
    #res = FitResults(model)
    #res.printResults()
    #plotResults(model)

    # Start from scratch
    model = makeModel()
    parkOptimize(model)
    res = FitResults(model)
    res.printResults()
    plotResults(model)


# End of file
