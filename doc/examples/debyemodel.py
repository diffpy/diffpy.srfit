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
"""Example of fitting the Debye model to experimental Debye-Waller factors.

This is an example of building a FitModel in order to fit experimental data.
It is assumed that the function we need cannot be modified by us (although we
define it below). This will help us demonstrate how to extend a function using
SrFit.

The makeModel function shows how to build a FitModel that will fit our model to
the data. The scipyOptimize and parkOptimize functions show two different ways
of refining the model, using scipy and park, respectively.

Once you understand this, move on to the intensitycalculator.py example.
"""


from functools import partial

import numpy

from diffpy.srfit.fitbase import Contribution, FitModel, Profile
from diffpy.srfit.park import FitnessAdapter

# Functions required for calculation of Debye curve. Feel free to skip these,
# as we treat them as if existing in some external library that we cannot
# modify.
def adps(m,thetaD,T):
    """Calculates atomic displacement factors within the Debye model

    <u^2> = (3h^2/4 pi^2 m kB thetaD)(phi(thetaD/T)/(ThetaD/T) + 1/4)

    arguments:
    m -- float -- mass of the ion in atomic mass units (e.g., C = 12)
    thetaD -- float -- Debye Temperature
    T -- float -- temperature.

    return:
    Uiso -- float -- the thermal factor from the Debye model at temp T

    """
    h = 6.6260755e-34   # Planck's constant. J.s of m^2.kg/s
    kB = 1.3806503e-23  # Boltzmann's constant. J/K
    amu = 1.66053886e-27 # Atomic mass unit. kg

    def __phi(x):
        """evaluates the phi integral needed in Debye calculation

        phi(x) = (1/x) int_0^x xi/(exp(xi)-1) dxi

        arguments:
        x -- float -- value of thetaD (Debye temperature)/T

        returns:
        phi -- float -- value of the phi function

        """
        def __debyeKernel(xi):
            """function needed by debye calculators

            """
            y = xi/(numpy.exp(xi)-1)
            return y

        import scipy.integrate

        int = scipy.integrate.quad(__debyeKernel, 0, x)
        phi = (1/x) * int[0]

        return phi


    m = m * amu
    u2 = (3*h**2 / (4 * numpy.pi**2 *m *kB *thetaD))*(__phi(thetaD/T)/(thetaD/T) + 1./4.)

    return u2*1e20

def debye(T, m, thetaD):
    """A wrapped version of 'adps' that can handle an array of T-values."""
    y = numpy.array([adps(m, thetaD, x) for x in T])
    return y

# The data
data = """\
012.0 0.00244
050.0 0.00366
100.0 0.00541
150.0 0.00779
200.0 0.01019
250.0 0.01200
300.0 0.01408
350.0 0.01572
400.0 0.01775
450.0 0.01946
500.0 0.02275
"""

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
    in this example since the signature of the fitting equation (the 'debye'
    function) is so simple. The contribution also defines the residual function
    to optimize for the data/equation pair. This can be modified, but we won't
    do that here.

    Once we define the FitModel, we can send it an optimizer to be optimized.
    See the scipyOptimize and parkOptimize functions.
    
    """

        
    ## The Profile
    # Create a Profile to hold the experimental and calculated signal.
    profile = Profile()

    # Load data and add it to the profile. It is our responsibility to get our
    # data into the profile.
    xy = numpy.array(data.split(), dtype=float)
    x = xy[0::2]
    y = xy[1::2]
    profile.setObservedProfile(x, y)

    ## The Contribution
    # The Contribution associates the profile with the Debye Calculator. 
    contribution = Contribution("pb")
    # Tell the contribution about the Profile. We will need to use the
    # independent variable (the temperature) from the data to calculate the
    # theoretical signal, so give it an informative name ('T') that we can use
    # later.
    contribution.setProfile(profile, xname="T")

    # We now want to tell the contribution to use the 'debye' function defined
    # above. The 'registerFunction' method will let us do this. Since we
    # haven't told it otherwise, 'registerFunction' will extract the name of
    # the function ('debye') and the names of the arguments ('T', 'm',
    # 'thetaD'). (Note that we could have given it other names.) Since we named
    # the x-variable 'T' above, the 'T' in the 'debye' equation will refer  to
    # this x-variable when it gets called.
    contribution.registerFunction(debye)

    # We have a function registered to the contribution, but we have yet to
    # define the fitting equation. On top of that, we need a vertical offset in
    # our equation that does not appear in 'debye'.  We could have written a
    # function that calls 'debye' that includes an offset, but this will not
    # always be an option. Thus, we will add the offset when we define the
    # equation.  We don't need to specify the parameters to the 'debye'
    # function since the contribution already knows what they are.  (In fact,
    # the '()' is optional.)
    contribution.setEquation("debye()+offset")

    # While we're at it, we should give some initial values to parameters.  We
    # know the mass of lead, so we'll set that here. Note that we can access
    # the fitting parameters by name as attributes of the contribution.  We're
    # varying the other two parameters ('thetaD' and 'offset'), so we'll give
    # them initial values below. Parameters have a 'setValue' and 'getValue'
    # method.
    contribution.m.setValue(207.2)

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
    model.addVar(contribution.offset, 0)

    # We will handle thetaD in a convoluted way to demonstrate how to use
    # constraints. We want thetaD to be positive, so we'll create a new fit
    # variable named "tvar" and constrain thetaD to be the absolute value of
    # this variable. In this way thetaD will be be given the value abs(tvar)
    # whenever tvar is adjusted.
    model.newVar("tvar", 300)
    model.constrain(contribution.thetaD, "abs(tvar)")

    # While we're at it, let's keep the offset positive. We could do the
    # constraint method above, but we'll use a bounds restraint in order to
    # demonstrate the syntax. This restraint will add infinity to model's cost
    # function if the offset becomes negative. (Any resonable optimizer will
    # therefore keep offset >= 0).
    from numpy import inf
    model.confine(contribution.offset, 0, inf)

    # Return the  model. See the scipyOptimize and parkOptimize functions to
    # see how it is used.
    return model

def scipyOptimize():
    """Optimize the model created above using scipy.

    The FitModel we created in makeModel has a 'residual' method that we can be
    minimized using a scipy optimizer. The details are described in the source.

    """

    model = makeModel()

    # We're going to use the least-squares (Levenberg-Marquardt) optimizer from
    # scipy. We simply have to give it the function to minimize
    # (model.residual) and the starting values of the variables
    # (model.getValues()). We defined offset and tvar as variables, so that is
    # what will be adjusted by the optimizer.
    from scipy.optimize.minpack import leastsq
    print "Fit using scipy's LM optimizer"
    leastsq(model.residual, model.getValues())


    # We'll look at the results in another function.
    displayResults(model)
    return

def parkOptimize():
    """Optimize the model created above using PARK."""

    model = makeModel()

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

    displayResults(model)

    return

def displayResults(model):
    """Display the results contained within a refined FitModel."""

    # For the basic info about the fit, we can use the FitModel directly
    chiv = model.residual()

    print "Chi^2 = ", numpy.dot(chiv, chiv)
    # Get the refined variable values, noting that we didn't refine thetaD
    # directly. If we want uncertainties, we have to go to the optimizer
    # directly.
    offset, tvar = model.getValues()

    print "tvar =", tvar
    print "offset =", offset
    
    # Plot this.
    # Note that since the contribution was given the name "pb", it is
    # accessible from the model with this name. This is a useful way to
    # organize multiple contributions to a fit.
    T = model.pb.profile.x
    U = model.pb.profile.y
    Ucalc = model.pb.profile.ycalc

    import pylab
    pylab.plot(T,U,'o',label="Pb $U_{iso}$ Data")
    lbl1 = "$T_d$=%3.1f K, off=%1.5f $\AA^2$"% (abs(tvar),offset)
    pylab.plot(T,Ucalc,label=lbl1)
    pylab.xlabel("T (K)")
    pylab.ylabel("$U_{iso} (\AA^2)$")
    pylab.legend(loc = (0.0,0.8))

    pylab.show()
    return

if __name__ == "__main__":

    scipyOptimize()
    parkOptimize()
    print \
"""\nNote that the solutions are equivalent (to several digits). We cannot assess
the parameter uncertainty without uncertainties on the data points.\
"""


# End of file
