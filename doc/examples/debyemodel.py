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
"""Example of a calculator for the debye model."""

from functools import partial

import numpy
import scipy.integrate

from diffpy.srfit.fitbase import Calculator, Contribution, FitModel, Profile
from diffpy.srfit.park import FitnessAdapter

class DebyeCalculator(Calculator):
    """A class for calculating adps from the Debye model."""

    def __init__(self):
        Calculator.__init__(self, "debye")
        self._newParameter("m", 12)
        self._newParameter("thetaD", 300)
        self._newParameter("offset", 0)
        return

    def __call__(self, x):

        m = self.m.getValue()
        tD = self.thetaD.getValue()
        offset = self.offset.getValue()

        y = [adps(m, tD, T) + offset for T in x]

        return y

# End class DebyeCalculator

class DebyeContribution(Contribution):
    """A Contribution that uses the DebyeCalculator."""

    def __init__(self, name):
        """Auto-initialize the calculator."""
        Contribution.__init__(self, name)
        self.setCalculator(DebyeCalculator())
        return

# End class DebyeContribution


# Functions required for calculation of Debye curve

def adps(m,thetaD,T):
    """calculates atomic displacement factors within the Debye model

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

    m = m * amu
    u2 = (3*h**2 / (4 * numpy.pi**2 *m *kB *thetaD))*(__phi(thetaD/T)/(thetaD/T) + 1./4.)

    return u2*1e20

def __debyeKernel(xi):
    """function needed by debye calculators

    """
    y = xi/(numpy.exp(xi)-1)
    return y

def __phi(x):
    """evaluates the phi integral needed in Debye calculation

    phi(x) = (1/x) int_0^x xi/(exp(xi)-1) dxi

    arguments:
    x -- float -- value of thetaD (Debye temperature)/T

    returns:
    phi -- float -- value of the phi function

    """
    int = scipy.integrate.quad(__debyeKernel, 0, x)
    phi = (1/x) * int[0]

    return phi

####### Example Code

def makeModel():
    """Create a model that encompasses a profile calculator and an experimental
    profile. The model can then be optimized using an external optimizer.

    """

    # Create a Profile. This will hold the experimental and calculated signal.
    profile = Profile()

    # Load data and add it to the profile
    x, y = numpy.loadtxt("data/PbADPs.dat", unpack=True)
    profile.setObservedProfile(x, y)


    # Create the Contribution, that will associate the profile with the Debye
    # Calculator. We have created a custom DebyeContribution above that
    # performs half of the association for us.
    contribution = DebyeContribution("pb")
    contribution.setProfile(profile)

    # Make a FitModel where we can create variables, constraints and
    # restraints. If we had multiple profiles to fit simultaneously, the
    # contribution from each could be added to the model.
    model = FitModel()
    model.addContribution(contribution)

    # Give the calculator parameters values.  The DebyeCalculator has
    # the name "debye", so we can access its variables from the "debye"
    # attribute of the contribution.
    #
    # We know the mass of lead, so we'll set that here. We're varying the other
    # two parameters, so we'll give them initial values below.
    contribution.debye.m.setValue(207.2)

    # Specify which parameters we want to refine. We can give them initial
    # values in the process. 
    model.addVar(contribution.debye.offset, 0)

    # We will handle the thetaD parameter in a convoluted, yet instructive way.
    # We want this to be positive, so we'll create a new variable to refine,
    # named "tvar" and constrain the thetaD to be the absolute of this variable.
    model.newVar("tvar", 300)
    model.constrain(contribution.debye.thetaD, "abs(tvar)")

    # While we're at it, let's keep the offset positive. We could do the
    # constraint method above, but we'll use a restraint instead. This
    # restraint will add infinity to the chi^2 if the offset goes negative.
    from numpy import inf
    model.restrain(contribution.debye.offset, lb=0, ub=inf)

    # Give the model away so it can be used!
    return model

def scipyOptimize():
    """Optimize the model created above using scipy."""

    model = makeModel()

    # We're going to use the least-squares (Levenberg-Marquardt) optimizer from
    # scipy.
    from scipy.optimize.minpack import leastsq
    from numpy import dot
    out = leastsq(model.residual, model.getValues())
    offset, tvar = out[0]

    chiv = model.residual()
    print "Fit using scipy's leastsq optimizer"
    print "Chi^2 = ", numpy.dot(chiv, chiv)
    # Get the refined variable values, noting that we didn't refine thetaD
    # directly. If we want uncertainties, we have to go to the optimizer
    # directly.
    offset, tvar = model.getValues()

    print "tvar =", tvar
    print "offset =", offset

    return

def parkOptimize():
    """Optimize the model created above using PARK."""

    model = makeModel()

    # We have to turn the model into something that PARK can use. In PARK, a
    # Fitness object is the equivalent of a SrFit Contribution. However, we
    # want a very lean interface to any optimizer, so we treat the entire
    # FitModel as a Fitness object. To do this, we have written a special
    # FitnessAdapter class in the diffpy.srfit.park package.
    f = FitnessAdapter(model)

    # Now we can fit this
    from park.fitting.fit import fit
    result = fit([f])

    # For the basic info about the fit, we can use the FitModel directly
    chiv = model.residual()

    print "Fit using the default PARK optimizer"
    print "Chi^2 = ", numpy.dot(chiv, chiv)
    # Get the refined variable values, noting that we didn't refine thetaD
    # directly. If we want uncertainties, we have to go to the optimizer
    # directly.
    offset, tvar = result.pvec

    print "tvar =", tvar
    print "offset =", offset

    return

if __name__ == "__main__":

    scipyOptimize()
    parkOptimize()
    print \
"""Note that the solutions are equivalent (to several digits). We cannot assess
the parameter uncertainty without uncertainties on the data points.\
"""


# End of file
