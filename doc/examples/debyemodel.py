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

from diffpy.srfit.fitbase.calculator import Calculator
from diffpy.srfit.fitbase.profile import Profile

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

if __name__ == "__main__":

    # Demonstrate the Debye calculation

    xy = numpy.loadtxt("data/PbADPs.dat")
    x = xy[:,0]
    y = xy[:,1]

    profile = Profile()
    profile.setObservedProfile(x, y)

    debye = DebyeCalculator()
    debye.setProfile(profile)

    debye.m.setValue(207.2)
    debye.thetaD.setValue(240)
    debye.offset.setValue(0.0035)

    debye.eval()
    ycalc = profile.ycalc

    from pylab import plot, show
    plot(x, y, x, ycalc)
    show()

