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
"""Example of a calculator for calculating diffraction intensity using the
debye equation and the diffpy.Structure adapter.

"""

import os

import numpy

from diffpy.srfit.fitbase import Calculator, Contribution, FitModel, Profile
from diffpy.srfit.park import FitnessAdapter
from diffpy.srfit.structure import StructureParSet

class IntensityCalculator(Calculator):
    """A class for calculating adps from the Debye model."""

    def __init__(self):
        Calculator.__init__(self, "intensity")
        # We need a handle to the structure for it for the calculation.
        self.stru = None
        # Count the calls
        self.count = 0
        return

    def setStructure(self, strufile):
        """Set the structure used in the calculation.

        This will create the refinement parameters using the Structure adapter
        from diffpy.srfit.structure. Any changes made by the optimizer to the
        parameters go through the adapter, which propagates the changes to the
        diffpy.Structure.Structure object.
        """
        from diffpy.Structure import Structure
        self.stru = Structure()
        self.stru.read(strufile)
        parset = StructureParSet(self.stru, "structure")
        self.addParameterSet(parset)
        return

    def __call__(self, q):
        """Calculate the intensity.

        By the time this function is called the structure has been updated by
        the optimizer. Thus, we need only call iofq with the internal structure
        object.

        """
        self.count += 1
        print self.count
        return iofq(self.stru, q)

# End class IntensityCalculator

class IntensityContribution(Contribution):
    """A Contribution that uses the IntensityCalculator."""

    def __init__(self, name):
        """Auto-initialize the calculator."""
        Contribution.__init__(self, name)
        # We're naming the calculator "I". This is the name that will be used
        # in equations. If we left of the "I", then the name would default to
        # the name defined in the IntensityCalculator class, which is
        # "intensity".
        self.setCalculator(IntensityCalculator(), "I")
        return

    def setStructure(self, strufile):
        """Set the structure used in the calculation.

        Calls the calculator method of the same name.
        
        """
        self.intensity.setStructure(strufile)
        return

# End class IntensityContribution

def iofq(S, q):
    """Calculate I(Q) (X-ray) using the Debye Equation.

    I(Q) = 2 sum(i,j) f_i(Q) f_j(Q) sinc(rij Q) exp(-0.5 ssij Q**2)
    (The exponential term is the Debye-Waller factor.)

    S   --  A diffpy.Structure.Structure instance. It is assumed that the
            structure is that of an isolated scatterer. Periodic boundary
            conditions are not applied.
    q   --  The q-points to calculate over.

    The calculator uses cctbx for the calculation of the f_i if it is
    available, otherwise f_i = 1.

    """
    # The functions we need
    sinc = numpy.sinc
    exp = numpy.exp
    pi = numpy.pi

    # The brute-force calculation is very slow. Thus we optimize a little bit.

    # The precision of distance measurements
    deltad = 1e-5
    dmult = int(1/deltad)

    pairdict = {}
    elcount = {}
    n = len(S)
    for i in xrange(n):

        # count the number of each element
        eli = S[i].element
        m = elcount.get(eli, 0)
        elcount[eli] = m + 1

        for j in xrange(i+1,n):

            elj = S[j].element

            # Get the pair
            els = [eli, elj]
            els.sort()

            # Get the distance to the desired precision
            d = S.distance(i, j)
            D = int(d*dmult)

            # Get the DW factor to the same precision
            ss = S[i].Uisoequiv + S[j].Uisoequiv
            SS = int(ss*dmult)

            # Record the multiplicity of this pair
            key = (els[0], els[1], D, SS)
            mult = pairdict.get(key, 0)
            pairdict[key] = mult + 1

    # Now we can calculate IofQ from the pair dictionary. Making the dictionary
    # first reduces the amount of calls to sinc and exp we have to make.

    # First we must cache the scattering factors
    fdict = {}
    for el in elcount:
        fdict[el] = getXScatteringFactor(el, q)

    # Now we can compute I(Q) for the i != j pairs
    y = 0
    x = q * deltad / pi
    for key, mult in pairdict.items():
        eli = key[0]
        elj = key[1]
        fi = fdict[eli]
        fj = fdict[elj]
        D = key[2]
        SS = key[3]

        # Note that numpy's sinc(x) = sin(x*pi)/(x*pi)
        y += mult * sinc(x*D) * exp(-0.5*x*SS)

    # We must multiply by 2 since we only counted j > i pairs.
    y *= 2

    # Now we must add in the i == j pairs.
    for el, f in fdict.items():
        y += f * elcount[el]

    # And that's it!

    return y

def getXScatteringFactor(el, q):
    """Get the x-ray scattering factor for an element over the q range.
    
    If cctbx is not available, f(q) = 1 is used.

    """
    try:
        import cctbx.eltbx.xray_scattering as xray
        wk1995 = xray.wk1995(el)
        g = wk1995.fetch()
        # at_stol - at sin(theta)/lambda = Q/(4*pi)
        f = numpy.asarray( map( g.at_stol, q/(4*numpy.pi) ) )
        return f
    except ImportError:
        return 1

####### Example Code

def makeData(strufile, q, datname):
    """Make some fake data and save it to file.

    Use the iofq function to make some data. This adds noise and a background
    to the data and saves it to the specified filename. It only generates the
    data if it does not exist.

    strufile--  A filename holding the sample structure
    q       --  The q-range to calculate over.
    datname --  The name of the file we're saving to.

    """

    import os.path
    if os.path.exists(datname): return

    # Expand the lattice by +/= 5%
    from diffpy.Structure import Structure
    S = Structure()
    S.read(strufile)
    import random
    expand = 10*(0.5 - random.random())
    a = (1 + expand/100) * S.lattice.a
    S.lattice.setLatPar(a, a, a)
    y = iofq(S, q)

    # Add a polynomial background.
    bkgd = (q + 1)**2 * (1.5*max(q) - q)**5
    bkgd *= 0.2 * max(y) / max(bkgd)

    y += bkgd

    # Now add uniform noise at +/-2% of the max intensity
    nrange = 0.04*max(y)
    noise = numpy.zeros_like(q)
    for i in xrange(len(q)):
        noise[i] = (random.random() - 0.5)* nrange

    y += noise

    # Multipy by an arbitrary scale factor
    y *= random.uniform(1, 15)

    # Calculate the uncertainty (uniform distribution)
    u = numpy.ones_like(q)* nrange / 12**0.5

    # Now save it
    numpy.savetxt(datname, zip(q, y, u))
    return

def makeModel(strufile, datname):
    """Create a model that encompasses a profile calculator and an experimental
    profile. The model can then be optimized using an external optimizer.

    We need to fit a background to I(Q), but our calculator doesn't have one.
    We can easily expand the calculator by adding the background we want to
    fit.

    """

    # Create a Profile. This will hold the experimental and calculated signal.
    profile = Profile()

    # Load data and add it to the profile
    x, y, u = numpy.loadtxt(datname, unpack=True)
    profile.setObservedProfile(x, y, u)


    # Create the Contribution, that will associate the profile with the
    # Calculator. We have created a custom IntensityContribution above that
    # performs half of the association for us. We also want to tell the
    # contribution to name the x-variable of the profile "q", so we can use it
    # in equations with this name.
    contribution = IntensityContribution("C60")
    contribution.setProfile(profile, xname = "q")

    # We need to add the Structure we want to fit to the contribution.
    contribution.setStructure(strufile)

    # We want to add an arbitrary background. We'll assume we need 10
    # coefficients.
    for i in xrange(10):
        contribution.newParameter("b%i"%i, 0)

    # We also need a scale factor
    contribution.newParameter("scale", 1)

    # Write the equation for the background. We will write a python function
    # for this and add it to the contribution. As long as we use the same names
    # as things we have already defined (q for the x-name, b0, b1, ... for the
    # coefficients), then the contribution will be able to use it.
    def bkgd(q, b0, b1, b2, b3, b4, b5, b6, b7, b8, b9):
        return b0 + b1*q + b2*q**2 + b3*q**3 + b4*q**4 + b5*q*5 +\
                b6*q**6 + b7*q**7 +b8*q**8 + b9*q**9

    contribution.registerFunction("bkgd", bkgd)

    # Now we can incorporate the scale and bkgd into our calculation
    contribution.setEquation("scale * I() + bkgd()")

    # Make a FitModel where we can create variables, constraints and
    # restraints. If we had multiple profiles to fit simultaneously, the
    # contribution from each could be added to the model.
    model = FitModel()
    model.addContribution(contribution)

    # Specify which parameters we want to refine. We can give them initial
    # values in the process. We want to refine the background variables that we
    # just defined in the contribution.
    model.addVar(contribution.b0, 0)
    model.addVar(contribution.b1, 0)
    model.addVar(contribution.b2, 0)
    model.addVar(contribution.b3, 0)
    model.addVar(contribution.b4, 0)
    model.addVar(contribution.b5, 0)
    model.addVar(contribution.b6, 0)
    model.addVar(contribution.b7, 0)
    model.addVar(contribution.b8, 0)
    model.addVar(contribution.b9, 0)

    # We also want to adjust the scale
    model.addVar(contribution.scale, 1)


    # We can also refine structural parameters. 
    structure = contribution.intensity.structure
    a = structure.lattice.a
    model.addVar(a)
    # We want to allow for isotropic expansion, so we'll make constraints for
    # that.
    model.constrain(structure.lattice.b, a)
    model.constrain(structure.lattice.c, a)

    # Give the model away so it can be used!
    return model

def scipyOptimize(strufile):
    """Optimize the model created above using scipy."""

    # Make the data and the model
    q = numpy.arange(1, 20, 0.05)
    makeData(strufile, q, "C60.iq")

    model = makeModel(strufile, "C60.iq")

    # We're going to use the least-squares (Levenberg-Marquardt) optimizer from
    # scipy.
    from scipy.optimize.minpack import leastsq
    from numpy import dot
    out = leastsq(model.residual, model.getValues(), full_output=1)

    # For the basic info about the fit, we can use the FitModel directly
    chiv = model.residual()
    names = model.getNames()
    vals = model.getValues()

    print "Fit using the scipy LM optimizer"
    chi2 = numpy.dot(chiv, chiv)
    rchi2 = chi2 / (len(q) - len(vals))
    print "Chi^2 =", chi2
    print "Red. Chi^2 =", rchi2
    
    for name, val in zip(names, vals):
        print "%s = %f" % (name, val)

    # Plot this for fun.
    q = model.C60.profile.x
    I = model.C60.profile.y
    Icalc = model.C60.profile.ycalc
    bkgd = model.C60.evaluateEquation("bkgd()")

    import pylab
    pylab.plot(q,I,'o',label="I(Q) Data")
    pylab.plot(q,Icalc,label="I(Q) Fit")
    pylab.plot(q,bkgd,label="Bkgd. Fit")
    pylab.xlabel("$Q (\AA^{-1})$")
    pylab.ylabel("Intensity (arb. units)")
    pylab.legend(loc=1)

    pylab.show()

    return

def parkOptimize(strufile):
    """Optimize the model created above using PARK."""

    # Make the data and the model
    q = numpy.arange(1, 20, 0.05)
    makeData(strufile, q, "C60.iq")

    model = makeModel(strufile, "C60.iq")

    # For the default optimizer, it is much easier to find a solution if we put
    # bounds on any parameter we know. We do this using our model, but we
    # could also do it using the park Fitness object that we create below.
    model.confine( model.a, 9, 11)
    model.confine( model.scale, 0.5, 20)
    for i in xrange(10):
        model.confine(getattr(model, "b%i"%i), -500, 500)

    # We have to turn the model into something that PARK can use. In PARK, a
    # Fitness object is the equivalent of a SrFit Contribution. However, we
    # want a very lean interface to any optimizer, so we treat the entire
    # FitModel as a Fitness object. To do this, we have written a special
    # FitnessAdapter class in the diffpy.srfit.park package.
    f = FitnessAdapter(model)

    # Now we can fit this
    from park.fitting.fitresult import ConsoleUpdate
    from park.optim.fitmc import FitMC
    handler = ConsoleUpdate(improvement_delta=0.1,progress_delta=1)
    fitter = FitMC(start_points=1)
    from park.fitting.fit import fit
    result = fit([f], handler=handler, fitter=fitter)

    # For the basic info about the fit, we can use the FitModel directly
    chiv = model.residual()
    names = model.getNames()
    vals = model.getValues()

    print "Fit using the default PARK optimizer"
    chi2 = numpy.dot(chiv, chiv)
    rchi2 = chi2 / (len(q) - len(vals))
    print "Chi^2 =", chi2
    print "Red. Chi^2 =", rchi2
    
    for name, val in zip(names, vals):
        print "%s = %f" % (name, val)

    # Plot this for fun.
    q = model.C60.profile.x
    I = model.C60.profile.y
    Icalc = model.C60.profile.ycalc
    bkgd = model.C60.evaluateEquation("bkgd()")

    import pylab
    pylab.plot(q,I,'o',label="I(Q) Data")
    pylab.plot(q,Icalc,label="I(Q) Fit")
    pylab.plot(q,bkgd,label="Bkgd. Fit")
    pylab.xlabel("$Q (\AA^{-1})")
    pylab.ylabel("Intensity (arb. units)")
    pylab.legend(loc=1)

    pylab.show()

    return

if __name__ == "__main__":

    print "Optimize with scipy"
    scipyOptimize("data/C60.stru")
    print "Optimize with park"
    parkOptimize("data/C60.stru")

# End of file
