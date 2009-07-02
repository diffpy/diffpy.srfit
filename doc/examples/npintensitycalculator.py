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
"""Example of using Calculators in FitModels.

This is an example of building a Calculator and using it in a FitModel in order
to fit theoretical intensity data.

The IntensityCalculator class is an example of a Calculator that can be used by
a Contribution to help generate a signal.

The makeModel function shows how to build a FitModel that uses the
IntensityCalculator.

"""

import os

import numpy

from diffpy.srfit.fitbase import Calculator, Contribution, FitModel, Profile
from diffpy.srfit.fitbase import FitResults
from diffpy.srfit.park import FitnessAdapter
from diffpy.srfit.structure.diffpystructure import StructureParSet

from gaussianmodel import scipyOptimize, parkOptimize

class IntensityCalculator(Calculator):
    """A class for calculating intensity using the Debye equation.

    Calculating intensity from a structure is difficult in general. This class
    takes a diffpy.Structure.Structure object and from that calculates a
    theoretical intensity signal. Unlike the example in gaussianmodel.py, the
    intensity calculator is not simple, so we must define this Calculator to
    help us interface with a FitModel.

    The purpose of a Calculator is to
    1) provide a function that calculates a profile signal
    2) organize the Parameters required for the calculation

    This calculator wraps the 'iofq' function defined below.
    
    """

    def __init__(self, name):
        """Define our calculator.

        Keep count of how many times the function has been called (self.count).

        """
        Calculator.__init__(self, name)
        # Count the calls
        self.count = 0
        return

    def setStructure(self, strufile):
        """Set the structure used in the calculation.

        This will create the refinement Parameters using the Structure adapter
        from diffpy.srfit.structure. The created Parameters are proxies for
        attributes of the Structure object that can be interfaced within SrFit.
        The Parameters will be accessible by name under the 'structure'
        attribute of this calculator, and are organized hierarchically:
        structure
          - lattice
          - atom1 (the name depends on the element)
            - x
            - y
            - z
            - occ
            - U11
            - U22
            - U33
            - U12
            - U13
            - U23
            - Uiso
          - ...
          - atomN

        See the makeModel code to see how these Parameters are accessed.
        
        """
        # Load the structure from file
        from diffpy.Structure import Structure
        stru = Structure()
        stru.read(strufile)

        # Create a ParameterSet designed to interface with
        # diffpy.Structure.Structure objects that organizes the Parameter
        # hierarchy. Note that the StructureParSet holds a handle to the loaded
        # structure that we use in the __call__ method below.
        parset = StructureParSet(stru, "structure")
        # Put this ParameterSet in the Calculator.
        self.addParameterSet(parset)
        return

    def __call__(self, q):
        """Calculate the intensity.

        This Calculator will be used in a contribution equation that will be
        optimized to fit some data.  By the time this function is evaluated,
        the structure has been updated by the optimizer via the ParameterSet
        defined in setStructure. Thus, we need only call iofq with the internal
        structure object.

        """
        self.count += 1
        print "iofq called", self.count
        return iofq(self.structure.stru, q)

# End class IntensityCalculator

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
    deltad = 1e-6
    dmult = int(1/deltad)
    deltau = deltad**2
    umult = int(1/deltau)

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
            SS = int(ss*umult)

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
        y += fi * fj * mult * sinc(x * D) * exp(-0.5 * SS * deltau * q**2)

    # We must multiply by 2 since we only counted j > i pairs.
    y *= 2

    # Now we must add in the i == j pairs.
    for el, f in fdict.items():
        y += f**2 * elcount[el]

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

def makeData(strufile, q, datname, scale, a, Uiso, sig, bkgc, nl = 1):
    """Make some fake data and save it to file.

    Make some data to fit. This uses iofq to calculate an intensity curve, and
    adds to it a background, broadens the peaks, and noise.

    strufile--  A filename holding the sample structure
    q       --  The q-range to calculate over.
    datname --  The name of the file we're saving to.
    scale   --  The scale factor
    a       --  The lattice constant to use
    Uiso    --  The thermal factor for all atoms
    sig     --  The broadening factor
    bkgc    --  A parameter that gives minor control of the background.
    nl      --  Noise level (0, inf), default 1, larger -> less noise.

    """

    from diffpy.Structure import Structure
    S = Structure()
    S.read(strufile)

    # Set the lattice parameters
    S.lattice.setLatPar(a, a, a)

    # Set a DW factor
    for a in S:
        a.Uisoequiv = Uiso
    y = iofq(S, q)

    # We want to broaden the peaks as well. This simulates instrument effects.
    q0 = q[len(q)/2]
    g = numpy.exp(-0.5*((q-q0)/sig)**2)
    y = numpy.convolve(y, g, mode='same')/sum(g)

    # Add a polynomial background.
    bkgd = (q + bkgc)**2 * (1.5*max(q) - q)**5
    bkgd *= 0.2 * max(y) / max(bkgd)

    y += bkgd

    # Multipy by a scale factor
    y *= scale

    # Calculate the uncertainty 
    u = (y/nl)**0.5

    # And apply the noise
    if nl > 0:
        y = numpy.random.poisson(y*nl) / nl

    # Now save it
    numpy.savetxt(datname, zip(q, y, u))
    return

####### Example Code

def makeModel(strufile, datname):
    """Create a model that uses the IntensityCalculator.

    This will create a Contribution that uses the IntensityCalculator,
    associate this with a Profile, and use this to define a FitModel.

    """

    ## The Profile
    # Create a Profile. This will hold the experimental and calculated signal.
    profile = Profile()

    # Load data and add it to the profile
    x, y, u = numpy.loadtxt(datname, unpack=True)
    profile.setObservedProfile(x, y, u)

    ## The Calculator
    # Create an IntensityCalculator named "I". This will be the name we use to
    # refer to the calculator from within the Contribution equation.  We also
    # need to load the model structure we're using.
    calculator = IntensityCalculator("I")
    calculator.setStructure(strufile)
    
    ## The Contribution
    # Create a Contribution, that will associate the Profile with the
    # Calculator.  The calculator will be accessible as an attribute of the
    # Contribution by its name ("I"), or simply by "calculator".  We also want
    # to tell the contribution to name the x-variable of the profile "q", so we
    # can use it in equations with this name.
    contribution = Contribution("bucky")
    contribution.setCalculator(calculator)
    contribution.setProfile(profile, xname = "q")

    # Now we're ready to define the contribution equation. We need to modify
    # the Calcultor, and we'll do that from within the Contribution eqation for
    # the sake of instruction. We want to modify the calculator in three ways.
    # We need a scale factor, a polynomial background, and we want to broaden
    # the peaks. 
    #
    # There is added benefit for defining these operations outside of the
    # IntensityCalculator. By combining the different parts of the calculation
    # within the contribution equation, the time-consuming iofq calculation is
    # only performed when a structural parameter is changed. If only
    # non-structural parameters are changed, such as the background and
    # broadening parameters, then then previously computed iofq value will be
    # used to compute the contribution equation.  The benefit in this is very
    # apparent when refining the model with the LM optimizer, which only
    # changes two variables at a time most of the time. Note in the refinement
    # output how many times the residual is calculated, versus how many times
    # iofq is called when using the scipyOptimize function.

    # We will define the background as a string.

    bkgdstr = "b0 + b1*q + b2*q**2 + b3*q**3 + b4*q**4 + b5*q*5 + b6*q**6 +\
               b7*q**7 + b8*q**8 + b9*q**9"

    # This creates a callable equation within the Contribution, and turns
    # the polynomial coefficients into Parameters.
    contribution.registerStringFunction(bkgdstr, "bkgd")

    # We will create the broadening function that we need by creating a
    # python function and registering it with the Contribution.
    pi = numpy.pi
    exp = numpy.exp
    def gaussian(q, q0, width):
        return 1/(2*pi*width**2)**0.5 * exp(-0.5 * ((q-q0)/width)**2)

    # This registers the python function and extracts the name and creates
    # Parameters from the arguments.
    contribution.registerFunction(gaussian)

    # Center the Gaussian.
    contribution.q0.setValue(x[len(x)/2])

    # Now we can incorporate the scale and bkgd into our calculation. We also
    # convolve the signal with the gaussian to broaden it. Recall that we don't
    # need to supply arguments to the registered functions unless we want to
    # make changes to their input values.
    contribution.setEquation("scale * convolve(I, gaussian) + bkgd")

    # Make the FitModel and add the Contribution.
    model = FitModel()
    model.addContribution(contribution)

    # Specify which parameters we want to refine.
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

    # We also want to adjust the scale and the convolution width
    model.addVar(contribution.scale, 1)
    model.addVar(contribution.width, 0.1)

    # We can also refine structural parameters. Here we extract the 'structure'
    # ParameterSet from the intensity calculator and use the parameters like we
    # would any others.
    structure = calculator.structure

    # We want to allow for isotropic expansion, so we'll constrain the lattice
    # parameters to the same value.
    a = structure.lattice.a
    model.addVar(a)
    model.constrain(structure.lattice.b, a)
    model.constrain(structure.lattice.c, a)
    # We want to refine the thermal paramters as well. We will add a new
    # variable that we call "Uiso" and constrain the atomic Uiso values to
    # this. The structure ParameterSet has an 'atoms' list that we can use to
    # make this easier.
    Uiso = model.newVar("Uiso", 0.01)
    for atom in structure.atoms:
        model.constrain(atom.Uiso, Uiso)

    # Give the model away so it can be used!
    return model

def plotResults(model):
    """Plot the results contained within a refined FitModel."""

    # All this should be pretty familiar by now.
    names = model.getNames()
    vals = model.getValues()

    q = model.bucky.profile.x

    I = model.bucky.profile.y
    Icalc = model.bucky.profile.ycalc
    bkgd = model.bucky.evaluateEquation("bkgd")
    diff = I - Icalc

    import pylab
    pylab.plot(q,I,'o',label="I(Q) Data")
    pylab.plot(q,Icalc,label="I(Q) Fit")
    pylab.plot(q,diff,label="I(Q) diff")
    pylab.plot(q,bkgd,label="Bkgd. Fit")
    pylab.xlabel("$Q (\AA^{-1})$")
    pylab.ylabel("Intensity (arb. units)")
    pylab.legend(loc=1)

    pylab.show()
    return

if __name__ == "__main__":

    # Make the data and the model
    strufile = "data/C60.stru"
    q = numpy.arange(1, 20, 0.05)
    makeData(strufile, q, "C60.iq", 1.0, 100.68, 0.005, 0.13, 2)

    # Make the model
    model = makeModel(strufile, "C60.iq")

    # Optimize
    scipyOptimize(model)
    #parkOptimize(model)

    # Generate and print the FitResults
    res = FitResults(model)
    # Get the number of calls to the residual function from the FitModel, and
    # the number of calls to 'iofq' from the IntensityCalculator.
    rescount = model.fithook.count
    calcount = model.bucky.calculator.count
    footer = "iofq called %i%% of the time"%int(100.0*calcount/rescount)
    res.printResults(footer = footer)

    # Plot!
    plotResults(model)

# End of file
