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

This is an example of building a FitModel in order to fit theoretical intensity
data. To understand the basics of FitModels, see the debyemodel.py.

The CCTBXIntensityCalculator class is an example of a Calculator that can be
used by a Contribution to help generate a signal.

The makeModel function shows how to build a FitModel that uses the
CCTBXIntensityCalculator. 
"""

import os

import numpy

from cctbx import xray
from cctbx import crystal
from cctbx.array_family import flex

from diffpy.srfit.fitbase import Calculator, Contribution, FitModel, Profile
from diffpy.srfit.fitbase import FitResults
from diffpy.srfit.park import FitnessAdapter
from diffpy.srfit.structure import CCTBXStructureParSet

from npintensitycalculator import makeData, getXScatteringFactor
from debyemodel import scipyOptimize, parkOptimize

class CCTBXIntensityCalculator(Calculator):
    """A class for calculating intensity using the Debye equation.

    Calculating intensity from a structure is difficult in general. This class
    takes a cctbx Structure object and from that calculates a
    theoretical intensity signal.
    
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

        This will create the refinement parameters using the CCTBXStructure
        adapter from diffpy.srfit.structure. Thus, the calculator will have its
        own parameters, each of which will be a proxy for some part of the
        structure. The parameters will be accessible by name under the
        'structure' attribute of this calculator.
        
        """
        # Load the structure from file
        from diffpy.Structure import Structure
        ds = Structure()
        ds.read(strufile)

        # Convert this to a cctbx xray scatterer
        cs = d2cStructure(ds)

        # Turn this into a parameterset
        parset = CCTBXStructureParSet(cs, "structure")

        # Use this parameter set for the calculator
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
        # Due to the design of cctbx structure, we have to create a new
        # unit_cell whenever a lattice parameter changes. By calling the update
        # function we can avoid remaking the structure multiple times between
        # uses of the structure.
        self.structure.update()
        self.count += 1
        print "iofq called", self.count
        return iofq(self.structure.stru, q)

# End class CCTBXIntensityCalculator

def d2cStructure(ds):
    """Turn a diffpy.Structure object into a cctbx.xray.structure object.

    This is used to simplify the input of structure files into a cctbx
    structure.
    
    """
    symm = crystal.symmetry( 
            unit_cell = ds.lattice.abcABG(),
            space_group_symbol = 1
            )
    scatt = []

    for a in ds:
        sc = xray.scatterer(
                label = a.element,
                site = a.xyz,
                u = a.Uisoequiv,
                occupancy = a.occupancy)
        scatt.append(sc)

    cs = xray.structure(
            crystal_symmetry = symm,
            scatterers = flex.xray_scatterer(scatt)
            )

    return cs

def iofq(cs, q):
    """Calculate I(Q) (X-ray) using the Debye Equation.

    I(Q) = 2 sum(i,j) f_i(Q) f_j(Q) sinc(rij Q) exp(-0.5 ssij Q**2)
    (The exponential term is the Debye-Waller factor.)

    cs  --  A cctbx.xray.scatterer object (P1 symmetry assumed).
    q   --  The q-points to calculate over.

    The calculator uses cctbx for the calculation of the f_i if it is
    available, otherwise f_i = 1.

    """
    # This will give us a distance generator
    mappings = cs.asu_mappings(buffer_thickness = 0)
    pairgen = crystal.neighbors_fast_pair_generator(
            mappings, distance_cutoff = 100)

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
    sites = set()
    elcount = {}
    scatterers = cs.scatterers()
    n = scatterers.size()
    for pair in pairgen:
        i = pair.i_seq
        j = pair.j_seq
        d = pair.dist_sq**0.5

        # count the number of each element
        si = scatterers[i]
        sj = scatterers[j]
        eli = si.element_symbol()
        elj = sj.element_symbol()

        if i not in sites:
            m = elcount.get(eli, 0)
            elcount[eli] = m + 1
            sites.add(i)

        # Get the pair
        els = [eli, elj]
        els.sort()

        # Get the distance to the desired precision
        D = int(d*dmult)

        # Get the DW factor to the same precision
        ss = si.u_iso + sj.u_iso
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

    y *= 2

    # Now we must add in the i == j pairs.
    for el, f in fdict.items():
        y += f**2 * elcount[el]

    # And that's it!

    return y

####### Example Code

def makeModel(strufile, datname):
    """Create a model that uses the CCTBXIntensityCalculator.

    This will create a Contribution that uses the CCTBXIntensityCalculator,
    associate this with a Profile, and use this to define a FitModel.

    """

    ## The Profile
    # Create a Profile. This will hold the experimental and calculated signal.
    profile = Profile()

    # Load data and add it to the profile
    x, y, u = numpy.loadtxt(datname, unpack=True)
    profile.setObservedProfile(x, y, u)

    ## The Calculator
    # Create a CCTBXIntensityCalculator named "I". This will be the name we use
    # to refer to the calculator from within the Contribution equation.  We
    # also need to load the model structure we're using.
    calculator = CCTBXIntensityCalculator("I")
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
    # CCTBXIntensityCalculator. By combining the different parts of the
    # calculation within the contribution equation, the time-consuming iofq
    # calculation is only performed when a structural parameter is changed. If
    # only non-structural parameters are changed, such as the background and
    # broadening parameters, then then previously computed iofq value will be
    # used to compute the contribution equation.  The benefit in this is very
    # apparent when refining the model with the LM optimizer, which only
    # changes two variables at a time in most cases. Note in the refinement
    # output how many times the residual is calculated, versus how many times
    # iofq is called when using the scipyOptimize function.

    # We will define the background as a string.
    bkgdstr = "b0 + b1*q + b2*q**2 + b3*q**3 + b4*q**4 + b5*q*5 + b6*q**6+\
               b7*q**7 +b8*q**8 + b9*q**9"

    contribution.registerStringFunction(bkgdstr, "bkgd")

    # We will create the broadening function by registering a python function.
    pi = numpy.pi
    exp = numpy.exp
    def gaussian(q, q0, width):
        return 1/(2*pi*width**2)**0.5 * exp(-0.5 * ((q-q0)/width)**2)

    contribution.registerFunction(gaussian)
    # Center the gaussian
    contribution.q0.setValue(x[len(x)/2])

    # Now we can incorporate the scale and bkgd into our calculation. We also
    # convolve the signal with the gaussian to broaden it.
    contribution.setEquation("scale * convolve(I, gaussian) + bkgd")

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

    # We also want to adjust the scale and the convolution width
    model.addVar(contribution.scale, 1)
    model.addVar(contribution.width, 0.1)

    # We can also refine structural parameters. 
    structure = calculator.structure
    a = structure.unitcell.a
    model.addVar(a)
    # We want to allow for isotropic expansion, so we'll make constraints for
    # that.
    model.constrain(structure.unitcell.b, a)
    model.constrain(structure.unitcell.c, a)
    # We want to refine the thermal paramters as well. We will add a new
    # variable that we call "Uiso" and constrain the atomic Uiso values to
    # this.
    uiso = model.newVar("uiso", 0.01)
    for s in structure.scatterers:
        model.constrain(s.uiso, uiso)

    # Give the model away so it can be used!
    return model

def plotResults(model):
    """Plot the results contained within a refined FitModel."""

    names = model.getNames()
    vals = model.getValues()

    q = model.bucky.profile.x

    # Plot this for fun.
    I = model.bucky.profile.y
    Icalc = model.bucky.profile.ycalc
    bkgd = model.bucky.evaluateEquation("bkgd()")
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
    makeData(strufile, q, "C60.iq", 8, 100.68, 0.005, 0.13, 2)

    model = makeModel(strufile, "C60.iq")
    scipyOptimize(model)
    rescount = model.fithook.count
    calcount = model.bucky.calculator.count
    footer = "iofq called %i%% of the time"%int(100.0*calcount/rescount)
    res = FitResults(model)
    res.printResults(footer = footer)
    plotResults(model)

    #model = makeModel(strufile, "C60.iq")
    #parkOptimize(model)
    #res = FitResults(model)
    #res.printResults()
    #plotResults(model)

# End of file
