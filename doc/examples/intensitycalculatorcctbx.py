#!/usr/bin/env python

from cctbx import crystal
from cctbx import xray
from cctbx import miller
from cctbx.array_family import flex

from diffpy.srfit.fitbase import Calculator, Contribution, FitModel, Profile
from diffpy.srfit.fitbase import FitResults
from diffpy.srfit.park import FitnessAdapter
from diffpy.srfit.structure import CCTBXStructureParSet

from debyemodel import scipyOptimize, parkOptimize

import numpy

class CCTBXIntensityCalculator(Calculator):
    """A class for calculating intensity using the Debye equation.

    Calculating intensity from a structure is difficult in general. This class
    takes a cctbx Structure object and from that calculates a theoretical
    intensity signal.
    
    """

    def __init__(self, name):
        """Define our calculator.

        Keep count of how many times the function has been called (self.count).
        """
        Calculator.__init__(self, name)
        # Count the calls
        self.count = 0
        self.wl = 0

        # Another parameter for zeroing the signal
        self.newParameter("zero", 0)

        # We need parameters for the Gaussian width of the peak (following
        # Rietveld, as described in the GSAS manual).
        self.newParameter("U", 0)
        self.newParameter("V", 0)
        self.newParameter("W", 0)
        self.newParameter("As", 0)
        self.newParameter("F1", 0)
        self.newParameter("F2", 0)

        return

    def setWavelength(self, wl):
        """Set the wavelength of the scattering probe.

        This is needed to calculate the scattering angle for various
        corrections. The wavelength is in Angstroms.

        """
        self.wl = wl
        return

    def setStructure(self, cs):
        """Set the structure used in the calculation.

        This will create the refinement parameters using the CCTBXStructure
        adapter from diffpy.srfit.structure. Thus, the calculator will have its
        own parameters, each of which will be a proxy for some part of the
        structure. The parameters will be accessible by name under the
        'structure' attribute of this calculator.
        
        """
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

        # Get the current parameter values and pass them along to the
        # calculation function. By virtue of the CCTBXStructureParSet, the
        # structure is already up to date.
        z = self.zero.getValue()
        u = self.U.getValue()
        v = self.V.getValue()
        w = self.W.getValue()
        asym = self.As.getValue()
        f1 = self.F1.getValue()
        f2 = self.F2.getValue()

        print "iofq called", self.count
        return iofq(self.structure.stru, q-z, self.wl, u, v, w, asym, f1, f2)

# End class CCTBXIntensityCalculator

# Calculation code

def iofq(cs, q, wl, u, v, w, asym, f1, f2):
    """Calculate the diffraction intensity of a structure over a q-range."""
    dmin = 1.9*numpy.pi/q[-1]

    millerset = miller.build_set( 
            crystal_symmetry=cs, 
            d_min=dmin, 
            anomalous_flag=False)

    sffs = xray.structure_factors.from_scatterers(
            miller_set = millerset,
            wing_cutoff = 1e-6,
            exp_table_one_over_step_size=0)

    fcalc = sffs(
            xray_structure = cs,
            miller_set = millerset,
            algorithm="direct").f_calc().data()
    fcalc = flex.pow(flex.abs(fcalc), 2)

    m = millerset.multiplicities().data()

    # Convert the miller set into q
    refq = 4*numpy.pi*millerset.unit_cell().stol(millerset.indices())
    # FIXME - figure out how to sort properly!
    qfm = zip(refq, fcalc, m)
    qfm.sort()
    refq = numpy.array([a[0] for a in qfm])
    f = numpy.array([a[1] for a in qfm])
    m = numpy.array([a[2] for a in qfm])

    y = shapePeaks(q, refq, f, m, wl, u, v, w, asym, f1, f2)

    return y

def shapePeaks(q, refq, f, m, wl, u, v, w, asym, f1, f2):
    """Broaden the peaks using a gaussian of width.

    This applies the Rietveld peak broadening as a function of theta. See the
    GSAS manual for a description of these terms.

    returns a new intensity array
    
    """
    exp = numpy.exp
    pi = numpy.pi
    tan = numpy.tan
    abs = numpy.abs
    arcsin = numpy.arcsin

    y = 0

    # Calculate theta (radians)
    th = arcsin(q*wl/(4*pi))

    for q0, amp, mult in zip(refq, f, m):

        # scattering angle in radians
        th0 = arcsin(q0*wl/(4*pi))
        tt0 = tan(th0)

        # The gaussian width and normalization
        ss = u*tt0**2 + v*tt0 + w + f2/(tt0**4)
        N  = (2*pi*ss)**-0.5

        # The displacement of theta
        dt = th - th0
        dtp = dt + f1 / tan(2*th0)

        # the peak shape function
        y0 = 1 - dtp * abs(dt) * asym / tt0
        y0 *= mult * N
        y0 *= exp(-0.5 * dtp**2 / ss)

        y += amp * y0

    return y

def makeStructure():
    """Make the nickel structure."""
    symm = crystal.symmetry( 
            unit_cell = (3.52, 3.52, 3.52, 90, 90, 90),
            space_group_symbol = 225
            )

    ni = xray.scatterer(
            label = "Ni",
            site = (0, 0, 0),
            u = 0.002,
            occupancy = 1)
    cs = xray.structure(
            crystal_symmetry = symm,
            scatterers = flex.xray_scatterer([ni])
            )

    return cs

####### Example Code

def makeModel(cs, datname):
    """Create a model that uses the CCTBXIntensityCalculator.

    This will create a Contribution that uses the CCTBXIntensityCalculator,
    associate this with a Profile, and use this to define a FitModel.

    """

    ## The Profile
    # Create a Profile. This will hold the experimental and calculated signal.
    profile = Profile()

    # Load data and add it to the profile
    x, y = numpy.loadtxt(datname, unpack=True)
    profile.setObservedProfile(x, y)

    # Set the range we need
    profile.setCalculationRange(2.5, 8)

    ## The Calculator
    # Create a CCTBXIntensityCalculator named "I". This will be the name we use
    # to refer to the calculator from within the Contribution equation.  We
    # also need to load the model structure we're using.
    calculator = CCTBXIntensityCalculator("I")
    calculator.setStructure(cs)
    calculator.setWavelength(0.2128)
    
    ## The Contribution
    # Create a Contribution, that will associate the Profile with the
    # Calculator.  The calculator will be accessible as an attribute of the
    # Contribution by its name ("I"), or simply by "calculator".  We also want
    # to tell the contribution to name the x-variable of the profile "q", so we
    # can use it in equations with this name.
    contribution = Contribution("nickel")
    contribution.setCalculator(calculator)
    contribution.setProfile(profile, xname = "q")

    # Now we're ready to define the contribution equation. We need to modify
    # the Calcultor, and we'll do that from within the Contribution eqation for
    # the sake of instruction. We want to modify the calculator in two ways.
    # We need a scale factor and a polynomial background.

    # We will define the background as a string.
    bkgdstr = "b0 + b1*q + b2*q**2 + b3*q**3 + b4*q**4 + b5*q*5"

    contribution.registerStringFunction(bkgdstr, "bkgd")

    # Now we can incorporate the scale and bkgd into our calculation.
    contribution.setEquation("scale * I + bkgd")

    # Make a FitModel where we can create variables, constraints and
    # restraints. If we had multiple profiles to fit simultaneously, the
    # contribution from each could be added to the model.
    model = FitModel()
    model.addContribution(contribution)

    # Specify which parameters we want to refine. We can give them initial
    # values in the process. We want to refine the background variables that we
    # just defined in the contribution. We only need about 5 of them.
    model.addVar(contribution.b0, 0)
    model.addVar(contribution.b1, 0)
    model.addVar(contribution.b2, 0)
    model.addVar(contribution.b3, 0)
    model.addVar(contribution.b4, 0)
    model.addVar(contribution.b5, 0)

    # We also want to adjust the scale
    model.addVar(contribution.scale, 1)

    # We can also refine structural parameters. 
    structure = calculator.structure
    a = structure.unitcell.a
    model.addVar(a)
    # We want to allow for isotropic expansion, so we'll make constraints for
    # that. This will be done automatically in the future.
    model.constrain(structure.unitcell.b, a)
    model.constrain(structure.unitcell.c, a)
    # We want to refine the thermal paramters as well. We will add a new
    # variable that we call "uiso" and constrain the atomic uiso values to
    # this.
    uiso = model.newVar("uiso", 0.00127)
    for s in structure.scatterers:
        model.constrain(s.uiso, uiso)

    # Don't forget about the width of the peaks
    model.addVar(contribution.calculator.V, 1e-6)
    model.addVar(contribution.calculator.W, 1e-6)

    # And the zero offset
    model.addVar(contribution.calculator.zero, 0)

    # Give the model away so it can be used!
    return model

def plotResults(model):
    """Plot the results contained within a refined FitModel."""

    names = model.getNames()
    vals = model.getValues()

    q = model.nickel.profile.x

    # Plot this for fun.
    I = model.nickel.profile.y
    Icalc = model.nickel.profile.ycalc
    bkgd = model.nickel.evaluateEquation("bkgd()")
    diff = I - Icalc - 0.1*max(I)

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

def fitModel(model):
    """Fit the model according to a recommended heuristic.

    It is unlikely that the Rietveld refinement will converge given the number
    of free parameters without some guidance. This function provides that
    guidance.
    
    """

    # Clean slate, nothing is refined.
    model.fixAll()

    # Start with the scale factor
    model.freeVar(model.scale)
    scipyOptimize(model)


    # Next the zero-shift of the q vector, a linear background and the lattice
    # parameter.
    model.freeVar(model.zero)
    model.freeVar(model.b0)
    model.freeVar(model.b1)
    model.freeVar(model.a)
    scipyOptimize(model)

    # now the isotropic thermal parameters
    model.freeVar(model.uiso)
    scipyOptimize(model)

    # Next the peak width terms
    model.freeVar(model.V)
    model.freeVar(model.W)
    scipyOptimize(model)

    # Get the rest of the background terms now
    model.freeAll()
    scipyOptimize(model)

    return

if __name__ == "__main__":

    # Make the data and the model
    cs = makeStructure()
    model = makeModel(cs, "data/ni.iq")
    fitModel(model)
    rescount = model.fithook.count
    calcount = model.nickel.calculator.count
    footer = "iofq called %i%% of the time"%int(100.0*calcount/rescount)
    res = FitResults(model)
    res.printResults(footer = footer)
    plotResults(model)

