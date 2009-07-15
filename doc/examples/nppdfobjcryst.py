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
"""Example of using more complex Calculators.

This is an example of using a more complex Calculator in a FitContribution.

The PDFCalculator class is an example of a Calculator that can be used by a
FitContribution to help generate a signal. It uses the ObjCrystParSet to hold
crystal and molecular information from a pyobjcryst Crystal.

The makeRecipe function shows how to build a FitRecipe that uses the
PDFCalculator.

"""

import numpy

from diffpy.srfit.fitbase import Calculator, FitContribution, FitRecipe, Profile
from diffpy.srfit.fitbase import FitResults
from diffpy.srfit.structure.objcryststructure import ObjCrystParSet

from gaussianrecipe import scipyOptimize, parkOptimize

class PDFCalculator(Calculator):
    """A class for calculating the PDF for an isolated scatterer.

    This is another example of using a Calculator class with a ParameterSet
    adapter to refine a structure recipe to data using a FitRecipe.
    
    """

    def __init__(self, name):
        """Initialize our calculator.

        Here we add a calculator-level Parameter called 'delta2' that
        describes the vibrational correlation between atoms. We also add some
        metadata required by the calculator.
        
        """
        Calculator.__init__(self, name)

        # Add any non-structural parameters here
        self.newParameter("delta2", 0)

        # Non-Parameters that are needed by the calculator.
        self.meta["qmin"] = 1
        self.meta["qmax"] = 20
        return

    def setCrystal(self, cryst):
        """Set the pyobjcryst.Crystal instance for the calculator.

        This converts the Crystal to a ObjCrystParSet that is used to organize
        the fitting parameters for the calculator.  The calculator will have
        its own parameters, each of which will be a proxy for some part of the
        crystal. The parameters will be accessible by name under the 'crystal'
        attribute of this calculator.
        
        """

        # Create a custom ParameterSet designed to interface with
        # pyobjcryst.Crystal
        parset = ObjCrystParSet(cryst, "crystal")
        # Put this ParameterSet in the Calculator.
        self.addParameterSet(parset)
        return

    def __call__(self, r):
        """Calculate the PDF.

        This Calculator will be used in a fit equation that will be optimized
        to fit some data.  By the time this function is evaluated, the crystal
        has been updated by the optimizer via the ObjCrystParSet created in
        setCrystal. Thus, we need only call pdf with the internal structure
        object.

        """
        qmin = self.meta["qmin"]
        qmax = self.meta["qmax"]
        return pdf(self.crystal.cryst, r, self.delta2.getValue(), qmin, qmax)

# End class PDFCalculator

def pdf(cryst, r, delta2 = 0, qmin = 1, qmax = 20):
    """Calculate the PDF from a diffpy.Structure

    cryst   --  A pyobjcryst.Crystal instance. It is assumed that the structure
                is that of an isolated scatterer. Periodic boundary conditions
                are not applied.
    r       --      The r-points to calculate over.
    delta2  --  The correlation term in the Debye-Waller factor (default 0).
    qmin    --  The minimum observed q-value (default 1).
    qmax    --  The maximum observed q-value (default 20).
    
    """
    pi = numpy.pi

    # First we must calculate F(Q). To do this, we need to create the q-points
    # corresponding to the r-points. It is assumed that the r points are
    # equidistant.
    rmax = r[-1]
    dr = r[1] - r[0]

    qmaxcalc = 2*pi/dr
    # Note that this should be 2*pi/rmax, but the fft calculates a signal for
    # both positive and negative values of frequency, and we're going to throw
    # away half of the transformed values.
    dq = pi/rmax

    # We must create the q-points that correspond to the r-points
    q = numpy.arange(0, qmaxcalc, dq)

    # Calculate F(Q)
    fq = fofq(cryst, q, delta2)

    # and enforce the qmin and qmax bounds
    fq[numpy.logical_or(q < qmin, q > qmax)] = 0


    # Now we transform this to gr
    from numpy.fft import ifft, ifftshift

    # sin FT = imaginary part of ifft
    grcalc = ifft(fq).imag

    # Properly normalize this. 
    grcalc *= (2 / pi) * qmaxcalc

    # Now we need to pick out the positive frequency terms. 
    rmaxidx = len(grcalc)/2 + 1
    # We also need to discard values for r < rmin.
    rminidx = int(numpy.ceil(r[0]/dr))

    gr = grcalc[rminidx:rmaxidx]

    assert(len(r) == len(gr))

    return gr

def fofq(cryst, q, delta2):
    """Calculate F(Q) (X-ray) using the Debye Equation.

    F(Q) = 2/(N<f>^2) 
           sum(i, j > i) f_i(Q) f_j(Q) sin(rij Q)/rij exp(-0.5 ssij Q**2)
    (The exponential term is the Debye-Waller factor.)

    cryst   --  A pyobjcryst.Crystal instance. It is assumed that the structure
                is that of an isolated scatterer. Periodic boundary conditions
                are not applied.
    q       --  The q-points to calculate over.
    delta2  --  The correlation term in the Debye-Waller factor.

    The calculator uses cctbx for the calculation of the f_i if it is
    available, otherwise f_i = 1.

    """
    # The functions we need
    sin = numpy.sin
    exp = numpy.exp
    pi = numpy.pi

    # The brute-force calculation is very slow. Thus we optimize a little bit.

    # The precision of distance measurements.
    deltad = 1e-6
    dmult = int(1/deltad)
    deltau = deltad**2
    umult = int(1/deltau)

    pairdict = {}
    elcount = {}

    # Loop over scattering components in the crystal.

    scl = cryst.GetScatteringComponentList()
    n = len(scl)


    for i in xrange(n):

        # Get the ScatteringPower
        scli = scl[i]
        spi = scli.mpScattPow

        # Check for dummy element
        if spi is None: continue

        # Count the element
        eli = spi.GetSymbol().title()
        m = elcount.get(eli, 0)
        elcount[eli] = m + 1

        xif = scli.X
        yif = scli.Y
        zif = scli.Z
        xi, yi, zi = cryst.FractionalToOrthonormalCoords(xif, yif, zif)

        for j in xrange(i+1,n):

            sclj = scl[j]
            # Get the ScatteringPower
            spj = sclj.mpScattPow

            # Check for dummy element
            if spj is None: continue

            elj = spj.GetSymbol().title()

            # Get the pair
            els = [eli, elj]
            els.sort()

            # Get the coordinates
            xjf = sclj.X
            yjf = sclj.Y
            zjf = sclj.Z
            xj, yj, zj = cryst.FractionalToOrthonormalCoords(xjf, yjf, zjf)

            # Get the distance to the desired precision
            d = ((xi - xj)**2 + (yi - yj)**2 + (zi - zj)**2)**0.5
            D = int(d*dmult)

            # Get the DW factor to the same precision (calculated from Biso,
            # not Uiso)
            ss = spi.Biso + spj.Biso
            corfact = 1 - delta2/d**2
            if corfact > 0:
                ss *= corfact
            SS = int(ss*umult)

            # Record the multiplicity of this pair
            key = (els[0], els[1], D, SS)
            mult = pairdict.get(key, 0)
            pairdict[key] = mult + 1


    # Now we can calculate G(r) from the pair dictionary. Making the dictionary
    # first reduces the amount of calls to sin and exp we have to make.

    # First we must cache the scattering factors
    fdict = {}
    nfavg = 0
    for el, m in elcount.iteritems():
        fdict[el] = f = getXScatteringFactor(el, q)
        nfavg += m * f

    # Now we can compute F(Q) for the i != j pairs
    y = 0
    BtoU = 1.0 / (8 * pi**2)
    for key, mult in pairdict.items():
        eli = key[0]
        elj = key[1]
        fi = fdict[eli]
        fj = fdict[elj]
        D = key[2]
        SS = key[3]

        d = D * deltad

        # Add in the fitcontribution
        y += fi * fj * mult * sin(q * d) / d *\
                exp(-0.5 * (SS * deltau) * BtoU * q**2)


    # We must multiply by 2 since we only counted j > i pairs.
    y *= 2

    # Now divide by N <f>^2

    y *= n / nfavg**2

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
        return 1.0


c60xyz = \
"""
3.451266498   0.685000000   0.000000000
3.451266498  -0.685000000   0.000000000
-3.451266498   0.685000000   0.000000000
-3.451266498  -0.685000000   0.000000000
0.685000000   0.000000000   3.451266498
-0.685000000   0.000000000   3.451266498
0.685000000   0.000000000  -3.451266498
-0.685000000   0.000000000  -3.451266498
0.000000000   3.451266498   0.685000000
0.000000000   3.451266498  -0.685000000
0.000000000  -3.451266498   0.685000000
0.000000000  -3.451266498  -0.685000000
3.003809890   1.409000000   1.171456608
3.003809890   1.409000000  -1.171456608
3.003809890  -1.409000000   1.171456608
3.003809890  -1.409000000  -1.171456608
-3.003809890   1.409000000   1.171456608
-3.003809890   1.409000000  -1.171456608
-3.003809890  -1.409000000   1.171456608
-3.003809890  -1.409000000  -1.171456608
1.409000000   1.171456608   3.003809890
1.409000000  -1.171456608   3.003809890
-1.409000000   1.171456608   3.003809890
-1.409000000  -1.171456608   3.003809890
1.409000000   1.171456608  -3.003809890
1.409000000  -1.171456608  -3.003809890
-1.409000000   1.171456608  -3.003809890
-1.409000000  -1.171456608  -3.003809890
1.171456608   3.003809890   1.409000000
-1.171456608   3.003809890   1.409000000
1.171456608   3.003809890  -1.409000000
-1.171456608   3.003809890  -1.409000000
1.171456608  -3.003809890   1.409000000
-1.171456608  -3.003809890   1.409000000
1.171456608  -3.003809890  -1.409000000
-1.171456608  -3.003809890  -1.409000000
2.580456608   0.724000000   2.279809890
2.580456608   0.724000000  -2.279809890
2.580456608  -0.724000000   2.279809890
2.580456608  -0.724000000  -2.279809890
-2.580456608   0.724000000   2.279809890
-2.580456608   0.724000000  -2.279809890
-2.580456608  -0.724000000   2.279809890
-2.580456608  -0.724000000  -2.279809890
0.724000000   2.279809890   2.580456608
0.724000000  -2.279809890   2.580456608
-0.724000000   2.279809890   2.580456608
-0.724000000  -2.279809890   2.580456608
0.724000000   2.279809890  -2.580456608
0.724000000  -2.279809890  -2.580456608
-0.724000000   2.279809890  -2.580456608
-0.724000000  -2.279809890  -2.580456608
2.279809890   2.580456608   0.724000000
-2.279809890   2.580456608   0.724000000
2.279809890   2.580456608  -0.724000000
-2.279809890   2.580456608  -0.724000000
2.279809890  -2.580456608   0.724000000
-2.279809890  -2.580456608   0.724000000
2.279809890  -2.580456608  -0.724000000
-2.279809890  -2.580456608  -0.724000000
"""

def makeC60():
    """Make the C60 molecule using pyobjcryst."""

    from pyobjcryst.crystal import Crystal
    from pyobjcryst.molecule import Molecule
    from pyobjcryst.scatteringpower import ScatteringPowerAtom

    pi = numpy.pi
    c = Crystal(100, 100, 100, "P1")
    c.SetName("c60frame")
    m = Molecule(c, "c60")

    c.AddScatterer(m)

    sp = ScatteringPowerAtom("C", "C")
    sp.SetBiso(0.25)

    # Create a dummy atom at the center.
    m.AddAtom(0, 0, 0, None, "center")

    # Add the other atoms. They will be named C0, C1, ..., C59.
    for i, l in enumerate(c60xyz.strip().splitlines()):
        x, y, z = map(float, l.split())
        m.AddAtom(x, y, z, sp, "C%i"%i)

    return c

####### Example Code

def makeRecipe(cryst, datname):
    """Create a recipe that uses the PDFCalculator.

    This will create a FitContribution that uses the PDFCalculator,
    associate this with a Profile, and use this to define a FitRecipe.

    """

    ## The Profile
    # Create a Profile. This will hold the experimental and calculated signal.
    profile = Profile()

    # Load data and add it to the profile
    x, y, dy = numpy.loadtxt(datname, unpack=True)
    profile.setObservedProfile(x, y, dy)
    profile.setCalculationRange(xmin=1.6, xmax=8)

    ## The Calculator
    # Create an PDFCalculator named "G". This will be the name we use to refer
    # to the calculator from within the FitContribution equation.  We also need
    # to load the pyobjcryst Crystal we're using.
    calculator = PDFCalculator("G")
    calculator.setCrystal(cryst)
    # These are metadata needed by the calculator
    calculator.meta["qmin"] = 0.68
    calculator.meta["qmax"] = 22
    
    ## The FitContribution
    # Create a FitContribution that will associate the Profile with the
    # Calculator.  The Calculator will be accessible as an attribute of the
    # FitContribution by its name ("G"), or simply by "calculator".  We also
    # want to tell the FitContribution to name the x-variable of the profile
    # "r" so we can use it in equations with this name.
    contribution = FitContribution("bucky")
    contribution.setCalculator(calculator)
    contribution.setProfile(profile, xname = "r")

    # Now we're ready to define the FitContribution equation. We need to modify
    # the Calcultor, and we'll do that from within the FitContribution equation
    # for the sake of instruction. Here we add a peak dampening factor and a
    # scale factor.
    contribution.setEquation("scale * exp(-0.5 * (qdamp * r)**2) * G")

    # Make a FitRecipe where we can create variables, constraints and
    # restraints. If we had multiple profiles to fit simultaneously, the
    # contribution from each could be added to the recipe.
    recipe = FitRecipe()
    recipe.addContribution(contribution)

    # Specify which parameters we want to refine. We'll be using the
    # MoleculeParSet within the calculator's ObjCrystParSet directly, so let's
    # get a handle to it. See the diffpy.srfit.structure.objcryststructure
    # module for more information about the ObjCrystParSet hierarchy.
    c60 = calculator.crystal.c60

    # First, the isotropic thermal displacement factor.
    biso = recipe.newVar("biso", 0.25)
    for atom in c60.atoms:
        # We have to check for reference atoms. Dummy atoms have no biso
        if not atom.isDummy():
            recipe.constrain(atom.biso, biso)

    # And the correlation term
    recipe.addVar(calculator.delta2, 2)

    # We need to let the molecule expand. If we were modeling it as a crystal,
    # we could let the unit cell expand. For instruction purposes, we use a
    # Molecule to model C60, and molecules have different modeling options than
    # crystals. To make the molecule expand from a central point, we will
    # constrain the distance from each atom to a dummy center atom that was
    # created with the molecule, and allow that distance to vary. (We could
    # also let the nearest-neighbor bond lengths vary, but that would be much
    # more difficult to set up.)
    center = c60.center
    radius = recipe.newVar("radius", 3.5)
    for i, atom in enumerate(c60.atoms[1:]):
        # This creates a Parameter that moves atoms according to the bond
        # length. Note that each Parameter needs a unique name.
        par = c60.addBondLengthParameter("rad%i"%i, center, atom)
        recipe.constrain(par, "abs(radius)")

    # We also want to adjust the scale and qdamp
    recipe.addVar(contribution.scale, 1.3e4)
    recipe.addVar(contribution.qdamp, 0.1)

    # Give the recipe away so it can be used!
    return recipe

def plotResults(recipe):
    """Plot the results contained within a refined FitRecipe."""

    names = recipe.getNames()
    vals = recipe.getValues()

    r = recipe.bucky.profile.x

    # Plot this.
    G = recipe.bucky.profile.y
    Gcalc = recipe.bucky.profile.ycalc
    diff = G - Gcalc - 10 * recipe.scale.getValue()

    import pylab
    pylab.plot(r,G,'ob',label="G(r) Data")
    pylab.plot(r,Gcalc,'-r',label="G(r) Fit")
    pylab.plot(r,diff,'-g',label="G(r) diff")
    pylab.xlabel("$r (\AA)$")
    pylab.ylabel("$G (\AA^{-2})$")
    pylab.legend(loc=1)

    pylab.show()
    return

if __name__ == "__main__":

    cryst = makeC60()
    # Make the data and the recipe
    recipe = makeRecipe(cryst, "data/C60.gr")
    # Tell the fithook that we want very verbose output.
    recipe.fithook.verbose = 3
    
    # Optimize
    scipyOptimize(recipe)
    #parkOptimize(recipe)

    # Print results
    res = FitResults(recipe)
    res.printResults()

    # Plot results
    plotResults(recipe)

# End of file
