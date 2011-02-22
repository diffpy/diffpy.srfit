#!/usr/bin/env python

from diffpy.srreal.pdfcalculator import PDFCalculator
from diffpy.Structure import Structure
from diffpy.srfit.fit.recipe import Recipe
import diffpy.srfit.pdf
from diffpy.srfit import loadProfile, residual, FitResults

def main():
    data = "data/ni-q27r100-neutron.gr"
    profile = loadProfile(data, "PDF")
    pdfc = PDFCalculator()
    pdfc.qmax = profile.meta["qmax"]
    pdfc.setScatteringFactorTableByType(profile.meta["stype"])

    stru = Structure(filename="./data/ni.cif")
    # XXX - made change here. Adapted objects must be registered with the
    # recipe. The reason for this is because these objects are held as adapters
    # inside the recipe. Changing this would require an overhaul of this
    # branch, which is premature at this point.
    recipe = Recipe(pdfc=pdfc,stru=stru,profile=profile)
    recipe.vars(a=stru.lattice.a, Uiso=0.003, scale=1, qdamp=0.01, delta2=5)

    # Link the fit range between the profile and the calculator
    recipe.pars(rmin = 0.1, rmax = 20, rstep=0.05)
    recipe.constrain(pdfc, recipe.rmin)('pdfc.rmin = rmin')
    recipe.constrain(pdfc, recipe.rmax)('pdfc.rmax = rmax')
    recipe.constrain(pdfc, recipe.rstep)('pdfc.rstep = rstep')
    # XXX pdfc is not a parameter, but it works very well here
    recipe.constrain(profile, pdfc)('profile.setPoints(pdfc.rgrid)')

    # Set of other constraints
    recipe.constrain(stru, recipe.a)('stru.lattice.setLatPar(a=a, b=a, c=a)')
    recipe.constrain(stru, recipe.Uiso)('stru.Uisoequiv = Uiso')
    recipe.constrain(pdfc, recipe.scale)('pdfc.scale = scale')
    recipe.constrain(pdfc, recipe.qdamp)('pdfc.qdamp = qdamp')
    recipe.constrain(pdfc, recipe.delta2)('pdfc.delta2 = delta2')

    # Make the equation to be fit
    gdiff = recipe.build('(profile.y - pdfc(stru)[1]) / profile.dy')
    res = residual(gdiff)

    from scipy.optimize import leastsq
    recipe.rstep.value=0.1
    leastsq(res.vec, res.values)
    recipe.rstep.value=0.01
    leastsq(res.vec, res.values)

    results = FitResults(res)
    results.show()

    # and plot 
    from pylab import plot, show
    r, gr, dgr = profile
    r, rcalc = pdfc(stru)
    plot(r, gr, 'bo')
    plot(r, rcalc, color='r', linestyle='-')
    diff = gr - rcalc
    plot(r, diff + 1.15*min(gr), color='g', linestyle='-')
    show()

if __name__ == "__main__":
    main()
