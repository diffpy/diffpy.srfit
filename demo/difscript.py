################################################################################
#
#   SrRietveld Refinement of RaPDF Ni data.
#
################################################################################
__id__ = "$Id$"

from SrRietveld     import *

def initialize(fit):
    # Set the number of refinement cycles
    fit.setNumCycles(-1)

    # set component
    comp = RietveldComponent()
    comp.setName("Ni")
    fit.addComponent(comp)

    ### Instantiate the pattern
    pat = CWXPattern()
    comp.setPattern(pat)
    pat.loadData("nickel_5f_300K_1500msec_2theta.dat", XYSigmaParser())
    pat.setPrimaryWavelength(0.125677)

    lp = StandardDebyeScherrerLPFactor()
    lp.setMonoPolarizationCorrection(0.91)
    pat.setLorentzPolarizationFactor(lp)

    # setup the excluded region
    comp.setFitRange(xmin = 2.50)

    # Set the zero parameters
    pat.setValue("zshift", 0)
    fit.mapVP("v_zshift", pat, "zshift")

    # Instantiate the background
    bkgd = InterpolatedBackground("Cubic")
    comp.setBackground(bkgd)
    # setup the background
    bkgdpointlist = [
     ( 2.506 , 1060.166 ), 
     ( 3.010 , 872.325 ),
     ( 5.176 , 801.970 ),
     ( 6.310 , 799.834 ),
     ( 7.871 , 663.524 ),
     ( 8.451 , 672.373 ),
     ( 9.685 , 627.120 ),
     ( 10.365 , 605.252 ),
     ( 11.247 , 532.866 ),
     ( 11.776 , 584.413 ),
     ( 12.733 , 545.038 ),
     ( 13.186 , 530.719 ),
     ( 13.942 , 498.045 ),
     ( 15.050 , 525.185 ),
     ( 16.234 , 482.769 ),
     ( 16.562 , 475.984 ),
     ( 17.166 , 483.853 ),
     ( 18.148 , 482.349 ),
     ( 19.080 , 470.769 ),
     ( 19.912 , 442.536 ),
     ( 20.743 , 466.518 ),
     ( 21.675 , 443.857 ),
     ( 23.060 , 422.947 ),
     ( 23.363 , 431.474 ),
     ( 23.841 , 432.576 ),
     ( 24.572 , 418.514 ),
     ( 25.227 , 407.441 ),
     ( 25.529 , 412.225 ),
     ( 25.932 , 399.844 ),
     ( 26.587 , 403.386 ),
     ( 27.468 , 396.797 ),
     ( 27.897 , 392.290 ),
     ( 28.451 , 381.667 ),
     ( 29.080 , 384.864 ),
     ( 29.685 , 378.912 ),
     ( 30.264 , 379.603 ),
     ( 30.793 , 365.727 )
    ]
    bkgd.addBackgroundPoints(bkgdpointlist)

    ### Instantiate the phase 1
    pha1 = CrystalPhase()
    pha1.setName("Ni")
    comp.addPhase(pha1)
    pha1.setATZ(540979.1880)
    # Set the space group
    pha1.setSpaceGroup("Fm-3m")

    # Create a lattice for the phase
    lat1 = Lattice()
    pha1.setLattice(lat1)
    lat1.setValue("a", 3.524)
    fit.mapVP("v_a", lat1, "a")
    # Add an atom
    a1 = configureAtom(pha1, "Ni",  "Ni", 0.0,  0.00000,  0.0,  0.38787,   1.0)
    fit.mapVP("v_biso", a1.getADP(), "Biso")

    # Instantiate the peak-shape profile function.
    prof1 = PseudoVoigtProfile()
    pha1.setPeakProfile(prof1)
    prof1.setValue("scale" , 0.34741e-7) 
    prof1.setValue("eta0", 0.8)
    prof1.setValue("V", -0.007618)
    prof1.setValue("W", 0.006255)
    prof1.setValue("X", 0.018961)
    # do some constraints
    fit.mapVP("v_pscale",  prof1, "scale")
    fit.mapVP("v_eta0",  prof1, "eta0")
    fit.mapVP("v_V",  prof1, "V")
    fit.mapVP("v_W",  prof1, "W")
    fit.mapVP("v_X",  prof1, "X")


def runfit(fit):

    comp = fit.getComponent(0)
    # Fitting
    stage1(fit)
    stage2(fit)
    stage3(fit)
    stage3b(fit)
    stage4(fit)
    stage5(fit)
    stage6(fit)
    stage7(fit)
    stage8(fit)
    stage9(fit)
    stage10(fit)
    stage11(fit)
    fit.printResults()
    fit.saveResults("ni.dif.res")
    comp.saveFitArrays("ni.dif.fit")

def stage1(fit):
    print "stage 1"
    fit.fixAll()
    fit.guessV("v_pscale")
    fit.refine()
    fit.printResults()
    return

def stage2(fit):
    print "stage 2"
    fit.fixAll()
    fit.guessV("v_pscale")
    fit.guessV("v_zshift")
    fit.refine()
    fit.printResults()
    return

def stage3(fit):
    print "stage 3"
    fit.fixAll()
    fit.guessV("v_pscale")
    fit.guessV("v_zshift")
    fit.guessV("v_biso")
    fit.refine()
    fit.printResults()
    return

def stage3b(fit):
    print "stage 3b"
    fit.fixAll()
    fit.guessV("v_eta0")
    fit.refine()
    fit.printResults()
    return

def stage4(fit):
    print "stage 4"
    fit.fixAll()
    fit.guessV("v_eta0")
    fit.guessV("v_pscale")
    fit.guessV("v_zshift")
    fit.guessV("v_biso")
    fit.refine()
    fit.printResults()
    return

def stage5(fit):
    print "stage 5"
    fit.fixAll()
    fit.guessV("v_pscale")
    fit.guessV("v_eta0")
    fit.guessV("v_zshift")
    fit.guessV("v_biso")
    fit.guessV("v_V")
    fit.refine()
    fit.printResults()
    return

def stage6(fit):
    print "stage 6"
    fit.fixAll()
    fit.guessV("v_eta0")
    fit.guessV("v_pscale")
    fit.guessV("v_zshift")
    fit.guessV("v_biso")
    fit.guessV("v_W")
    fit.refine()
    fit.printResults()
    return

def stage7(fit):
    print "stage 7"
    fit.fixAll()
    fit.guessV("v_eta0")
    fit.guessV("v_pscale")
    fit.guessV("v_zshift")
    fit.guessV("v_biso")
    fit.guessV("v_X")
    fit.refine()
    fit.printResults()
    return

def stage8(fit):
    print "stage 8"
    fit.fixAll()
    fit.guessV("v_eta0")
    fit.guessV("v_pscale")
    fit.guessV("v_zshift")
    fit.guessV("v_biso")
    fit.guessV("v_W")
    fit.guessV("v_X")
    fit.refine()
    fit.printResults()
    return

def stage9(fit):
    print "stage 9"
    fit.fixAll()
    fit.guessV("v_eta0")
    fit.guessV("v_pscale")
    fit.guessV("v_zshift")
    fit.guessV("v_biso")
    #fit.guessV("v_U")
    fit.guessV("v_V")
    #fit.guessV("v_W")
    fit.guessV("v_X")
    fit.refine()
    fit.printResults()
    return

def stage10(fit):
    print "stage 10"
    fit.fixAll()
    fit.guessV("v_eta0")
    fit.guessV("v_pscale")
    fit.guessV("v_zshift")
    fit.guessV("v_biso")
    #fit.guessV("v_U")
    fit.guessV("v_V")
    fit.guessV("v_W")
    fit.guessV("v_X")
    fit.refine()
    fit.printResults()
    return

def stage11(fit):
    print "stage 11"
    fit.guessV("v_zshift")
    fit.guessV("v_pscale")
    fit.guessV("v_V")
    fit.guessV("v_W")
    fit.guessV("v_X")
    fit.guessV("v_eta0")
    fit.guessV("v_biso")
    fit.guessV("v_a")
    fit.refine()
    fit.printResults()
    return
    
if __name__ == "__main__":
    fit   = Fit()
    initialize(fit)
    runfit(fit)
