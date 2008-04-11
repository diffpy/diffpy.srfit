################################################################################
#
#   PDF and Rietveld Corefinement of Ni data.
#
################################################################################
__id__ = "$Id$"

def runfit(fit):

    import difscript

    # Do the rietveld refinement first (it's more tempermental)
    difscript.runfit(fit)
    difcomp = fit.getComponent(0)

    ###############################################################################
    # PDF                                                                         #
    ###############################################################################
    from PDFAPI import PDFComponent, PDFData, PDFParser
    from PDFAPI import CrystalPhase as PDFPhase
    from PDFAPI import Lattice as PDFLattice
    from PDFAPI import Atom as PDFAtom
    from PDFAPI import IsotropicAtomicDisplacementFactor as PDFIADP

    # set component
    pdfcomp = PDFComponent()
    pdfcomp.setName("PDFNi")
    fit.addComponent(pdfcomp)

    ### Instantiate the pattern
    pdfpat = PDFData()
    pdfcomp.setData(pdfpat)
    pdfpat.loadData("ni.dat", PDFParser())
    pdfpat.qdamp = 0.06
    fit.mapVP("v_qdamp", pdfpat, "qdamp")

    # setup the excluded region
    pdfcomp.setFitRange(1.5, 20.0, 0.05)

    ### Instantiate the phase
    pdfpha = PDFPhase()
    pdfpha.setName("PDFNi")
    pdfcomp.addPhase(pdfpha)
    pdfpha.weight = 0.72
    fit.mapVP("v_weight", pdfpha, "weight")
    pdfpha.delta2 = 0.07
    fit.mapVP("v_delta2", pdfpha, "delta2")

    # Create a lattice for the phase
    pdflat = PDFLattice()
    pdfpha.setLattice(pdflat)
    pdflat.a = pdflat.b = pdflat.c = 3.524
    fit.mapVP("v_a", pdflat, "a")
    fit.mapVP("v_a", pdflat, "b")
    fit.mapVP("v_a", pdflat, "c")
    # Add some atom
    for i in range(4):
        a = PDFAtom("Ni")
        a.x = a.y = a.z = 0
        a.setADP( PDFIADP() )
        fit.mapVP("v_biso", a.getADP(), "Biso")
        fit.guessV("v_biso", 0.60)
        pdfpha.addAtom(a)
    pdfpha.getAtom(1).y = pdfpha.getAtom(1).z = 0.5
    pdfpha.getAtom(2).x = pdfpha.getAtom(2).z = 0.5
    pdfpha.getAtom(3).x = pdfpha.getAtom(3).y = 0.5

    # Fitting
    fit.setWeight(pdfcomp, 0.5)
    fit.setWeight(difcomp, 0.5)
    fit.refine()
    fit.printResults()
    fit.saveResults("ni.co.res")
    pdfcomp.saveFitArrays("ni.co.pdf.fit")
    difcomp.saveFitArrays("ni.co.dif.fit")

if __name__ == "__main__":
    from SrRietveld import Fit
    fit   = Fit()
    runfit(fit)
