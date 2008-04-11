################################################################################
#
#   PDF and Rietveld Corefinement of Ni data.
#
################################################################################
__id__ = "$Id$"

def runfit(fit):

    import difscript

    # Do the rietveld refinement first (it's more tempermental)
    difscript.initialize(fit)
    difcomp = fit.getComponent(0)
    # Let's see where we're starting
    if 1:
        calc = difcomp.getCalculator()
        calc._prepareSrRietUI()
        calc.calculate()
        from pylab import plot, show, title, xlabel, ylabel, legend
        x0,y0,u0 = difcomp.getFitDataArrays()
        x1,y1 = difcomp.getFitArrays()
        dy = y0-y1
        dy -= 1.5*( abs(min(dy)) + abs(min(y0)))
        plot(x1, y0, 'bo', x1, y1, 'r-', x1, dy, 'g--')
        title("Initial Rietveld configuration")
        xlabel("Q (A^-1)")
        ylabel("I (arb.)")
        legend(("data", "fit (R = %f)"%difcomp.getR(), "difference"))
        show()
    print "Rietveld fit"
    difscript.runfit(fit)

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

    # Check where we're at
    if 1:
        calc = pdfcomp.getCalculator()
        calc.calculate()
        from pylab import plot, show, title, xlabel, ylabel, legend
        x0,y0,u0 = pdfcomp.getFitDataArrays()
        x1,y1 = pdfcomp.getFitArrays()
        dy = y0-y1
        dy -= 1.5*( abs(min(dy)) + abs(min(y0)))
        plot(x1, y0, 'bo', x1, y1, 'r-', x1, dy, 'g--')
        title("Initial PDF configuration")
        xlabel("r (A)")
        ylabel("G (A^-2)")
        legend(("data", "fit (R = %f)"%pdfcomp.getR(), "difference"))
        show()
        raw_input("Press any key to continue ... ")

    # Fitting
    fit.setWeight(pdfcomp, 0.5)
    fit.setWeight(difcomp, 0.5)
    print "Rietveld+PDF fit"
    fit.refine()
    fit.printResults()
    fit.saveResults("ni.co.res")
    pdfcomp.saveFitArrays("ni.co.pdf.fit")
    difcomp.saveFitArrays("ni.co.dif.fit")

    # Show final Rietveld fit
    if 1:
        from pylab import plot, show, title, xlabel, ylabel, legend
        x0,y0,u0 = difcomp.getFitDataArrays()
        x1,y1 = difcomp.getFitArrays()
        dy = y0-y1
        dy -= 1.5*( abs(min(dy)) + abs(min(y0)))
        plot(x1, y0, 'bo', x1, y1, 'r-', x1, dy, 'g--')
        title("Rietveld fit from co-refinement")
        xlabel("Q (A^-1)")
        ylabel("I (arb.)")
        legend(("data", "fit (R = %f)"%difcomp.getR(), "difference"))
        show()
        raw_input("Press any key to continue ... ")

    
    # Show final PDF fit
    if 1:
        calc = pdfcomp.getCalculator()
        calc.calculate()
        from pylab import plot, show, title, xlabel, ylabel, legend
        x0,y0,u0 = pdfcomp.getFitDataArrays()
        x1,y1 = pdfcomp.getFitArrays()
        dy = y0-y1
        dy -= 1.5*( abs(min(dy)) + abs(min(y0)))
        plot(x1, y0, 'bo', x1, y1, 'r-', x1, dy, 'g--')
        title("PDF fit from co-refinement")
        xlabel("r (A)")
        ylabel("G (A^-2)")
        legend(("data", "fit (R = %f)"%pdfcomp.getR(), "difference"))
        show()
        raw_input("Press any key to continue ... ")

if __name__ == "__main__":
    from SrRietveld import Fit
    fit   = Fit()
    runfit(fit)
