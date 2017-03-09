#!/usr/bin/env python
##############################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2008 The Trustees of Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################

"""PDF profile generator using the Debye equation.

The DebyePDFGenerator class can take a diffpy.Structure,
pyobjcryst.crystal.Crystal or pyobjcryst.molecule.Molecule object and calculate
the PDF from it. This generator is especially appropriate for isolated
scatterers, such as nanoparticles and molecules.
"""

__all__ = ["DebyePDFGenerator"]

from diffpy.srreal.pdfcalculator import DebyePDFCalculator
from diffpy.srfit.pdf.basepdfgenerator import BasePDFGenerator

class DebyePDFGenerator(BasePDFGenerator):
    """A class for calculating the PDF from an isolated scatterer.

    This works with diffpy.Structure.Structure, pyobjcryst.crystal.Crystal and
    pyobjcryst.molecule.Molecule instances. Note that the managed Parameters
    are not created until the structure is added.

    Attributes:
    _calc   --  DebyePDFCalculator instance for calculating the PDF
    _phase  --  The structure ParameterSets used to calculate the profile.
    stru    --  The structure objected adapted by _phase.
    _lastr  --  The last value of r over which the PDF was calculated. This is
                used to configure the calculator when r changes.

    Managed Parameters:
    scale   --  Scale factor
    delta1  --  Linear peak broadening term
    delta2  --  Quadratic peak broadening term
    qbroad  --  Resolution peak broadening term
    qdamp   --  Resolution peak dampening term

    Managed ParameterSets:
    The structure ParameterSet (SrRealStructure instance) used to calculate the
    profile is named by the user.

    Usable Metadata:
    stype   --  The scattering type "X" for x-ray, "N" for neutron (see
                'setScatteringType').
    qmax    --  The maximum scattering vector used to generate the PDF (see
                setQmax).
    qmin    --  The minimum scattering vector used to generate the PDF (see
                setQmin).
    scale   --  See Managed Parameters.
    delta1  --  See Managed Parameters.
    delta2  --  See Managed Parameters.
    qbroad  --  See Managed Parameters.
    qdamp   --  See Managed Parameters.

    """

    def setStructure(self, stru, name = "phase", periodic = False):
        """Set the structure that will be used to calculate the PDF.

        This creates a DiffpyStructureParSet, ObjCrystCrystalParSet or
        ObjCrystMoleculeParSet that adapts stru to a ParameterSet interface.
        See those classes (located in diffpy.srfit.structure) for how they are
        used. The resulting ParameterSet will be managed by this generator.

        stru    --  diffpy.Structure.Structure, pyobjcryst.crystal.Crystal or
                    pyobjcryst.molecule.Molecule instance.  Default None.
        name    --  A name to give to the managed ParameterSet that adapts stru
                    (default "phase").
        periodic -- The structure should be treated as periodic (default
                    False). Note that some structures do not support
                    periodicity, in which case this will have no effect on the
                    PDF calculation.

        """
        return BasePDFGenerator.setStructure(self, stru, name, periodic)


    def setPhase(self, parset, periodic = False):
        """Set the phase that will be used to calculate the PDF.

        Set the phase directly with a DiffpyStructureParSet,
        ObjCrystCrystalParSet or ObjCrystMoleculeParSet that adapts a structure
        object (from diffpy or pyobjcryst).  The passed ParameterSet will be
        managed by this generator.

        parset  --  A SrRealParSet that holds the structural information.
                    This can be used to share the phase between multiple
                    BasePDFGenerators, and have the changes in one reflect in
                    another.
        periodic -- The structure should be treated as periodic (default True).
                    Note that some structures do not support periodicity, in
                    which case this will be ignored.

        """
        return BasePDFGenerator.setPhase(self, parset, periodic)


    def __init__(self, name = "pdf"):
        """Initialize the generator.

        """
        BasePDFGenerator.__init__(self, name)
        self._setCalculator(DebyePDFCalculator())
        return

# End class DebyePDFGenerator

# End of file
