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
"""PDF profile generator using the Debye equation.

The DebyePDFGenerator class can take a diffpy.Structure,
pyobjcryst.crystal.Crystal or pyobjcryst.molecule.Molecule object and calculate
the PDF from it. This generator is especially appropriate for isolated
scatterers, such as nanoparticles and molecules.

"""
__all__ = ["PDFDebyeGenerator"]

from diffpy.srreal.pdfcalculator import DebyePDFCalculator
from diffpy.srfit.pdf.basepdfgenerator import BasePDFGenerator

class DebyePDFGenerator(BasePDFGenerator):
    """A class for calculating the PDF from an isolated scatterer.

    This works with diffpy.Structure.Structure, pyobjcryst.crystal.Crystal and
    pyobjcryst.molecule.Molecule instances. Note that the managed Parameters
    are not created until the structure is added.

    Attributes:
    _calc   --  DebyePDFCalculator_ext instance for calculating the PDF
    _phase  --  The structure ParameterSets used to calculate the profile.
    _lastr  --  The last value of r over which the PDF was calculated. This is
                used to configure the calculator when r changes.

    Managed Parameters:
    scale   --  Scale factor
    delta1  --  Linear peak broadening term
    delta2  --  Quadratic peak broadening term
    qbroad  --  Resolution peak broadening term
    qdamp   --  Resolution peak dampening term

    Managed ParameterSets:
    The structure ParameterSet (BaseStructure instance) used to calculate the
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

    def setPhase(self, stru = None, name = None, parset = None, periodic =
            False):
        """Set the phase that will be used to calculate the PDF.

        This creates a StructureParSet or ObjCrystParSet that adapts stru to a
        ParameterSet interface. See those classes (located in
        diffpy.srfit.structure) for how they are used. The resulting
        ParameterSet will be managed by this generator.

        stru    --  diffpy.Structure.Structure, pyobjcryst.crystal.Crystal or
                    pyobjcryst.molecule.Molecule instance . Default None.
        name    --  A name to give the structure. If name is None (default),
                    then the name will be set as "phase".
        parset  --  A ParameterSet that holds the structural information. This
                    can be used to share the phase between multiple
                    PDFGenerators, and have the changes in one reflect in
                    another. If both stru and parset are specified, only parset
                    is used. Default None. 
        periodic -- The structure should be treated as periodic (default
                    False). Note that some structures do not support
                    periodicity, in which case this will be ignored.

        Raises ValueError if neither stru nor parset is specified.

        """
        BasePDFGenerator.setPhase(self, stru, name, parset, periodic)

    def __init__(self, name = "pdf"):
        """Initialize the generator.
        
        """
        BasePDFGenerator.__init__(self, name)
        self._calc = DebyePDFCalculator()
        return

# End class DebyePDFGenerator

# End of file
