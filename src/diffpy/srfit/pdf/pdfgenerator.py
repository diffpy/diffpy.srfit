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
"""PDF profile generator.

The PDFGenerator class can take a diffpy.structure,
pyobjcryst.crystal.Crystal or pyobjcryst.molecule.Molecule object and
calculate the crystal PDF from it. The passed structure object is
wrapped in a StructureParameter set, which makes its attributes
refinable. See the class definition for more details and the examples
for its use.
"""

__all__ = ["PDFGenerator"]

from diffpy.srfit.pdf.basepdfgenerator import BasePDFGenerator


class PDFGenerator(BasePDFGenerator):
    """A class for calculating the PDF from a single crystal structure.

    This works with diffpy.structure.Structure, pyobjcryst.crystal.Crystal and
    pyobjcryst.molecule.Molecule instances. Note that the managed Parameters
    are not created until the structure is added.

    Attributes:
    _calc   --  PDFCalculator instance for calculating the PDF
    _phase  --  The structure ParameterSet used to calculate the profile.
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

    def __init__(self, name="pdf"):
        """Initialize the generator."""
        from diffpy.srreal.pdfcalculator import PDFCalculator

        BasePDFGenerator.__init__(self, name)
        self._setCalculator(PDFCalculator())
        return


# End class PDFGenerator
