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
"""PDF calculation tools."""

__all__ = ["PDFGenerator", "DebyePDFGenerator", "PDFContribution", "PDFParser"]

import warnings

from diffpy.srfit.pdf.debyepdfgenerator import DebyePDFGenerator
from diffpy.srfit.pdf.pdfcontribution import PDFContribution
from diffpy.srfit.pdf.pdfgenerator import PDFGenerator
from diffpy.srfit.pdf.pdfparser import PDFParser
from diffpy.utils._deprecator import build_deprecation_message

removal_version = "4.0.0"

warnings.warn(
    build_deprecation_message(
        "diffpy.srfit",
        "pdf",
        "cmipdf",
        removal_version,
        new_base="diffpy",
    ),
    DeprecationWarning,
    stacklevel=2,
)

# End of file
