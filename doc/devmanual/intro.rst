.. _developers-manual-intro:

===================
Introduction
===================

SrFit modular configurable structure modeling code. Traditional structure
refinement and modeling codes are monolithic applications for a particular
task: refinement of structure models for crystals from powder data (Rietveld),
refinement of structure models from single-crystal data, structure solution
from single-crystal or powder data, and so on.

These codes are not flexible enough for modern materials structure problems.
The types of data that must be used, and the classes of models employed are
becoming more diverse. Increasingly we need a model that is specific for a
molecular system that is best expressed as a Z-matrix, or for a discrete
nanoparticle, and so on. Data may include time-of-flight neutron powder
diffraction data that are analyzed in reciprocal-space, or real-space or both,
x-ray diffraction data, EXAFS data and so on. We may want to employ a fast
local-search regression algorithm such as Levenberg-Marquardt, or increase
convergence using a global Monte-Carlo or genetic algorithm, or in general
combine both approaches in a basin-hopping scheme.

This diversity is not served by current codes. SrFit is an application that
will allow a modeling code to be built on the fly from components such as
function calculators (that calculate different data spectra), regression
algorithms and structure models. The target function being optimized can be
specified by the user by adding terms for more than one function calculator
depending on the data available. SrFit will be extensible: as new opportunities
and methods arise new modules can be written in the future.

SrFit will support real-space (PDF) and small-angle reciprocal space (SAS)
refinement. SrFit will contain support for parallel optimizers. Tools for post
analysis will be provided, such as parametric plotting, calculation of derived
quantities such as bond-lengths, bond-angles, bond valence sums etc., automatic
table generation in LaTex, and so on. Refinement results will be stored in a
database structure to aid post analysis. We will also support CIF input and
output of structure models. SrFit will support user-specified analytic math
functions and numerical algorithms for rapid prototyping of new ideas. A
contingent item will be to support experimental feedback to the
data-acquisition system. 

Features
-----------

 * Supports simultaneous PDF refinement of one or more crystal phases from one
   or more data sets.
 * Refine diffpy.Structure or pyobjcryst.crystal.Crystal objects to data.
 * Fast PDF calculation provided by SrReal.
   (http://danse.us/trac/diffraction/wiki/SrReal)
 * Powerful molecule constraints and restraints from ObjCryst++ via PyObjCryst
   (http://vincefn.net/ObjCryst/)
 * Automatic and user-defined space group constraints.
 * Data parser compatible with PDFGetX2 and PDFGetN output.
 * Automatic configuration from parsed metadata.
 * Small-angle scattering models leveraging the DANSE SANS software
   (http://danse.chem.utk.edu/).
 * Various nanoparticle form factors for nanoparticle PDF modeling - extended
   by small-angle scattering inversion.

Release Definition
--------------------

The SrFit Alpha-9 release definition can be found at
http://danse.us/trac/diffraction/wiki/SrFitReleaseAlpha9.

Major changes since last release
----------------------------------------

 * Support for parallel PDF calculators from SrReal.
 * Various interface changes.

   * pdfnpgenerator module renamed debyepdfgenerator.
   * nanoformfactors module renamed characteristicfunctions.  
   * Form factors renamed from xyzFF to xyzCF.
   * Structure adapters renamed for clarification. The monikers "Structure" and
     "stru" are reserved for structure representations from outside of SrFit
     (e.g. diffpy.Structure.Structure).  The monikers "ParSet" and "Phase" are
     used to represent the adapted structures (e.g. DiffpyStructureParSet).

 * PDFContribution class introduced to simpliy PDF fits. See :doc:`simplepdf`
   and :doc:`simplepdftwophase` examples.
 * Classes with parameters now support list-like parameter access.
 * Support for sequence of FitHooks.
 * Parameter values can be assigned with an ``=`` sign from containers. For
   example, recipe.a = 3.4 now works. Update operations, such as ``+=`` are not
   supported.
 * FitResults now reports on fixed variables.
 * New ``show`` method prints a summary of the fit configuration.  The display
   (and implementatation) are a work-in-progress.
 * Constraints and restraints can be cleared recursively.
