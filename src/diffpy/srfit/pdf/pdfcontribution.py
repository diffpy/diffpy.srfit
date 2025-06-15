#!/usr/bin/env python
##############################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2010 The Trustees of Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################
"""PDFContribution class.

This is a custom FitContribution that simplifies the creation of PDF
fits.
"""

__all__ = ["PDFContribution"]

from diffpy.srfit.fitbase import FitContribution, Profile


class PDFContribution(FitContribution):
    """PDFContribution class.

    PDFContribution is a FitContribution that is customized for PDF fits. Data
    and phases can be added directly to the PDFContribution. Setup of
    constraints and restraints requires direct interaction with the generator
    attributes (see setPhase).

    Attributes
    name            --  A name for this FitContribution.
    profile         --  A Profile that holds the measured (and calculated)
                        signal.
    _meta           --  Metadata dictionary. This is specific to this object,
                        and not shared with the profile. This is used to record
                        configuration options, like qmax.
    _calculators    --  A managed dictionary of Calculators, indexed by name.
    _constraints    --  A set of constrained Parameters. Constraints can be
                        added using the 'constrain' methods.
    _generators     --  A managed dictionary of ProfileGenerators.
    _parameters     --  A managed OrderedDict of parameters.
    _restraints     --  A set of Restraints. Restraints can be added using the
                        'restrain' or 'confine' methods.
    _parsets        --  A managed dictionary of ParameterSets.
    _eqfactory      --  A diffpy.srfit.equation.builder.EquationFactory
                        instance that is used to create constraints and
                        restraints from string
    _eq             --  The FitContribution equation that will be optimized.
    _reseq          --  The residual equation.
    _xname          --  Name of the x-variable
    _yname          --  Name of the y-variable
    _dyname         --  Name of the dy-variable

    Managed Parameters:
    scale   --  Scale factor
    qbroad  --  Resolution peak broadening term
    qdamp   --  Resolution peak dampening term
    """

    def __init__(self, name):
        """Create the PDFContribution.

        name        --  The name of the contribution.
        """
        FitContribution.__init__(self, name)
        self._meta = {}
        # Add the profile
        profile = Profile()
        self.setProfile(profile, xname="r")

        # Need a parameter for the overall scale, in the case that this is a
        # multi-phase fit.
        self.newParameter("scale", 1.0)
        # Profile-related parameters that will be shared between the generators
        self.newParameter("qdamp", 0)
        self.newParameter("qbroad", 0)
        return

    # Data methods

    def loadData(self, data):
        """Load the data in various formats.

        This uses the PDFParser to load the data and then passes it to the
        built-in profile with loadParsedData.

        data    --  An open file-like object, name of a file that contains data
                    or a string containing the data.
        """
        # Get the data into a string
        from diffpy.srfit.util.inpututils import inputToString

        datstr = inputToString(data)

        # Load data with a PDFParser
        from diffpy.srfit.pdf.pdfparser import PDFParser

        parser = PDFParser()
        parser.parseString(datstr)

        # Pass it to the profile
        self.profile.loadParsedData(parser)
        return

    def setCalculationRange(self, xmin=None, xmax=None, dx=None):
        """Set epsilon-inclusive calculation range.

        Adhere to the observed ``xobs`` points when ``dx`` is the same
        as in the data.  ``xmin`` and ``xmax`` are clipped at the bounds
        of the observed data.

        Parameters
        ----------

        xmin : float or "obs", optional
            The minimum value of the independent variable.  Keep the
            current minimum when not specified.  If specified as "obs"
            reset to the minimum observed value.
        xmax : float or "obs", optional
            The maximum value of the independent variable.  Keep the
            current maximum when not specified.  If specified as "obs"
            reset to the maximum observed value.
        dx : float or "obs", optional
            The sample spacing in the independent variable.  When different
            from the data, resample the ``x`` as anchored at ``xmin``.

        Note that xmin is always inclusive (unless clipped). xmax is inclusive
        if it is within the bounds of the observed data.

        Raises
        ------
        AttributeError
            If there is no observed data.
        ValueError
            When xmin > xmax or if dx <= 0.  Also if dx > xmax - xmin.
        """
        return self.profile.setCalculationRange(xmin, xmax, dx)

    def savetxt(self, fname, **kwargs):
        """Call numpy.savetxt with x, ycalc, y, dy.

        This calls on the built-in Profile.

        Arguments are passed to numpy.savetxt.
        """
        return self.profile.savetxt(fname, **kwargs)

    # Phase methods

    def addStructure(self, name, stru, periodic=True):
        """Add a phase that goes into the PDF calculation.

        name    --  A name to give the generator that will manage the PDF
                    calculation from the passed structure. The adapted
                    structure will be accessible via the name "phase" as an
                    attribute of the generator, e.g.
                    contribution.name.phase, where 'contribution' is this
                    contribution and 'name' is passed name.
                    (default), then the name will be set as "phase".
        stru    --  diffpy.structure.Structure, pyobjcryst.crystal.Crystal or
                    pyobjcryst.molecule.Molecule instance.  Default None.
        periodic -- The structure should be treated as periodic.  If this is
                    True (default), then a PDFGenerator will be used to
                    calculate the PDF from the phase. Otherwise, a
                    DebyePDFGenerator will be used. Note that some structures
                    do not support periodicity, in which case this may be
                    ignored.

        Returns the new phase (ParameterSet appropriate for what was passed in
        stru.)
        """
        # Based on periodic, create the proper generator.
        if periodic:
            from diffpy.srfit.pdf.pdfgenerator import PDFGenerator

            gen = PDFGenerator(name)
        else:
            from diffpy.srfit.pdf.debyepdfgenerator import DebyePDFGenerator

            gen = DebyePDFGenerator(name)

        # Set up the generator
        gen.setStructure(stru, "phase", periodic)
        self._setupGenerator(gen)

        return gen.phase

    def addPhase(self, name, parset, periodic=True):
        """Add a phase that goes into the PDF calculation.

        name    --  A name to give the generator that will manage the PDF
                    calculation from the passed parameter phase. The parset
                    will be accessible via the name "phase" as an attribute
                    of the generator, e.g., contribution.name.phase, where
                    'contribution' is this contribution and 'name' is passed
                    name.
        parset  --  A SrRealParSet that holds the structural information.
                    This can be used to share the phase between multiple
                    BasePDFGenerators, and have the changes in one reflect in
                    another.
        periodic -- The structure should be treated as periodic.  If this is
                    True (default), then a PDFGenerator will be used to
                    calculate the PDF from the phase. Otherwise, a
                    DebyePDFGenerator will be used. Note that some structures
                    do not support periodicity, in which case this may be
                    ignored.

        Returns the new phase (ParameterSet appropriate for what was passed in
        stru.)
        """
        # Based on periodic, create the proper generator.
        if periodic:
            from diffpy.srfit.pdf.pdfgenerator import PDFGenerator

            gen = PDFGenerator(name)
        else:
            from diffpy.srfit.pdf.debyepdfgenerator import DebyePDFGenerator

            gen = DebyePDFGenerator(name)

        # Set up the generator
        gen.setPhase(parset, periodic)
        self._setupGenerator(gen)

        return gen.phase

    def _setupGenerator(self, gen):
        """Setup a generator.

        The generator must already have a managed SrRealParSet, added
        with setStructure or setPhase.
        """
        # Add the generator to this FitContribution
        self.addProfileGenerator(gen)

        # Set the proper equation for the fit, depending on the number of
        # phases we have.
        gnames = self._generators.keys()
        eqstr = " + ".join(gnames)
        eqstr = "scale * (%s)" % eqstr
        self.setEquation(eqstr)

        # Update with our metadata
        gen.meta.update(self._meta)
        gen.processMetaData()

        # Constrain the shared parameters
        self.constrain(gen.qdamp, self.qdamp)
        self.constrain(gen.qbroad, self.qbroad)
        return

    # Calculation setup methods

    def _getMetaValue(self, kwd):
        """Get metadata according to object hierarchy."""
        # Check self, then generators then profile
        if kwd in self._meta:
            return self._meta[kwd]
        for gen in self._generators.values():
            if kwd in gen.meta:
                return gen.meta[kwd]
        val = self.profile.meta.get(kwd)
        return val

    def setScatteringType(self, type="X"):
        """Set the scattering type.

        type    --   "X" for x-ray or "N" for neutron

        Raises ValueError if type is not "X" or "N"
        """
        self._meta["stype"] = type
        for gen in self._generators.values():
            gen.setScatteringType(type)
        return

    def getScatteringType(self):
        """Get the scattering type.

        See 'setScatteringType'.
        """
        return self._getMetaValue("stype")

    def setQmax(self, qmax):
        """Set the qmax value."""
        self._meta["qmax"] = qmax
        for gen in self._generators.values():
            gen.setQmax(qmax)
        return

    def getQmax(self):
        """Get the qmax value."""
        return self._getMetaValue("qmax")

    def setQmin(self, qmin):
        """Set the qmin value."""
        self._meta["qmin"] = qmin
        for gen in self._generators.values():
            gen.setQmin(qmin)
        return

    def getQmin(self):
        """Get the qmin value."""
        return self._getMetaValue("qmin")


# End of file
