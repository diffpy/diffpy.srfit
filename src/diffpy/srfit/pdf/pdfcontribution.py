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
from diffpy.utils._deprecator import build_deprecation_message, deprecated

removal_version = "4.0.0"
contribution_base = "diffpy.srfit.pdf.pdfcontribution.PDFContribution"
contribution_new_base = "diffpy.cmipdf.pdfcontribution.PDFContribution"


def _build_contribution_dep_msg(old_name, new_name):
    message = build_deprecation_message(
        contribution_base,
        old_name,
        new_name,
        removal_version,
        new_base=contribution_new_base,
    )
    return message


loadData_dep_msg = _build_contribution_dep_msg("loadData", "load_data")
setCalculationRange_dep_msg = _build_contribution_dep_msg(
    "setCalculationRange", "set_calculation_range"
)
addStructure_dep_msg = (
    _build_contribution_dep_msg("addStructure", "add_structure")
    + " Additionally, the signature has changed. Please pass the structure "
    "before the name, as add_structure(structure, name)."
)
addPhase_dep_msg = (
    _build_contribution_dep_msg("addPhase", "add_structure_from_parset")
    + " Additionally, the signature has changed. Please pass the parset "
    "before the name, as add_structure_from_parset(parset, name)."
)
setScatteringType_dep_msg = _build_contribution_dep_msg(
    "setScatteringType", "set_scattering_type"
)
getScatteringType_dep_msg = _build_contribution_dep_msg(
    "getScatteringType", "get_scattering_type"
)
setQmax_dep_msg = _build_contribution_dep_msg("setQmax", "set_qmax")
getQmax_dep_msg = _build_contribution_dep_msg("getQmax", "get_qmax")
setQmin_dep_msg = _build_contribution_dep_msg("setQmin", "set_qmin")
getQmin_dep_msg = _build_contribution_dep_msg("getQmin", "get_qmin")


class PDFContribution(FitContribution):
    """A FitContribution that is customized for PDF fits.

    Data and phases can be added directly to the PDFContribution. Setup of
    constraints and restraints requires direct interaction with the generator
    attributes (see add_structure_from_parset).

    Attributes
    ----------
    name
        A name for this FitContribution.
    profile
        A Profile that holds the measured (and calculated)
        signal.
    _meta
        Metadata dictionary. This is specific to this object,
        and not shared with the profile. This is used to record
        configuration options, like qmax.
    _calculators
        A managed dictionary of Calculators, indexed by name.
    _constraints
        A set of constrained Parameters. Constraints can be
        added using the 'constrain' methods.
    _generators
        A managed dictionary of ProfileGenerators.
    _parameters
        A managed OrderedDict of parameters.
    _restraints
        A set of Restraints. Restraints can be added using the
        'restrain' or 'confine' methods.
    _parsets
        A managed dictionary of ParameterSets.
    _eqfactory
        A diffpy.srfit.equation.builder.EquationFactory
        instance that is used to create constraints and
        restraints from string
    _eq
        The FitContribution equation that will be optimized.
    _reseq
        The residual equation.
    _xname
        Name of the x-variable
    _yname
        Name of the y-variable
    _dyname
        Name of the dy-variable

    Managed Parameters
    ------------------
    scale
        Scale factor
    qbroad
        Resolution peak broadening term
    qdamp
        Resolution peak dampening term
    """

    def __init__(self, name):
        """Create the PDFContribution.

        Parameters
        ----------
        name : str
            The name of the contribution.
        """
        FitContribution.__init__(self, name)
        self._meta = {}
        # Add the profile
        profile = Profile()
        self.set_profile(profile, xname="r")

        # Need a parameter for the overall scale, in the case that this is a
        # multi-phase fit.
        self.new_parameter("scale", 1.0)
        # Profile-related parameters that will be shared between the generators
        self.new_parameter("qdamp", 0)
        self.new_parameter("qbroad", 0)
        return

    # Data methods

    def load_data(self, datafile):
        """Load the data from a data file.

        This uses the PDFParser to load the data and then passes it to
        the built-in profile with load_parsed_data.

        Parameters
        ----------
        datafile : str or Path
            The path to the file that contains the data.
        """
        from diffpy.srfit.pdf.pdfparser import PDFParser

        parser = PDFParser()
        parser.parse_file(datafile)

        # Pass it to the profile
        self.profile.load_parsed_data(parser)
        return

    @deprecated(loadData_dep_msg)
    def loadData(self, datafile):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use
        diffpy.cmipdf.pdfcontribution.PDFContribution.load_data
        instead.
        """
        return self.load_data(datafile)

    def set_calculation_range(self, xmin=None, xmax=None, dx=None):
        """Set epsilon-inclusive calculation range.

        Adhere to the observed ``xobs`` points when ``dx`` is the same
        as in the data.  ``xmin`` and ``xmax`` are clipped at the bounds
        of the observed data.

        Parameters
        ----------
        xmin : float or ``obs``, optional
            The minimum value of the independent variable.  Keep the
            current minimum when not specified.  If specified as "obs"
            reset to the minimum observed value.
        xmax : float or ``obs``, optional
            The maximum value of the independent variable.  Keep the
            current maximum when not specified.  If specified as "obs"
            reset to the maximum observed value.
        dx : float or ``obs``, optional
            The sample spacing in the independent variable.  When different
            from the data, resample the ``x`` as anchored at ``xmin``.

        Note that ``xmin`` is always inclusive (unless clipped).
        ``xmax`` is inclusive if it is within the bounds of the observed data.

        Raises
        ------
        AttributeError
            If there is no observed data.
        ValueError
            When xmin > xmax or if dx <= 0.  Also if dx > xmax - xmin.
        """
        return self.profile.set_calculation_range(xmin, xmax, dx)

    @deprecated(setCalculationRange_dep_msg)
    def setCalculationRange(self, xmin=None, xmax=None, dx=None):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use
        diffpy.cmipdf.pdfcontribution.PDFContribution.set_calculation_range
        instead.
        """
        return self.set_calculation_range(xmin, xmax, dx)

    def savetxt(self, fname, **kwargs):
        """Call numpy.savetxt with x, ycalc, y, dy.

        This calls on the built-in Profile.

        Parameters
        ----------
        fname : str or Path
            The file or filename to which the data is saved.
        **kwargs
            Additional arguments passed to numpy.savetxt.
        """
        return self.profile.savetxt(fname, **kwargs)

    # Phase methods

    def add_structure(self, structure, name="phase", periodic=True):
        """Add a phase that goes into the PDF calculation.

        Parameters
        ----------
        structure
            The diffpy.structure.Structure, pyobjcryst.crystal.Crystal or
            pyobjcryst.molecule.Molecule instance to adapt.
        name : str, optional
            A name to give the generator that will manage the PDF
            calculation from the passed structure. The adapted
            structure will be accessible via the name "phase" as an
            attribute of the generator, e.g.
            contribution.name.phase, where 'contribution' is this
            contribution and 'name' is passed name (default "phase").
        periodic : bool, optional
            The structure should be treated as periodic.  If this is
            True (default), then a PDFGenerator will be used to
            calculate the PDF from the phase. Otherwise, a
            DebyePDFGenerator will be used. Note that some structures
            do not support periodicity, in which case this may be
            ignored.

        Returns
        -------
        ParameterSet
            The new phase, of the type appropriate for what was passed
            in structure.
        """
        # Based on periodic, create the proper generator.
        if periodic:
            from diffpy.srfit.pdf.pdfgenerator import PDFGenerator

            gen = PDFGenerator(name)
        else:
            from diffpy.srfit.pdf.debyepdfgenerator import DebyePDFGenerator

            gen = DebyePDFGenerator(name)

        # Set up the generator
        gen.set_structure(structure, "phase", periodic)
        self._setup_generator(gen)

        return gen.phase

    @deprecated(addStructure_dep_msg)
    def addStructure(self, name, stru, periodic=True):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use
        diffpy.cmipdf.pdfcontribution.PDFContribution.add_structure
        instead.
        """
        return self.add_structure(stru, name, periodic)

    def add_structure_from_parset(self, parset, name, periodic=True):
        """Add a phase that goes into the PDF calculation.

        Parameters
        ----------
        parset
            The SrRealParSet that holds the structural information.
            This can be used to share the phase between multiple
            BasePDFGenerators, and have the changes in one reflect in
            another.
        name : str
            A name to give the generator that will manage the PDF
            calculation from the passed parameter phase. The parset
            will be accessible via the name "phase" as an attribute
            of the generator, e.g., contribution.name.phase, where
            'contribution' is this contribution and 'name' is passed
            name.
        periodic : bool, optional
            The structure should be treated as periodic.  If this is
            True (default), then a PDFGenerator will be used to
            calculate the PDF from the phase. Otherwise, a
            DebyePDFGenerator will be used. Note that some structures
            do not support periodicity, in which case this may be
            ignored.

        Returns
        -------
        ParameterSet
            The new phase, of the type appropriate for what was passed
            in parset.
        """
        # Based on periodic, create the proper generator.
        if periodic:
            from diffpy.srfit.pdf.pdfgenerator import PDFGenerator

            gen = PDFGenerator(name)
        else:
            from diffpy.srfit.pdf.debyepdfgenerator import DebyePDFGenerator

            gen = DebyePDFGenerator(name)

        # Set up the generator
        gen.set_structure_from_parset(parset, periodic)
        self._setup_generator(gen)

        return gen.phase

    @deprecated(addPhase_dep_msg)
    def addPhase(self, name, parset, periodic=True):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use
        diffpy.cmipdf.pdfcontribution.PDFContribution.add_structure_from_parset
        instead.
        """
        return self.add_structure_from_parset(parset, name, periodic)

    def _setup_generator(self, gen):
        """Setup a generator.

        The generator must already have a managed SrRealParSet, added
        with set_structure or set_structure_from_parset.
        """
        # Add the generator to this FitContribution
        self.add_profile_generator(gen)

        # Set the proper equation for the fit, depending on the number of
        # phases we have.
        gnames = self._generators.keys()
        eqstr = " + ".join(gnames)
        eqstr = "scale * (%s)" % eqstr
        self.set_equation(eqstr)

        # Update with our metadata
        gen.meta.update(self._meta)
        gen._process_metadata()

        # Constrain the shared parameters
        self.add_constraint(gen.qdamp, self.qdamp)
        self.add_constraint(gen.qbroad, self.qbroad)
        return

    # Calculation setup methods

    def _get_meta_value(self, kwd):
        """Get metadata according to object hierarchy."""
        # Check self, then generators then profile
        if kwd in self._meta:
            return self._meta[kwd]
        for gen in self._generators.values():
            if kwd in gen.meta:
                return gen.meta[kwd]
        val = self.profile.meta.get(kwd)
        return val

    def set_scattering_type(self, type="X"):
        """Set the scattering type.

        Parameters
        ----------
        type : str, optional
            "X" for x-ray or "N" for neutron (default "X").

        Raises
        ------
        ValueError
            If type is not "X" or "N".
        """
        self._meta["stype"] = type
        for gen in self._generators.values():
            gen.set_scattering_type(type)
        return

    @deprecated(setScatteringType_dep_msg)
    def setScatteringType(self, type="X"):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use
        diffpy.cmipdf.pdfcontribution.PDFContribution.set_scattering_type
        instead.
        """
        return self.set_scattering_type(type)

    def get_scattering_type(self):
        """Get the scattering type.

        See ``set_scattering_type``.

        Returns
        -------
        str
            The scattering type used to calculate the PDF.
        """
        scattering_type = self._get_meta_value("stype")
        return scattering_type

    @deprecated(getScatteringType_dep_msg)
    def getScatteringType(self):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use
        diffpy.cmipdf.pdfcontribution.PDFContribution.get_scattering_type
        instead.
        """
        return self.get_scattering_type()

    def set_qmax(self, qmax):
        """Set the qmax value.

        Parameters
        ----------
        qmax : float
            The maximum scattering vector used to generate the PDF.
        """
        self._meta["qmax"] = qmax
        for gen in self._generators.values():
            gen.set_qmax(qmax)
        return

    @deprecated(setQmax_dep_msg)
    def setQmax(self, qmax):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use
        diffpy.cmipdf.pdfcontribution.PDFContribution.set_qmax
        instead.
        """
        return self.set_qmax(qmax)

    def get_qmax(self):
        """Get the qmax value.

        Returns
        -------
        float
            The maximum scattering vector used to generate the PDF.
        """
        qmax = self._get_meta_value("qmax")
        return qmax

    @deprecated(getQmax_dep_msg)
    def getQmax(self):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use
        diffpy.cmipdf.pdfcontribution.PDFContribution.get_qmax
        instead.
        """
        return self.get_qmax()

    def set_qmin(self, qmin):
        """Set the qmin value.

        Parameters
        ----------
        qmin : float
            The minimum scattering vector used to generate the PDF.
        """
        self._meta["qmin"] = qmin
        for gen in self._generators.values():
            gen.set_qmin(qmin)
        return

    @deprecated(setQmin_dep_msg)
    def setQmin(self, qmin):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use
        diffpy.cmipdf.pdfcontribution.PDFContribution.set_qmin
        instead.
        """
        return self.set_qmin(qmin)

    def get_qmin(self):
        """Get the qmin value.

        Returns
        -------
        float
            The minimum scattering vector used to generate the PDF.
        """
        qmin = self._get_meta_value("qmin")
        return qmin

    @deprecated(getQmin_dep_msg)
    def getQmin(self):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use
        diffpy.cmipdf.pdfcontribution.PDFContribution.get_qmin
        instead.
        """
        return self.get_qmin()


# End of file
