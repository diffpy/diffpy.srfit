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
"""PDF profile generator base class.

The BasePDFGenerator class interfaces with SrReal PDF calculators and is
used as a base for the PDFGenerator and DebyePDFGenerator classes.
"""

__all__ = ["BasePDFGenerator"]

import numpy

from diffpy.srfit.exceptions import SrFitError
from diffpy.srfit.fitbase import ProfileGenerator
from diffpy.srfit.fitbase.parameter import ParameterAdapter
from diffpy.srfit.structure import struToParameterSet
from diffpy.utils._deprecator import build_deprecation_message, deprecated

removal_version = "4.0.0"
generator_base = "diffpy.srfit.pdf.basepdfgenerator.BasePDFGenerator"
generator_new_base = "diffpy.cmipdf.basepdfgenerator.BasePDFGenerator"


def _build_generator_dep_msg(old_name, new_name):
    message = build_deprecation_message(
        generator_base,
        old_name,
        new_name,
        removal_version,
        new_base=generator_new_base,
    )
    return message


setScatteringType_dep_msg = _build_generator_dep_msg(
    "setScatteringType", "set_scattering_type"
)
getScatteringType_dep_msg = _build_generator_dep_msg(
    "getScatteringType", "get_scattering_type"
)
setQmax_dep_msg = _build_generator_dep_msg("setQmax", "set_qmax")
getQmax_dep_msg = _build_generator_dep_msg("getQmax", "get_qmax")
setQmin_dep_msg = _build_generator_dep_msg("setQmin", "set_qmin")
getQmin_dep_msg = _build_generator_dep_msg("getQmin", "get_qmin")
setStructure_dep_msg = (
    _build_generator_dep_msg("setStructure", "set_structure")
    + " Additionally, the signature has changed. Please pass the parameter "
    "'stru' as 'structure'."
)
setPhase_dep_msg = _build_generator_dep_msg(
    "setPhase", "set_structure_from_parset"
)

# FIXME - Parameter creation will have to be smarter once deeper calculator
# configuration is enabled.


class BasePDFGenerator(ProfileGenerator):
    """Base class for calculating PDF profiles using SrReal.

    This works with diffpy.structure.Structure, pyobjcryst.crystal.Crystal and
    pyobjcryst.molecule.Molecule instances. Note that the managed Parameters
    are not created until the structure is added.

    Attributes
    ----------
    _calc
        PDFCalculator or DebyePDFCalculator instance for calculating
        the PDF.
    _phase
        The structure ParameterSet used to calculate the profile.
    stru
        The structure objected adapted by _phase.
    _lastr
        The last value of r over which the PDF was calculated. This is
        used to configure the calculator when r changes.
    _pool
        A multiprocessing.Pool for managing parallel computation.

    Managed Parameters
    ------------------
    scale
        Scale factor
    delta1
        Linear peak broadening term
    delta2
        Quadratic peak broadening term
    qbroad
        Resolution peak broadening term
    qdamp
        Resolution peak dampening term
    Managed ParameterSets
        The structure ParameterSet (SrRealParSet instance) used to
        calculate the profile is named by the user.

    Usable Metadata
    ---------------
    stype
        The scattering type "X" for x-ray, "N" for neutron (see
        'set_scattering_type').
    qmax
        The maximum scattering vector used to generate the PDF (see
        set_qmax).
    qmin
        The minimum scattering vector used to generate the PDF (see
        set_qmin).
    scale
        See Managed Parameters.
    delta1
        See Managed Parameters.
    delta2
        See Managed Parameters.
    qbroad
        See Managed Parameters.
    qdamp
        See Managed Parameters.
    """

    def __init__(self, name="pdf"):
        """Initialize the generator.

        Parameters
        ----------
        name : str, optional
            The name for this generator (default "pdf").
        """
        ProfileGenerator.__init__(self, name)

        self._phase = None
        self.stru = None
        self.meta = {}
        self._lastr = numpy.empty(0)
        self._calc = None

        self._pool = None

        return

    _parnames = ["delta1", "delta2", "qbroad", "scale", "qdamp"]

    def _set_calculator(self, calc):
        """Set the SrReal calculator instance.

        Setting the calculator creates Parameters from the variable
        attributes of the SrReal calculator.
        """
        self._calc = calc
        for pname in self.__class__._parnames:
            self.add_parameter(ParameterAdapter(pname, self._calc, attr=pname))
        self._process_metadata()
        return

    def parallel(self, ncpu, mapfunc=None):
        """Run calculation in parallel.

        Parameters
        ----------
        ncpu : int
            The number of parallel processes. Revert to serial mode
            when 1.
        mapfunc : callable, optional
            The mapping function to use. If this is None (default),
            ``multiprocessing.Pool.imap_unordered`` will be used.
        """
        from diffpy.srreal.parallel import createParallelCalculator

        calc_serial = self._calc
        if hasattr(calc_serial, "pqobj"):
            calc_serial = calc_serial.pqobj
        # revert to serial calculator for ncpu <= 1
        if ncpu <= 1:
            self._calc = calc_serial
            self._pool = None
            return
        # Why don't we let the user shoot his foot or test on single CPU?
        # ncpu = min(ncpu, multiprocessing.cpu_count())
        if mapfunc is None:
            import multiprocessing

            self._pool = multiprocessing.Pool(ncpu)
            mapfunc = self._pool.imap_unordered

        self._calc = createParallelCalculator(calc_serial, ncpu, mapfunc)
        return

    def _process_metadata(self):
        """Process the metadata once it gets set."""
        ProfileGenerator._process_metadata(self)

        stype = self.meta.get("stype")
        if stype is not None:
            self.set_scattering_type(stype)

        qmax = self.meta.get("qmax")
        if qmax is not None:
            self.set_qmax(qmax)

        qmin = self.meta.get("qmin")
        if qmin is not None:
            self.set_qmin(qmin)

        for name in self.__class__._parnames:
            val = self.meta.get(name)
            if val is not None:
                par = self.get(name)
                par.set_value(val)

        return

    def set_scattering_type(self, stype="X"):
        """Set the scattering type.

        Parameters
        ----------
        stype : str, optional
            "X" for x-ray, "N" for neutron, "E" for electrons,
            or any registered type from diffpy.srreal from
            ScatteringFactorTable.getRegisteredTypes() (default "X").

        Raises
        ------
        ValueError
            If stype is not a recognized scattering type.
        """
        self._calc.scatteringfactortable = stype
        # update the meta dictionary only if there was no exception
        self.meta["stype"] = self.get_scattering_type()
        return

    @deprecated(setScatteringType_dep_msg)
    def setScatteringType(self, stype="X"):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use
        diffpy.cmipdf.basepdfgenerator.BasePDFGenerator.set_scattering_type
        instead.
        """
        return self.set_scattering_type(stype)

    def get_scattering_type(self):
        """Get the scattering type.

        See ``set_scattering_type``.

        Returns
        -------
        str
            The scattering type used to calculate the PDF.
        """
        scattering_type = self._calc.getRadiationType()
        return scattering_type

    @deprecated(getScatteringType_dep_msg)
    def getScatteringType(self):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use
        diffpy.cmipdf.basepdfgenerator.BasePDFGenerator.get_scattering_type
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
        self._calc.qmax = qmax
        self.meta["qmax"] = self.get_qmax()
        return

    @deprecated(setQmax_dep_msg)
    def setQmax(self, qmax):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use
        diffpy.cmipdf.basepdfgenerator.BasePDFGenerator.set_qmax
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
        qmax = self._calc.qmax
        return qmax

    @deprecated(getQmax_dep_msg)
    def getQmax(self):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use
        diffpy.cmipdf.basepdfgenerator.BasePDFGenerator.get_qmax
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
        self._calc.qmin = qmin
        self.meta["qmin"] = self.get_qmin()
        return

    @deprecated(setQmin_dep_msg)
    def setQmin(self, qmin):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use
        diffpy.cmipdf.basepdfgenerator.BasePDFGenerator.set_qmin
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
        qmin = self._calc.qmin
        return qmin

    @deprecated(getQmin_dep_msg)
    def getQmin(self):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use
        diffpy.cmipdf.basepdfgenerator.BasePDFGenerator.get_qmin
        instead.
        """
        return self.get_qmin()

    def set_structure(self, structure, name="phase", periodic=True):
        """Set the structure that will be used to calculate the PDF.

        This creates a DiffpyStructureParSet, ObjCrystCrystalParSet or
        ObjCrystMoleculeParSet that adapts the structure to a ParameterSet
        interface. See those classes (located in diffpy.srfit.structure) for
        how they are used. The resulting ParameterSet will be managed by this
        generator.

        Parameters
        ----------
        structure
            The diffpy.structure.Structure, pyobjcryst.crystal.Crystal or
            pyobjcryst.molecule.Molecule instance to adapt.
        name : str, optional
            The name to give to the managed ParameterSet that adapts the
            structure (default "phase").
        periodic : bool, optional
            The structure should be treated as periodic (default
            True). Note that some structures do not support
            periodicity, in which case this will have no effect on the
            PDF calculation.
        """
        parset = struToParameterSet(name, structure)
        self.set_structure_from_parset(parset, periodic)
        return

    @deprecated(setStructure_dep_msg)
    def setStructure(self, stru, name="phase", periodic=True):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use
        diffpy.cmipdf.basepdfgenerator.BasePDFGenerator.set_structure
        instead.
        """
        return self.set_structure(stru, name, periodic)

    def set_structure_from_parset(self, parset, periodic=True):
        """Set the phase that will be used to calculate the PDF.

        Set the phase directly with a DiffpyStructureParSet,
        ObjCrystCrystalParSet or ObjCrystMoleculeParSet that adapts a structure
        object (from diffpy or pyobjcryst).  The passed ParameterSet will be
        managed by this generator.

        Parameters
        ----------
        parset
            The SrRealParSet that holds the structural information.
            This can be used to share the phase between multiple
            BasePDFGenerators, and have the changes in one reflect in
            another.
        periodic : bool, optional
            The structure should be treated as periodic (default True).
            Note that some structures do not support periodicity, in
            which case this will be ignored.
        """
        # Store the ParameterSet for easy access
        self._phase = parset
        self.stru = self._phase.stru

        # Put this ParameterSet in the ProfileGenerator.
        self.add_parameter_set(parset)

        # Set periodicity
        self._phase.useSymmetry(periodic)
        return

    @deprecated(setPhase_dep_msg)
    def setPhase(self, parset, periodic=True):
        """This function has been deprecated and will be removed in version
        4.0.0.

        Please use
        diffpy.cmipdf.basepdfgenerator.BasePDFGenerator.set_structure_from_parset
        instead.
        """
        return self.set_structure_from_parset(parset, periodic)

    def _prepare(self, r):
        """Prepare the calculator when a new r-value is passed."""
        self._lastr = r.copy()
        lo, hi = r.min(), r.max()
        ndiv = max(len(r) - 1, 1)
        self._calc.rstep = (hi - lo) / ndiv
        self._calc.rmin = lo
        self._calc.rmax = hi + 0.5 * self._calc.rstep
        return

    def _validate(self):
        """Validate my state.

        This validates that the phase is not None. This performs
        ProfileGenerator validations.

        Raises
        ------
        SrFitError
            If validation fails.
        """
        if self._calc is None:
            raise SrFitError("_calc is None")
        if self._phase is None:
            raise SrFitError("_phase is None")
        ProfileGenerator._validate(self)
        return

    def __call__(self, r):
        """Calculate the PDF.

        This ProfileGenerator will be used in a fit equation that will
        be optimized to fit some data.  By the time this function is
        evaluated, the crystal has been updated by the optimizer via the
        ObjCrystParSet created in setCrystal. Thus, we need only call
        pdf with the internal structure object.

        Parameters
        ----------
        r : np.ndarray
            The independent variable over which to calculate the PDF.

        Returns
        -------
        np.ndarray
            The calculated PDF, evaluated over r.
        """
        if not numpy.array_equal(r, self._lastr):
            self._prepare(r)

        rcalc, y = self._calc(self._phase._get_srreal_structure())

        if numpy.isnan(y).any():
            y = numpy.zeros_like(r)
        else:
            y = numpy.interp(r, rcalc, y)
        return y


# End class BasePDFGenerator
