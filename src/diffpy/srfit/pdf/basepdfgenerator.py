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

# FIXME - Parameter creation will have to be smarter once deeper calculator
# configuration is enabled.


class BasePDFGenerator(ProfileGenerator):
    """Base class for calculating PDF profiles using SrReal.

    This works with diffpy.structure.Structure, pyobjcryst.crystal.Crystal and
    pyobjcryst.molecule.Molecule instances. Note that the managed Parameters
    are not created until the structure is added.

    Attributes:
    _calc   --  PDFCalculator or DebyePDFCalculator instance for calculating
                the PDF.
    _phase  --  The structure ParameterSet used to calculate the profile.
    stru    --  The structure objected adapted by _phase.
    _lastr  --  The last value of r over which the PDF was calculated. This is
                used to configure the calculator when r changes.
    _pool   --  A multiprocessing.Pool for managing parallel computation.

    Managed Parameters:
    scale   --  Scale factor
    delta1  --  Linear peak broadening term
    delta2  --  Quadratic peak broadening term
    qbroad  --  Resolution peak broadening term
    qdamp   --  Resolution peak dampening term

    Managed ParameterSets:
    The structure ParameterSet (SrRealParSet instance) used to calculate the
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
        ProfileGenerator.__init__(self, name)

        self._phase = None
        self.stru = None
        self.meta = {}
        self._lastr = numpy.empty(0)
        self._calc = None

        self._pool = None

        return

    _parnames = ["delta1", "delta2", "qbroad", "scale", "qdamp"]

    def _setCalculator(self, calc):
        """Set the SrReal calculator instance.

        Setting the calculator creates Parameters from the variable
        attributes of the SrReal calculator.
        """
        self._calc = calc
        for pname in self.__class__._parnames:
            self.addParameter(ParameterAdapter(pname, self._calc, attr=pname))
        self.processMetaData()
        return

    def parallel(self, ncpu, mapfunc=None):
        """Run calculation in parallel.

        ncpu    -- Number of parallel processes.  Revert to serial mode when 1.
        mapfunc -- A mapping function to use. If this is None (default),
                   multiprocessing.Pool.imap_unordered will be used.

        No return value.
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

    def processMetaData(self):
        """Process the metadata once it gets set."""
        ProfileGenerator.processMetaData(self)

        stype = self.meta.get("stype")
        if stype is not None:
            self.setScatteringType(stype)

        qmax = self.meta.get("qmax")
        if qmax is not None:
            self.setQmax(qmax)

        qmin = self.meta.get("qmin")
        if qmin is not None:
            self.setQmin(qmin)

        for name in self.__class__._parnames:
            val = self.meta.get(name)
            if val is not None:
                par = self.get(name)
                par.setValue(val)

        return

    def setScatteringType(self, stype="X"):
        """Set the scattering type.

        stype   --   "X" for x-ray, "N" for neutron, "E" for electrons,
                     or any registered type from diffpy.srreal from
                     ScatteringFactorTable.getRegisteredTypes().

        Raises ValueError for unknown scattering type.
        """
        self._calc.setScatteringFactorTableByType(stype)
        # update the meta dictionary only if there was no exception
        self.meta["stype"] = self.getScatteringType()
        return

    def getScatteringType(self):
        """Get the scattering type.

        See 'setScatteringType'.
        """
        return self._calc.getRadiationType()

    def setQmax(self, qmax):
        """Set the qmax value."""
        self._calc.qmax = qmax
        self.meta["qmax"] = self.getQmax()
        return

    def getQmax(self):
        """Get the qmax value."""
        return self._calc.qmax

    def setQmin(self, qmin):
        """Set the qmin value."""
        self._calc.qmin = qmin
        self.meta["qmin"] = self.getQmin()
        return

    def getQmin(self):
        """Get the qmin value."""
        return self._calc.qmin

    def setStructure(self, stru, name="phase", periodic=True):
        """Set the structure that will be used to calculate the PDF.

        This creates a DiffpyStructureParSet, ObjCrystCrystalParSet or
        ObjCrystMoleculeParSet that adapts stru to a ParameterSet interface.
        See those classes (located in diffpy.srfit.structure) for how they are
        used. The resulting ParameterSet will be managed by this generator.

        stru    --  diffpy.structure.Structure, pyobjcryst.crystal.Crystal or
                    pyobjcryst.molecule.Molecule instance.  Default None.
        name    --  A name to give to the managed ParameterSet that adapts stru
                    (default "phase").
        periodic -- The structure should be treated as periodic (default
                    True). Note that some structures do not support
                    periodicity, in which case this will have no effect on the
                    PDF calculation.
        """

        # Create the ParameterSet
        parset = struToParameterSet(name, stru)

        # Set the phase
        self.setPhase(parset, periodic)
        return

    def setPhase(self, parset, periodic=True):
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
        # Store the ParameterSet for easy access
        self._phase = parset
        self.stru = self._phase.stru

        # Put this ParameterSet in the ProfileGenerator.
        self.addParameterSet(parset)

        # Set periodicity
        self._phase.useSymmetry(periodic)
        return

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

        Raises SrFitError if validation fails.
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
        """
        if not numpy.array_equal(r, self._lastr):
            self._prepare(r)

        rcalc, y = self._calc(self._phase._getSrRealStructure())

        if numpy.isnan(y).any():
            y = numpy.zeros_like(r)
        else:
            y = numpy.interp(r, rcalc, y)
        return y


# End class BasePDFGenerator
