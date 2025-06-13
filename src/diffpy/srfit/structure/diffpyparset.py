#!/usr/bin/env python
##############################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2009 The Trustees of Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################
"""Adapters for interfacing a diffpy.structure.Structure with SrFit.

A diffpy.structure.Structure object is meant to be passed to a
DiffpyStructureParSet object from this module, which can then be used as a
ParameterSet. (It has other methods for interfacing with SrReal calculator
adapters.) Any change to the lattice or existing atoms will be registered with
the Structure. Changes in the number of atoms will not be recognized.  Thus,
the diffpy.structure.Structure object should be fully configured before passing
it to DiffpyStructureParSet.

DiffpyStructureParSet --  Adapter for diffpy.structure.Structure
DiffpyLatticeParSet   --  Adapter for diffpy.structure.Lattice
DiffpyAtomParSet      --  Adapter for diffpy.structure.Atom
"""

__all__ = ["DiffpyStructureParSet"]

from diffpy.srfit.fitbase.parameter import ParameterAdapter, ParameterProxy
from diffpy.srfit.fitbase.parameterset import ParameterSet
from diffpy.srfit.structure.srrealparset import SrRealParSet
from diffpy.srfit.util.argbinders import bind2nd


# Accessor for xyz of atoms
class _xyzgetter(object):

    def __init__(self, i):
        self.i = i

    def __call__(self, atom):
        return atom.xyz[self.i]


class _xyzsetter(object):

    def __init__(self, i):
        self.i = i

    def __call__(self, atom, value):
        atom.xyz[self.i] = value


class DiffpyAtomParSet(ParameterSet):
    """A wrapper for diffpy.structure.Atom.

    This class derives from diffpy.srfit.fitbase.parameterset.ParameterSet. See
    this class for base attributes.

    Attributes:
    atom        --  The diffpy.structure.Atom this is adapting
    element     --  The element name (property).

    Managed Parameters:
    x (y, z)    --  Atom position in crystal coordinates (ParameterAdapter)
    occupancy   --  Occupancy of the atom on its crystal location
                    (ParameterAdapter)
    occ         --  Proxy for occupancy (ParameterProxy).
    U11, U22, U33, U12, U21, U23, U32, U13, U31
                --  Anisotropic displacement factor for atom (ParameterAdapter
                    or ParameterProxy). Note that the Uij and Uji parameters
                    are the same.
    Uiso        --  Isotropic ADP (ParameterAdapter).
    B11, B22, B33, B12, B21, B23, B32, B13, B31
                --  Anisotropic displacement factor for atom (ParameterAdapter
                    or ParameterProxy). Note that the Bij and Bji parameters
                    are the same. (Bij = 8*pi**2*Uij)
    Biso        --  Isotropic ADP (ParameterAdapter).
    """

    def __init__(self, name, atom):
        """Initialize.

        atom    --  A diffpy.structure.Atom instance
        """
        ParameterSet.__init__(self, name)
        self.atom = atom
        a = atom
        # x, y, z, occupancy
        self.addParameter(
            ParameterAdapter("x", a, _xyzgetter(0), _xyzsetter(0))
        )
        self.addParameter(
            ParameterAdapter("y", a, _xyzgetter(1), _xyzsetter(1))
        )
        self.addParameter(
            ParameterAdapter("z", a, _xyzgetter(2), _xyzsetter(2))
        )
        occupancy = ParameterAdapter("occupancy", a, attr="occupancy")
        self.addParameter(occupancy)
        self.addParameter(ParameterProxy("occ", occupancy))
        # U
        self.addParameter(ParameterAdapter("U11", a, attr="U11"))
        self.addParameter(ParameterAdapter("U22", a, attr="U22"))
        self.addParameter(ParameterAdapter("U33", a, attr="U33"))
        U12 = ParameterAdapter("U12", a, attr="U12")
        U21 = ParameterProxy("U21", U12)
        U13 = ParameterAdapter("U13", a, attr="U13")
        U31 = ParameterProxy("U31", U13)
        U23 = ParameterAdapter("U23", a, attr="U23")
        U32 = ParameterProxy("U32", U23)
        self.addParameter(U12)
        self.addParameter(U21)
        self.addParameter(U13)
        self.addParameter(U31)
        self.addParameter(U23)
        self.addParameter(U32)
        self.addParameter(ParameterAdapter("Uiso", a, attr="Uisoequiv"))
        # B
        self.addParameter(ParameterAdapter("B11", a, attr="B11"))
        self.addParameter(ParameterAdapter("B22", a, attr="B22"))
        self.addParameter(ParameterAdapter("B33", a, attr="B33"))
        B12 = ParameterAdapter("B12", a, attr="B12")
        B21 = ParameterProxy("B21", B12)
        B13 = ParameterAdapter("B13", a, attr="B13")
        B31 = ParameterProxy("B31", B13)
        B23 = ParameterAdapter("B23", a, attr="B23")
        B32 = ParameterProxy("B32", B23)
        self.addParameter(B12)
        self.addParameter(B21)
        self.addParameter(B13)
        self.addParameter(B31)
        self.addParameter(B23)
        self.addParameter(B32)
        self.addParameter(ParameterAdapter("Biso", a, attr="Bisoequiv"))
        return

    def __repr__(self):
        return repr(self.atom)

    def _getElem(self):
        return self.atom.element

    def _setElem(self, el):
        self.atom.element = el

    element = property(_getElem, _setElem, "type of atom")


# End class DiffpyAtomParSet


def _latgetter(par):
    return bind2nd(getattr, par)


def _latsetter(par):
    return bind2nd(setattr, par)


class DiffpyLatticeParSet(ParameterSet):
    """A wrapper for diffpy.structure.Lattice.

    This class derives from diffpy.srfit.fitbase.parameterset.ParameterSet. See
    this class for base attributes.

    Attributes
    lattice     --  The diffpy.structure.Lattice this is adapting
    name        --  Always "lattice"
    angunits    --  "deg", the units of angle

    Managed Parameters:
    a, b, c, alpha, beta, gamma --  The lattice parameters (ParameterAdapter).
    """

    def __init__(self, lattice):
        """Initialize.

        lattice --  A diffpy.structure.Lattice instance
        """
        ParameterSet.__init__(self, "lattice")
        self.angunits = "deg"
        self.lattice = lattice
        lat = lattice
        self.addParameter(
            ParameterAdapter("a", lat, _latgetter("a"), _latsetter("a"))
        )
        self.addParameter(
            ParameterAdapter("b", lat, _latgetter("b"), _latsetter("b"))
        )
        self.addParameter(
            ParameterAdapter("c", lat, _latgetter("c"), _latsetter("c"))
        )
        self.addParameter(
            ParameterAdapter(
                "alpha", lat, _latgetter("alpha"), _latsetter("alpha")
            )
        )
        self.addParameter(
            ParameterAdapter(
                "beta", lat, _latgetter("beta"), _latsetter("beta")
            )
        )
        self.addParameter(
            ParameterAdapter(
                "gamma", lat, _latgetter("gamma"), _latsetter("gamma")
            )
        )
        return

    def __repr__(self):
        return repr(self.lattice)


# End class DiffpyLatticeParSet


class DiffpyStructureParSet(SrRealParSet):
    """A wrapper for diffpy.structure.Structure.

    This class derives from diffpy.srfit.fitbase.parameterset.ParameterSet. See
    this class for base attributes.

    Attributes:
    atoms   --  The list of DiffpyAtomParSets, provided for convenience.
    stru    --  The diffpy.structure.Structure this is adapting

    Managed ParameterSets:
    lattice     --  The managed DiffpyLatticeParSet
    <el><idx>   --  A managed DiffpyAtomParSets. <el> is the atomic element and
                    <idx> is the index of that element in the structure,
                    starting from zero. Thus, for nickel in P1 symmetry, the
                    managed DiffpyAtomParSets will be named "Ni0", "Ni1", "Ni2"
                    and "Ni3".
    """

    def __init__(self, name, stru):
        """Initialize.

        name    --  A name for the structure
        stru    --  A diffpy.structure.Structure instance
        """
        SrRealParSet.__init__(self, name)
        self.stru = stru
        self.addParameterSet(DiffpyLatticeParSet(stru.lattice))
        self.atoms = []

        cdict = {}
        for a in stru:
            el = a.element.title()
            # Try to sanitize the name.
            el = el.replace("+", "p")
            el = el.replace("-", "m")
            i = cdict.get(el, 0)
            aname = "%s%i" % (el, i)
            cdict[el] = i + 1
            atom = DiffpyAtomParSet(aname, a)
            self.addParameterSet(atom)
            self.atoms.append(atom)

        return

    def __repr__(self):
        return repr(self.stru)

    def getLattice(self):
        """Get the ParameterSet containing the lattice Parameters."""
        return self.lattice

    @classmethod
    def canAdapt(self, stru):
        """Return whether the structure can be adapted by this class."""
        from diffpy.structure import Structure

        return isinstance(stru, Structure)

    def getScatterers(self):
        """Get a list of ParameterSets that represents the scatterers.

        The site positions must be accessible from the list entries via
        the names "x", "y", and "z". The ADPs must be accessible as
        well, but the name and nature of the ADPs (U-factors, B-factors,
        isotropic, anisotropic) depends on the adapted structure.
        """
        return self.atoms

    def _getSrRealStructure(self):
        """Get the structure object for use with SrReal calculators.

        If this is periodic, then return the structure, otherwise, pass
        it inside of a nosymmetry wrapper. This takes the extra step of
        wrapping the structure in a nometa wrapper.
        """
        from diffpy.srreal.structureadapter import nometa

        stru = SrRealParSet._getSrRealStructure(self)
        return nometa(stru)


# End class DiffpyStructureParSet
