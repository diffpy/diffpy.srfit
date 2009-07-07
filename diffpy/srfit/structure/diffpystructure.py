#!/usr/bin/env python
"""Wrappers for interfacing a diffpy.Structure.Structure with SrFit.

A diffpy.Structure.Structure object is meant to be passed to a StrucureParSet
object from this module, which can then be used as a ParameterSet. Any change
to the lattice or existing atoms will be registered with the Structure. Changes
in the number of atoms will not be recognized. Thus, the
diffpy.Structure.Structure object should be fully configured before passing it
to Structure.

StructureParSet --  Wrapper for diffpy.Structure.Structure
LatticeParSet   --  Wrapper for diffpy.Structure.Lattice
AtomParSet      --  Wrapper for diffpy.Structure.Atom

"""
__id__ = "$Id$"

from diffpy.srfit.fitbase.parameter import Parameter, ParameterProxy
from diffpy.srfit.fitbase.parameter import ParameterWrapper
from diffpy.srfit.fitbase.parameterset import ParameterSet


# Accessor for xyz of atoms
def _xyzgetter(i):

    def f(atom):
        return atom.xyz[i]

    return f

def _xyzsetter(i):

    def f(atom, value):
        atom.xyz[i] = value
        return

    return f


class AtomParSet(ParameterSet):
    """A wrapper for diffpy.Structure.Atom.

    This class derives from ParameterSet.

    Attributes:
    x (y, z)    --  Atom position in crystal coordinates (ParameterWrapper)
    occupancy   --  Occupancy of the atom on its crystal location
                    (ParameterWrapper)
    occ         --  Proxy for occupancy (ParameterProxy).
    U11, U22, U33, U12, U21, U23, U32, U13, U31
                --  Anisotropic displacement factor for atom (ParameterWrapper
                    or ParameterProxy). Note that the Uij and Uji parameters
                    are the same.
    Uiso        --  Isotropic ADP (ParameterWrapper).
    B11, B22, B33, B12, B21, B23, B32, B13, B31
                --  Anisotropic displacement factor for atom (ParameterWrapper
                    or ParameterProxy). Note that the Bij and Bji parameters
                    are the same. (Bij = 8*pi**2*Uij)
    Biso        --  Isotropic ADP (ParameterWrapper).
    element     --  The element name.
    
    """

    def __init__(self, atom, name):
        """Initialize

        atom    --  A diffpy.Structure.Atom instance

        """
        ParameterSet.__init__(self, name)
        self.atom = atom
        a = atom
        # x, y, z, occupancy
        self.addParameter(ParameterWrapper(a, "x", _xyzgetter(0),
            _xyzsetter(0)))
        self.addParameter(ParameterWrapper(a, "y", _xyzgetter(1),
            _xyzsetter(1)))
        self.addParameter(ParameterWrapper(a, "z", _xyzgetter(2),
            _xyzsetter(2)))
        occupancy = ParameterWrapper(a, "occupancy", attr = "occupancy")
        self.addParameter(occupancy)
        self.addParameter(ParameterProxy("occ", occupancy))
        # U
        self.addParameter(ParameterWrapper(a, "U11", attr = "U11"))
        self.addParameter(ParameterWrapper(a, "U22", attr = "U22"))
        self.addParameter(ParameterWrapper(a, "U33", attr = "U33"))
        U12 = ParameterWrapper(a, "U12", attr = "U12")
        U21 = ParameterProxy("U21", U12)
        U13 = ParameterWrapper(a, "U13", attr = "U13")
        U31 = ParameterProxy("U31", U13)
        U23 = ParameterWrapper(a, "U23", attr = "U23")
        U32 = ParameterProxy("U32", U23)
        self.addParameter(U12)
        self.addParameter(U21)
        self.addParameter(U13)
        self.addParameter(U31)
        self.addParameter(U23)
        self.addParameter(U32)
        self.addParameter(ParameterWrapper(a, "Uiso", attr = "Uisoequiv"))
        # B
        self.addParameter(ParameterWrapper(a, "B11", attr = "B11"))
        self.addParameter(ParameterWrapper(a, "B22", attr = "B22"))
        self.addParameter(ParameterWrapper(a, "B33", attr = "B33"))
        B12 = ParameterWrapper(a, "B12", attr = "B12")
        B21 = ParameterProxy("B21", B12)
        B13 = ParameterWrapper(a, "B13", attr = "B13")
        B31 = ParameterProxy("B31", B13)
        B23 = ParameterWrapper(a, "B23", attr = "B23")
        B32 = ParameterProxy("B32", B23)
        self.addParameter(B12)
        self.addParameter(B21)
        self.addParameter(B13)
        self.addParameter(B31)
        self.addParameter(B23)
        self.addParameter(B32)
        self.addParameter(ParameterWrapper(a, "Biso", attr = "Bisoequiv"))

        # Other setup
        self.__repr__ = a.__repr__
        return

    def _getElem(self):
        return self.atom.element

    def _setElem(self, el):
        self.atom.element = el

    element = property(_getElem, _setElem, "type of atom")

# End class AtomParSet


def _latgetter(par):

    def f(lat):
        return getattr(lat, par)

    return f

def _latsetter(par):

    def f(lat, value):
        setattr(lat, par, value)
        lat.setLatPar()
        return

    return f


class LatticeParSet(ParameterSet):
    """A wrapper for diffpy.Structure.Lattice.

    Attributes
    name    --  Always "lattice"
    a, b, c, alpha, beta, gamma --  The lattice parameters (ParameterWrapper).
    
    """

    def __init__(self, lattice):
        """Initialize

        lattice --  A diffpy.Structure.Lattice instance

        """
        ParameterSet.__init__(self, "lattice")
        self.lattice = lattice
        l = lattice
        self.addParameter(ParameterWrapper(l, "a", _latgetter("a"),
            _latsetter("a")))
        self.addParameter(ParameterWrapper(l, "b", _latgetter("b"),
            _latsetter("b")))
        self.addParameter(ParameterWrapper(l, "c", _latgetter("c"),
            _latsetter("c")))
        self.addParameter(ParameterWrapper(l, "alpha", _latgetter("alpha"),
            _latsetter("alpha")))
        self.addParameter(ParameterWrapper(l, "beta", _latgetter("beta"),
            _latsetter("beta")))
        self.addParameter(ParameterWrapper(l, "gamma", _latgetter("gamma"),
            _latsetter("gamma")))

        # Other setup
        self.__repr__ = l.__repr__
        return

# End class LatticeParSet

class StructureParSet(ParameterSet):
    """A wrapper for diffpy.Structure.Structure.

    atoms   --  The list of atoms.
    
    """

    def __init__(self, stru, name):
        """Initialize

        stru    --  A diffpy.Structure.Lattice instance

        """
        ParameterSet.__init__(self, name)
        self.stru = stru
        self.addParameterSet(LatticeParSet(stru.lattice))
        self.atoms = []

        cdict = {}
        for a in stru:
            el = a.element
            i = cdict.get(el, 0)
            aname = "%s%i"%(el,i)
            cdict[el] = i+1
            atom = AtomParSet(a, aname)
            self.addParameterSet(atom)
            self.atoms.append(atom)


        # other setup
        self.__repr__ = stru.__repr__
        return

# End class StructureParSet

