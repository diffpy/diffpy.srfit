#!/usr/bin/env python
"""Wrappers for interfacing a diffpy.Structure.Structure as a ParameterSet
with the same hierarchy.

Structure   --  Name required. Contains a Lattice ParameterSet and several
                Atom parameter sets.
                Other Attributes:
Lattice     --  Named "lattice". Contains Parameters "a", "b", "c", "alpha",
                "beta", "gamma".
                Other Attributes:
Atom        --  Named "%s%i" % (element, number). Contains Parameters "x", "y",
                "z", "occupancy", "B11", "B22", "B33", "B12", "B23", "B13",
                "B11", "B22", "B33", "B12", "B23", "B13". The asymmetric
                parameters also have proxies with inverted indices.
                Other Attributes:
                element

"""
__id__ = "$Id$"

from diffpy.srfit.fitbase.parameter import Parameter, ParameterProxy
from diffpy.srfit.fitbase.parameterset import ParameterSet

from .wrappers import ParameterWrapper

# Accessor for xyz of atoms
def _getter(i):
    def _f(atom):
        return atom.xyz[i]
    return _f

def _setter(i):
    def _f(atom, value):
        atom.xyz[i] = value
        return
    return _f

class Atom(ParameterSet):
    """A wrapper for diffpy.Structure.Atom."""

    def __init__(self, atom, name):
        """Initialize

        atom    --  A diffpy.Structure.Atom instance
        """
        ParameterSet.__init__(self, name)
        self.atom = atom
        a = atom
        # x, y, z, occupancy
        self.addParameter(ParameterWrapper(a, "x", _getter(0), _setter(0)))
        self.addParameter(ParameterWrapper(a, "y", _getter(1), _setter(1)))
        self.addParameter(ParameterWrapper(a, "z", _getter(2), _setter(2)))
        self.addParameter(ParameterWrapper(a, "occupancy", attr = "occupancy"))
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
        self.addParameter(ParameterWrapper(a, "B11", attr = "B11"))
        self.addParameter(ParameterWrapper(a, "B22", attr = "B22"))
        self.addParameter(ParameterWrapper(a, "B33", attr = "B33"))
        self.addParameter(U12)
        self.addParameter(U21)
        self.addParameter(U13)
        self.addParameter(U31)
        self.addParameter(U23)
        self.addParameter(U32)
        # B
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

        # Other setup
        self.__repr__ = a.__repr__
        return

    def _getElem(self):
        return self.atom.element

    def _setElem(self, el):
        self.atom.element = el

    element = property(self._getElem, self._setElem, "type of atom")

# End class Atom

class Lattice(ParameterSet):
    """A wrapper for diffpy.Structure.Lattice."""

    def __init__(self, lattice):
        """Initialize

        lattice --  A diffpy.Structure.Lattice instance
        """
        ParameterSet.__init__(self, "lattice")
        self.lattice = lattice
        l = lattice
        self.addParameter(ParameterWrapper(l, "a", attr = "a"))
        self.addParameter(ParameterWrapper(l, "b", attr = "b"))
        self.addParameter(ParameterWrapper(l, "c", attr = "c"))
        self.addParameter(ParameterWrapper(l, "alpha", attr = "alpha"))
        self.addParameter(ParameterWrapper(l, "beta", attr = "beta"))
        self.addParameter(ParameterWrapper(l, "gamma", attr = "gamma"))

        # Other setup
        self.__repr__ = l.__repr__
        return

# End class Lattice

class Structure(ParameterSet):
    """A wrapper for diffpy.Structure.Structure."""

    def __init__(self, stru, name):
        """Initialize

        stru    --  A diffpy.Structure.Lattice instance
        """
        ParameterSet.__init__(self, name)
        self.stru = stru
        self.addParameterSet(Lattice(stru.lattice))

        cdict = {}
        for a in stru:
            el = a.element
            i = cdict.get(el, 0)
            aname = "%s%i"%(el,i)
            cdict[el] = i+1
            self.addParameterSet(Atom(a, aname))

        # other setup
        self.__stru__ = stru.__stru__
        return

# End class Structure

