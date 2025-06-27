#!/usr/bin/env python
##############################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2010 The Trustees of Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################
"""Tests for diffpy.srfit.structure package."""

import pickle
import unittest

import numpy as np
import pytest

from diffpy.srfit.structure.diffpyparset import DiffpyStructureParSet


def testDiffpyStructureParSet(diffpy_structure_available):
    """Test the structure conversion."""
    if not diffpy_structure_available:
        pytest.skip("diffpy.structure package not available")
    from diffpy.structure import Atom, Lattice, Structure

    a1 = Atom("Cu", xyz=np.array([0.0, 0.1, 0.2]), Uisoequiv=0.003)
    a2 = Atom("Ag", xyz=np.array([0.3, 0.4, 0.5]), Uisoequiv=0.002)
    lattice = Lattice(2.5, 2.5, 2.5, 90, 90, 90)

    dsstru = Structure([a1, a2], lattice)
    # Structure makes copies
    a1 = dsstru[0]
    a2 = dsstru[1]

    s = DiffpyStructureParSet("CuAg", dsstru)

    assert s.name == "CuAg"

    def _testAtoms():
        # Check the atoms thoroughly
        assert a1.element == s.Cu0.element
        assert a2.element == s.Ag0.element
        assert a1.Uisoequiv == s.Cu0.Uiso.getValue()
        assert a2.Uisoequiv == s.Ag0.Uiso.getValue()
        assert a1.Bisoequiv == s.Cu0.Biso.getValue()
        assert a2.Bisoequiv == s.Ag0.Biso.getValue()
        for i in range(1, 4):
            for j in range(i, 4):
                uijstru = getattr(a1, "U%i%i" % (i, j))
                uij = getattr(s.Cu0, "U%i%i" % (i, j)).getValue()
                uji = getattr(s.Cu0, "U%i%i" % (j, i)).getValue()
                assert uijstru == uij
                assert uijstru == uji
                bijstru = getattr(a1, "B%i%i" % (i, j))
                bij = getattr(s.Cu0, "B%i%i" % (i, j)).getValue()
                bji = getattr(s.Cu0, "B%i%i" % (j, i)).getValue()
                assert bijstru == bij
                assert bijstru == bji

        assert a1.xyz[0] == s.Cu0.x.getValue()
        assert a1.xyz[1] == s.Cu0.y.getValue()
        assert a1.xyz[2] == s.Cu0.z.getValue()
        return

    def _testLattice():

        # Test the lattice
        assert dsstru.lattice.a == s.lattice.a.getValue()
        assert dsstru.lattice.b == s.lattice.b.getValue()
        assert dsstru.lattice.c == s.lattice.c.getValue()
        assert dsstru.lattice.alpha == s.lattice.alpha.getValue()
        assert dsstru.lattice.beta == s.lattice.beta.getValue()
        assert dsstru.lattice.gamma == s.lattice.gamma.getValue()

    _testAtoms()
    _testLattice()

    # Now change some values from the diffpy Structure
    a1.xyz[1] = 0.123
    a1.U11 = 0.321
    a1.B32 = 0.111
    dsstru.lattice.setLatPar(a=3.0, gamma=121)
    _testAtoms()
    _testLattice()

    # Now change values from the srfit DiffpyStructureParSet
    s.Cu0.x.setValue(0.456)
    s.Cu0.U22.setValue(0.441)
    s.Cu0.B13.setValue(0.550)
    d = dsstru.lattice.dist(a1.xyz, a2.xyz)
    s.lattice.b.setValue(4.6)
    s.lattice.alpha.setValue(91.3)
    _testAtoms()
    _testLattice()
    # Make sure the distance changed
    assert d != dsstru.lattice.dist(a1.xyz, a2.xyz)
    return


def test___repr__(diffpy_structure_available):
    """Test representation of DiffpyStructureParSet objects."""
    if not diffpy_structure_available:
        pytest.skip("diffpy.structure package not available")
    from diffpy.structure import Atom, Lattice, Structure

    lat = Lattice(3, 3, 2, 90, 90, 90)
    atom = Atom("C", [0, 0.2, 0.5])
    stru = Structure([atom], lattice=lat)
    dsps = DiffpyStructureParSet("dsps", stru)
    assert repr(stru) == repr(dsps)
    assert repr(lat) == repr(dsps.lattice)
    assert repr(atom) == repr(dsps.atoms[0])
    return


def test_pickling(diffpy_structure_available):
    """Test pickling of DiffpyStructureParSet."""
    if not diffpy_structure_available:
        pytest.skip("diffpy.structure package not available")
    from diffpy.structure import Atom, Structure

    stru = Structure([Atom("C", [0, 0.2, 0.5])])
    dsps = DiffpyStructureParSet("dsps", stru)
    data = pickle.dumps(dsps)
    dsps2 = pickle.loads(data)
    assert 1 == len(dsps2.atoms)
    assert 0.2 == dsps2.atoms[0].y.value
    return


# End of class TestParameterAdapter

if __name__ == "__main__":
    unittest.main()
