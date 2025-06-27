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
"""Tests space group constraints."""

import unittest

import numpy
import pytest

# ----------------------------------------------------------------------------


def test_ObjCryst_constrainSpaceGroup(pyobjcryst_available):
    """Make sure that all Parameters are constrained properly.

    This tests constrainSpaceGroup from
    diffpy.srfit.structure.sgconstraints, which is performed
    automatically when an ObjCrystCrystalParSet is created.
    """
    if not pyobjcryst_available:
        pytest.skip("pyobjcrysta package not available")

    from diffpy.srfit.structure.objcrystparset import ObjCrystCrystalParSet

    pi = numpy.pi

    occryst = makeLaMnO3()
    stru = ObjCrystCrystalParSet(occryst.GetName(), occryst)
    # Make sure we actually create the constraints
    stru._constrainSpaceGroup()
    # Make the space group parameters individually
    stru.sgpars.latpars
    stru.sgpars.xyzpars
    stru.sgpars.adppars

    # Check the orthorhombic lattice
    lattice = stru.getLattice()
    assert lattice.alpha.const
    assert lattice.beta.const
    assert lattice.gamma.const
    assert pi / 2 == lattice.alpha.getValue()
    assert pi / 2 == lattice.beta.getValue()
    assert pi / 2 == lattice.gamma.getValue()

    assert not lattice.a.const
    assert not lattice.b.const
    assert not lattice.c.const
    assert 0 == len(lattice._constraints)

    # Now make sure the scatterers are constrained properly
    scatterers = stru.getScatterers()
    la = scatterers[0]
    assert not la.x.const
    assert not la.y.const
    assert la.z.const
    assert 0 == len(la._constraints)

    mn = scatterers[1]
    assert mn.x.const
    assert mn.y.const
    assert mn.z.const
    assert 0 == len(mn._constraints)

    o1 = scatterers[2]
    assert not o1.x.const
    assert not o1.y.const
    assert o1.z.const
    assert 0 == len(o1._constraints)

    o2 = scatterers[3]
    assert not o2.x.const
    assert not o2.y.const
    assert not o2.z.const
    assert 0 == len(o2._constraints)

    # Make sure we can't constrain these
    with pytest.raises(ValueError):
        mn.constrain(mn.x, "y")

    with pytest.raises(ValueError):
        mn.constrain(mn.y, "z")

    with pytest.raises(ValueError):
        mn.constrain(mn.z, "x")

    # Nor can we make them into variables
    from diffpy.srfit.fitbase.fitrecipe import FitRecipe

    f = FitRecipe()
    with pytest.raises(ValueError):
        f.addVar(mn.x)

    return


def test_DiffPy_constrainAsSpaceGroup(datafile, pyobjcryst_available):
    """Test the constrainAsSpaceGroup function."""
    if not pyobjcryst_available:
        pytest.skip("pyobjcrysta package not available")

    from diffpy.srfit.structure.diffpyparset import DiffpyStructureParSet
    from diffpy.srfit.structure.sgconstraints import constrainAsSpaceGroup

    stru = makeLaMnO3_P1(datafile)
    parset = DiffpyStructureParSet("LaMnO3", stru)

    sgpars = constrainAsSpaceGroup(
        parset,
        "P b n m",
        scatterers=parset.getScatterers()[::2],
        constrainadps=True,
    )

    # Make sure that the new parameters were created
    for par in sgpars:
        assert par is not None
        assert par.getValue() is not None

    # Test the unconstrained atoms
    for scatterer in parset.getScatterers()[1::2]:
        assert not scatterer.x.const
        assert not scatterer.y.const
        assert not scatterer.z.const
        assert not scatterer.U11.const
        assert not scatterer.U22.const
        assert not scatterer.U33.const
        assert not scatterer.U12.const
        assert not scatterer.U13.const
        assert not scatterer.U23.const
        assert 0 == len(scatterer._constraints)

    proxied = [p.par for p in sgpars]

    def _consttest(par):
        return par.const

    def _constrainedtest(par):
        return par.constrained

    def _proxytest(par):
        return par in proxied

    def _alltests(par):
        return _consttest(par) or _constrainedtest(par) or _proxytest(par)

    for idx, scatterer in enumerate(parset.getScatterers()[::2]):
        # Under this scheme, atom 6 is free to vary
        test = False
        for par in [scatterer.x, scatterer.y, scatterer.z]:
            test |= _alltests(par)
        assert test

        test = False
        for par in [
            scatterer.U11,
            scatterer.U22,
            scatterer.U33,
            scatterer.U12,
            scatterer.U13,
            scatterer.U23,
        ]:
            test |= _alltests(par)

        assert test

    return


def test_ConstrainAsSpaceGroup_args(pyobjcryst_available, datafile):
    """Test the arguments processing of constrainAsSpaceGroup function."""
    if not pyobjcryst_available:
        pytest.skip("pyobjcrysta package not available")

    from diffpy.srfit.structure.diffpyparset import DiffpyStructureParSet
    from diffpy.srfit.structure.sgconstraints import constrainAsSpaceGroup
    from diffpy.structure.spacegroups import GetSpaceGroup

    stru = makeLaMnO3_P1(datafile)
    parset = DiffpyStructureParSet("LaMnO3", stru)
    sgpars = constrainAsSpaceGroup(parset, "P b n m")
    sg = GetSpaceGroup("P b n m")
    parset2 = DiffpyStructureParSet("LMO", makeLaMnO3_P1(datafile))
    sgpars2 = constrainAsSpaceGroup(parset2, sg)
    list(sgpars)
    list(sgpars2)
    assert sgpars.names == sgpars2.names
    return


def makeLaMnO3_P1(datafile):
    from diffpy.structure import Structure

    stru = Structure()
    stru.read(datafile("LaMnO3.stru"))
    return stru


def makeLaMnO3():
    from pyobjcryst.atom import Atom
    from pyobjcryst.crystal import Crystal
    from pyobjcryst.scatteringpower import ScatteringPowerAtom

    pi = numpy.pi
    # It appears that ObjCryst only supports standard symbols
    crystal = Crystal(5.486341, 5.619215, 7.628206, "P b n m")
    crystal.SetName("LaMnO3")
    # La1
    sp = ScatteringPowerAtom("La1", "La")
    sp.SetBiso(8 * pi * pi * 0.003)
    atom = Atom(0.996096, 0.0321494, 0.25, "La1", sp)
    crystal.AddScatteringPower(sp)
    crystal.AddScatterer(atom)
    # Mn1
    sp = ScatteringPowerAtom("Mn1", "Mn")
    sp.SetBiso(8 * pi * pi * 0.003)
    atom = Atom(0, 0.5, 0, "Mn1", sp)
    crystal.AddScatteringPower(sp)
    crystal.AddScatterer(atom)
    # O1
    sp = ScatteringPowerAtom("O1", "O")
    sp.SetBiso(8 * pi * pi * 0.003)
    atom = Atom(0.0595746, 0.496164, 0.25, "O1", sp)
    crystal.AddScatteringPower(sp)
    crystal.AddScatterer(atom)
    # O2
    sp = ScatteringPowerAtom("O2", "O")
    sp.SetBiso(8 * pi * pi * 0.003)
    atom = Atom(0.720052, 0.289387, 0.0311126, "O2", sp)
    crystal.AddScatteringPower(sp)
    crystal.AddScatterer(atom)

    return crystal


# ----------------------------------------------------------------------------

if __name__ == "__main__":
    unittest.main()
