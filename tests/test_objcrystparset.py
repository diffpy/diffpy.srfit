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

import unittest

import numpy
import pytest

# Global variables to be assigned in setUp
ObjCrystCrystalParSet = spacegroups = None
Crystal = Atom = Molecule = ScatteringPowerAtom = None


c60xyz = """\
3.451266498   0.685000000   0.000000000
3.451266498  -0.685000000   0.000000000
-3.451266498   0.685000000   0.000000000
-3.451266498  -0.685000000   0.000000000
0.685000000   0.000000000   3.451266498
-0.685000000   0.000000000   3.451266498
0.685000000   0.000000000  -3.451266498
-0.685000000   0.000000000  -3.451266498
0.000000000   3.451266498   0.685000000
0.000000000   3.451266498  -0.685000000
0.000000000  -3.451266498   0.685000000
0.000000000  -3.451266498  -0.685000000
3.003809890   1.409000000   1.171456608
3.003809890   1.409000000  -1.171456608
3.003809890  -1.409000000   1.171456608
3.003809890  -1.409000000  -1.171456608
-3.003809890   1.409000000   1.171456608
-3.003809890   1.409000000  -1.171456608
-3.003809890  -1.409000000   1.171456608
-3.003809890  -1.409000000  -1.171456608
1.409000000   1.171456608   3.003809890
1.409000000  -1.171456608   3.003809890
-1.409000000   1.171456608   3.003809890
-1.409000000  -1.171456608   3.003809890
1.409000000   1.171456608  -3.003809890
1.409000000  -1.171456608  -3.003809890
-1.409000000   1.171456608  -3.003809890
-1.409000000  -1.171456608  -3.003809890
1.171456608   3.003809890   1.409000000
-1.171456608   3.003809890   1.409000000
1.171456608   3.003809890  -1.409000000
-1.171456608   3.003809890  -1.409000000
1.171456608  -3.003809890   1.409000000
-1.171456608  -3.003809890   1.409000000
1.171456608  -3.003809890  -1.409000000
-1.171456608  -3.003809890  -1.409000000
2.580456608   0.724000000   2.279809890
2.580456608   0.724000000  -2.279809890
2.580456608  -0.724000000   2.279809890
2.580456608  -0.724000000  -2.279809890
-2.580456608   0.724000000   2.279809890
-2.580456608   0.724000000  -2.279809890
-2.580456608  -0.724000000   2.279809890
-2.580456608  -0.724000000  -2.279809890
0.724000000   2.279809890   2.580456608
0.724000000  -2.279809890   2.580456608
-0.724000000   2.279809890   2.580456608
-0.724000000  -2.279809890   2.580456608
0.724000000   2.279809890  -2.580456608
0.724000000  -2.279809890  -2.580456608
-0.724000000   2.279809890  -2.580456608
-0.724000000  -2.279809890  -2.580456608
2.279809890   2.580456608   0.724000000
-2.279809890   2.580456608   0.724000000
2.279809890   2.580456608  -0.724000000
-2.279809890   2.580456608  -0.724000000
2.279809890  -2.580456608   0.724000000
-2.279809890  -2.580456608   0.724000000
2.279809890  -2.580456608  -0.724000000
-2.279809890  -2.580456608  -0.724000000
"""


def makeC60():
    """Make a crystal containing the C60 molecule using pyobjcryst."""
    pi = numpy.pi
    c = Crystal(100, 100, 100, "P1")
    c.SetName("c60frame")
    m = Molecule(c, "c60")

    c.AddScatterer(m)

    sp = ScatteringPowerAtom("C", "C")
    sp.SetBiso(8 * pi * pi * 0.003)
    # c.AddScatteringPower(sp)

    for i, l in enumerate(c60xyz.strip().splitlines()):
        x, y, z = map(float, l.split())
        m.AddAtom(x, y, z, sp, "C%i" % i)

    return c


# ----------------------------------------------------------------------------


class TestParameterAdapter:
    @pytest.fixture(autouse=True)
    def setup(self, pyobjcryst_available):
        # shared setup
        if not pyobjcryst_available:
            pytest.skip("pyobjcryst package not available")

        global ObjCrystCrystalParSet, Crystal, Atom, Molecule
        global ScatteringPowerAtom
        from pyobjcryst.atom import Atom
        from pyobjcryst.crystal import Crystal
        from pyobjcryst.molecule import Molecule
        from pyobjcryst.scatteringpower import ScatteringPowerAtom

        from diffpy.srfit.structure.objcrystparset import ObjCrystCrystalParSet

        self.occryst = makeC60()
        self.ocmol = self.occryst.GetScatterer("c60")
        return

    def tearDown(self):
        del self.occryst
        del self.ocmol
        return

    def testImplicitBondAngleRestraints(self):
        """Test the structure with implicit bond angles."""
        occryst = self.occryst
        ocmol = self.ocmol

        # Add some bond angles to the molecule
        ocmol.AddBondAngle(ocmol[0], ocmol[5], ocmol[8], 1.1, 0.1, 0.1)
        ocmol.AddBondAngle(ocmol[0], ocmol[7], ocmol[44], 1.3, 0.1, 0.1)

        # make our crystal
        cryst = ObjCrystCrystalParSet("bucky", occryst)
        m = cryst.c60
        m.wrapRestraints()

        # make sure that we have some restraints in the molecule
        assert 2 == len(m._restraints)

        # make sure these evaluate to whatver we get from objcryst
        res0, res1 = m._restraints
        p0 = set([res0.penalty(), res1.penalty()])
        angles = ocmol.GetBondAngleList()
        p1 = set([angles[0].GetLogLikelihood(), angles[1].GetLogLikelihood()])
        assert p0 == p1

        return

    def testObjCrystParSet(self):
        """Test the structure conversion."""

        occryst = self.occryst
        ocmol = self.ocmol
        cryst = ObjCrystCrystalParSet("bucky", occryst)
        m = cryst.c60

        assert cryst.name == "bucky"

        def _testCrystal():
            # Test the lattice
            assert occryst.a == pytest.approx(cryst.a.value)
            assert occryst.b == pytest.approx(cryst.b.getValue())
            assert occryst.c == pytest.approx(cryst.c.getValue())
            assert occryst.alpha == pytest.approx(cryst.alpha.getValue())
            assert occryst.beta == pytest.approx(cryst.beta.getValue())
            assert occryst.gamma == pytest.approx(cryst.gamma.getValue())
            return

        def _testMolecule():

            # Test position / occupancy
            assert ocmol.X == pytest.approx(m.x.getValue())
            assert ocmol.Y == pytest.approx(m.y.getValue())
            assert ocmol.Z == pytest.approx(m.z.getValue())
            assert ocmol.Occupancy == pytest.approx(m.occ.getValue())

            # Test orientation
            assert ocmol.Q0 == pytest.approx(m.q0.getValue())
            assert ocmol.Q1 == pytest.approx(m.q1.getValue())
            assert ocmol.Q2 == pytest.approx(m.q2.getValue())
            assert ocmol.Q3 == pytest.approx(m.q3.getValue())

            # Check the atoms thoroughly
            for i in range(len(ocmol)):
                oca = ocmol[i]
                ocsp = oca.GetScatteringPower()
                a = m.atoms[i]
                assert ocsp.GetSymbol() == a.element
                assert oca.X == pytest.approx(a.x.getValue())
                assert oca.Y == pytest.approx(a.y.getValue())
                assert oca.Z == pytest.approx(a.z.getValue())
                assert oca.Occupancy == pytest.approx(a.occ.getValue())
                assert ocsp.Biso == pytest.approx(a.Biso.getValue())
            return

        _testCrystal()
        _testMolecule()

        # Now change some values from ObjCryst
        ocmol[0].X *= 1.1
        ocmol[0].Occupancy *= 1.1
        ocmol[0].GetScatteringPower().Biso *= 1.1
        ocmol.Q0 *= 1.1
        occryst.a *= 1.1

        _testCrystal()
        _testMolecule()

        # Now change values from the srfit StructureParSet
        cryst.c60.C44.x.setValue(1.1)
        cryst.c60.C44.occ.setValue(1.1)
        cryst.c60.C44.Biso.setValue(1.1)
        cryst.c60.q3.setValue(1.1)
        cryst.a.setValue(1.1)

        _testCrystal()
        _testMolecule()
        return

    def testImplicitBondLengthRestraints(self):
        """Test the structure with implicit bond lengths."""
        occryst = self.occryst
        ocmol = self.ocmol

        # Add some bonds to the molecule
        ocmol.AddBond(ocmol[0], ocmol[5], 3.3, 0.1, 0.1)
        ocmol.AddBond(ocmol[0], ocmol[7], 3.3, 0.1, 0.1)

        # make our crystal
        cryst = ObjCrystCrystalParSet("bucky", occryst)
        m = cryst.c60
        m.wrapRestraints()

        # make sure that we have some restraints in the molecule
        assert 2 == len(m._restraints)

        # make sure these evaluate to whatver we get from objcryst
        res0, res1 = m._restraints
        p0 = set([res0.penalty(), res1.penalty()])
        bonds = ocmol.GetBondList()
        p1 = set([bonds[0].GetLogLikelihood(), bonds[1].GetLogLikelihood()])
        assert p0 == p1

        return

    def testImplicitDihedralAngleRestraints(self):
        """Test the structure with implicit dihedral angles."""
        occryst = self.occryst
        ocmol = self.ocmol

        # Add some bond angles to the molecule
        ocmol.AddDihedralAngle(
            ocmol[0], ocmol[5], ocmol[8], ocmol[41], 1.1, 0.1, 0.1
        )
        ocmol.AddDihedralAngle(
            ocmol[0], ocmol[7], ocmol[44], ocmol[2], 1.3, 0.1, 0.1
        )

        # make our crystal
        cryst = ObjCrystCrystalParSet("bucky", occryst)
        m = cryst.c60
        m.wrapRestraints()

        # make sure that we have some restraints in the molecule
        assert 2 == len(m._restraints)

        # make sure these evaluate to whatver we get from objcryst
        res0, res1 = m._restraints
        p0 = set([res0.penalty(), res1.penalty()])
        angles = ocmol.GetDihedralAngleList()
        p1 = set([angles[0].GetLogLikelihood(), angles[1].GetLogLikelihood()])
        assert p0 == p1

        return

    def testImplicitStretchModes(self):
        """Test the molecule with implicit stretch modes."""
        # Not sure how to make this happen.
        pass

    def testExplicitBondLengthRestraints(self):
        """Test the structure with explicit bond lengths."""
        occryst = self.occryst
        ocmol = self.ocmol

        # make our crystal
        cryst = ObjCrystCrystalParSet("bucky", occryst)
        m = cryst.c60

        # make some bond angle restraints
        res0 = m.restrainBondLength(m.atoms[0], m.atoms[5], 3.3, 0.1, 0.1)
        res1 = m.restrainBondLength(m.atoms[0], m.atoms[7], 3.3, 0.1, 0.1)

        # make sure that we have some restraints in the molecule
        assert 2 == len(m._restraints)

        # make sure these evaluate to whatver we get from objcryst
        p0 = [res0.penalty(), res1.penalty()]
        bonds = ocmol.GetBondList()
        assert 2 == len(bonds)
        p1 = [b.GetLogLikelihood() for b in bonds]
        assert p0 == p1

        return

    def testExplicitBondAngleRestraints(self):
        """Test the structure with explicit bond angles.

        Note that this cannot work with co-linear points as the
        direction of rotation cannot be defined in this case.
        """
        occryst = self.occryst
        ocmol = self.ocmol

        # make our crystal
        cryst = ObjCrystCrystalParSet("bucky", occryst)
        m = cryst.c60

        # restrain some bond angles
        res0 = m.restrainBondAngle(
            m.atoms[0], m.atoms[5], m.atoms[8], 3.3, 0.1, 0.1
        )
        res1 = m.restrainBondAngle(
            m.atoms[0], m.atoms[7], m.atoms[44], 3.3, 0.1, 0.1
        )

        # make sure that we have some restraints in the molecule
        assert 2 == len(m._restraints)

        # make sure these evaluate to whatver we get from objcryst
        p0 = set([res0.penalty(), res1.penalty()])
        angles = ocmol.GetBondAngleList()
        p1 = set([angles[0].GetLogLikelihood(), angles[1].GetLogLikelihood()])
        assert p0 == p1

        return

    def testExplicitDihedralAngleRestraints(self):
        """Test the structure with explicit dihedral angles."""
        occryst = self.occryst
        ocmol = self.ocmol

        # make our crystal
        cryst = ObjCrystCrystalParSet("bucky", occryst)
        m = cryst.c60

        # Restrain some dihedral angles.
        res0 = m.restrainDihedralAngle(
            m.atoms[0], m.atoms[5], m.atoms[8], m.atoms[41], 1.1, 0.1, 0.1
        )
        res1 = m.restrainDihedralAngle(
            m.atoms[0], m.atoms[7], m.atoms[44], m.atoms[2], 1.1, 0.1, 0.1
        )

        # make sure that we have some restraints in the molecule
        assert 2 == len(m._restraints)

        # make sure these evaluate to whatver we get from objcryst
        p0 = set([res0.penalty(), res1.penalty()])
        angles = ocmol.GetDihedralAngleList()
        p1 = set([angles[0].GetLogLikelihood(), angles[1].GetLogLikelihood()])
        assert p0 == p1

        return

    def testExplicitBondLengthParameter(self):
        """Test adding bond length parameters to the molecule."""
        occryst = self.occryst

        # make our crystal
        cryst = ObjCrystCrystalParSet("bucky", occryst)
        m = cryst.c60

        a0 = m.atoms[0]
        a7 = m.atoms[7]
        a20 = m.atoms[20]

        # Add a parameter
        p1 = m.addBondLengthParameter("C07", a0, a7)
        # Have another atom tag along for the ride
        p1.addAtoms([a20])

        xyz0 = numpy.array([a0.x.getValue(), a0.y.getValue(), a0.z.getValue()])
        xyz7 = numpy.array([a7.x.getValue(), a7.y.getValue(), a7.z.getValue()])
        xyz20 = numpy.array(
            [a20.x.getValue(), a20.y.getValue(), a20.z.getValue()]
        )

        dd = xyz0 - xyz7
        d0 = numpy.dot(dd, dd) ** 0.5
        assert d0 == pytest.approx(p1.getValue(), abs=1e-6)

        # Record the unit direction of change for later
        u = dd / d0

        # Change the value
        scale = 1.05
        p1.setValue(scale * d0)

        # Verify that it has changed.
        assert scale * d0 == pytest.approx(p1.getValue(), abs=1e-6)

        xyz0a = numpy.array(
            [a0.x.getValue(), a0.y.getValue(), a0.z.getValue()]
        )
        xyz7a = numpy.array(
            [a7.x.getValue(), a7.y.getValue(), a7.z.getValue()]
        )
        xyz20a = numpy.array(
            [a20.x.getValue(), a20.y.getValue(), a20.z.getValue()]
        )

        dda = xyz0a - xyz7a
        d1 = numpy.dot(dda, dda) ** 0.5

        assert scale * d0 == pytest.approx(d1, abs=1e-6)

        # Verify that only the second and third atoms have moved.

        assert numpy.array_equal(xyz0, xyz0a)

        xyz7calc = xyz7 + (1 - scale) * d0 * u
        for i in range(3):
            assert xyz7a[i] == pytest.approx(xyz7calc[i], abs=1e-5)

        xyz20calc = xyz20 + (1 - scale) * d0 * u
        for i in range(3):
            assert xyz20a[i] == pytest.approx(xyz20calc[i], abs=1e-6)

        return

    def testExplicitBondAngleParameter(self):
        """Test adding bond angle parameters to the molecule."""
        occryst = self.occryst

        # make our crystal
        cryst = ObjCrystCrystalParSet("bucky", occryst)
        m = cryst.c60

        a0 = m.atoms[0]
        a7 = m.atoms[7]
        a20 = m.atoms[20]
        a25 = m.atoms[25]

        xyz0 = numpy.array([a0.x.getValue(), a0.y.getValue(), a0.z.getValue()])
        xyz7 = numpy.array([a7.x.getValue(), a7.y.getValue(), a7.z.getValue()])
        xyz20 = numpy.array(
            [a20.x.getValue(), a20.y.getValue(), a20.z.getValue()]
        )
        xyz25 = numpy.array(
            [a25.x.getValue(), a25.y.getValue(), a25.z.getValue()]
        )

        v1 = xyz7 - xyz0
        d1 = numpy.dot(v1, v1) ** 0.5
        v2 = xyz7 - xyz20
        d2 = numpy.dot(v2, v2) ** 0.5

        angle0 = numpy.arccos(numpy.dot(v1, v2) / (d1 * d2))

        # Add a parameter
        p1 = m.addBondAngleParameter("C0720", a0, a7, a20)
        # Have another atom tag along for the ride
        p1.addAtoms([a25])

        assert angle0 == pytest.approx(p1.getValue(), abs=1e-6)

        # Change the value
        scale = 1.05
        p1.setValue(scale * angle0)

        # Verify that it has changed.
        assert scale * angle0 == pytest.approx(p1.getValue(), abs=1e-6)

        xyz0a = numpy.array(
            [a0.x.getValue(), a0.y.getValue(), a0.z.getValue()]
        )
        xyz7a = numpy.array(
            [a7.x.getValue(), a7.y.getValue(), a7.z.getValue()]
        )
        xyz20a = numpy.array(
            [a20.x.getValue(), a20.y.getValue(), a20.z.getValue()]
        )
        xyz25a = numpy.array(
            [a25.x.getValue(), a25.y.getValue(), a25.z.getValue()]
        )

        v1a = xyz7a - xyz0a
        d1a = numpy.dot(v1a, v1a) ** 0.5
        v2a = xyz7a - xyz20a
        d2a = numpy.dot(v2a, v2a) ** 0.5

        angle1 = numpy.arccos(numpy.dot(v1a, v2a) / (d1a * d2a))

        assert scale * angle0 == pytest.approx(angle1, abs=1e-6)

        # Verify that only the last two atoms have moved.

        assert numpy.array_equal(xyz0, xyz0a)
        assert numpy.array_equal(xyz7, xyz7a)
        assert not numpy.array_equal(xyz20, xyz20a)
        assert not numpy.array_equal(xyz25, xyz25a)

        return

    def testExplicitDihedralAngleParameter(self):
        """Test adding dihedral angle parameters to the molecule."""
        occryst = self.occryst

        # make our crystal
        cryst = ObjCrystCrystalParSet("bucky", occryst)
        m = cryst.c60

        a0 = m.atoms[0]
        a7 = m.atoms[7]
        a20 = m.atoms[20]
        a25 = m.atoms[25]
        a33 = m.atoms[33]

        xyz0 = numpy.array([a0.x.getValue(), a0.y.getValue(), a0.z.getValue()])
        xyz7 = numpy.array([a7.x.getValue(), a7.y.getValue(), a7.z.getValue()])
        xyz20 = numpy.array(
            [a20.x.getValue(), a20.y.getValue(), a20.z.getValue()]
        )
        xyz25 = numpy.array(
            [a25.x.getValue(), a25.y.getValue(), a25.z.getValue()]
        )
        xyz33 = numpy.array(
            [a33.x.getValue(), a33.y.getValue(), a33.z.getValue()]
        )

        v12 = xyz0 - xyz7
        v23 = xyz7 - xyz20
        v34 = xyz20 - xyz25
        v123 = numpy.cross(v12, v23)
        v234 = numpy.cross(v23, v34)

        d123 = numpy.dot(v123, v123) ** 0.5
        d234 = numpy.dot(v234, v234) ** 0.5
        angle0 = -numpy.arccos(numpy.dot(v123, v234) / (d123 * d234))

        # Add a parameter
        p1 = m.addDihedralAngleParameter("C072025", a0, a7, a20, a25)
        # Have another atom tag along for the ride
        p1.addAtoms([a33])

        assert angle0 == pytest.approx(p1.getValue(), abs=1e-6)

        # Change the value
        scale = 1.05
        p1.setValue(scale * angle0)

        # Verify that it has changed.
        assert scale * angle0 == pytest.approx(p1.getValue(), abs=1e-6)

        xyz0a = numpy.array(
            [a0.x.getValue(), a0.y.getValue(), a0.z.getValue()]
        )
        xyz7a = numpy.array(
            [a7.x.getValue(), a7.y.getValue(), a7.z.getValue()]
        )
        xyz20a = numpy.array(
            [a20.x.getValue(), a20.y.getValue(), a20.z.getValue()]
        )
        xyz25a = numpy.array(
            [a25.x.getValue(), a25.y.getValue(), a25.z.getValue()]
        )
        xyz33a = numpy.array(
            [a33.x.getValue(), a33.y.getValue(), a33.z.getValue()]
        )

        v12a = xyz0a - xyz7a
        v23a = xyz7a - xyz20a
        v34a = xyz20a - xyz25a
        v123a = numpy.cross(v12a, v23a)
        v234a = numpy.cross(v23a, v34a)

        d123a = numpy.dot(v123a, v123a) ** 0.5
        d234a = numpy.dot(v234a, v234a) ** 0.5
        angle1 = -numpy.arccos(numpy.dot(v123a, v234a) / (d123a * d234a))
        assert scale * angle0 == pytest.approx(angle1, abs=1e-6)

        # Verify that only the last two atoms have moved.

        assert numpy.array_equal(xyz0, xyz0a)
        assert numpy.array_equal(xyz7, xyz7a)
        assert numpy.array_equal(xyz20, xyz20a)
        assert not numpy.array_equal(xyz25, xyz25a)
        assert not numpy.array_equal(xyz33, xyz33a)

        return


class TestCreateSpaceGroup:
    """Test space group creation from pyobjcryst structures.

    This makes sure that the space groups created by the structure
    parameter set are correct.
    """

    @pytest.fixture(autouse=True)
    def setup(self, diffpy_structure_available, pyobjcryst_available):
        # shared setup
        if not diffpy_structure_available:
            pytest.skip("diffpy.structure package not available")
        if not pyobjcryst_available:
            pytest.skip("pyobjcryst package not available")

        global ObjCrystCrystalParSet, spacegroups
        from diffpy.srfit.structure.objcrystparset import ObjCrystCrystalParSet
        from diffpy.structure import spacegroups

    @staticmethod
    def getObjCrystParSetSpaceGroup(sg):
        """Make an ObjCrystCrystalParSet with the proper space group."""
        from pyobjcryst.spacegroup import SpaceGroup

        sgobjcryst = SpaceGroup(sg.short_name)
        sgnew = ObjCrystCrystalParSet._createSpaceGroup(sgobjcryst)
        return sgnew

    @staticmethod
    def hashDiffPySpaceGroup(sg):
        lines = [str(sg.number % 1000)] + sorted(map(str, sg.iter_symops()))
        s = "\n".join(lines)
        return s

    def sgsEquivalent(self, sg1, sg2):
        """Check to see if two space group objects are the same."""
        hash1 = self.hashDiffPySpaceGroup(sg1)
        hash2 = self.hashDiffPySpaceGroup(sg2)
        return hash1 == hash2

    # FIXME: only about 50% of the spacegroups pass the assertion
    # test disabled even if cctbx is installed
    def xtestCreateSpaceGroup(self):
        """Check all sgtbx space groups for proper conversion to SpaceGroup."""

        try:
            from cctbx import sgtbx
        except ImportError:
            return

        for smbls in sgtbx.space_group_symbol_iterator():
            shn = smbls.hermann_mauguin()
            short_name = shn.replace(" ", "")
            if spacegroups.IsSpaceGroupIdentifier(short_name):
                sg = spacegroups.GetSpaceGroup(shn)
                sgnew = self.getObjCrystParSetSpaceGroup(sg)
                # print("dbsg: " + repr(self.sgsEquivalent(sg, sgnew)))
                assert self.sgsEquivalent(sg, sgnew)
        return


# End of class TestCreateSpaceGroup

if __name__ == "__main__":
    unittest.main()
