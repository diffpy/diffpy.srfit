#!/usr/bin/env python
"""Tests for diffpy.srfit.structure package."""

import unittest

import numpy

import diffpy.srfit.structure.objcrystadapters as oca
from diffpy.srfit.adapters.adaptersmod import selfgetter, nosetter
from diffpy.srfit.adapters.adaptersmod import adapt

from pyobjcryst.crystal import Crystal 
from pyobjcryst.atom import Atom 
from pyobjcryst.molecule import Molecule
from pyobjcryst.scatteringpower import ScatteringPowerAtom

def makeC60():
    """Make a crystal containing the C60 molecule using pyobjcryst."""
    pi = numpy.pi
    c = Crystal(100, 100, 100, "P1")
    c.SetName("c60frame")
    m = Molecule(c, "c60")

    c.AddScatterer(m)

    sp = ScatteringPowerAtom("C", "C")
    sp.SetBiso(0.25)
    #c.AddScatteringPower(sp)

    for i, l in enumerate(c60xyz.strip().splitlines()):
        x, y, z = map(float, l.split())
        m.AddAtom(x, y, z, sp, "C%i"%i)

    # Add an atom while we're at it
    a = Atom(0.1, 0.1, 0.1, "atom", sp)
    c.AddScatterer(a)

    return c

class TestObjCrystAdapters(unittest.TestCase):

    def testAdapt(self):
        """Test the adapt method."""
        cryst = makeC60()
        acryst = adapt(cryst)
        self.assertTrue(isinstance(acryst, oca.ObjCrystCrystalAdapter))
        self._testCrystalAdapter(acryst)

        mol = cryst.GetScatt(0)
        amol = adapt(mol)
        self.assertTrue(isinstance(amol, oca.ObjCrystMoleculeAdapter))
        self._testMoleculeAdapter(amol)

        molatom = cryst.GetScatt(0)[0]
        amolatom = adapt(molatom)
        self.assertTrue(isinstance(amolatom, oca.ObjCrystAtomAdapter))
        self._testMolAtomAdapter(amolatom)

        a = cryst.GetScatt(1)
        aa = adapt(a)
        self.assertTrue(isinstance(aa, oca.ObjCrystAtomAdapter))
        self._testAtomAdapter(aa)

        sp = cryst.GetScatt(1).GetScatteringPower()
        asp = adapt(sp)
        self.assertTrue(isinstance(asp, oca.ObjCrystScatteringPowerAdapter))
        self._testScatteringPowerAdapter(asp)
        return

    def _testScatteringPowerAdapter(self, asp):
        """Test adapted ScatteringPower."""
        self.assertEqual("C", asp.GetSymbol())
        self.assertAlmostEqual(0.25, asp.Biso.value)
        self.assertAlmostEqual(0, asp.B11.value)
        self.assertAlmostEqual(0, asp.B22.value)
        self.assertAlmostEqual(0, asp.B33.value)
        self.assertAlmostEqual(0, asp.B12.value)
        self.assertAlmostEqual(0, asp.B13.value)
        self.assertAlmostEqual(0, asp.B23.value)
        asp.Biso.value = 0.5
        asp.B11.value = 0.5
        asp.B22.value = 0.5
        asp.B33.value = 0.5
        asp.B12.value = 0.5
        asp.B13.value = 0.5
        asp.B23.value = 0.5
        self.assertAlmostEqual(0.5, asp.Biso.value)
        self.assertAlmostEqual(0.5, asp.B11.value)
        self.assertAlmostEqual(0.5, asp.B22.value)
        self.assertAlmostEqual(0.5, asp.B33.value)
        self.assertAlmostEqual(0.5, asp.B12.value)
        self.assertAlmostEqual(0.5, asp.B13.value)
        self.assertAlmostEqual(0.5, asp.B23.value)
        asp.Biso.value = 0.25
        asp.B11.value = 0.0
        asp.B22.value = 0.0
        asp.B33.value = 0.0
        asp.B12.value = 0.0
        asp.B13.value = 0.0
        asp.B23.value = 0.0
        return

    def _testMolAtomAdapter(self, aa):
        """Test molatom adapter."""
        self.assertAlmostEqual(3.451266489, aa.X.value)
        self.assertAlmostEqual(0.685, aa.Y.value)
        self.assertAlmostEqual(0.0, aa.Z.value)
        aa.X.value = 3.4
        aa.Y.value = 0
        aa.Z.value = 0.1
        self.assertAlmostEqual(3.4, aa.X.value)
        self.assertAlmostEqual(0.0, aa.Y.value)
        self.assertAlmostEqual(0.1, aa.Z.value)
        aa.X.value = 3.451266489
        aa.Y.value = 0.685
        aa.Z.value = 0.0
        asp = aa.GetScatteringPower()
        self._testScatteringPowerAdapter(asp)
        xyz = aa.getXYZ()
        self.assertTrue(aa.X is xyz[0])
        self.assertTrue(aa.Y is xyz[1])
        self.assertTrue(aa.Z is xyz[2])
        adps = aa.getADPs()
        self.assertTrue(asp.B11 is adps[0])
        self.assertTrue(asp.B22 is adps[1])
        self.assertTrue(asp.B33 is adps[2])
        self.assertTrue(asp.B12 is adps[3])
        self.assertTrue(asp.B13 is adps[4])
        self.assertTrue(asp.B23 is adps[5])
        self.assertTrue(asp.Biso is adps[6])
        return

    def _testAtomAdapter(self, aa):
        """Test adapted atoms."""
        #aa = oca.ObjCrystAtomAdapter("dummy", a, selfgetter, nosetter)
        self.assertAlmostEqual(0.1, aa.X.value)
        self.assertAlmostEqual(0.1, aa.Y.value)
        self.assertAlmostEqual(0.1, aa.Z.value)
        aa.X.value = 0.0
        aa.Y.value = 0.0
        aa.Z.value = 0.0
        self.assertAlmostEqual(0.0, aa.X.value)
        self.assertAlmostEqual(0.0, aa.Y.value)
        self.assertAlmostEqual(0.0, aa.Z.value)
        aa.X.value = 0.1
        aa.Y.value = 0.1
        aa.Z.value = 0.1
        asp = aa.GetScatteringPower()
        self._testScatteringPowerAdapter(asp)
        xyz = aa.getXYZ()
        self.assertTrue(aa.X is xyz[0])
        self.assertTrue(aa.Y is xyz[1])
        self.assertTrue(aa.Z is xyz[2])
        adps = aa.getADPs()
        self.assertTrue(asp.B11 is adps[0])
        self.assertTrue(asp.B22 is adps[1])
        self.assertTrue(asp.B33 is adps[2])
        self.assertTrue(asp.B12 is adps[3])
        self.assertTrue(asp.B13 is adps[4])
        self.assertTrue(asp.B23 is adps[5])
        self.assertTrue(asp.Biso is adps[6])
        return

    def _testCrystalAdapter(self, acryst):
        """Test crystal adapter."""
        self.assertEqual(2, acryst.GetNbScatterer())
        amol = acryst.GetScatt(0)
        self._testMoleculeAdapter(amol)
        aa = acryst.GetScatt(1)
        self._testAtomAdapter(aa)
        alat = acryst.getLattice()
        self.assertTrue(alat is acryst)
        self.assertEqual("rad", alat.angunits)
        self.assertEqual(alat.getLatPars(), (alat.a, alat.b, alat.c,
            alat.alpha, alat.beta, alat.gamma))
        self.assertAlmostEqual(100, alat.a.value)
        self.assertAlmostEqual(100, alat.b.value)
        self.assertAlmostEqual(100, alat.c.value)
        alat.a.value = 99
        alat.b.value = 99
        alat.c.value = 99
        self.assertAlmostEqual(99, alat.a.value)
        self.assertAlmostEqual(99, alat.b.value)
        self.assertAlmostEqual(99, alat.c.value)
        alat.a.value = 100
        alat.b.value = 100
        alat.c.value = 100
        from numpy import pi
        self.assertAlmostEqual(pi/2, alat.alpha.value)
        self.assertAlmostEqual(pi/2, alat.beta.value)
        self.assertAlmostEqual(pi/2, alat.gamma.value)
        alat.alpha.value = pi/3
        alat.beta.value = pi/3
        alat.gamma.value = pi/3
        self.assertAlmostEqual(pi/3, alat.alpha.value)
        self.assertAlmostEqual(pi/3, alat.beta.value)
        self.assertAlmostEqual(pi/3, alat.gamma.value)
        alat.alpha.value = pi/2
        alat.beta.value = pi/2
        alat.gamma.value = pi/2
        scatt = acryst.getScatterers()
        self.assertTrue(scatt[0] is amol)
        self.assertTrue(scatt[1] is aa)
        return

    def _testMoleculeAdapter(self, amol):
        """Test adapted molecule."""
        aa = amol[0]
        self._testMolAtomAdapter(aa)
        self.assertAlmostEqual(0, amol.X.value)
        self.assertAlmostEqual(0, amol.Y.value)
        self.assertAlmostEqual(0, amol.Z.value)
        amol.X.value = 0.1
        amol.Y.value = 0.2
        amol.Z.value = 0.3
        self.assertAlmostEqual(0.1, amol.X.value)
        self.assertAlmostEqual(0.2, amol.Y.value)
        self.assertAlmostEqual(0.3, amol.Z.value)
        amol.X.value = 0
        amol.Y.value = 0
        amol.Z.value = 0
        self.assertAlmostEqual(1, amol.Q0.value)
        self.assertAlmostEqual(0, amol.Q1.value)
        self.assertAlmostEqual(0, amol.Q2.value)
        self.assertAlmostEqual(0, amol.Q3.value)
        amol.Q0.value = 0.95
        amol.Q1.value = 0.01
        amol.Q2.value = 0.02
        amol.Q3.value = 0.03
        self.assertAlmostEqual(0.95, amol.Q0.value)
        self.assertAlmostEqual(0.01, amol.Q1.value)
        self.assertAlmostEqual(0.02, amol.Q2.value)
        self.assertAlmostEqual(0.03, amol.Q3.value)
        amol.Q0.value = 1
        amol.Q1.value = 0
        amol.Q2.value = 0
        amol.Q3.value = 0

        self._testBondLengthPar(amol)
        self._testBondAnglePar(amol)
        self._testDiAnglePar(amol)
        return

    def _testBondLengthPar(self, amol):
        """Test bond length parameters."""
        a0 = amol[0]
        a7 = amol[7]
        a20 = amol[20]

        # Add a parameter
        p1 = amol.GetBond(0, 7)
        p1.addAtoms(a20)

        xyz0 = numpy.array([a0.X.value, a0.Y.value, a0.Z.value])
        xyz7 = numpy.array([a7.X.value, a7.Y.value, a7.Z.value])
        xyz20 = numpy.array([a20.X.value, a20.Y.value, a20.Z.value])

        dd = xyz0 - xyz7
        d0 = numpy.dot(dd, dd)**0.5
        self.assertAlmostEquals(d0, p1.value, 6)
        
        # Record the unit direction of change for later
        u = dd/d0

        # Change the value 
        scale = 1.05
        p1.value = scale*d0

        # Verify that it has changed.
        self.assertAlmostEquals(scale*d0, p1.value)

        xyz0a = numpy.array([a0.X.value, a0.Y.value, a0.Z.value])
        xyz7a = numpy.array([a7.X.value, a7.Y.value, a7.Z.value])
        xyz20a = numpy.array([a20.X.value, a20.Y.value, a20.Z.value])

        dda = xyz0a - xyz7a
        d1 = numpy.dot(dda, dda)**0.5

        self.assertAlmostEquals(scale*d0, d1)

        # Verify that only the second and third atoms have moved.

        self.assertTrue(numpy.array_equal(xyz0, xyz0a))

        xyz7calc = xyz7 + (1-scale)*d0*u
        for i in range(3):
            self.assertAlmostEqual(xyz7a[i], xyz7calc[i], 6)

        xyz20calc = xyz20 + (1-scale)*d0*u
        for i in range(3):
            self.assertAlmostEqual(xyz20a[i], xyz20calc[i], 6)

        p1.value = d0

        return

    def _testBondAnglePar(self, amol):
        """Test bond angle parameters."""

        a0 = amol[0]
        a7 = amol[7]
        a20 = amol[20]
        a25 = amol[25]

        xyz0 = numpy.array([a0.X.value, a0.Y.value, a0.Z.value])
        xyz7 = numpy.array([a7.X.value, a7.Y.value, a7.Z.value])
        xyz20 = numpy.array([a20.X.value, a20.Y.value,
            a20.Z.value])
        xyz25 = numpy.array([a25.X.value, a25.Y.value,
            a25.Z.value])


        v1 = xyz7 - xyz0
        d1 = numpy.dot(v1, v1)**0.5
        v2 = xyz7 - xyz20
        d2 = numpy.dot(v2, v2)**0.5

        angle0 = numpy.arccos(numpy.dot(v1, v2)/(d1*d2))

        # Add a parameter
        p1 = amol.GetBondAngle(0, 7, 20)
        # Have another atom tag along for the ride
        p1.addAtoms(a25)

        self.assertAlmostEqual(angle0, p1.value, 6)

        # Change the value 
        scale = 1.05
        p1.value = scale*angle0

        # Verify that it has changed.
        self.assertAlmostEqual(scale*angle0, p1.value, 6)

        xyz0a = numpy.array([a0.X.value, a0.Y.value, a0.Z.value])
        xyz7a = numpy.array([a7.X.value, a7.Y.value, a7.Z.value])
        xyz20a = numpy.array([a20.X.value, a20.Y.value, a20.Z.value])
        xyz25a = numpy.array([a25.X.value, a25.Y.value, a25.Z.value])

        v1a = xyz7a - xyz0a
        d1a = numpy.dot(v1a, v1a)**0.5
        v2a = xyz7a - xyz20a
        d2a = numpy.dot(v2a, v2a)**0.5

        angle1 = numpy.arccos(numpy.dot(v1a, v2a)/(d1a*d2a))

        self.assertAlmostEquals(scale*angle0, angle1)

        # Verify that only the last two atoms have moved.

        self.assertTrue(numpy.array_equal(xyz0, xyz0a))
        self.assertTrue(numpy.array_equal(xyz7, xyz7a))
        self.assertFalse(numpy.array_equal(xyz20, xyz20a))
        self.assertFalse(numpy.array_equal(xyz25, xyz25a))

        p1.value = angle0
        
        return

    def _testDiAnglePar(self, amol):
        """Test dihedral angle parameters."""
        a0 = amol[0]
        a7 = amol[7]
        a20 = amol[20]
        a25 = amol[25]
        a33 = amol[33]

        xyz0 = numpy.array([a0.X.value, a0.Y.value, a0.Z.value])
        xyz7 = numpy.array([a7.X.value, a7.Y.value, a7.Z.value])
        xyz20 = numpy.array([a20.X.value, a20.Y.value, a20.Z.value])
        xyz25 = numpy.array([a25.X.value, a25.Y.value, a25.Z.value])
        xyz33 = numpy.array([a33.X.value, a33.Y.value, a33.Z.value])

        v12 = xyz0 - xyz7
        v23 = xyz7 - xyz20
        v34 = xyz20 - xyz25
        v123 = numpy.cross(v12, v23)
        v234 = numpy.cross(v23, v34)

        d123 = numpy.dot(v123, v123)**0.5
        d234 = numpy.dot(v234, v234)**0.5
        angle0 = -numpy.arccos(numpy.dot(v123, v234)/(d123*d234))

        # Add a parameter
        p1 = amol.GetDihedralAngle(0, 7, 20, 25)
        # Have another atom tag along for the ride
        p1.addAtoms(a33)

        self.assertAlmostEqual(angle0, p1.value, 6)

        # Change the value 
        scale = 1.05
        p1.value = scale*angle0

        # Verify that it has changed.
        self.assertAlmostEqual(scale*angle0, p1.value, 6)

        xyz0a = numpy.array([a0.X.value, a0.Y.value, a0.Z.value])
        xyz7a = numpy.array([a7.X.value, a7.Y.value, a7.Z.value])
        xyz20a = numpy.array([a20.X.value, a20.Y.value, a20.Z.value])
        xyz25a = numpy.array([a25.X.value, a25.Y.value, a25.Z.value])
        xyz33a = numpy.array([a33.X.value, a33.Y.value, a33.Z.value])

        v12a = xyz0a - xyz7a
        v23a = xyz7a - xyz20a
        v34a = xyz20a - xyz25a
        v123a = numpy.cross(v12a, v23a)
        v234a = numpy.cross(v23a, v34a)

        d123a = numpy.dot(v123a, v123a)**0.5
        d234a = numpy.dot(v234a, v234a)**0.5
        angle1 = -numpy.arccos(numpy.dot(v123a, v234a)/(d123a*d234a))

        self.assertAlmostEquals(scale*angle0, angle1)

        # Verify that only the last two atoms have moved.

        self.assertTrue(numpy.array_equal(xyz0, xyz0a))
        self.assertTrue(numpy.array_equal(xyz7, xyz7a))
        self.assertTrue(numpy.array_equal(xyz20, xyz20a))
        self.assertFalse(numpy.array_equal(xyz25, xyz25a))
        self.assertFalse(numpy.array_equal(xyz33, xyz33a))

        p1.value = angle0
        
        return

from diffpy.Structure import SpaceGroups

class TestCreateSpaceGroup(unittest.TestCase):
    """Test space group creation from pyobjcryst structures.

    This makes sure that the space groups created by the structure parameter
    set are correct.

    """

    @staticmethod
    def getObjCrystSpaceGroup(sg):
        """Make an ObjCrystCrystalAdapter with the proper space group."""
        from pyobjcryst.spacegroup import SpaceGroup
        sgobjcryst = SpaceGroup(sg.short_name)
        sgnew = oca.ObjCrystCrystalAdapter._createSpaceGroup(sgobjcryst)
        return sgnew

    @staticmethod
    def hashDiffPySpaceGroup(sg):
        lines = [str(sg.number % 1000)] + sorted(map(str, sg.iter_symops()))
        s = '\n'.join(lines)
        return s

    def sgsEquivalent(self, sg1, sg2):
        """Check to see if two space group objects are the same."""
        hash1 = self.hashDiffPySpaceGroup(sg1)
        hash2 = self.hashDiffPySpaceGroup(sg2)
        return hash1 == hash2

    def testCreateSpaceGroup(self):
        """Check all sgtbx space groups for proper conversion to SpaceGroup."""

        try:
            from cctbx import sgtbx
        except ImportError:
            return

        for smbls in sgtbx.space_group_symbol_iterator():
            shn = smbls.hermann_mauguin()
            short_name = shn.replace(' ', '')
            if SpaceGroups.IsSpaceGroupIdentifier(short_name):
                sg = SpaceGroups.GetSpaceGroup(shn)
                sgnew = self.getObjCrystSpaceGroup(sg)
                equiv = self.sgsEquivalent(sg, sgnew)
                self.assertTrue( self.sgsEquivalent(sg, sgnew) )
        return

c60xyz = \
"""
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

if __name__ == "__main__":

    unittest.main()

