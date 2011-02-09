#!/usr/bin/env python
import unittest
import os.path

from diffpy.Structure import Structure

import diffpy.srfit.structure.diffpyadapters as dpsa
from diffpy.srfit.adapters.adaptersmod import selfgetter, nosetter
from diffpy.srfit.adapters.adaptersmod import adapt 

thisfile = locals().get('__file__', 'testdiffpystructure.py')
tests_dir = os.path.dirname(os.path.abspath(thisfile))
testdata_dir = os.path.join(tests_dir, 'testdata')

class TestDiffpyStructureAdapters(unittest.TestCase):
    """Test diffpy structure adapters nodes."""

    def testAdapt(self):
        """Make sure adapt works properly."""
        fpath = os.path.join(testdata_dir, "ni.cif")
        stru = Structure(filename = fpath)

        astru = adapt(stru)
        self.assertTrue(isinstance(astru, dpsa.DiffpyStructureAdapter))
        self._testStructureAdapter(astru)

        a1 = adapt(stru[0])
        self.assertTrue(isinstance(a1, dpsa.DiffpyAtomAdapter))
        self._testAtomAdapter(a1)

        alattice = adapt(stru.lattice)
        self.assertTrue(isinstance(alattice, dpsa.DiffpyLatticeAdapter))
        self._testLatticeAdapter(alattice)
        return

    def _testAtomAdapter(self, a1):
        """Test the atom adapter."""
        self.assertAlmostEqual(0, a1.x.value)
        self.assertAlmostEqual(0, a1.y.value)
        self.assertAlmostEqual(0, a1.z.value)
        self.assertAlmostEqual(0, a1.U11.value)
        self.assertAlmostEqual(0, a1.U22.value)
        self.assertAlmostEqual(0, a1.U33.value)
        self.assertAlmostEqual(0, a1.U12.value)
        self.assertAlmostEqual(0, a1.U13.value)
        self.assertAlmostEqual(0, a1.U23.value)
        self.assertAlmostEqual(0, a1.Uisoequiv.value)
        self.assertAlmostEqual(0, a1.B11.value)
        self.assertAlmostEqual(0, a1.B22.value)
        self.assertAlmostEqual(0, a1.B33.value)
        self.assertAlmostEqual(0, a1.B12.value)
        self.assertAlmostEqual(0, a1.B13.value)
        self.assertAlmostEqual(0, a1.B23.value)
        self.assertAlmostEqual(0, a1.Bisoequiv.value)
        a1.x.value = 0.1
        a1.y.value = 0.2
        a1.z.value = 0.3
        self.assertAlmostEqual(0.1, a1.x.value)
        self.assertAlmostEqual(0.2, a1.y.value)
        self.assertAlmostEqual(0.3, a1.z.value)
        a1.U11.value = 0.1
        a1.U22.value = 0.2
        a1.U33.value = 0.3
        a1.U12.value = 0.4
        a1.U13.value = 0.5
        a1.U23.value = 0.6
        self.assertAlmostEqual(0.1, a1.U11.value)
        self.assertAlmostEqual(0.2, a1.U22.value)
        self.assertAlmostEqual(0.3, a1.U33.value)
        self.assertAlmostEqual(0.4, a1.U12.value)
        self.assertAlmostEqual(0.5, a1.U13.value)
        self.assertAlmostEqual(0.6, a1.U23.value)
        a1.B11.value = 0.1
        a1.B22.value = 0.2
        a1.B33.value = 0.3
        a1.B12.value = 0.4
        a1.B13.value = 0.5
        a1.B23.value = 0.6
        self.assertAlmostEqual(0.1, a1.B11.value)
        self.assertAlmostEqual(0.2, a1.B22.value)
        self.assertAlmostEqual(0.3, a1.B33.value)
        self.assertAlmostEqual(0.4, a1.B12.value)
        self.assertAlmostEqual(0.5, a1.B13.value)
        self.assertAlmostEqual(0.6, a1.B23.value)
        a1.Uisoequiv.value = 0.1
        self.assertAlmostEqual(0.1, a1.Uisoequiv.value)
        a1.Bisoequiv.value = 0.2
        self.assertAlmostEqual(0.2, a1.Bisoequiv.value)
        a1.x.value = 0
        a1.y.value = 0
        a1.z.value = 0
        a1.U11.value = 0
        a1.U22.value = 0
        a1.U33.value = 0
        a1.U12.value = 0
        a1.U13.value = 0
        a1.U23.value = 0
        a1.Uisoequiv.value = 0
        xyz = a1.getXYZ()
        self.assertTrue(a1.x is xyz[0])
        self.assertTrue(a1.y is xyz[1])
        self.assertTrue(a1.z is xyz[2])
        adps = a1.getADPs()
        self.assertTrue(a1.U11 is adps[0])
        self.assertTrue(a1.U22 is adps[1])
        self.assertTrue(a1.U33 is adps[2])
        self.assertTrue(a1.U12 is adps[3])
        self.assertTrue(a1.U13 is adps[4])
        self.assertTrue(a1.U23 is adps[5])
        self.assertTrue(a1.Uisoequiv is adps[6])
        return

    def _testLatticeAdapter(self, alattice):
        """Test lattice adapter."""
        self.assertEqual("deg", alattice.angunits)
        self.assertEqual(alattice.getLatPars(),
                (alattice.a, alattice.b, alattice.c, alattice.alpha,
                    alattice.beta, alattice.gamma))
        self.assertEqual(3.52387, alattice.a.value)
        self.assertEqual(3.52387, alattice.b.value)
        self.assertEqual(3.52387, alattice.c.value)
        self.assertEqual(90, alattice.alpha.value)
        self.assertEqual(90, alattice.beta.value)
        self.assertEqual(90, alattice.gamma.value)
        alattice.a.value = 3.5
        alattice.b.value = 3.5
        alattice.c.value = 3.5
        alattice.alpha.value = 91
        alattice.beta.value = 91
        alattice.gamma.value = 91
        self.assertEqual(3.5, alattice.a.value)
        self.assertEqual(3.5, alattice.b.value)
        self.assertEqual(3.5, alattice.c.value)
        self.assertEqual(91, alattice.alpha.value)
        self.assertEqual(91, alattice.beta.value)
        self.assertEqual(91, alattice.gamma.value)
        alattice.a.value = 3.52387
        alattice.b.value = 3.52387
        alattice.c.value = 3.52387
        alattice.alpha.value = 90
        alattice.beta.value = 90
        alattice.gamma.value = 90
        return


    def _testStructureAdapter(self, astru):
        """Test structure adapter."""
        # Test atom access
        stru = astru.get()
        a1 = astru.getAtom(0)
        self.assertTrue(a1.get() is stru[0])
        self.assertTrue(isinstance(a1, dpsa.DiffpyAtomAdapter))
        self._testAtomAdapter(a1)

        self.assertEqual(4, len(astru))
        for a1, a2, a3 in zip(astru, astru.getScatterers(), 
                map(astru.getAtom, range(4))):
            self.assertTrue(a1 is a2)
            self.assertTrue(a1 is a3)

        alattice = astru.getLattice()
        self.assertTrue(isinstance(alattice, dpsa.DiffpyLatticeAdapter))
        self.assertTrue(alattice is astru.lattice)
        self._testLatticeAdapter(alattice)

        # Test distance accessor
        dpar = astru.distance(0, 1)
        self.assertAlmostEqual(3.52387 / 2**0.5, dpar.value)
        alattice.a.value = 3.5
        alattice.b.value = 3.5
        alattice.c.value = 3.5
        self.assertAlmostEqual(3.5 / 2**0.5, dpar.value)
        alattice.a.value = 3.52387
        alattice.b.value = 3.52387
        alattice.c.value = 3.52387

        # Test angle accessor
        dpar = astru.angle(0, 1, 2)
        self.assertAlmostEqual(60, dpar.value)
        return

    def testPickle(self):
        """Test pickling diffpy Structure adapters."""
        fpath = os.path.join(testdata_dir, "ni.cif")
        stru = Structure(filename = fpath)
        astru = dpsa.DiffpyStructureAdapter("nickel", stru, selfgetter,
                nosetter)
        # Test the structure first so attributes get made
        self._testStructureAdapter(astru)
        # Now pickle and reconstitute and test again
        from pickle import dumps, loads
        pstr = dumps(astru)
        astru2 = loads(pstr)
        self._testStructureAdapter(astru2)
        return

if __name__ == "__main__":

    unittest.main()

