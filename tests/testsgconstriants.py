#!/usr/bin/env python
"""Tests space group constraints."""

import unittest
import os

import numpy

from diffpy.Structure import Structure
from diffpy.srfit import adapt
from diffpy.srfit.util import isFixed
from diffpy.srfit.structure import constrainAsSpaceGroup

thisfile = locals().get('__file__', 'testdiffpystructure.py')
tests_dir = os.path.dirname(os.path.abspath(thisfile))
testdata_dir = os.path.join(tests_dir, 'testdata')

class TestSGConstraints(unittest.TestCase):
    """Test creation of space group constraints."""

    def testConstrainAsSpaceGroup1(self):
        """Test the constrainAsSpaceGroup function."""
        fpath = os.path.join(testdata_dir, "ni.cif")
        stru = Structure(filename = fpath)
        phase = adapt(stru, "nickel")
        sgpars = constrainAsSpaceGroup(phase, 225)

        lattice = phase.getLattice()
        self.assertTrue(lattice.a.isVaried())
        self.assertTrue(lattice.b.isConstrained())
        self.assertTrue(lattice.b._constraint is lattice.a)
        self.assertTrue(lattice.c.isConstrained())
        self.assertTrue(lattice.c._constraint is lattice.a)
        self.assertTrue(isFixed(lattice.alpha))
        self.assertTrue(isFixed(lattice.beta))
        self.assertTrue(isFixed(lattice.gamma))

        self.assertEqual(2, len(sgpars))
        scatterers = phase.getScatterers()

        # Ni1
        Ni1 = scatterers[0]
        self.assertTrue(isFixed(Ni1.x))
        self.assertTrue(isFixed(Ni1.y))
        self.assertTrue(isFixed(Ni1.z))
        self.assertTrue(isFixed(Ni1.U11))
        self.assertTrue(isFixed(Ni1.U22))
        self.assertTrue(isFixed(Ni1.U33))
        self.assertTrue(isFixed(Ni1.U12))
        self.assertTrue(isFixed(Ni1.U13))
        self.assertTrue(isFixed(Ni1.U23))
        self.assertTrue(Ni1.Uisoequiv.isVaried())

        # Ni2
        Ni2 = scatterers[1]
        self.assertTrue(isFixed(Ni2.x))
        self.assertTrue(isFixed(Ni2.y))
        self.assertTrue(isFixed(Ni2.z))
        self.assertTrue(isFixed(Ni2.U11))
        self.assertTrue(isFixed(Ni2.U22))
        self.assertTrue(isFixed(Ni2.U33))
        self.assertTrue(isFixed(Ni2.U12))
        self.assertTrue(isFixed(Ni2.U13))
        self.assertTrue(isFixed(Ni2.U23))
        self.assertTrue(Ni2.Uisoequiv.isConstrained())
        self.assertTrue(Ni2.Uisoequiv._constraint is Ni1.Uisoequiv)

        # Ni3
        Ni3 = scatterers[2]
        self.assertTrue(isFixed(Ni3.x))
        self.assertTrue(isFixed(Ni3.y))
        self.assertTrue(isFixed(Ni3.z))
        self.assertTrue(isFixed(Ni3.U11))
        self.assertTrue(isFixed(Ni3.U22))
        self.assertTrue(isFixed(Ni3.U33))
        self.assertTrue(isFixed(Ni3.U12))
        self.assertTrue(isFixed(Ni3.U13))
        self.assertTrue(isFixed(Ni3.U23))
        self.assertTrue(Ni3.Uisoequiv.isConstrained())
        self.assertTrue(Ni3.Uisoequiv._constraint is Ni1.Uisoequiv)

        # Ni4
        Ni4 = scatterers[3]
        self.assertTrue(isFixed(Ni4.x))
        self.assertTrue(isFixed(Ni4.y))
        self.assertTrue(isFixed(Ni4.z))
        self.assertTrue(isFixed(Ni4.U11))
        self.assertTrue(isFixed(Ni4.U22))
        self.assertTrue(isFixed(Ni4.U33))
        self.assertTrue(isFixed(Ni4.U12))
        self.assertTrue(isFixed(Ni4.U13))
        self.assertTrue(isFixed(Ni4.U23))
        self.assertTrue(Ni4.Uisoequiv.isConstrained())
        self.assertTrue(Ni4.Uisoequiv._constraint is Ni1.Uisoequiv)

        return

    def testConstrainAsSpaceGroup2(self):
        """Test the constrainAsSpaceGroup function."""
        stru = makeLaMnO3_P1()
        phase = adapt(stru, "LaMnO3")
        sgpars = constrainAsSpaceGroup(phase, "P b n m",
                scatterers = phase.getScatterers()[::2])

        lattice = phase.getLattice()
        self.assertTrue(lattice.a.isVaried())
        self.assertTrue(lattice.b.isVaried())
        self.assertTrue(lattice.c.isVaried())
        self.assertTrue(isFixed(lattice.alpha))
        self.assertTrue(isFixed(lattice.beta))
        self.assertTrue(isFixed(lattice.gamma))

        # Test the unconstrained atoms
        for scatterer in phase.getScatterers()[1::2]:
            self.assertFalse(scatterer.x.isVaried())
            self.assertFalse(scatterer.y.isVaried())
            self.assertFalse(scatterer.z.isVaried())
            self.assertFalse(scatterer.U11.isVaried())
            self.assertFalse(scatterer.U22.isVaried())
            self.assertFalse(scatterer.U33.isVaried())
            self.assertFalse(scatterer.U12.isVaried())
            self.assertFalse(scatterer.U13.isVaried())
            self.assertFalse(scatterer.U23.isVaried())
            self.assertFalse(scatterer.x.isConstrained())
            self.assertFalse(scatterer.y.isConstrained())
            self.assertFalse(scatterer.z.isConstrained())
            self.assertFalse(scatterer.U11.isConstrained())
            self.assertFalse(scatterer.U22.isConstrained())
            self.assertFalse(scatterer.U33.isConstrained())
            self.assertFalse(scatterer.U12.isConstrained())
            self.assertFalse(scatterer.U13.isConstrained())
            self.assertFalse(scatterer.U23.isConstrained())

        # Now test the constrained parameters
        self.assertEqual(30, len(sgpars))
        scatterers = phase.getScatterers()[::2]
        # La1
        La1 = scatterers[0]
        self.assertTrue(La1.x.isVaried())
        self.assertTrue(La1.x in sgpars)
        self.assertTrue(La1.y.isVaried())
        self.assertTrue(La1.y in sgpars)
        self.assertTrue(isFixed(La1.z))
        self.assertTrue(La1.U11.isVaried())
        self.assertTrue(La1.U11 in sgpars)
        self.assertTrue(La1.U22.isVaried())
        self.assertTrue(La1.U22 in sgpars)
        self.assertTrue(La1.U33.isVaried())
        self.assertTrue(La1.U33 in sgpars)
        self.assertTrue(La1.U12.isVaried())
        self.assertTrue(La1.U12 in sgpars)
        self.assertTrue(isFixed(La1.U13))
        self.assertTrue(isFixed(La1.U23))
        self.assertTrue(isFixed(La1.Uisoequiv))

        # La3
        La3 = scatterers[1]
        self.assertTrue(La3.x.isConstrained())
        self.assertTrue(La3.y.isConstrained())
        self.assertTrue(isFixed(La3.z))
        self.assertTrue(La3.U11.isConstrained())
        self.assertTrue(La3.U22.isConstrained())
        self.assertTrue(La3.U33.isConstrained())
        self.assertTrue(La3.U12.isConstrained())
        self.assertTrue(isFixed(La3.U13))
        self.assertTrue(isFixed(La3.U23))
        self.assertTrue(isFixed(La3.Uisoequiv))

        # Mn1
        Mn1 = scatterers[2]
        self.assertTrue(isFixed(Mn1.x))
        self.assertTrue(isFixed(Mn1.y))
        self.assertTrue(isFixed(Mn1.z))
        self.assertTrue(Mn1.U11.isVaried())
        self.assertTrue(Mn1.U11 in sgpars)
        self.assertTrue(Mn1.U22.isVaried())
        self.assertTrue(Mn1.U22 in sgpars)
        self.assertTrue(Mn1.U33.isVaried())
        self.assertTrue(Mn1.U33 in sgpars)
        self.assertTrue(Mn1.U12.isVaried())
        self.assertTrue(Mn1.U12 in sgpars)
        self.assertTrue(Mn1.U13.isVaried())
        self.assertTrue(Mn1.U13 in sgpars)
        self.assertTrue(Mn1.U23.isVaried())
        self.assertTrue(Mn1.U23 in sgpars)
        self.assertTrue(isFixed(Mn1.Uisoequiv))

        # Mn3
        Mn3 = scatterers[3]
        self.assertTrue(isFixed(Mn3.x))
        self.assertTrue(isFixed(Mn3.y))
        self.assertTrue(isFixed(Mn3.z))
        self.assertTrue(Mn3.U11.isConstrained())
        self.assertTrue(Mn3.U22.isConstrained())
        self.assertTrue(Mn3.U33.isConstrained())
        self.assertTrue(Mn3.U12.isConstrained())
        self.assertTrue(Mn3.U13.isConstrained())
        self.assertTrue(Mn3.U23.isConstrained())
        self.assertTrue(isFixed(Mn3.Uisoequiv))

        # O1
        O1 = scatterers[4]
        self.assertTrue(O1.x.isVaried())
        self.assertTrue(O1.x in sgpars)
        self.assertTrue(O1.y.isVaried())
        self.assertTrue(O1.y in sgpars)
        self.assertTrue(isFixed(O1.z))
        self.assertTrue(O1.U11.isVaried())
        self.assertTrue(O1.U11 in sgpars)
        self.assertTrue(O1.U22.isVaried())
        self.assertTrue(O1.U22 in sgpars)
        self.assertTrue(O1.U33.isVaried())
        self.assertTrue(O1.U33 in sgpars)
        self.assertTrue(O1.U12.isVaried())
        self.assertTrue(O1.U12 in sgpars)
        self.assertTrue(isFixed(O1.U13))
        self.assertTrue(isFixed(O1.U23))
        self.assertTrue(isFixed(O1.Uisoequiv))

        # O3
        O3 = scatterers[5]
        self.assertTrue(O3.x.isConstrained())
        self.assertTrue(O3.y.isConstrained())
        self.assertTrue(isFixed(O3.z))
        self.assertTrue(O3.U11.isConstrained())
        self.assertTrue(O3.U22.isConstrained())
        self.assertTrue(O3.U33.isConstrained())
        self.assertTrue(O3.U12.isConstrained())
        self.assertTrue(isFixed(O3.U13))
        self.assertTrue(isFixed(O3.U23))
        self.assertTrue(isFixed(O3.Uisoequiv))

        # O5
        O5 = scatterers[6]
        self.assertTrue(O5.x.isVaried())
        self.assertTrue(O5.x in sgpars)
        self.assertTrue(O5.y.isVaried())
        self.assertTrue(O5.y in sgpars)
        self.assertTrue(O5.z.isVaried())
        self.assertTrue(O5.z in sgpars)
        self.assertTrue(O5.U11.isVaried())
        self.assertTrue(O5.U11 in sgpars)
        self.assertTrue(O5.U22.isVaried())
        self.assertTrue(O5.U22 in sgpars)
        self.assertTrue(O5.U33.isVaried())
        self.assertTrue(O5.U33 in sgpars)
        self.assertTrue(O5.U12.isVaried())
        self.assertTrue(O5.U12 in sgpars)
        self.assertTrue(O5.U13.isVaried())
        self.assertTrue(O5.U13 in sgpars)
        self.assertTrue(O5.U23.isVaried())
        self.assertTrue(O5.U23 in sgpars)
        self.assertTrue(isFixed(O5.Uisoequiv))

        # O7
        O7 = scatterers[7]
        self.assertTrue(O7.x.isConstrained())
        self.assertTrue(O7.y.isConstrained())
        self.assertTrue(O7.z.isConstrained())
        self.assertTrue(O7.U11.isConstrained())
        self.assertTrue(O7.U22.isConstrained())
        self.assertTrue(O7.U33.isConstrained())
        self.assertTrue(O7.U12.isConstrained())
        self.assertTrue(O7.U13.isConstrained())
        self.assertTrue(O7.U23.isConstrained())
        self.assertTrue(isFixed(O7.Uisoequiv))

        # O9
        O9 = scatterers[8]
        self.assertTrue(O9.x.isConstrained())
        self.assertTrue(O9.y.isConstrained())
        self.assertTrue(O9.z.isConstrained())
        self.assertTrue(O9.U11.isConstrained())
        self.assertTrue(O9.U22.isConstrained())
        self.assertTrue(O9.U33.isConstrained())
        self.assertTrue(O9.U12.isConstrained())
        self.assertTrue(O9.U13.isConstrained())
        self.assertTrue(O9.U23.isConstrained())
        self.assertTrue(isFixed(O9.Uisoequiv))

        # O11
        O11 = scatterers[9]
        self.assertTrue(O11.x.isConstrained())
        self.assertTrue(O11.y.isConstrained())
        self.assertTrue(O11.z.isConstrained())
        self.assertTrue(O11.U11.isConstrained())
        self.assertTrue(O11.U22.isConstrained())
        self.assertTrue(O11.U33.isConstrained())
        self.assertTrue(O11.U12.isConstrained())
        self.assertTrue(O11.U13.isConstrained())
        self.assertTrue(O11.U23.isConstrained())
        self.assertTrue(isFixed(O11.Uisoequiv))

        return


def makeLaMnO3_P1():
    stru = Structure()
    stru.readStr(lamno3stru)
    return stru


lamno3stru =\
"""\
title  Cell structure file of LaMnO3.0
format pdffit
scale   1.000000
sharp   0.000000,  0.000000,  1.000000,  3.500000
spcgr   Pbnm
cell    5.486341,  5.619215,  7.628206, 90.000000, 90.000000, 90.000000
dcell   0.000118,  0.000156,  0.000118,  0.000000,  0.000000,  0.000000
ncell          1,         1,         1,        20
atoms
LA          0.99609631        0.03214940        0.25000000       1.0000
            0.00003041        0.00000852        0.00000000       0.0000
            0.00253993        0.00253993        0.00253993
            0.00000214        0.00000214        0.00000214
            0.00000000        0.00000000        0.00000000
            0.00000000        0.00000000        0.00000000
LA          0.49609631        0.46785060        0.75000000       1.0000
            0.00003041        0.00000852        0.00000000       0.0000
            0.00253993        0.00253993        0.00253993
            0.00000214        0.00000214        0.00000214
            0.00000000        0.00000000        0.00000000
            0.00000000        0.00000000        0.00000000
LA          0.00390369        0.96785063        0.75000000       1.0000
            0.00003041        0.00000852        0.00000000       0.0000
            0.00253993        0.00253993        0.00253993
            0.00000214        0.00000214        0.00000214
            0.00000000        0.00000000        0.00000000
            0.00000000        0.00000000        0.00000000
LA          0.50390369        0.53214937        0.25000000       1.0000
            0.00003041        0.00000852        0.00000000       0.0000
            0.00253993        0.00253993        0.00253993
            0.00000214        0.00000214        0.00000214
            0.00000000        0.00000000        0.00000000
            0.00000000        0.00000000        0.00000000
MN          0.00000000        0.50000000        0.00000000       1.0000
            0.00000000        0.00000000        0.00000000       0.0000
            0.00065337        0.00065337        0.00065337
            0.00000165        0.00000165        0.00000165
            0.00000000        0.00000000        0.00000000
            0.00000000        0.00000000        0.00000000
MN          0.50000000        0.00000000        0.00000000       1.0000
            0.00000000        0.00000000        0.00000000       0.0000
            0.00065337        0.00065337        0.00065337
            0.00000165        0.00000165        0.00000165
            0.00000000        0.00000000        0.00000000
            0.00000000        0.00000000        0.00000000
MN          0.00000000        0.50000000        0.50000000       1.0000
            0.00000000        0.00000000        0.00000000       0.0000
            0.00065337        0.00065337        0.00065337
            0.00000165        0.00000165        0.00000165
            0.00000000        0.00000000        0.00000000
            0.00000000        0.00000000        0.00000000
MN          0.50000000        0.00000000        0.50000000       1.0000
            0.00000000        0.00000000        0.00000000       0.0000
            0.00065337        0.00065337        0.00065337
            0.00000165        0.00000165        0.00000165
            0.00000000        0.00000000        0.00000000
            0.00000000        0.00000000        0.00000000
O           0.05957463        0.49616399        0.25000000       1.0000
            0.00001546        0.00001610        0.00000000       0.0000
            0.00082010        0.00082010        0.00082010
            0.00000137        0.00000137        0.00000137
            0.00000000        0.00000000        0.00000000
            0.00000000        0.00000000        0.00000000
O           0.55957460        0.00383601        0.75000000       1.0000
            0.00001546        0.00001610        0.00000000       0.0000
            0.00082010        0.00082010        0.00082010
            0.00000137        0.00000137        0.00000137
            0.00000000        0.00000000        0.00000000
            0.00000000        0.00000000        0.00000000
O           0.94042540        0.50383604        0.75000000       1.0000
            0.00001546        0.00001610        0.00000000       0.0000
            0.00082010        0.00082010        0.00082010
            0.00000137        0.00000137        0.00000137
            0.00000000        0.00000000        0.00000000
            0.00000000        0.00000000        0.00000000
O           0.44042537        0.99616396        0.25000000       1.0000
            0.00001546        0.00001610        0.00000000       0.0000
            0.00082010        0.00082010        0.00082010
            0.00000137        0.00000137        0.00000137
            0.00000000        0.00000000        0.00000000
            0.00000000        0.00000000        0.00000000
O           0.72005206        0.28938726        0.03111255       1.0000
            0.00001528        0.00001560        0.00002506       0.0000
            0.00512371        0.00512371        0.00512371
            0.00000153        0.00000153        0.00000153
            0.00000000        0.00000000        0.00000000
            0.00000000        0.00000000        0.00000000
O           0.22005206        0.21061274        0.96888745       1.0000
            0.00001528        0.00001560        0.00002506       0.0000
            0.00512371        0.00512371        0.00512371
            0.00000153        0.00000153        0.00000153
            0.00000000        0.00000000        0.00000000
            0.00000000        0.00000000        0.00000000
O           0.27994794        0.71061277        0.53111255       1.0000
            0.00001528        0.00001560        0.00002506       0.0000
            0.00512371        0.00512371        0.00512371
            0.00000153        0.00000153        0.00000153
            0.00000000        0.00000000        0.00000000
            0.00000000        0.00000000        0.00000000
O           0.77994794        0.78938723        0.46888745       1.0000
            0.00001528        0.00001560        0.00002506       0.0000
            0.00512371        0.00512371        0.00512371
            0.00000153        0.00000153        0.00000153
            0.00000000        0.00000000        0.00000000
            0.00000000        0.00000000        0.00000000
O           0.27994794        0.71061277        0.96888745       1.0000
            0.00001528        0.00001560        0.00002506       0.0000
            0.00512371        0.00512371        0.00512371
            0.00000153        0.00000153        0.00000153
            0.00000000        0.00000000        0.00000000
            0.00000000        0.00000000        0.00000000
O           0.77994794        0.78938723        0.03111255       1.0000
            0.00001528        0.00001560        0.00002506       0.0000
            0.00512371        0.00512371        0.00512371
            0.00000153        0.00000153        0.00000153
            0.00000000        0.00000000        0.00000000
            0.00000000        0.00000000        0.00000000
O           0.72005206        0.28938726        0.46888745       1.0000
            0.00001528        0.00001560        0.00002506       0.0000
            0.00512371        0.00512371        0.00512371
            0.00000153        0.00000153        0.00000153
            0.00000000        0.00000000        0.00000000
            0.00000000        0.00000000        0.00000000
O           0.22005206        0.21061274        0.53111255       1.0000
            0.00001528        0.00001560        0.00002506       0.0000
            0.00512371        0.00512371        0.00512371
            0.00000153        0.00000153        0.00000153
            0.00000000        0.00000000        0.00000000
            0.00000000        0.00000000        0.00000000
"""

if __name__ == "__main__":

    unittest.main()

