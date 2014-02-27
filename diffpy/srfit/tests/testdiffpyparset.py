#!/usr/bin/env python
##############################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2010 Trustees of the Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
##############################################################################
"""Tests for diffpy.srfit.structure package."""

import unittest

import numpy

from utils import testcase, TestCaseStructure

# Global variables to be assigned in setUp
Atom = Lattice = Structure = DiffpyStructureParSet = None


class TestParameterAdapter(testcase(TestCaseStructure)):

    def setUp(self):
        global Atom, Lattice, Structure, DiffpyStructureParSet
        from diffpy.Structure import Atom, Lattice, Structure
        from diffpy.srfit.structure.diffpyparset import DiffpyStructureParSet

    def testDiffpyStructureParSet(self):
        """Test the structure conversion."""

        a1 = Atom("Cu", xyz = numpy.array([.0, .1, .2]), Uisoequiv = 0.003)
        a2 = Atom("Ag", xyz = numpy.array([.3, .4, .5]), Uisoequiv = 0.002)
        l = Lattice(2.5, 2.5, 2.5, 90, 90, 90)

        dsstru = Structure([a1,a2], l)
        # Structure makes copies
        a1 = dsstru[0]
        a2 = dsstru[1]

        s = DiffpyStructureParSet("CuAg", dsstru)

        self.assertEquals(s.name, "CuAg")

        def _testAtoms():
            # Check the atoms thoroughly
            self.assertEquals(a1.element, s.Cu0.element)
            self.assertEquals(a2.element, s.Ag0.element)
            self.assertEquals(a1.Uisoequiv, s.Cu0.Uiso.getValue())
            self.assertEquals(a2.Uisoequiv, s.Ag0.Uiso.getValue())
            self.assertEquals(a1.Bisoequiv, s.Cu0.Biso.getValue())
            self.assertEquals(a2.Bisoequiv, s.Ag0.Biso.getValue())
            for i in xrange(1,4):
                for j in xrange(i,4):
                    uijstru = getattr(a1, "U%i%i"%(i,j))
                    uij = getattr(s.Cu0, "U%i%i"%(i,j)).getValue()
                    uji = getattr(s.Cu0, "U%i%i"%(j,i)).getValue()
                    self.assertEquals(uijstru, uij)
                    self.assertEquals(uijstru, uji)
                    bijstru = getattr(a1, "B%i%i"%(i,j))
                    bij = getattr(s.Cu0, "B%i%i"%(i,j)).getValue()
                    bji = getattr(s.Cu0, "B%i%i"%(j,i)).getValue()
                    self.assertEquals(bijstru, bij)
                    self.assertEquals(bijstru, bji)

            self.assertEquals(a1.xyz[0], s.Cu0.x.getValue())
            self.assertEquals(a1.xyz[1], s.Cu0.y.getValue())
            self.assertEquals(a1.xyz[2], s.Cu0.z.getValue())
            return

        def _testLattice():

            # Test the lattice
            self.assertEquals(dsstru.lattice.a, s.lattice.a.getValue())
            self.assertEquals(dsstru.lattice.b, s.lattice.b.getValue())
            self.assertEquals(dsstru.lattice.c, s.lattice.c.getValue())
            self.assertEquals(dsstru.lattice.alpha, s.lattice.alpha.getValue())
            self.assertEquals(dsstru.lattice.beta, s.lattice.beta.getValue())
            self.assertEquals(dsstru.lattice.gamma, s.lattice.gamma.getValue())

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
        self.assertNotEquals(d, dsstru.lattice.dist(a1.xyz, a2.xyz))
        return



if __name__ == "__main__":
    unittest.main()
