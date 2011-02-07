#!/usr/bin/env python
import unittest

from numpy import ones, array, array_equal, dot, concatenate

from diffpy.srfit.fit import residualmod
from diffpy.srfit.fit import functions
from diffpy.srfit.fit.parameters import Var, Par, Pars

class TestFunctions(unittest.TestCase):
    """Test the chi and Rw functions."""

    def testChi(self):
        """Test unbound operator."""
        a, b, c, d = Pars(("a", 3), ("b", 4), ("c", 7), ("d", 0.5))

        res = residualmod.chi(a + b, c)
        self.assertEqual(0, res.value)

        res = residualmod.chi(a + b, c.set(8), d)
        self.assertEqual(-2, res.value)

        y = ones(2) * 7
        dy = ones(2) * 0.5
        calc = ones(2) * 8
        chi = (calc - y)/dy
        chi2 = dot(chi, chi)

        res = residualmod.chi(calc, y, dy)
        retval = dot(res.value, res.value)
        self.assertAlmostEqual(retval, chi2)

        # Make sure we can accept Parameters as well
        res = residualmod.chi(*Pars(("calc",calc), ("y",y), ("dy",dy)))
        retval = dot(res.value, res.value)
        self.assertAlmostEqual(retval, chi2)
        return

    def testRw(self):
        """Test unbound operator."""
        a, b, c, d = Pars(("a", 3), ("b", 4), ("c", 7), ("d", 0.5))

        res = residualmod.Rw(a + b, c)
        self.assertEqual(0, res.value)

        res = residualmod.Rw(a + b, c.set(8), d)
        self.assertEqual(-0.125, res.value)

        y = ones(2) * 7
        w = ones(2) * 0.5
        calc = ones(2) * 8
        diff = calc-y
        Rw2 = dot(w * diff, diff) / dot(w * y, y)

        res = residualmod.Rw(calc, y, w)
        retval = dot(res.value, res.value)
        self.assertAlmostEqual(retval, Rw2)

        # Make sure we can accept Parameters as well
        res = residualmod.Rw(*Pars(("calc",calc), ("y",y), ("w",w)))
        retval = dot(res.value, res.value)
        self.assertAlmostEqual(retval, Rw2)
        return


class TestResidual(unittest.TestCase):
    """Test the chi and Rw functions."""

    def setUp(self):
        y = ones(2) * 8
        dy = ones(2) * 0.5
        p1 = Var("p1", ones(2) * 0.9)
        p2 = Var("p2", ones(2) * 7)
        p3 = Var("p3", ones(2) * 9)
        res1 = residualmod.chi(p2 + p1, y, dy)
        res2 = residualmod.chi(p3 - p1, y, dy)

        self.residual = residualmod.residual(res1, res2)
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.y = y
        self.dy = dy
        return

    def testResiduals(self):
        residual = self.residual
        residual.fithooks[0].verbose = 0
        p1 = self.p1
        p2 = self.p2
        p3 = self.p3
        
        r1 = (p2.value + p1.value - self.y) / self.dy
        r2 = (p3.value - p1.value - self.y) / self.dy

        chiv = concatenate([r1, r2])
        self.assertTrue(array_equal(chiv, residual.vec()))
        self.assertAlmostEqual(dot(chiv, chiv), residual())
        return

    def testProperties(self):
        residual = self.residual
        p1 = self.p1
        p2 = self.p2
        p3 = self.p3
        self.assertEqual(residual.variables, [p1, p2, p3])
        self.assertEqual(residual.names, ["p1", "p2", "p3"])
        self.assertEqual(residual.values, [p1.value, p2.value, p3.value])
        p3.value = 4
        self.assertEqual(residual.values, [p1.value, p2.value, p3.value])
        p3.fix(3)
        self.assertEqual(residual.variables, [p1, p2])
        self.assertEqual(residual.names, ["p1", "p2"])
        self.assertEqual(residual.values, [p1.value, p2.value])
        return

    def testAppendRemove(self):
        residual = self.residual
        p1 = self.p1
        p2 = self.p2
        p3 = self.p3
        p4 = Var("p4", 3)
        residual.append(p4)
        self.assertEqual([p1, p2, p3, p4], residual.variables)
        residual.remove(p4)
        self.assertEqual([p1, p2, p3], residual.variables)
        return

if __name__ == "__main__":

    unittest.main()

