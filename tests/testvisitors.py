#!/usr/bin/env python
"""Tests for refinableobj module."""

import diffpy.srfit.equation.visitors as visitors
import diffpy.srfit.equation.literals as literals
import unittest

from utils import _makeArgs


class TestEvaluator(unittest.TestCase):

    def testSimpleFunction(self):
        """Test a simple function."""

        # Make some variables
        v1, v2, v3, v4 = _makeArgs(4)

        # Make some operations
        mult = literals.MultiplicationOperator()
        plus = literals.AdditionOperator()
        minus = literals.SubtractionOperator()

        # Create the equation (v1+v3)*(v4-v2)
        plus.addLiteral(v1)
        plus.addLiteral(v3)
        minus.addLiteral(v4)
        minus.addLiteral(v2)
        mult.addLiteral(plus)
        mult.addLiteral(minus)

        # Set the values of the variables.
        # The equation should evaluate to (1+3)*(4-2) = 8
        v1.setValue(1)
        v2.setValue(2)
        v3.setValue(3)
        v4.setValue(4)

        # Evaluate this
        evaluator = visitors.Evaluator()
        mult.identify(evaluator)
        evaluator.click()
        self.assertEqual(8, evaluator.value)
        self.assertTrue(evaluator._clicker > plus.clicker)
        self.assertTrue(evaluator._clicker > minus.clicker)
        self.assertTrue(evaluator._clicker > mult.clicker)
        self.assertTrue(evaluator._clicker > v1.clicker)
        self.assertTrue(evaluator._clicker > v2.clicker)
        self.assertTrue(evaluator._clicker > v3.clicker)
        self.assertTrue(evaluator._clicker > v4.clicker)

        # Change one of the variables
        v1.setValue(7)
        self.assertTrue(v1.clicker > evaluator._clicker)
        self.assertTrue(plus.clicker > evaluator._clicker)
        self.assertTrue(mult.clicker > evaluator._clicker)
        self.assertTrue(v1.clicker == plus.clicker)
        self.assertTrue(v1.clicker == mult.clicker)
        mult.identify(evaluator)
        evaluator.click()
        self.assertEqual(20, evaluator.value)
        return
        
    def testCustomFunction(self):
        """Test a custom function."""
        evaluator = visitors.Evaluator()

        # Make some variables
        v1, v2, v3 = _makeArgs(3)

        # Make some operations
        mult = literals.MultiplicationOperator()
        plus = literals.AdditionOperator()
        minus = literals.SubtractionOperator()
        import numpy
        sin = literals.UFuncOperator(numpy.sin)

        # Create the equation v1*sin(v2) + v3
        sin.addLiteral(v2)
        mult.addLiteral(sin)
        mult.addLiteral(v1)
        plus.addLiteral(v3)
        plus.addLiteral(mult)

        # Give the variables values. The equation should evaluate to
        # 2*sin(pi/6)+3
        v1.setValue(2)
        v2.setValue(numpy.pi/6)
        v3.setValue(3)

        # Evaluate this
        plus.identify(evaluator)
        evaluator.click()
        self.assertAlmostEqual(4, evaluator.value)
        return

    def testArray(self):
        """Test functions operating on a numpy array."""
        evaluator = visitors.Evaluator()

        # Make some variables
        v1, v2 = _makeArgs(2)

        # Make some operations
        mult = literals.MultiplicationOperator()
        import numpy
        sum = literals.Operator()
        sum.name = sum.symbol = "sum"
        sum.nin = 1
        sum.operation = numpy.sum

        # Create the equation sum(v1*v2)
        sum.addLiteral(mult)
        mult.addLiteral(v1)
        mult.addLiteral(v2)

        # Give the variables values. 
        # v1 = 5
        # v2 = array([0, 1, 2, 3, 4])
        v1.setValue(5)
        v2.setValue(numpy.arange(5))

        # Evaluate this. It should give 50.
        sum.identify(evaluator)
        evaluator.click()
        self.assertAlmostEqual(50, evaluator.value)
        return

    def testTrivialEquation(self):
        v1 = literals.Argument(value = 1)
        evaluator = visitors.Evaluator()
        v1.identify(evaluator)
        evaluator.click()
        self.assertEqual(1, evaluator.value)

        v1.setValue(2)
        v1.identify(evaluator)
        evaluator.click()
        self.assertEqual(2, evaluator.value)

        # can we do it again?
        v1.identify(evaluator)
        evaluator.click()
        self.assertEqual(2, evaluator.value)
        return

    def testTrivialPartition(self):
        """Test a function with a partition in it."""
        p1 = literals.Partition("p1")

        # Make the partition 1|2
        v1, v2, v3 = _makeArgs(3)
        p1.addArgument(v1)
        p1.addArgument(v2)
        p1.addArgument(v3)
        evaluator = visitors.Evaluator()
        p1.identify(evaluator)
        evaluator.click()
        self.assertEqual(6, evaluator.value)

        # Change a value
        v3.setValue(4)
        p1.identify(evaluator)
        evaluator.click()
        self.assertEqual(7, evaluator.value)

        # Can we do it again?
        p1.identify(evaluator)
        evaluator.click()
        self.assertEqual(7, evaluator.value)
        return

    def testSimplePartition(self):
        """Test a function with a partition in it."""
        p1 = literals.Partition("p1")

        # Make the partition 1|2
        v1, v2, v3 = _makeArgs(3)
        p1.addArgument(v1)
        p1.addArgument(v2)

        # make the equation sin(3*partition).
        mult = literals.MultiplicationOperator()
        import numpy
        sin = literals.UFuncOperator(numpy.sin)

        mult.addLiteral(p1)
        mult.addLiteral(v3)
        sin.addLiteral(mult)

        # With no operations able to combine, the equation will be combined
        # after the root operation.
        # sin(3*1) + sin(3*2)
        evaluator = visitors.Evaluator()
        sin.identify(evaluator)
        evaluator.click()
        self.assertEqual(numpy.sin(3*1)+numpy.sin(3*2), evaluator.value)

        # Change the value of an argument.
        v1.setValue(2)
        sin.identify(evaluator)
        evaluator.click()
        self.assertEqual(numpy.sin(3*2)+numpy.sin(3*2), evaluator.value)
        v1.setValue(1)

        # Now let the '*' operator combine. This should give
        # sin(3*1 + 3*2)
        mult.setCombine()
        sin.identify(evaluator)
        evaluator.click()
        self.assertEqual(numpy.sin(3*1 + 3*2), evaluator.value)

        # Do that again to prove that it works.
        sin.identify(evaluator)
        evaluator.click()
        self.assertEqual(numpy.sin(3*1 + 3*2), evaluator.value)
        return

    def testCombinePartition(self):
        """Test a function with a partition in it and the CombineOperator."""

        p1 = literals.Partition("p1")
        # Make the partition 1|2
        v1, v2, v3 = _makeArgs(3)
        p1.addArgument(v1)
        p1.addArgument(v2)

        # make the equation sin(3*partition).
        mult = literals.MultiplicationOperator()
        import numpy
        sin = literals.UFuncOperator(numpy.sin)
        comb = literals.CombineOperator()

        # This combines after the '*' operator. This should give
        # sin(3*1 + 3*2)
        mult.addLiteral(p1)
        mult.addLiteral(v3)
        comb.addLiteral(mult)
        sin.addLiteral(comb)

        evaluator = visitors.Evaluator()
        sin.identify(evaluator)
        evaluator.click()
        self.assertEqual(numpy.sin(3*1+3*2), evaluator.value)

        evaluator = visitors.Evaluator()
        sin.identify(evaluator)
        self.assertEqual(numpy.sin(3*1 + 3*2), evaluator.value)
        return

    def testTwoPartitions(self):
        """Test a function with two partitions."""

        p1 = literals.Partition("p1")
        p2 = literals.Partition("p2")
        v1, v2, v3, v4, v5 = _makeArgs(5)

        # make the equation (4|5) + ( (1|2) + 3 )
        p1.addArgument(v1)
        p1.addArgument(v2)

        p2.addArgument(v4)
        p2.addArgument(v5)

        addroot = literals.AdditionOperator()
        add = literals.AdditionOperator()

        add.addLiteral(p1)
        add.addLiteral(v3)
        addroot.addLiteral(p2)
        addroot.addLiteral(add)

        # The partitions colide in the top '+' and are therefore collected
        # before then. This should give
        # (4|5) + (4|5) = 9 + 9 = 18
        evaluator = visitors.Evaluator()
        addroot.identify(evaluator)
        evaluator.click()
        self.assertEqual(18, evaluator.value)
        return

    def testTaggedPartition(self):
        """Test a function with a partition in it."""
        p1 = literals.Partition("p1")

        # Make the partition 1|2
        v1, v2, v3 = _makeArgs(3)
        p1.addArgument(v1, "v1")
        p1.addArgument(v2, "v2")

        # make the equation sin(3*partition).
        mult = literals.MultiplicationOperator()
        import numpy
        sin = literals.UFuncOperator(numpy.sin)

        mult.addLiteral(p1)
        mult.addLiteral(v3)
        sin.addLiteral(mult)

        # Make it so that the multiplication argument only works on arguments
        # with that "v1" tag.
        mult.clearTags()
        mult.addTags("v1")
        #This should then give
        # sin(3*1) + sin(2)
        evaluator = visitors.Evaluator()
        sin.identify(evaluator)
        evaluator.click()
        self.assertEqual(numpy.sin(3*1)+numpy.sin(2), evaluator.value)

        # Now let it act on all tags
        mult.clearTags()
        mult.addTags("v1", "v2")
        #This should then give
        # sin(3*1) + sin(3*2)
        evaluator = visitors.Evaluator()
        sin.identify(evaluator)
        evaluator.click()
        self.assertEqual(numpy.sin(3*1)+numpy.sin(3*2), evaluator.value)

        # Now let mult work on "v1" and sin work on "v2"
        mult.clearTags()
        mult.addTags("v1")
        sin.clearTags()
        sin.addTags("v2")
        # This should give
        # 3*1 + sin(2)
        evaluator = visitors.Evaluator()
        sin.identify(evaluator)
        evaluator.click()
        self.assertEqual(3*1+numpy.sin(2), evaluator.value)

        # Don't let mult operate, give it a tag that is not in the partition.
        mult.clearTags()
        mult.addTags("xyz")
        sin.clearTags()
        #This should then give
        # sin(1) + sin(2)
        evaluator = visitors.Evaluator()
        sin.identify(evaluator)
        evaluator.click()
        self.assertEqual(numpy.sin(1)+numpy.sin(2), evaluator.value)

        # Don't let mult either operate
        mult.clearTags()
        mult.addTags("xyz")
        sin.clearTags()
        sin.addTags("xyz")
        #This should then give
        # 1 + 2
        evaluator = visitors.Evaluator()
        sin.identify(evaluator)
        evaluator.click()
        self.assertEqual(3, evaluator.value)
        return

    def testArrayPartition(self):
        """Test an array partition."""

        import numpy

        _x = numpy.arange(0, 10, 0.05)
        g1 = numpy.exp(-0.5*((_x - 3.0)/0.2)**2)
        g2 = numpy.exp(-0.5*((_x - 5.0)/0.2)**2)
        g3 = numpy.exp(-0.5*((_x - 7.0)/0.2)**2)

        # Three gaussians with different centers
        a1 = literals.Argument(name = "a1", value = g1)
        a2 = literals.Argument(name = "a2", value = g2)
        a3 = literals.Argument(name = "a3", value = g3)

        class ThreePeak(literals.Partition):

            def __init__(self):
                literals.Partition.__init__(self, "threepeak")

                self.addArgument(a1, "A-A", "A")
                self.addArgument(a2, "A-B", "B-A", "A", "B")
                self.addArgument(a3, "B-B", "B")
                return

        amp = literals.Argument(name = "amp", value = 1.0)
        mult = literals.MultiplicationOperator()
        mult.addLiteral(amp)
        mult.addLiteral(ThreePeak())

        # Test this.
        evaluator = visitors.Evaluator()
        mult.identify(evaluator)
        evaluator.click()
        self.assertTrue( numpy.array_equal(g1+g2+g3, evaluator.value) )

        # Now we want to selectively increase the peak height of all peaks
        # involving "A"
        amp.setValue(2.0)
        mult.addTags("A")
        mult.identify(evaluator)
        evaluator.click()
        self.assertTrue( numpy.array_equal(2*g1+2*g2+g3, evaluator.value) )

        # Back to normal
        amp.setValue(1.0)
        mult.clearTags()
        mult.identify(evaluator)
        evaluator.click()
        self.assertTrue( numpy.array_equal(g1+g2+g3, evaluator.value) )

        # Let't try back-to-back amplitude changes
        amp.setValue(3.0)
        amp2 = literals.Argument(name = "amp2", value = 5.0)
        mult2 = literals.MultiplicationOperator()

        mult2.addLiteral(amp2)
        mult2.addLiteral(mult)

        mult.addTags("A")
        mult2.addTags("B")

        evaluator = visitors.Evaluator()
        mult2.identify(evaluator)
        evaluator.click()
        r = evaluator.value - (3*g1 + 3*5*g2 + 5*g3)
        self.assertTrue( numpy.dot(r, r)**0.5 < 1e-6)
        return

    def testSimpleGenerator(self):
        """Test a simple Generator."""

        import numpy

        _x = numpy.arange(0, 10, 0.05)
        g1 = numpy.exp(-0.5*((_x - 3.0)/2)**2)
        g2 = numpy.exp(-0.5*((_x - 5.0)/2)**2)
        g3 = numpy.exp(-0.5*((_x - 7.0)/2)**2)

        # Three gaussians with different centers
        a1 = literals.Argument(name = "a1", value = g1)
        a2 = literals.Argument(name = "a2", value = g2)
        a3 = literals.Argument(name = "a3", value = g3)

        class ThreePeak(literals.Partition):

            def __init__(self):
                literals.Partition.__init__(self, "threepeak")

                self.addArgument(a1, "A-A", "A")
                self.addArgument(a2, "A-B", "B-A", "A", "B")
                self.addArgument(a3, "B-B", "B")
                return

        gen1 = literals.Generator("threepeak")
        gen1.literal = ThreePeak()

        amp = literals.Argument(name = "amp", value = 1.0)
        mult = literals.MultiplicationOperator()
        mult.addLiteral(amp)
        mult.addLiteral(gen1)

        # Test this.
        evaluator = visitors.Evaluator()
        mult.identify(evaluator)
        evaluator.click()
        self.assertTrue( numpy.array_equal(g1+g2+g3, evaluator.value) )

        # Now we want to selectively increase the peak height of all peaks
        # involving "A"
        amp.setValue(2.0)
        mult.addTags("A")
        mult.identify(evaluator)
        evaluator.click()
        self.assertTrue( numpy.array_equal(2*g1+2*g2+g3, evaluator.value) )

        # Back to normal
        amp.setValue(1.0)
        mult.clearTags()
        mult.identify(evaluator)
        evaluator.click()
        self.assertTrue( numpy.array_equal(g1+g2+g3, evaluator.value) )
        return


    def testGenerator(self):
        """Test a real partitioning of a PDF using a generator."""
        
        import numpy

        _x = numpy.arange(0, 10, 0.05)

        def _gaussian(sigma, x0, x):
            """Gaussian function."""
            return 1.0/(numpy.sqrt(2*numpy.pi)*sigma) \
                    * numpy.exp( -0.5 * ((x-x0)/sigma)**2 )

        from diffpy.srfit.equation import builder
        g = builder.wrapFunction("g", _gaussian, 3)

        s1 = builder.ArgumentBuilder(name = "s1", value = 0.2)
        s2 = builder.ArgumentBuilder(name = "s2", value = 0.2)
        s3 = builder.ArgumentBuilder(name = "s3", value = 0.2)
        x1 = builder.ArgumentBuilder(name = "x1", value = 3.0)
        x2 = builder.ArgumentBuilder(name = "x2", value = 5.0)
        x3 = builder.ArgumentBuilder(name = "x3", value = 7.0)
        A1 = builder.ArgumentBuilder(name = "A1", value = 1.0)
        A2 = builder.ArgumentBuilder(name = "A2", value = 1.0)
        A3 = builder.ArgumentBuilder(name = "A3", value = 1.0)
        x = builder.ArgumentBuilder(name = "x", value = _x, const=True)
        ggen1 = A1*g(s1, x1, x)
        ggen2 = A2*g(s2, x2, x)
        ggen3 = A3*g(s3, x3, x)
        g1 = ggen1.getEquation()
        g2 = ggen2.getEquation()
        g3 = ggen3.getEquation()

        class ThreePeak(literals.Generator):
            """Generator of three gaussians."""

            def __init__(self):
                literals.Generator.__init__(self, "threepeak")
                self.literal = None
                self.__setup()
                return

            def __setup(self):
                self.g1 = g1
                self.g2 = g2
                self.g3 = g3

                self.a1 = literals.Argument(value = self.g1())
                self.a2 = literals.Argument(value = self.g2())
                self.a3 = literals.Argument(value = self.g3())

                self.addLiteral(self.g1.root)
                self.addLiteral(self.g2.root)
                self.addLiteral(self.g3.root)

                # Make the initial partition
                self.literal = literals.Partition("threepeak")
                self.literal.addArgument(self.a1, "A-A", "A")
                self.literal.addArgument(self.a2, "A-B", "B-A", "A", "B")
                self.literal.addArgument(self.a3, "B-B", "B")
                return

            def generate(self, clicker):

                self.a1.setValue( self.g1() )
                self.a2.setValue( self.g2() )
                self.a3.setValue( self.g3() )

                return

        # Now our three peaks are generated dynamically
        gen1 = ThreePeak()

        amp = literals.Argument(name = "amp", value = 1.0)
        mult = literals.MultiplicationOperator()
        mult.addLiteral(amp)
        mult.addLiteral(gen1)

        # Test this.
        evaluator = visitors.Evaluator()
        mult.identify(evaluator)
        evaluator.click()
        self.assertTrue( numpy.array_equal(g1()+g2()+g3(), evaluator.value) )

        # Now we want to selectively increase the peak height of all peaks
        # involving "A"
        amp.setValue(2.0)
        mult.addTags("A")
        mult.identify(evaluator)
        evaluator.click()
        self.assertTrue( numpy.array_equal(2*g1()+2*g2()+g3(), evaluator.value) )
        # Back to normal
        amp.setValue(1.0)
        mult.clearTags()
        mult.identify(evaluator)
        evaluator.click()
        self.assertTrue( numpy.array_equal(g1()+g2()+g3(), evaluator.value) )

        # Now change some values
        finder = visitors.ArgFinder() 
        gen1.identify(finder)

        for arg in finder.args:
            arg.setValue(0.12345)
            mult.identify(evaluator)
            evaluator.click()
            self.assertTrue( 
                    numpy.array_equal(g1()+g2()+g3(), evaluator.value) )

        return




class TestValidator(unittest.TestCase):

    def testSimpleFunction(self):
        """Test a simple function."""

        # Make some variables
        v1, v2, v3, v4 = _makeArgs(4)

        # Make some operations
        mult = literals.MultiplicationOperator()
        plus = literals.AdditionOperator()
        minus = literals.SubtractionOperator()

        # Let's hobble the plus operator
        plus.name = None
        plus.symbol = None
        plus.operation = None

        # Partially define the equation (v1+v3)*(v4-v2). Let's only give one
        # variable to the '-' operation.
        plus.addLiteral(v1)
        plus.addLiteral(v3)
        minus.addLiteral(v4)
        mult.addLiteral(plus)
        mult.addLiteral(minus)

        # Now validate
        validator = visitors.Validator()
        mult.identify(validator)
        self.assertEqual(4, len(validator.errors))

        # Fix the equation
        minus.addLiteral(v3)
        validator.reset()
        mult.identify(validator)
        self.assertEqual(3, len(validator.errors))

        # Fix the name of plus
        plus.name = "add"
        validator.reset()
        mult.identify(validator)
        self.assertEqual(2, len(validator.errors))

        # Fix the symbol of plus
        plus.symbol = "+"
        validator.reset()
        mult.identify(validator)
        self.assertEqual(1, len(validator.errors))

        # Fix the operation of plus
        import numpy
        plus.operation = numpy.add
        validator.reset()
        mult.identify(validator)
        self.assertEqual(0, len(validator.errors))

        # Add another literal to minus
        minus.addLiteral(v1)
        validator.reset()
        mult.identify(validator)
        self.assertEqual(1, len(validator.errors))

        return

class TestArgFinder(unittest.TestCase):

    def testSimpleFunction(self):
        """Test a simple function."""

        # Make some variables
        v1, v2, v3, v4 = _makeArgs(4)

        # Make some operations
        mult = literals.MultiplicationOperator()
        plus = literals.AdditionOperator()
        minus = literals.SubtractionOperator()

        # Create the equation (v1+v3)*(v4-v2)
        plus.addLiteral(v1)
        plus.addLiteral(v3)
        minus.addLiteral(v4)
        minus.addLiteral(v2)
        mult.addLiteral(plus)
        mult.addLiteral(minus)

        # Set the values of the variables.
        # The equation should evaluate to (1+3)*(4-2) = 8
        v1.setValue(1)
        v2.setValue(2)
        v3.setValue(3)
        v4.setValue(4)

        # now get the args
        argfinder = visitors.ArgFinder()
        mult.identify(argfinder)
        args = argfinder.args
        self.assertEqual(4, len(args))
        self.assertTrue(v1 in args)
        self.assertTrue(v2 in args)
        self.assertTrue(v3 in args)
        self.assertTrue(v4 in args)

        return

    def testSimplePartition(self):
        """Test an equation with a partition."""
        p1 = literals.Partition("p1")

        # Make the partition 1|2
        v1, v2, v3 = _makeArgs(3)
        p1.addArgument(v1)
        p1.addArgument(v2)

        # make the equation sin(3*partition).
        mult = literals.MultiplicationOperator()
        import numpy
        sin = literals.UFuncOperator(numpy.sin)

        mult.addLiteral(p1)
        mult.addLiteral(v3)
        sin.addLiteral(mult)

        argfinder = visitors.ArgFinder()
        sin.identify(argfinder)
        args = argfinder.args
        self.assertEqual(3, len(args))
        self.assertTrue(v1 in args)
        self.assertTrue(v2 in args)
        self.assertTrue(v3 in args)
        return

if __name__ == "__main__":

    unittest.main()

