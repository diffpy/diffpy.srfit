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
"""Tests for refinableobj module."""

import unittest

import pytest

import diffpy.srfit.equation.literals as literals
import diffpy.srfit.equation.visitors as visitors


class TestValidator:
    @pytest.fixture(autouse=True)
    def setup(self, make_args):
        self.make_args = make_args

    def testSimpleFunction(self):
        """Test a simple function."""

        # Make some variables
        v1, v2, v3, v4 = self.make_args(4)

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
        assert 4 == len(validator.errors)

        # Fix the equation
        minus.addLiteral(v3)
        validator.reset()
        mult.identify(validator)
        assert 3 == len(validator.errors)

        # Fix the name of plus
        plus.name = "add"
        validator.reset()
        mult.identify(validator)
        assert 2 == len(validator.errors)

        # Fix the symbol of plus
        plus.symbol = "+"
        validator.reset()
        mult.identify(validator)
        assert 1 == len(validator.errors)

        # Fix the operation of plus
        import numpy

        plus.operation = numpy.add
        validator.reset()
        mult.identify(validator)
        assert 0 == len(validator.errors)

        # Add another literal to minus
        minus.addLiteral(v1)
        validator.reset()
        mult.identify(validator)
        assert 1 == len(validator.errors)

        return


class TestArgFinder:
    @pytest.fixture(autouse=True)
    def setup(self, make_args):
        self.make_args = make_args

    def testSimpleFunction(self):
        """Test a simple function."""

        # Make some variables
        v1, v2, v3, v4 = self.make_args(4)

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
        args = visitors.getArgs(mult)
        assert 4 == len(args)
        assert v1 in args
        assert v2 in args
        assert v3 in args
        assert v4 in args

        return

    def testArg(self):
        """Test just an Argument equation."""
        # Make some variables
        v1 = self.make_args(1)[0]

        args = visitors.getArgs(v1)

        assert 1 == len(args)
        assert args[0] == v1
        return


class TestSwapper:
    @pytest.fixture(autouse=True)
    def setup(self, make_args):
        self.make_args = make_args

    def testSimpleFunction(self):
        """Test a simple function."""

        # Make some variables
        v1, v2, v3, v4, v5 = self.make_args(5)

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
        v5.setValue(5)

        # Evaluate
        assert 8 == mult.value

        # Now swap an argument
        visitors.swap(mult, v2, v5)

        # Check that the operator value is invalidated
        assert mult._value is None
        assert not v2.hasObserver(minus._flush)
        assert v5.hasObserver(minus._flush)

        # now get the args
        args = visitors.getArgs(mult)
        assert 4 == len(args)
        assert v1 in args
        assert v2 not in args
        assert v3 in args
        assert v4 in args
        assert v5 in args

        # Re-evaluate (1+3)*(4-5) = -4
        assert -4 == mult.value

        # Swap out the "-" operator
        plus2 = literals.AdditionOperator()
        visitors.swap(mult, minus, plus2)
        assert mult._value is None
        assert not minus.hasObserver(mult._flush)
        assert plus2.hasObserver(mult._flush)

        # plus2 has no arguments yet. Verify this.
        with pytest.raises(TypeError):
            mult.getValue()
        # Add the arguments to plus2.
        plus2.addLiteral(v4)
        plus2.addLiteral(v5)

        # Re-evaluate (1+3)*(4+5) = 36
        assert 36 == mult.value

        return


if __name__ == "__main__":
    unittest.main()
