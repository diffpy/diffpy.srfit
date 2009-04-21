#!/usr/bin/env python
"""Tests for refinableobj module."""

import unittest

from diffpy.srfit.fitbase.fitmodel import FitModel

class TestFitModel(unittest.TestCase):

    def setUp(self):
        self.model = FitModel()
        return

if __name__ == "__main__":

    unittest.main()

