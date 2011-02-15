#!/usr/bin/env python
import unittest

import diffpy.srfit.adapters.adaptersmod as adapters
from diffpy.srfit.structure import *
adapt = adapters.adapt

class TestAdapt(unittest.TestCase):
    """Test ContainerAdapter nodes."""

    def testAdaptBuiltin(self):
        """Test the adapt method for builtin types
        
        Other types must be tested elsewhere.
        
        """

        # Test adapter types
        def test(a, b):
            return a+b
        UnboundOperator = adapters.UnboundOperator
        self.assertTrue(isinstance(adapt(test), UnboundOperator))
        self.assertTrue(isinstance(adapt(lambda z: z), UnboundOperator))
        from numpy import sin, vectorize, arange, array
        self.assertTrue(isinstance(adapt(sin), UnboundOperator))
        testv = vectorize(test)
        self.assertTrue(isinstance(adapt(testv), UnboundOperator))
        class Test(object):
            def test1(self, a, b):
                return a+b
        MethodAdapter = adapters.MethodAdapter
        self.assertTrue(isinstance(adapt(Test.test1), MethodAdapter))
        ParameterAdapter = adapters.ParameterAdapter
        self.assertTrue(isinstance(adapt(True), ParameterAdapter))
        self.assertTrue(isinstance(adapt(3.0), ParameterAdapter))
        self.assertTrue(isinstance(adapt(3), ParameterAdapter))
        self.assertTrue(isinstance(adapt(long(3)), ParameterAdapter))
        from numpy import float, float32, float64
        self.assertTrue(isinstance(adapt(float(3)), ParameterAdapter))
        self.assertTrue(isinstance(adapt(float32(3)), ParameterAdapter))
        self.assertTrue(isinstance(adapt(float64(3)), ParameterAdapter))
        ContainerAdapter = adapters.ContainerAdapter
        self.assertTrue(isinstance(adapt({}), ContainerAdapter))
        self.assertTrue(isinstance(adapt([]), ContainerAdapter))
        self.assertTrue(isinstance(adapt(Test()), ContainerAdapter))
        self.assertTrue(isinstance(adapt(1+1j), ContainerAdapter))
        self.assertTrue(isinstance(adapt(arange(2)), ContainerAdapter))
        self.assertTrue(isinstance(adapt(None), ContainerAdapter))
        testgetter = lambda obj: "testgetter"
        def testsetter(obj, val):
            raise RuntimeError
        adapter = adapt(Test(), "test", testgetter, testsetter, ["test1"])
        self.assertTrue(isinstance(adapter, ContainerAdapter))
        self.assertEqual("test", adapter.name)
        self.assertEqual("testgetter", adapter.get())
        self.assertRaises(RuntimeError, adapter.set, 0)
        self.assertTrue("test1" in adapter.ignore)
        return

if __name__ == "__main__":

    unittest.main()

