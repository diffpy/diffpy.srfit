#!/usr/bin/env python
import unittest

from diffpy.srfit.adapters.adaptersmod import ObjectAdapter
from diffpy.srfit.adapters.adaptersmod import selfgetter, nosetter
from diffpy.srfit.adapters.adaptersmod import adapt
from diffpy.srfit.adapters.nodes import Parameter
import utils

class TestObjectAdapter(unittest.TestCase):
    """Test ObjectAdapter nodes."""

    def setUp(self):
        self.tlist = utils.TestListContainer()
        self.alist = ObjectAdapter("list", self.tlist, selfgetter, nosetter)
        self.alist.accessors = ["getc", "getpar", "getwhatever"]
        self.alist.ignore = ["configure"]

        self.tdict = utils.TestDictContainer()
        self.adict = ObjectAdapter("dict", self.tdict, selfgetter, nosetter)
        self.adict.accessors = ["getc", "getpar"]
        self.adict.ignore = ["configure"]
        return

    def test__getattr__(self):
        """Test __getattr__ for attributes.

        Behavior: Adapt contained object. Return same object if called for
        again.

        """
        tlist = self.tlist
        alist = self.alist

        aa = alist.a
        self.assertTrue(aa in alist.adapters)
        self.assertTrue(aa in alist._viewers)
        self.assertTrue(alist in aa._viewers)
        self.assertTrue(aa._container is alist)
        # Make sure all is adapted properly
        self.assertTrue(aa.name is "a")
        self.assertEqual(aa.value, tlist.a)
        self.assertTrue(aa is alist.a)
        return

    def testAccessor(self):
        """Test __getattr__ for accessors.

        Behavior: Adapt contained object using accessor as setter. Return same
        object if called for again.

        """
        tlist = self.tlist
        alist = self.alist

        c = alist.getc().rename("c")
        self.assertEqual(c.value, tlist.getc())
        self.assertRaises(AttributeError, c.set, 1.234)
        tlist._c = 4
        self.assertEqual(c._get(), tlist.getc())
        self.assertTrue(c in alist.adapters)
        self.assertTrue(c in alist._viewers)
        self.assertTrue(alist in c._viewers)
        self.assertTrue(c._container is alist)
        self.assertTrue(c is alist.getc())

        b = alist.getpar("b").rename("b")
        self.assertEqual(b.value, tlist.getpar("b"))
        self.assertRaises(AttributeError, b.set, 1.234)
        tlist.b = 4
        self.assertEqual(b._get(), tlist.getpar("b"))
        self.assertTrue(b in alist.adapters)
        self.assertTrue(b in alist._viewers)
        self.assertTrue(alist in b._viewers)
        self.assertTrue(b._container is alist)
        self.assertTrue(b is alist.getpar("b"))

        w = alist.getwhatever(None)
        self.assertRaises(AttributeError, w.set, 1.234)
        self.assertTrue(w in alist.adapters)
        self.assertTrue(w in alist._viewers)
        self.assertTrue(alist in w._viewers)
        self.assertTrue(w._container is alist)
        self.assertTrue(w is alist.getwhatever(None))
        return

    def testIgnore(self):
        """Test __getattr__ when attribute is in the ignore list."""
        tlist = self.tlist
        alist = self.alist

        self.assertFalse(hasattr(tlist, "configured"))
        self.assertTrue(None is alist.configure())
        self.assertTrue(tlist.configured)
        return

    def testMethod(self):
        """Test __getattr__ for methods not ignored or accessors."""
        alist = self.alist
        d = Parameter("d", 4)
        a = alist.a
        func = alist.calc(d)
        self.assertTrue(func in alist._viewers)
        self.assertTrue(alist in func._viewers)
        self.assertTrue(func._container is alist)
        self.assertEquals(func.value, 1+2+3+4)
        d.value = 5
        self.assertEquals(func.value, 1+2+3+5)
        b = Parameter("b", 2)
        a.constrain(b)
        self.assertEquals(func.value, 2+2+3+5)
        d.value = 4
        self.assertEquals(func.value, 2+2+3+4)
        b.value = 1
        self.assertEquals(func.value, 1+2+3+4)
        return

    def test__call__(self):
        """Test __call__ """
        alist = self.alist
        v = Parameter("v", 4)
        f = alist(v)
        self.assertTrue(alist in f._viewers)
        self.assertTrue(f in alist._viewers)
        self.assertTrue(v in f._args)
        self.assertEqual({}, f._kw)
        self.assertTrue(f._isfunction)
        self.assertEqual((1,(1,2,2)), f.value)
        # Try to reference a value from the output
        p1 = f[0]
        self.assertTrue(f in p1._viewers)
        self.assertTrue(p1 in f._viewers)
        self.assertEqual(1, p1.value)
        alist.a.value = 3
        self.assertEqual((3,(1,2,2)), f.value)
        self.assertEqual(3, p1.value)
        p2 = f[1][2]
        self.assertEqual(2, p2.value)
        alist.b.value = 9
        self.assertEqual(9, p2.value)
        return

    # Test for refining functions of containers where the constraints are
    # hiding underneath. This is analogous to refining structure parameters but
    # only passing the structure itself to the calculator.
    def testHiddenConstraint(self):
        """Test a hidden constraint."""
        class TestClass1(object):
            def __init__(self):
                self.a = TestClass2()
                self.b = TestClass2()
                self.b.x = 2
        class TestClass2(object):
            def __init__(self):
                self.x = 1

        t = TestClass1()
        tt = ObjectAdapter("t", t, selfgetter, nosetter)

        def func(obj):
            return t.a.x + t.b.x
        afunc = adapt(func)

        eq = afunc(tt)
        self.assertEqual(3, eq.value)
        tt.a.x.constrain(tt.b.x)
        self.assertEqual(4, eq.value)
        tt.b.x.value = 4
        self.assertEqual(8, eq.value)
        return

    def test__getitem__list(self):
        """Test __getitem__ for list-like containers."""
        tlist = self.tlist
        alist = self.alist

        v1 = alist[0].rename("v1")
        v2 = alist[1].rename("v2")

        getit = lambda : alist[2]
        self.assertRaises(IndexError, getit)

        self.assertEquals(4, v1.value)
        self.assertEquals(5, v2.value)
        v1.value = 6
        self.assertEquals(6, v1.value)
        self.assertEquals(6, tlist[0])
        self.assertTrue(v1 in alist.adapters)
        self.assertTrue(v1 in alist._viewers)
        self.assertTrue(alist in v1._viewers)
        self.assertTrue(v1._container is alist)
        return

    def testlenlist(self):
        """Test len for list-like containers."""
        self.assertTrue(len(self.alist), 2)
        return

    def testkeyslist(self):
        """Test keys for list-like containers."""
        self.assertRaises(AttributeError, self.alist.keys)
        return

    def testiterkeyslist(self):
        """Test iterkeys for list-like containers."""
        self.assertRaises(AttributeError, self.alist.iterkeys)
        return

    def testvalueslist(self):
        """Test values for list-like containers."""
        self.assertRaises(AttributeError, self.alist.values)
        return

    def testitervalueslist(self):
        """Test itervalues for list-like containers."""
        self.assertRaises(AttributeError, self.alist.itervalues)
        return

    def testitemslist(self):
        """Test items for list-like containers."""
        self.assertRaises(AttributeError, self.alist.items)
        return

    def testiteritemslist(self):
        """Test iteritems for list-like containers."""
        self.assertRaises(AttributeError, self.alist.iteritems)
        return

    def test__getitem__dict(self):
        """Test __getitem__ for dict-like containers."""
        tdict = self.tdict
        adict = self.adict

        v1 = adict["d"]
        v2 = adict["e"]

        self.assertTrue("d" is v1.name)
        self.assertTrue("e" is v2.name)

        getit = lambda : adict["f"]
        self.assertRaises(KeyError, getit)

        self._testDictItems(v1, v2)
        return

    def _testDictItems(self, v1, v2):
        adict = self.adict
        tdict = self.tdict
        self.assertEquals(4, v1.value)
        self.assertEquals(5, v2.value)
        v1.value = 6
        self.assertEquals(6, v1.value)
        self.assertEquals(6, tdict["d"])
        self.assertTrue(v1 in adict.adapters)
        self.assertTrue(v1 in adict._viewers)
        self.assertTrue(adict in v1._viewers)
        self.assertTrue(v1._container is adict)
        return

    def testlendict(self):
        """Test len for dict-like containers."""
        self.assertTrue(len(self.adict), 2)
        return

    def testkeysdict(self):
        """Test keys for dict-like containers."""
        self.assertEquals(set(["d", "e"]), set(self.adict.keys()))
        return

    def testiterkeysdict(self):
        """Test iterkeys for dict-like containers."""
        self.assertEquals(set(["d", "e"]), set(self.adict.iterkeys()))
        return

    def testvaluesdict(self):
        """Test values for dict-like containers."""
        v1, v2 = self.adict.values()
        self.assertTrue(v1.name in ("d", "e"))
        self.assertTrue(v2.name in ("d", "e") and v1.name is not v2.name)
        if v1.name is "e": v2, v1 = v1, v2
        self._testDictItems(v1, v2)
        self.assertTrue(v1 is self.adict["d"])
        self.assertTrue(v2 is self.adict["e"])
        return

    def testitervaluesdict(self):
        """Test itervalues for dict-like containers."""
        v1, v2 = self.adict.itervalues()
        self.assertTrue(v1.name in ("d", "e"))
        self.assertTrue(v2.name in ("d", "e") and v1.name is not v2.name)
        if v1.name is "e": v2, v1 = v1, v2
        self._testDictItems(v1, v2)
        self.assertTrue(v1 is self.adict["d"])
        self.assertTrue(v2 is self.adict["e"])
        return

    def testitemsdict(self):
        """Test items for dict-like containers."""
        items = self.adict.items()
        d = dict(items)
        self.assertTrue("d" is d["d"].name)
        self.assertTrue("e" is d["e"].name)
        self._testDictItems(d["d"], d["e"])
        self.assertTrue(d["d"] is self.adict["d"])
        self.assertTrue(d["e"] is self.adict["e"])
        return

    def testiteritemsdict(self):
        """Test iteritems for dict-like containers."""
        items = self.adict.iteritems()
        d = dict(items)
        self.assertTrue("d" is d["d"].name)
        self.assertTrue("e" is d["e"].name)
        self._testDictItems(d["d"], d["e"])
        self.assertTrue(d["d"] is self.adict["d"])
        self.assertTrue(d["e"] is self.adict["e"])
        return

    def testGet(self):
        """Test the get method.

        Behavior: Get's the adapted container after calling get for all
        sub-objects.
        
        """
        tlist = self.tlist
        alist = self.alist

        self.assertTrue(alist.get() is tlist)

        # Check get of sub-objects
        p = Parameter("p", 4)
        alist.a.constrain(p)
        self.assertEqual(alist.a.value, p.value)
        p.value = 3
        self.assertEqual(alist.a.value, p.value)

        return

    def testSet(self):
        """Test the set method.

        Behavior: Sets the adapted object. If the change is invalid, this is
        detected by depended adapters.

        """
        tlist = self.tlist
        alist = self.alist

        # Create embedded containers
        tlist.tlist = utils.TestListContainer()
        tlist.tlist.tlist = utils.TestListContainer()

        alist2 = alist.tlist
        alist3 = alist2.tlist
        a3 = alist3.a

        self.assertTrue(ObjectAdapter is type(alist2))
        self.assertTrue(ObjectAdapter is type(alist3))

        self.assertTrue(alist2 in alist.adapters)
        self.assertTrue(alist3 in alist2.adapters)
        self.assertTrue(a3 in alist3.adapters)

        alist2.set(None)
        self.assertTrue(alist2 in alist.adapters)
        self.assertTrue(alist3 in alist2.adapters)
        self.assertTrue(a3 in alist3.adapters)

        # Make sure the objects are indeed abandoned.
        self.assertRaises(AttributeError, alist3.get)
        self.assertRaises(AttributeError, a3.get)

        # Set to something we can actually use and recheck the adapters.
        tlist4 = utils.TestListContainer()
        tlist4.tlist = utils.TestListContainer()
        tlist4.tlist.a = 9
        alist2.set(tlist4)
        self.assertTrue(alist2.get() is tlist4)
        self.assertEqual(9, a3.value)
        a3.value = 10
        self.assertEqual(10, tlist4.tlist.a)
        return

    def testPickle(self):
        """Test pickling of containers."""
        from pickle import dumps, loads

        # Try the list.
        alist = self.alist
        # Access lots of stuff before pickling
        a = alist.a
        b = alist.getpar("b").rename("b")
        c = alist.getc().rename("c")
        d = Parameter("d", 4)
        # Use a key that will have a different representation after pickling.
        w = alist.getwhatever(testkey)
        func = alist.calc(d)
        v1 = alist[0].rename("v1")
        v2 = alist[1].rename("v2")
        dobj = dumps((alist, a, b, c, d, func, v1, v2, w, testkey))
        alist2, a2, b2, c2, d2, func2, v12, v22, w2, testkey2, = loads(dobj)
        self.assertEqual(alist2.name, alist.name)
        self.assertEqual(alist2.value, alist.value)
        self.assertTrue(a2 is alist2.a)
        self.assertEqual(a2.name, a.name)
        self.assertEqual(a2.value, a.value)
        self.assertTrue(w2 is alist2.getwhatever(testkey2))
        self.assertTrue(b2 is alist2.getpar("b"))
        self.assertEqual(b2.name, b.name)
        self.assertEqual(b2.value, b.value)
        self.assertTrue(c2 is alist2.getc())
        self.assertEqual(c2.name, c.name)
        self.assertEqual(c2.value, c.value)
        self.assertTrue(v12 is alist2[0])
        self.assertEqual(v12.name, v1.name)
        self.assertEqual(v12.value, v1.value)
        self.assertTrue(v22 is alist2[1])
        self.assertEqual(v22.name, v2.name)
        self.assertEqual(v22.value, v2.value)
        self.assertTrue(d2 in func2._args)
        # We don't have to preserve function identity, we just need to make
        # sure that they compute to the same values
        self.assertEqual(func2.value, alist2.calc(d2).value)
        return

class Test(object): pass
testkey = Test()

if __name__ == "__main__":

    unittest.main()

