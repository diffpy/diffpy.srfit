#!/usr/bin/env python

"""Unit tests for tagmanager.py
"""

# version
__id__ = '$Id$'

import os
import unittest

from diffpy.srfit.util.tagmanager import TagManager

# useful variables
thisfile = locals().get('__file__', 'file.py')
tests_dir = os.path.dirname(os.path.abspath(thisfile))
# testdata_dir = os.path.join(tests_dir, 'testdata')

##############################################################################
class TestTagManager(unittest.TestCase):

    def setUp(self):
        self.m = TagManager()
        return

    def tearDown(self):
        return

    def test_tag(self):
        """check TagManager.tag()
        """
        m = self.m
        obj = 3
        m.tag(obj, "3", "three")
        m.tag(obj, "tri")
        tags = set(m.tags(obj))
        self.assertEquals(set(["3", "three", "tri"]), tags)

        # Try an unhashable object
        obj = set()
        self.assertRaises(TypeError, obj, "unhashable")
        return

    def test_untag(self):
        """check TagManager.untag()
        """
        m = self.m
        obj = 3
        m.tag(obj, "3", "three", "tri", "tres", "trois")
        tags = set(m.tags(obj))
        self.assertEquals(set(["3", "three", "tri", "tres", "trois"]), tags)
        self.assertRaises(KeyError, m.untag, obj, "4")
        self.assertRaises(KeyError, m.untag, 4, "3")
        self.assertRaises(KeyError, m.untag, 4, "4")
        m.untag(obj, "3")
        tags = set(m.tags(obj))
        self.assertEquals(set(["three", "tri", "tres", "trois"]), tags)
        m.untag(obj, "three", "tri")
        tags = set(m.tags(obj))
        self.assertEquals(set(["tres", "trois"]), tags)
        m.untag(obj)
        tags = set(m.tags(obj))
        self.assertEquals(set(), tags)
        return

    def test_objects(self):
        """check TagManager.objects()
        """
        m = self.m
        m.tag(3, "3", "number")
        m.tag(4, "4", "number")
        objs = set([3,4])
        self.assertEqual(m.objects(), set())
        self.assertEqual(m.objects("number"), objs)
        self.assertEqual(m.objects("3"), set([3]))
        self.assertRaises(KeyError, m.objects, "fail")
        return

# End of class TestTagManager

if __name__ == '__main__':
    unittest.main()

# End of file
