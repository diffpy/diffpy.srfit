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

"""Unit tests for tagmanager.py
"""

import unittest

from diffpy.srfit.util.tagmanager import TagManager

##############################################################################
class TestTagManager(unittest.TestCase):

    def setUp(self):
        self.m = TagManager()
        self.m.silent = False
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
        self.assertEqual(set(["3", "three", "tri"]), tags)

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
        self.assertEqual(set(["3", "three", "tri", "tres", "trois"]), tags)
        self.assertRaises(KeyError, m.untag, obj, "4")
        self.assertRaises(KeyError, m.untag, 4, "3")
        self.assertRaises(KeyError, m.untag, 4, "4")
        m.untag(obj, "3")
        tags = set(m.tags(obj))
        self.assertEqual(set(["three", "tri", "tres", "trois"]), tags)
        m.untag(obj, "three", "tri")
        tags = set(m.tags(obj))
        self.assertEqual(set(["tres", "trois"]), tags)
        m.untag(obj)
        tags = set(m.tags(obj))
        self.assertEqual(set(), tags)
        return

    def test_union_and_intersection(self):
        """check TagManager.union() and TagManager.intersection()
        """
        m = self.m
        m.tag(3, "3", "number")
        m.tag(4, "4", "number")
        objs = set([3,4])
        self.assertEqual(m.union(), set())
        self.assertEqual(m.union("number"), objs)
        self.assertEqual(m.union("3"), set([3]))
        self.assertEqual(m.union("3", "4"), objs)
        self.assertRaises(KeyError, m.union, "fail")
        m.silent = True
        self.assertEqual(set(), m.union("fail"))
        self.assertEqual(set([3]), m.union("fail", "3"))
        m.silent = False
        self.assertEqual(m.intersection(), set())
        self.assertEqual(m.intersection("number"), objs)
        self.assertEqual(m.intersection("3"), set([3]))
        self.assertEqual(m.intersection("3", "4"), set())
        self.assertRaises(KeyError, m.intersection, "fail")
        m.silent = True
        self.assertEqual(set(), m.intersection("fail"))
        return

    def test_hasTags(self):
        """check TagManager.hasTags()
        """
        m = self.m
        m.tag(3, "3", "number")
        m.tag(4, "4", "number")
        self.assertTrue( m.hasTags(3, "3") )
        self.assertTrue( m.hasTags(3, "3", "number") )
        self.assertFalse( m.hasTags(3, "3", "4") )
        self.assertRaises(KeyError, m.hasTags, 3, "fail")
        m.silent = True
        self.assertFalse(m.hasTags(3, "fail"))
        return

# End of class TestTagManager

if __name__ == '__main__':
    unittest.main()

# End of file
