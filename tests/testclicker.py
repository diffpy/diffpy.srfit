#!/usr/bin/env python
"""Tests for refinableobj module."""

from diffpy.srfit.equation.clicker import clickerFactory
import unittest


class TestClicker(unittest.TestCase):

    def testClicker(self):
        """Test all aspects of the Clicker."""
        Clicker = clickerFactory()

        c1 = Clicker()
        c2 = Clicker()
        observer = Clicker()
        ref = Clicker()
        ref2 = Clicker()

        # See if these things are at the same spot
        self.assertTrue( c1 == c2 )

        # Click one to make greater than other
        c1.click()
        self.assertTrue( c1 > c2 )
        self.assertTrue( c2 < c1 )
        self.assertFalse( c2 > c1)
        self.assertFalse( c1 < c2)
        self.assertTrue( c1 >= c2 )
        self.assertTrue( c2 <= c1 )
        self.assertFalse( c2 >= c1)
        self.assertFalse( c1 <= c2)

        # Now do the other
        c2.click()
        self.assertTrue( c2 > c1 )
        self.assertTrue( c1 < c2 )
        self.assertFalse( c1 > c2)
        self.assertFalse( c2 < c1)
        self.assertTrue( c2 >= c1 )
        self.assertTrue( c1 <= c2 )
        self.assertFalse( c1 >= c2)
        self.assertFalse( c2 <= c1)

        # Observe these two. Note that this does not change the state of the
        # observer or the subjects.
        observer.addSubject(c1)
        observer.addSubject(c2)
        self.assertTrue( observer < c1 )
        self.assertTrue( observer < c2 )
        self.assertTrue( c2 > c1 )

        # Check relations
        self.assertTrue( c1.hasObserver(observer) )
        self.assertFalse( c1.hasObserver(c2) )
        self.assertTrue( observer.hasSubject(c1) )
        self.assertFalse( c1.hasSubject(c2) )

        # Click a subject
        c1.click()
        self.assertTrue( c1 > c2 )
        self.assertTrue( observer >= c1 )
        self.assertTrue( observer > c2 )

        # Click the observer
        observer.click()
        self.assertTrue( observer > c1 )
        self.assertTrue( observer > c2 )

        # Remove subjects
        observer.removeSubject(c1)
        self.assertFalse( c1.hasObserver(observer) )
        self.assertFalse( observer.hasSubject(c1) )
        c2.removeObserver(observer)
        self.assertFalse( c2.hasObserver(observer) )
        self.assertFalse( observer.hasSubject(c2) )
        return

    def testTwoClikerTypes(self):
        """Test using to clicker of different types."""
        Clicker1 = clickerFactory()
        Clicker2 = clickerFactory()
        c1 = Clicker1()
        c2 = Clicker2()
        c1.click()
        self.assertEqual(0, Clicker2._numclicks)
        self.assertRaises(TypeError, c1.__cmp__, c2)
        self.assertRaises(TypeError, c2.__cmp__, c1)
        return


if __name__ == "__main__":

    unittest.main()

