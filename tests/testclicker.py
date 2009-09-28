#!/usr/bin/env python
"""Tests for refinableobj module."""

from diffpy.srfit.util.clicker import Clicker
import unittest


class TestClicker(unittest.TestCase):

    def testClicker(self):
        """Test all aspects of the Clicker."""

        c1 = Clicker()
        c2 = Clicker()
        c3 = Clicker()
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

        # Observe these two. Note that this sets the state of the subject and
        # all observers equal to the global state.
        observer.addSubject(c1)
        observer.addSubject(c2)
        self.assertTrue( observer == c1 )
        self.assertTrue( observer == c2 )
        self.assertTrue( c2 == c1 )

        # Check relations
        self.assertTrue( c1.hasObserver(observer) )
        self.assertFalse( c1.hasObserver(c2) )
        self.assertTrue( observer.hasSubject(c1) )
        self.assertFalse( c1.hasSubject(c2) )

        # Add one more
        c2.addSubject(c3)
        self.assertTrue( observer >= c1 )
        self.assertTrue( c2 >= c1 )
        self.assertTrue( c3 >= c1 )
        self.assertTrue( observer == c2 == c3 )

        # Check relations
        self.assertTrue( c1.hasObserver(observer) )
        self.assertTrue( c2.hasObserver(observer) )
        self.assertTrue( observer.hasSubject(c1) )
        self.assertTrue( observer.hasSubject(c2) )
        self.assertTrue( c2.hasSubject(c3) )
        self.assertTrue( c3.hasObserver(c2) )

        self.assertFalse( c3.hasObserver(c1) )
        self.assertFalse( c3.hasObserver(observer) )
        self.assertFalse( c1.hasSubject(c3) )
        self.assertFalse( observer.hasSubject(c3) )

        self.assertFalse( c1.hasSubject(observer) )
        self.assertFalse( c2.hasSubject(observer) )
        self.assertFalse( observer.hasObserver(c1) )
        self.assertFalse( observer.hasObserver(c2) )
        self.assertFalse( c2.hasObserver(c3) )
        self.assertFalse( c3.hasSubject(c2) )

        # Try to create a loop
        self.assertRaises(ValueError, observer.addObserver, c1)
        self.assertRaises(ValueError, observer.addObserver, c2)
        self.assertRaises(ValueError, c1.addObserver, c1)

        # Try to add None
        self.assertRaises(ValueError, observer.addObserver, None)

        # Click c2
        c2.click()
        self.assertTrue( c2 == observer)
        self.assertTrue( c2 > c1 )
        self.assertTrue( c2 > c3 )

        # Click c3
        c3.click()
        self.assertTrue( observer > c1 )
        self.assertTrue( c2 > c1 )
        self.assertTrue( c3 > c1 )
        self.assertTrue( observer == c2 == c3 )

        # Click c1
        c1.click()
        self.assertTrue( c1 > c2 )
        self.assertTrue( observer >= c1 )
        self.assertTrue( observer > c2 )

        # Click the observer
        observer.click()
        self.assertTrue( observer > c1 )
        self.assertTrue( observer > c2 )

        # Link c3 directly to observer as well
        observer.addSubject(c3)

        # Check relations
        self.assertTrue( c1.hasObserver(observer) )
        self.assertTrue( c2.hasObserver(observer) )
        self.assertTrue( c3.hasObserver(observer) )
        self.assertTrue( observer.hasSubject(c1) )
        self.assertTrue( observer.hasSubject(c2) )
        self.assertTrue( observer.hasSubject(c3) )
        self.assertTrue( c2.hasSubject(c3) )
        self.assertTrue( c3.hasObserver(c2) )
        self.assertTrue( c3.hasObserver(observer) )

        self.assertFalse( c3.hasObserver(c1) )
        self.assertFalse( c1.hasSubject(c3) )

        self.assertFalse( c1.hasSubject(observer) )
        self.assertFalse( c2.hasSubject(observer) )
        self.assertFalse( c3.hasSubject(observer) )
        self.assertFalse( observer.hasObserver(c1) )
        self.assertFalse( observer.hasObserver(c2) )
        self.assertFalse( observer.hasObserver(c3) )
        self.assertFalse( c2.hasObserver(c3) )
        self.assertFalse( c3.hasSubject(c2) )
        self.assertFalse( c3.hasSubject(observer) )

        # Click c2
        c2.click()
        self.assertTrue( c2 == observer)
        self.assertTrue( c2 > c1 )
        self.assertTrue( c2 > c3 )

        # Click c3
        c3.click()
        self.assertTrue( observer > c1 )
        self.assertTrue( c2 > c1 )
        self.assertTrue( c3 > c1 )
        self.assertTrue( observer == c2 == c3 )

        # Click c1
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
        self.assertTrue( observer.hasSubject(c2) )
        self.assertTrue( c2.hasObserver(observer) )
        self.assertTrue( observer.hasSubject(c3) )
        self.assertTrue( c3.hasObserver(observer) )
        c2.removeObserver(observer)
        self.assertFalse( c2.hasObserver(observer) )
        self.assertFalse( observer.hasSubject(c2) )
        self.assertTrue( observer.hasSubject(c3) )
        self.assertTrue( c3.hasObserver(observer) )

        return

if __name__ == "__main__":

    unittest.main()

