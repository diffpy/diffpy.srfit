#!/usr/bin/env python
##############################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2008 The Trustees of Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################
"""TagManager class.

The TagManager class takes hashable objects and assigns tags to them.
Objects can then be easily referenced via their assigned tags.
"""

__all__ = ["TagManager"]

import functools


class TagManager(object):
    """TagManager class.

    Manage tags on hashable objects. Tags are strings that carry metadata.

    silent          --  Flag indicating whether to silently pass by when a tag
                        cannot be found (bool, True). If this is False, then a
                        KeyError will be thrown when a tag cannot be found.
    _tagdict        --  A dictionary of tags to sets of tagged objects.
    """

    def __init__(self):
        """Initialization."""
        self._tagdict = {}
        self.silent = True
        return

    def alltags(self):
        """Get all tags managed by the TagManager."""
        return self._tagdict.keys()

    def tag(self, obj, *tags):
        """Tag an object.

        Tags are stored as strings.

        obj     --  Any hashable object to be untagged.
        *tags   --  Tags to apply to obj.

        Raises TypeError if obj is not hashable.
        """
        for tag in tags:
            oset = self._tagdict.setdefault(str(tag), set())
            oset.add(obj)
        return

    def untag(self, obj, *tags):
        """Remove tags from an object.

        obj     --  Any hashable object to be untagged.
        *tags   --  Tags to remove from obj. If this is empty, then all
                    tags will be removed from obj.

        Raises KeyError if a passed tag does not apply to obj and self.silent
        is False
        """
        if not tags:
            tags = self.tags(obj)

        for tag in tags:
            oset = self.__getObjectSet(tag)
            if obj not in oset and not self.silent:
                raise KeyError("Tag '%s' does not apply" % tag)
            oset.discard(obj)

        return

    def tags(self, obj):
        """Get all tags on an object.

        Returns list
        """
        tags = [k for (k, v) in self._tagdict.items() if obj in v]
        return tags

    def hasTags(self, obj, *tags):
        """Determine if an object has all passed tags.

        Returns bool
        """
        setgen = (self.__getObjectSet(t) for t in tags)
        result = all(obj in s for s in setgen)
        return result

    def union(self, *tags):
        """Get all objects that have any of the passed tags.

        Returns set
        """
        if not tags:
            return set()
        setgen = (self.__getObjectSet(t) for t in tags)
        objs = functools.reduce(set.union, setgen)
        return objs

    def intersection(self, *tags):
        """Get all objects that have all of the passed tags.

        Returns set
        """
        if not tags:
            return set()
        setgen = (self.__getObjectSet(t) for t in tags)
        objs = functools.reduce(set.intersection, setgen)
        return objs

    def verifyTags(self, *tags):
        """Check that tags are all extant.

        Raises KeyError if a passed tag does not exist. This ignores
        self.silent.
        """
        keys = self._tagdict.keys()
        for tag in tags:
            if tag not in keys:
                raise KeyError("Tag '%s' does not exist" % tag)
        return True

    def __getObjectSet(self, tag):
        """Helper function for getting an object set with given tag.

        Raises KeyError if a passed tag does not exist and self.silent
        is False
        """
        oset = self._tagdict.get(str(tag))
        if oset is None:
            if not self.silent:
                raise KeyError("Tag '%s' does not exist" % tag)
            oset = set()
        return oset


# End class TagManager

# End of file
