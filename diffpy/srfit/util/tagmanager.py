#!/usr/bin/env python
########################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2008 Trustees of the Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
########################################################################
"""TagManager class.

The TagManager class takes hashable objects and assigns tags to them. Objects
can then be easily referenced via their assigned tags.

"""
__all__ = ["TagManager"]

import re

class TagManager(object):
    """TagManager class.

    Manage tags on hashable objects. Tags are strings that carry metadata.

    _tagdict        --  A dictionary of tags to sets of tagged objects.

    """

    def __init__(self):
        """Initialization."""
        self._tagdict = {}
        return

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

        Raises KeyError if a passed tag does not apply to obj.
        Raises KeyError if a passed tag does not exist.

        """
        if not tags:
            tags = self.tags(obj)

        for tag in tags:
            oset = self._tagdict.get(str(tag))
            if oset is None:
                raise KeyError("Tag '%s' does not exist" % tag)
            if obj not in oset:
                raise KeyError("Tag '%s' does not apply" % tag)
            oset.discard(obj)

        return

    def tags(self, obj):
        """Get all tags on an object.

        Returns list
        
        """
        tags = [k for (k, v) in self._tagdict.iteritems() if obj in v]
        return tags

    def objects(self, *tags):
        """Get all objects with passed tags.

        Raises KeyError if a passed tag does not exist.

        Returns set

        """
        objs = set()

        for tag in tags:
            oset = self._tagdict.get(str(tag))
            if oset is None:
                raise KeyError("Tag '%s' does not exist" % tag)
            objs.update(oset)

        return objs

# End class TagManager

# version
__id__ = "$Id$"

#
# End of file
