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

from diffpy.utils._deprecator import build_deprecation_message, deprecated

tagmanager_base = "diffpy.srfit.util.tagmanager.TagManager"
removal_version = "4.0.0"
hasTags_dep_msg = build_deprecation_message(
    tagmanager_base, "hasTags", "has_tags", removal_version
)
verifyTags_dep_msg = build_deprecation_message(
    tagmanager_base, "verifyTags", "verify_tags", removal_version
)


class TagManager(object):
    """TagManager class.

    Manage tags on hashable objects. Tags are strings that carry metadata.

    Attributes
    ----------
    silent : bool
        The flag that is True (default) to treat an unknown tag as having
        no objects. If False, an unknown tag raises a KeyError.
    _tagdict : dict
        The mapping of each tag to the set of objects it is applied to.
    """

    def __init__(self):
        self._tagdict = {}
        self.silent = True
        return

    def alltags(self):
        """Return all tags managed by the TagManager.

        Returns
        -------
        dict_keys
            The tags that have been applied to any object.
        """
        return self._tagdict.keys()

    def tag(self, obj, *tags):
        """Tag an object.

        Tags are stored as strings.

        Parameters
        ----------
        obj : hashable
            The object to tag.
        *tags
            The tags to apply to ``obj``.

        Raises
        ------
        TypeError
            If ``obj`` is not hashable.
        """
        for tag in tags:
            oset = self._tagdict.setdefault(str(tag), set())
            oset.add(obj)
        return

    def untag(self, obj, *tags):
        """Remove tags from an object.

        Parameters
        ----------
        obj : hashable
            The object to untag.
        *tags
            The tags to remove from ``obj``. If none are given, all tags
            are removed from ``obj``.

        Raises
        ------
        KeyError
            If a given tag does not apply to ``obj`` and ``silent`` is
            False.
        """
        if not tags:
            tags = self.tags(obj)

        for tag in tags:
            oset = self.__get_object_set(tag)
            if obj not in oset and not self.silent:
                raise KeyError("Tag '%s' does not apply" % tag)
            oset.discard(obj)

        return

    def tags(self, obj):
        """Return all tags on an object.

        Parameters
        ----------
        obj : hashable
            The object to look up.

        Returns
        -------
        list of str
            The tags applied to ``obj``.
        """
        tags = [k for (k, v) in self._tagdict.items() if obj in v]
        return tags

    def has_tags(self, obj, *tags):
        """Check whether an object has all of the given tags.

        Parameters
        ----------
        obj : hashable
            The object to check.
        *tags
            The tags to look for.

        Returns
        -------
        bool
            The flag that is True if ``obj`` has every tag in ``tags``.

        Raises
        ------
        KeyError
            If a tag does not exist and ``silent`` is False.
        """
        setgen = (self.__get_object_set(t) for t in tags)
        result = all(obj in s for s in setgen)
        return result

    def union(self, *tags):
        """Return all objects that have any of the given tags.

        Parameters
        ----------
        *tags
            The tags to look for.

        Returns
        -------
        set
            The objects that have at least one tag in ``tags``. Empty if
            no tags are given.

        Raises
        ------
        KeyError
            If a tag does not exist and ``silent`` is False.
        """
        if not tags:
            return set()
        setgen = (self.__get_object_set(t) for t in tags)
        objs = functools.reduce(set.union, setgen)
        return objs

    def intersection(self, *tags):
        """Return all objects that have all of the given tags.

        Parameters
        ----------
        *tags
            The tags to look for.

        Returns
        -------
        set
            The objects that have every tag in ``tags``. Empty if no tags
            are given.

        Raises
        ------
        KeyError
            If a tag does not exist and ``silent`` is False.
        """
        if not tags:
            return set()
        setgen = (self.__get_object_set(t) for t in tags)
        objs = functools.reduce(set.intersection, setgen)
        return objs

    def verify_tags(self, *tags):
        """Check that all of the given tags exist.

        This ignores ``silent``.

        Parameters
        ----------
        *tags
            The tags to check.

        Returns
        -------
        bool
            True when every tag exists.

        Raises
        ------
        KeyError
            If a tag does not exist.
        """
        keys = self._tagdict.keys()
        for tag in tags:
            if tag not in keys:
                raise KeyError("Tag '%s' does not exist" % tag)
        return True

    @deprecated(hasTags_dep_msg)
    def hasTags(self, obj, *tags):
        """This function has been deprecated and will be removed in
        version 4.0.0.

        Please use diffpy.srfit.util.tagmanager.TagManager.has_tags
        instead.
        """
        return self.has_tags(obj, *tags)

    @deprecated(verifyTags_dep_msg)
    def verifyTags(self, *tags):
        """This function has been deprecated and will be removed in
        version 4.0.0.

        Please use diffpy.srfit.util.tagmanager.TagManager.verify_tags
        instead.
        """
        return self.verify_tags(*tags)

    def __get_object_set(self, tag):
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
