#!/usr/bin/env python
########################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2011 Trustees of the Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
########################################################################


from diffpy.srfit.adapters.adaptersmod import registry, attrgetter
from diffpy.srfit.adapters.adaptersmod import ContainerAdapter
from diffpy.srfit.structure.baseadapters import\
        BaseScattererAdapter, BaseStructureAdapter, BaseLatticeAdapter

__all__ = []

class DiffpyAtomAdapter(BaseScattererAdapter):
    """Adapter for diffpy.Structure.Atom instances.

    Ignored:
    element, anisotropy
    
    """
    ignore = set(["element", "anisotropy"])

    def getXYZ(self):
        """Get tuple of adapted fractional coordinates.

        The coordinates must be supplied in order: x, y, z.
        
        """
        return (self.x, self.y, self.z)

    def getADPs(self):
        """Get tuple of adapted ADPs. 
        
        The ADPs must be supplied in order: 11, 22, 33, 12, 13, 23, iso.
        Whether the ADPs are U- or B-factors is dependent on the adapted
        scatterer. If any of these cannot be supplied, supply None instead. If
        none of these can be supplied, return None.

        """
        return (self.U11, self.U22, self.U33, self.U12, self.U13, self.U23,
                self.Uisoequiv)

    def _labelself(self):
        """Get a label for self."""
        atom = self._get()
        return atom.name or atom.element.title()

# End class DiffpyAtomAdapter

class DiffpyLatticeAdapter(BaseLatticeAdapter):
    """Adapter for diffpy Lattice instances."""
    angunits = "deg"

    def getLatPars(self):
        """Get a tuple of the adapted lattice parameters."""
        return (self.a, self.b, self.c, self.alpha, self.beta, self.gamma)

# End class DiffpyLatticeAdapter

class DiffpyStructureAdapter(BaseStructureAdapter):
    """Adapter for diffpy.Structure.Structure instances.

    getAtom is overloaded to return the approprate parameter.

    Ignored:
    getLabels
    
    """
    ignore = set(["getLabels"])

    def __init__(self, name, obj, getter, setter):
        """See the adapt function."""
        BaseStructureAdapter.__init__(self, name, obj, getter, setter)
        stru = self.get()
        # FIXME: is lookup by labels necessary here?
        hasuniquelabels = (len(stru) == len(set(stru.label)))
        if not hasuniquelabels:
            stru.assignUniqueLabels()
        self.__labels = stru.label
        # Reverse lookup for atom index based on label
        self.__lindex = dict((k, v) for (v, k) in enumerate(self.__labels))
        return

    def getAtom(self, id):
        """See Structure class.
        
        This is overloaded to return the same adapters as those accessed with
        list indexing.

        """
        try:
            if type(id) is not int:
                id = self.__lindex[id]
            rv = self[id]
            return rv
        except (IndexError, KeyError):
            emsg = "Invalid atom identifier %r." % id
            raise ValueError(emsg)
        return

    def getLattice(self):
        """Get the adapted lattice."""
        return self.lattice

    def getScatterers(self):
        """Get a tuple of adapted scatterers."""
        return tuple(self)

    def _labelitem(self, idx):
        """Label adapted atoms

        This uses getLabels to create the atom labels.
        
        """
        return self.__labels[idx]

# End class DiffpyStructureAdapter

import diffpy.Structure
registry[diffpy.Structure.Atom] = DiffpyAtomAdapter
registry[diffpy.Structure.Lattice] = DiffpyLatticeAdapter
registry[diffpy.Structure.Structure] = DiffpyStructureAdapter
