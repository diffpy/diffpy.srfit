#!/usr/bin/env python
########################################################################
#
# diffpy.srfit      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2009 Trustees of the Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
########################################################################
"""Adapters for pyobjcryst structure objects."""

__all__ = []

import numpy

from pyobjcryst.molecule import GetBondLength, GetBondAngle, GetDihedralAngle
from pyobjcryst.molecule import StretchModeBondLength, StretchModeBondAngle
from pyobjcryst.molecule import StretchModeTorsion

from diffpy.srfit.adapters.nodes import Parameter
from diffpy.srfit.adapters.adaptersmod import ContainerAdapter
from diffpy.srfit.structure.baseadapters import *

class ObjCrystScatteringPowerAdapter(ContainerAdapter):
    """Adapter for ScatteringPower objects.

    Accessors:
    None

    Ignored:
    GetSymbol
    
    """
    ignore = set(["GetSymbol"])

    def _labelself(self):
        obj = self._get()
        return obj.GetName()

# End class ObjCrystScatteringPowerAdapter

class ObjCrystScattererAdapter(BaseScattererAdapter):
    """A base adapter for an ObjCryst Scatterer.

    Base class for objcryst scatterer adapters. In addition to the
    BaseScattererAdapter interface, it provides a hasScatterers method to
    determin if the adapted object is a complex scatterer, such as a molecule.

    """

    def getXYZ(self):
        """Get tuple of adapted fractional coordinates.

        The coordinates must be supplied in order: x, y, z.
        
        """
        return (self.X, self.Y, self.Z)

    def getADPs(self):
        """Get tuple of adapted ADPs. 

        This returns None.
        """
        return None

    def hasScatterers(self):
        """Indicate if this scatterer has its own scatterers."""
        return hasattr(self, "getScatterers")

    def _labelself(self):
        obj = self._get()
        return obj.GetName()

# End class ObjCrystScattererAdapter

class ObjCrystAtomAdapter(ObjCrystScattererAdapter):
    """Adapter for Atom and MolAtom.

    Accessors:
    GetScatteringPower

    Ignored:
    None

    """
    accessors = set(["GetScatteringPower"])

    def getADPs(self):
        """Get tuple of adapted ADPs (B-factors)."""
        sp = self.GetScatteringPower()
        return (sp.B11, sp.B22, sp.B33, sp.B12, sp.B13, sp.B23, sp.Biso)

# End class ObjCrystAtomAdapter

class ObjCrystMoleculeAdapter(ObjCrystScattererAdapter):
    """A adapter for a pyobjcryst.Molecule.

    This provides getScatterers and getLattice, which are required of
    BaseStructureAdapter objects. GetBond, GetBondAngle and GetDihedralAngle
    have been overloaded to return StretchMode parameters.

    Accessors:
    GetScatteringPower
    GetAtom

    Ignored:
    GetNbAtoms, GetNbBonds, GetNbBondAngles, GetNbDihedralAngles,
    GetNbRigidGroups

    """
    accessors = ObjCrystScattererAdapter.accessors
    accessors.update(["GetAtom"])
    ignore = ObjCrystScattererAdapter.ignore
    ignore.update(["GetNbAtoms", "GetNbBonds", "GetNbBondAngles",
    "GetNbDihedralAngles", "GetNbRigidGroups"])

    def getScatterers(self):
        """Get a tuple of adapted scatterers."""
        return tuple(self)

    def getLattice(self):
        """Get the adapted lattice.

        Molecules don't have a lattice. This returns None.

        """
        return None

    def GetBond(self, idx1, idx2):
        """Overloaded to return a parameter representing the bond.

        idx1    --  Index of first atom in the bond.
        idx2    --  Index of second (mutated) atom in the bond.

        Returns ObjCrystBondLengthParameter instance.

        """
        atom1 = self.GetAtom(idx1)
        atom2 = self.GetAtom(idx2)
        name = "%s-%s" % (atom1.name, atom2.name)
        bond = ObjCrystBondLengthParameter(name, atom1, atom2)
        return bond

    def GetBondAngle(self, idx1, idx2, idx3):
        """Overloaded to return a parameter representing the bond angle.

        idx1    --  Index of first atom in the bond angle.
        idx2    --  Index of second atom in the bond angle.
        idx3    --  Index of third (mutated) atom in the bond angle.

        Returns ObjCrystBondAngleParameter instance.

        """
        atom1 = self.GetAtom(idx1)
        atom2 = self.GetAtom(idx2)
        atom3 = self.GetAtom(idx3)
        name = "%s-%s-%s" % (atom1.name, atom2.name, atom3.name)
        bondangle = ObjCrystBondAngleParameter(name, atom1, atom2, atom3)
        return bondangle

    def GetDihedralAngle(self, idx1, idx2, idx3, idx4):
        """Overloaded to return a parameter representing the dihedral angle.

        idx1    --  Index of first atom in the dihedral angle.
        idx2    --  Index of second atom in the dihedral angle.
        idx3    --  Index of third atom in the dihedral angle.
        idx4    --  Index of fourth (mutated) atom in the dihedral angle.

        Returns ObjCrystDihedralAngleParameter instance.

        """
        atom1 = self.GetAtom(idx1)
        atom2 = self.GetAtom(idx2)
        atom3 = self.GetAtom(idx3)
        atom4 = self.GetAtom(idx4)
        name = "%s-%s-%s-%s" % (atom1.name, atom2.name, atom3.name, atom4.name)
        diangle = ObjCrystDihedralAngleParameter(name, atom1, atom2, atom3,
                atom4)
        return diangle

# End class ObjCrystMoleculeAdapter

class ObjCrystCrystalAdapter(BaseStructureAdapter):
    """A adapter for Crystals.

    Accessors:
    GetScatt, GetBondValenceCost

    Ignored:
    GetNbScatterer
   
    """
    accessors = set(["GetScatt", "GetBondValenceCost"])
    ignore = set(["GetNbScatterer"])

    def __init__(self, name, obj, getter, setter):
        """See the adapt function."""
        BaseStructureAdapter.__init__(self, name, obj, getter, setter)
        self.angunits = "rad"
        self._constrainSpaceGroup()
        return

    def _labelself(self):
        obj = self._get()
        return obj.GetName()

    def getLattice(self):
        """Get the adapted lattice.

        Returns self.

        """
        return self

    def getLatPars(self):
        """Get a tuple of the adapted lattice parameters."""
        return (self.a, self.b, self.c, self.alpha, self.beta, self.gamma)

    def getScatterers(self):
        """Get a tuple of adapted scatterers."""
        nb = self.GetNbScatterer()
        return tuple(self.GetScatt(i) for i in range(nb))

    def _constrainSpaceGroup(self):
        """Constrain the space group."""
        stru = self.get()
        sg = self._createSpaceGroup(stru.GetSpaceGroup())
        from diffpy.srfit.structure.sgconstraints import \
                _constrainAsSpaceGroup
        _constrainAsSpaceGroup(self, sg)
        return

    @staticmethod
    def _createSpaceGroup(sgobjcryst):
        """Create a diffpy.Structure.SpaceGroup object from pyobjcryst.

        sgobjcryst  --  A pyobjcryst.spacegroup.SpaceGroup instance.

        This uses the actual space group operations from the
        pyobjcryst.spacegroup.SpaceGroup instance so there is no abiguity about
        the actual space group.

        """
        from diffpy.Structure.SpaceGroups import GetSpaceGroup, SymOp
        name = sgobjcryst.GetName()
        extnstr = ":%s" % sgobjcryst.GetExtension()
        if name.endswith(extnstr):
            name = name[:-len(extnstr)]

        # Get whatever spacegroup we can get by name. This will set the proper
        # crystal system.
        sg = GetSpaceGroup(name)

        # Replace the symmetry operations to guarantee that we get it right.
        symops = sgobjcryst.GetSymmetryOperations()
        tranops = sgobjcryst.GetTranslationVectors()
        sg.symop_list = []

        for trans in tranops:
            for shift, rot in symops:
                tv = trans + shift
                tv -= numpy.floor(tv)
                sg.symop_list.append( SymOp(rot, tv) )

        if sgobjcryst.IsCentrosymmetric():
            center = sgobjcryst.GetInversionCenter()
            for trans in tranops:
                for shift, rot in symops:
                    tv = center - trans - shift
                    tv -= numpy.floor(tv)
                    sg.symop_list.append( SymOp(-rot , tv) )

        return sg

# End class ObjCrystCrystalAdapter

class StretchModeParameter(Parameter):
    """Partial Parameter class encapsulating pyobjcryst stretch modes.

    This class relies upon attributes that do not belong to it. Do not
    instantiate this class.

    Required attributes:
    mode        --  The pyobjcryst.molecule.StretchMode used to change atomic
                    positions.
    keepcenter  --  Flag indicating whether to keep the center of mass of the
                    molecule stationary within the crystal when changing the
                    value of the parameter (bool, default True).

    """

    def __init__(self, name, value = None):
        """Initialization.
        
        name    --  Name of the parameter
        value   --  Value of the parameter

        """
        self.keepcenter = True
        Parameter.__init__(self, name, value)
        return

    def _view(self, *atoms):
        """Mutually view the added atoms."""
        for a in atoms:
            self._addViewer(a.X)
            self._addViewer(a.Y)
            self._addViewer(a.Z)
            a._addViewer(self)
        return

    def _set(self, val):
        """Set the parameter's value."""
        if val is None: return
        delta = val - self._get()
        if delta == 0: return
        self.mode.Stretch(delta, self.keepcenter)
        return

    def addAtoms(self, *atoms):
        """Associate adapted atoms with the Parameter.

        Associated atoms will be mutated in the exact same way as the primary
        mutated atoms.  This is useful when a group of atoms should move
        rigidly in response to a change in a bond property.

        """
        # Make sure we're viewing these atoms and they are viewing us.
        self._view(*atoms)
        scatlist = [a.get() for a in atoms]
        self.mode.AddAtoms(scatlist)
        return

# End class StretchModeParameter

class ObjCrystBondLengthParameter(StretchModeParameter):
    """Class for abstracting a bond length in a Molecule to a Parameter.

    This wraps up a pyobjcryst.molecule.StretchModeBondLength object so that
    the distance between two MolAtoms in a Molecule can be used as an
    adjustable Parameter. When a bond length is adjusted, the second MolAtom is
    moved, and the absolute position of the Molecule is altered to preserve the
    location of the center of mass within the Crystal. Thus, the x, y and z
    Parameters of the MolAtom and its parent Molecule are altered. This can be
    changed by setting the 'keepcenter' attribute of the parameter to False.

    This Parameter makes it possible to mutate a MolAtom multiple times in a
    single refinement step. If these mutations are not orthogonal, then this
    could lead to nonconvergence of a fit, depending on the optimizer. Consider
    mutating atom2 of a bond directly, and via a ObjCrystBondLengthParameter.
    The two mutations of atom2 may be determined independently by the
    optimizer, in which case the composed mutation will have an unexpected
    effect on the residual. It is best practice to either modify MolAtom
    positions directly, or thorough BondLengthParameters, BondAngleParameters
    and DihedralAngleParameters (which are mutually orthogonal).
    
    Attributes:
    atom1   --  The first ObjCrystAtomAdapter in the bond
    atom2   --  The second (mutated) ObjCrystAtomAdapter in the bond
    mode    --  The pyobjcryst.molecule.StretchModeBondLength for the bond

    Inherited Attributes
    name    --  A name for this parameter.
    value   --  The parameter value, property.

    """

    def __init__(self, name, atom1, atom2, value = None):
        """Create a ObjCrystBondLengthParameter.

        name    --  The name of the ObjCrystBondLengthParameter
        atom1   --  The first atom in the bond
        atom2   --  The second (mutated) atom in the bond
        value   --  An initial value for the bond length. If this is None
                    (default), then the current distance between the atoms will
                    be used.

        """

        # Create the mode
        self.mode = StretchModeBondLength(atom1.get(), atom2.get(), None)
        # We only add the last atom. This is the one that will move
        self.atom1 = atom1
        self.atom2 = atom2

        # We do this last so the atoms are defined before we set any values.
        StretchModeParameter.__init__(self, name, value)
        atom1._addViewer(self)
        self.addAtoms(atom2)
        return

    def _get(self):
        """Calculate the parameter value."""
        val = GetBondLength(self.atom1.get(), self.atom2.get())
        return val

# End class ObjCrystBondLengthParameter
                     
class ObjCrystBondAngleParameter(StretchModeParameter):
    """Class for abstracting a bond angle in a Molecule to a Parameter.

    This wraps up a pyobjcryst.molecule.StretchModeBondAngle object so that the
    angle defined by three MolAtoms in a Molecule can be used as an adjustable
    Parameter. When a bond angle is adjusted, the third MolAtom is moved, and
    the absolute position of the Molecule is altered to preserve the location
    of the center of mass within the crystal. This can be changed by setting
    the 'keepcenter' attribute of the parameter to False.

    See precautions in the ObjCrystBondLengthParameter class.

    Attributes
    atom1   --  The first ObjCrystAtomAdapter in the bond angle
    atom2   --  The second (central) ObjCrystAtomAdapter in the bond angle
    atom3   --  The third (mutated) ObjCrystAtomAdapter in the bond angle
    mode    --  The pyobjcryst.molecule.StretchModeBondAngle for the bond angle

    Inherited Attributes
    name    --  A name for this parameter.
    value   --  The parameter value, property.

    """

    def __init__(self, name, atom1, atom2, atom3, value = None):
        """Create a ObjCrystBondAngleParameter.

        name    --  The name of the ObjCrystBondAngleParameter.
        atom1   --  The first atom (ObjCrystAtomAdapter) in the bond angle
        atom2   --  The second (central) atom (ObjCrystAtomAdapter) in the
                    bond angle
        atom3   --  The third (mutated) atom (ObjCrystAtomAdapter) in the
                    bond angle
        value   --  An initial value for the bond length. If this is None
                    (default), then the current bond angle between the atoms
                    will be used.
        """

        # Create the stretch mode
        self.mode = StretchModeBondAngle(atom1.get(), atom2.get(), atom3.get(),
                None)
        # We only add the last atom. This is the one that will move
        self.mode.AddAtom(atom3.get())
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3

        # We do this last so the atoms are defined before we set any values.
        StretchModeParameter.__init__(self, name, value)
        atom1._addViewer(self)
        atom2._addViewer(self)
        self.addAtoms(atom3)
        return

    def _get(self):
        """Calculate the parameter value."""
        val = GetBondAngle(self.atom1.get(), self.atom2.get(),
                self.atom3.get())
        return val

# End class ObjCrystBondAngleParameter
                     
class ObjCrystDihedralAngleParameter(StretchModeParameter):
    """Class for abstracting a dihedral angle in a Molecule to a Parameter.

    This wraps up a pyobjcryst.molecule.StretchModeTorsion object so that the
    angle defined by four MolAtoms ([a1-a2].[a3-a4]) in a Molecule can be used
    as an adjustable parameter. When a dihedral angle is adjusted, the fourth
    MolAtom is moved, and the absolute position of the Molecule is altered to
    preserve the location of the center of mass within the crystal.  This can
    be changed by setting the 'keepcenter' attribute of the parameter to False.

    See precautions in the ObjCrystBondLengthParameter class.

    Attributes
    atom1   --  The first ObjCrystAtomAdapter in the dihedral angle
    atom2   --  The second (central) ObjCrystAtomAdapter in the dihedral
                angle
    atom3   --  The third (central) ObjCrystAtomAdapter in the dihedral angle
    atom4   --  The fourth (mutated) ObjCrystAtomAdapter in the dihedral
                angle

    Inherited Attributes
    name    --  A name for this parameter.
    value   --  The parameter value, property.

    """

    def __init__(self, name, atom1, atom2, atom3, atom4, value = None):
        """Create a ObjCrystDihedralAngleParameter.

        name    --  The name of the ObjCrystDihedralAngleParameter
        atom1   --  The first atom (ObjCrystAtomAdapter) in the dihderal
                    angle
        atom2   --  The second (central) atom (ObjCrystAtomAdapter) in the
                    dihderal angle
        atom3   --  The third (central) atom (ObjCrystAtomAdapter) in the
                    dihderal angle
        atom4   --  The fourth (mutated) atom (ObjCrystAtomAdapter) in the
                    dihderal angle
        value   --  An initial value for the bond length. If this is None
                    (default), then the current dihedral angle between atoms
                    will be used.

        """

        # Create the stretch mode
        self.mode = StretchModeTorsion(atom2.get(), atom3.get(), None)
        # We only add the last atom. This is the one that will move
        self.mode.AddAtom(atom4.get())
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4

        # We do this last so the atoms are defined before we set any values.
        StretchModeParameter.__init__(self, name, value)
        atom1._addViewer(self)
        atom2._addViewer(self)
        atom3._addViewer(self)
        self.addAtoms(atom4)
        return

    def _get(self):
        """Calculate the parameter value."""
        val = GetDihedralAngle(self.atom1.get(), self.atom2.get(),
                self.atom3.get(), self.atom4.get())
        return val

# End class ObjCrystDihedralAngleParameter

try:
    from diffpy.srfit.adapters.adaptersmod import registry
    import pyobjcryst
    registry[pyobjcryst.scatteringpower.ScatteringPower] =\
        ObjCrystScatteringPowerAdapter
    registry[pyobjcryst.atom.Atom] = ObjCrystAtomAdapter
    registry[pyobjcryst.molecule.Molecule] = ObjCrystMoleculeAdapter
    registry[pyobjcryst.molecule.MolAtom] = ObjCrystAtomAdapter
    registry[pyobjcryst.crystal.Crystal] = ObjCrystCrystalAdapter
except ImportError:
    pass


__id__ = "$Id$"
