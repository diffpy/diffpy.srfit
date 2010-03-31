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
"""Wrappers for adapting pyobjcryst.crystal.Crystal to a srfit ParameterSet.

This will adapt a pyobjcryst.Crystal into the ParameterSet interface. The
following classes are adapted.

ObjCrystParSet  --  Adapter for pyobjcryst.crystal.Crystal 
AtomParSet      --  Adapter for pyobjcryst.atom.Atom    
MoleculeParSet  --  Adapter for pyobjcryst.molecule.Molecule
MolAtomParSet   --  Adapter for pyobjcryst.molecule.MolAtom 

Related to the adaptation of Molecule and MolAtom, there are adaptors for
specifying molecule restraints.
BondLengthRestraint
BondAngleRestraint
DihedralAngleRestraint

There are also Parameters for encapsulating and modifying atoms via their
relative positions. These Parameters can also act like constraints, and can
modify the positions of multiple MolAtoms.
BondLengthParameter
BondAngleParameter
DihedralAngleParameter

"""

__all__ = ["ScattererParSet", "AtomParSet", "MoleculeParSet", "MolAtomParSet"
        "MoleculeRestraint", "BondLengthRestraint", "BondAngleRestraint",
        "DihedralAngleRestraint", "StretchModeParameter",
        "BondLengthParameter", "BondAngleParameter", "DihedralAngleParameter",
        "ObjCrystParSet"]

from pyobjcryst.molecule import GetBondLength, GetBondAngle, GetDihedralAngle
from pyobjcryst.molecule import StretchModeBondLength, StretchModeBondAngle
from pyobjcryst.molecule import StretchModeTorsion

from diffpy.srfit.fitbase.parameter import Parameter, ParameterAdapter
from diffpy.srfit.fitbase.parameter import ParameterProxy
from diffpy.srfit.fitbase.parameterset import ParameterSet
from diffpy.srfit.fitbase.restraint import Restraint
from diffpy.srfit.structure.basestructure import BaseStructure

class ScattererParSet(ParameterSet):
    """A base adaptor for an Objcryst Scatterer.

    This class derives from diffpy.srfit.fitbase.parameterset.ParameterSet and
    adapts pyobjcryst.scatterer.Scatterer derivatives (Molecule, Atom) and
    objects with a similar interface (MolAtom).  See the ParameterSet class for
    base attributes.

    Attributes:
    scat        --  The adapted pyobjcryst object.
    parent      --  The ParameterSet this belongs to

    Managed Parameters:
    x, y, z     --  Scatterer position in crystal coordinates (ParameterWraper)
    occ         --  Occupancy of the scatterer on its crystal site
                    (ParameterWraper)
    
    """

    def __init__(self, name, scat, parent):
        """Initialize

        name    --  The name of the scatterer
        scat    --  The pyobjcryst.Scatterer instance
        parent  --  The ParameterSet this belongs to

        """
        ParameterSet.__init__(self, name)
        self.scat = scat
        self.parent = parent

        # x, y, z, occ
        self.addParameter(ParameterAdapter("x", self.scat, attr = "X"))
        self.addParameter(ParameterAdapter("y", self.scat, attr = "Y"))
        self.addParameter(ParameterAdapter("z", self.scat, attr = "Z"))
        self.addParameter(ParameterAdapter("occ", self.scat, attr =
            "Occupancy"))
        return

    def isDummy(self):
        """Indicate whether this scatterer is a dummy atom."""
        return False

    def hasScatterers(self):
        """Indicate if this scatterer has its own scatterers."""
        return hasattr(self, "getScatterers")


# End class ScattererParSet

class AtomParSet(ScattererParSet):
    """A adaptor for a pyobjcryst.Atom.

    This class derives from ScattererParSet.

    Attributes:
    scat        --  The adapted pyobjcryst.atom.Atom.
    element     --  Non-refinable name of the element (property).
    parent      --  The ObjCrystParSet this belongs to.

    Managed Parameters:
    x, y, z     --  Atom position in crystal coordinates (ParameterAdapter)
    occ         --  Occupancy of the atom on its crystal location
                    (ParameterAdapter)
    Biso        --  Isotropic scattering factor (ParameterAdapter).
    B11, B22, B33, B12, B21, B23, B32, B13, B31
                --  Anisotropic displacement factor for scatterer
                (ParameterAdapter or ParameterProxy). Note that the Bij and Bji
                parameters are the same.
    
    """

    def __init__(self, name, atom, parent):
        """Initialize

        name    --  The name of the scatterer
        scat    --  The Scatterer instance
        parent  --  The ObjCrystParSet this belongs to

        """
        ScattererParSet.__init__(self, name, atom, parent)
        sp = atom.GetScatteringPower()

        # The B-parameters
        self.addParameter(ParameterAdapter("Biso", sp, attr = "Biso"))
        self.addParameter(ParameterAdapter("B11", sp, attr = "B11"))
        self.addParameter(ParameterAdapter("B22", sp, attr = "B22"))
        self.addParameter(ParameterAdapter("B33", sp, attr = "B33"))
        B12 = ParameterAdapter("B12", sp, attr = "B12")
        B21 = ParameterProxy("B21", B12)
        B13 = ParameterAdapter("B13", sp, attr = "B13")
        B31 = ParameterProxy("B31", B13)
        B23 = ParameterAdapter("B23", sp, attr = "B23")
        B32 = ParameterProxy("B32", B23)
        self.addParameter(B12)
        self.addParameter(B21)
        self.addParameter(B13)
        self.addParameter(B31)
        self.addParameter(B23)
        self.addParameter(B32)
        return

    def _getElem(self):
        """Getter for the element type."""
        return self.scat.GetScatteringPower().GetSymbol()

    element = property(_getElem)

# End class AtomParSet

class MoleculeParSet(ScattererParSet):
    """A adaptor for a pyobjcryst.Molecule.

    This class derives from ScattererParSet.

    Attributes:
    scat        --  The adapted pyobjcryst.molecule.Molecule.
    parent      --  The ObjCrystParSet this belongs to.

    Managed Parameters:
    x, y, z     --  Molecule position in crystal coordinates (ParameterAdapter)
    occ         --  Occupancy of the molecule on its crystal location
                    (ParameterAdapter)
    q0, q1, q2, q3  --  Orientational quaternion (ParameterAdapter)
    
    Other attributes are inherited from
    diffpy.srfit.fitbase.parameterset.ParameterSet

    """

    def __init__(self, name, molecule, parent):
        """Initialize

        name    --  The name of the scatterer
        molecule    --  The pyobjcryst.Molecule instance
        parent  --  The ObjCrystParSet this belongs to

        """
        ScattererParSet.__init__(self, name, molecule, parent)

        # Add orientiation quaternion
        self.addParameter(ParameterAdapter("q0", self.scat, attr = "Q0"))
        self.addParameter(ParameterAdapter("q1", self.scat, attr = "Q1"))
        self.addParameter(ParameterAdapter("q2", self.scat, attr = "Q2"))
        self.addParameter(ParameterAdapter("q3", self.scat, attr = "Q3"))

        # Wrap the MolAtoms within the molecule
        self.atoms = []
        anames = []

        for a in molecule:

            name = a.GetName()
            if not name:
                raise AttributeError("Each MolAtom must have a name")
            if name in anames:
                raise AttributeError("MolAtom name '%s' is duplicated"%name)

            atom = MolAtomParSet(name, a, self)
            atom.molecule = self
            self.addParameterSet(atom)
            self.atoms.append(atom)
            anames.append(name)

        return

    def getScatterers(self):
        """Get a list of ParameterSets that represents the scatterers.

        The site positions must be accessible from the list entries via the
        names "x", "y", and "z". The ADPs must be accessible as well, but the
        name and nature of the ADPs (U-factors, B-factors, isotropic,
        anisotropic) depends on the adapted structure.

        """
        return self.atoms

    def wrapRestraints(self):
        """Wrap the restraints implicit to the molecule.

        This will wrap MolBonds, MolBondAngles and MolDihedralAngles of the
        Molecule as MoleculeRestraint objects.

        """
        # Wrap restraints. Restraints wrapped in this way cannot be modified
        # from within this class.
        for b in self.scat.GetBondList():
            res = MoleculeRestraint(b)
            self._restraints.add(res)

        for ba in self.scat.GetBondAngleList():
            res = MoleculeRestraint(ba)
            self._restraints.add(res)

        for da in self.scat.GetDihedralAngleList():
            res = MoleculeRestraint(da)
            self._restraints.add(res)

        return

    def wrapStretchModeParameters(self):
        """Wrap the stretch modes implicit to the Molecule as Parameters.

        This will wrap StretchModeBondLengths and StretchModeBondAngles of the
        Molecule as Parameters. Note that this requires that the MolBondAtoms
        in the Molecule came in with unique names.  Torsion angles are not
        wrapped, as there is not enough information to determine each MolAtom
        in the angle.

        The Parameters will be given the concatenated name of its constituents.
        bond lengths: "bl_aname1_aname2"
        bond angles: "ba_aname1_aname2_aname3"

        """
        for mode in self.scat.GetStretchModeBondLengthList():
            name1 = mode.mpAtom0.GetName()
            name2 = mode.mpAtom1.GetName()

            name = "bl_" + "_".join((name1, name2))

            atom1 = getattr(self, name1)
            atom2 = getattr(self, name2)

            par = BondLengthParameter(name, atom1, atom2, mode = mode)

            atoms = []
            for a in mode.GetAtoms():
                name = a.GetName()
                atoms.append( getattr(self, name) )

            par.AddAtoms(atoms)

            self.addParameter(par)


        for mode in self.scat.GetStretchModeBondAngleList():
            name1 = mode.mpAtom0.GetName()
            name2 = mode.mpAtom1.GetName()
            name3 = mode.mpAtom2.GetName()

            name = "ba_" + "_".join((name1, name2, name3))

            atom1 = getattr(self, name1)
            atom2 = getattr(self, name2)
            atom3 = getattr(self, name3)

            par = BondAngleParameter(name, atom1, atom2, atom3, mode = mode)

            atoms = []
            for a in mode.GetAtoms():
                name = a.GetName()
                atoms.append( getattr(self, name) )
            par.AddAtoms(atoms)

            self.addParameter(par)

        return

    def restrainBondLength(self, atom1, atom2, length, sigma, delta, scaled =
            False):
        """Add a bond length restraint.

        This creates an instance of BondLengthRestraint and adds it to the
        MoleculeParSet.

        atom1   --  First atom (MolAtomParSet) in the bond
        atom2   --  Second atom (MolAtomParSet) in the bond
        length  --  The length of the bond (Angstroms)
        sigma   --  The uncertainty of the bond length (Angstroms)
        delta   --  The width of the bond (Angstroms)
        scaled  --  A flag indicating if the restraint is scaled (multiplied)
                    by the unrestrained point-average chi^2 (chi^2/numpoints)
                    (default False)

        Returns the BondLengthRestraint object for use with the 'unrestrain'
        method.

        """
        res = BondLengthRestraint(atom1, atom2, length, sigma, delta, scaled)
        self._restraints.add(res)

        return res

    def restrainBondLengthParameter(self, par, length, sigma, delta, scaled =
            False):
        """Add a bond length restraint.

        This creates an instance of BondLengthRestraint and adds it to the
        MoleculeParSet.

        par     --  A BondLengthParameter (see addBondLengthParameter)
        length  --  The length of the bond (Angstroms)
        sigma   --  The uncertainty of the bond length (Angstroms)
        delta   --  The width of the bond (Angstroms)
        scaled  --  A flag indicating if the restraint is scaled (multiplied)
                    by the unrestrained point-average chi^2 (chi^2/numpoints)
                    (default False)

        Returns the BondLengthRestraint object for use with the 'unrestrain'
        method.

        """
        return self.restrainBondLength(par.atom1, par.atom2, length, sigma,
                delta, scaled)

    def restrainBondAngle(self, atom1, atom2, atom3, angle, sigma, delta,
            scaled = False):
        """Add a bond angle restraint.

        This creates an instance of BondAngleRestraint and adds it to the
        MoleculeParSet.

        atom1   --  First atom (MolAtomParSet) in the bond angle
        atom2   --  Second (central) atom (MolAtomParSet) in the bond angle
        atom3   --  Third atom (MolAtomParSet) in the bond angle
        angle   --  The bond angle (radians)
        sigma   --  The uncertainty of the bond angle (radians)
        delta   --  The width of the bond angle (radians)
        scaled  --  A flag indicating if the restraint is scaled (multiplied)
                    by the unrestrained point-average chi^2 (chi^2/numpoints)
                    (default False).

        Returns the BondAngleRestraint object for use with the 'unrestrain'
        method.

        """
        res = BondAngleRestraint(atom1, atom2, atom3, angle, sigma, delta,
                scaled)
        self._restraints.add(res)

        return res

    def restrainBondAngleParameter(self, par, angle, sigma, delta,
            scaled = False):
        """Add a bond angle restraint.

        This creates an instance of BondAngleRestraint and adds it to the
        MoleculeParSet.

        par     --  A BondAngleParameter (see addBondAngleParameter)
        angle   --  The bond angle (radians)
        sigma   --  The uncertainty of the bond angle (radians)
        delta   --  The width of the bond angle (radians)
        scaled  --  A flag indicating if the restraint is scaled (multiplied)
                    by the unrestrained point-average chi^2 (chi^2/numpoints)
                    (default False).

        Returns the BondAngleRestraint object for use with the 'unrestrain'
        method.

        """
        return self.restrainBondAngle(par.atom1, par.atom2, par.atom3, angle,
                sigma, delta, scaled)

    def restrainDihedralAngle(self, atom1, atom2, atom3, atom4, angle, sigma,
            delta, scaled = False):
        """Add a dihedral angle restraint.

        This creates an instance of DihedralAngleRestraint and adds it to the
        MoleculeParSet.

        atom1   --  First atom (MolAtomParSet) in the angle
        atom2   --  Second (central) atom (MolAtomParSet) in the angle
        atom3   --  Third (central) atom (MolAtomParSet) in the angle
        atom4   --  Fourth atom in the angle (MolAtomParSet)
        angle   --  The dihedral angle (radians)
        sigma   --  The uncertainty of the dihedral angle (radians)
        delta   --  The width of the dihedral angle (radians)
        scaled  --  A flag indicating if the restraint is scaled (multiplied)
                    by the unrestrained point-average chi^2 (chi^2/numpoints)
                    (default False).

        Returns the DihedralAngleRestraint object for use with the 'unrestrain'
        method.

        """
        res = DihedralAngleRestraint(atom1, atom2, atom3, atom4, angle, sigma,
                delta, scaled)
        self._restraints.add(res)

        return res

    def restrainDihedralAngleParameter(self, par, angle, sigma, delta,
            scaled = False):
        """Add a dihedral angle restraint.

        This creates an instance of DihedralAngleRestraint and adds it to the
        MoleculeParSet.

        par     --  A DihedralAngleParameter (see addDihedralAngleParameter)
        angle   --  The dihedral angle (radians)
        sigma   --  The uncertainty of the dihedral angle (radians)
        delta   --  The width of the dihedral angle (radians)
        scaled  --  A flag indicating if the restraint is scaled (multiplied)
                    by the unrestrained point-average chi^2 (chi^2/numpoints)
                    (default False).

        Returns the DihedralAngleRestraint object for use with the 'unrestrain'
        method.

        """
        return self.restrainDihedralAngle(par.atom1, par.atom2, par.atom3,
                par.atom4, angle, sigma, delta, scaled)

    def addBondLengthParameter(self, name, atom1, atom2, value = None, const =
            False):
        """Add a bond length to the Molecule.

        This creates a BondLengthParameter to the MoleculeParSet that can be
        adjusted during the fit.

        name    --  The name of the BondLengthParameter
        atom1   --  The first atom (MolAtomParSet) in the bond
        atom2   --  The second (mutated) atom (MolAtomParSet) in the bond
        value   --  An initial value for the bond length. If this is None
                    (default), then the current distance between the atoms will
                    be used.
        const   --  A flag indicating whether the Parameter is constant
                    (default False)

        Returns the new BondLengthParameter.

        """
        par = BondLengthParameter(name, atom1, atom2, value, const)
        self.addParameter(par)

        return par

    def addBondAngleParameter(self, name, atom1, atom2, atom3, value = None,
            const = False):
        """Add a bond angle to the Molecule.

        This creates a BondAngleParameter to the MoleculeParSet that can be
        adjusted during the fit.

        name    --  The name of the BondAngleParameter
        atom1   --  The first atom (MolAtomParSet) in the bond angle
        atom2   --  The second (central) atom (MolAtomParSet) in the bond angle
        atom3   --  The third (mutated) atom (MolAtomParSet) in the bond angle
        value   --  An initial value for the bond angle. If this is None
                    (default), then the current bond angle between the atoms
                    will be used.
        const   --  A flag indicating whether the Parameter is constant
                    (default False).

        Returns the new BondAngleParameter.

        """
        par = BondAngleParameter(name, atom1, atom2, atom3, value, const)
        self.addParameter(par)

        return par

    def addDihedralAngleParameter(self, name, atom1, atom2, atom3, atom4, value
            = None, const = False):
        """Add a dihedral angle to the Molecule.

        This creates a DihedralAngleParameter to the MoleculeParSet that can be
        adjusted during the fit.

        name    --  The name of the DihedralAngleParameter.
        atom1   --  The first atom (MolAtomParSet) in the dihderal angle
        atom2   --  The second (central) atom (MolAtomParSet) in the dihderal
                    angle
        atom3   --  The third (central) atom (MolAtomParSet) in the dihderal
                    angle
        atom4   --  The fourth (mutated) atom (MolAtomParSet) in the dihderal
                    angle
        value   --  An initial value for the dihedral angle. If this is None
                    (default), then the current dihedral angle between atoms
                    will be used.
        const   --  A flag indicating whether the Parameter is constant
                    (default False).

        Returns the new DihedralAngleParameter.

        """
        par = DihedralAngleParameter(name, atom1, atom2, atom3, atom4, value,
                const)
        self.addParameter(par)

        return par

# End class MoleculeParSet

class MolAtomParSet(ScattererParSet):
    """A adaptor for an pyobjcryst.molecule.MolAtom.

    This class derives from srfit.fitbase.parameterset.ParameterSet. Note that
    MolAtom does not derive from Scatterer, but the relevant interface is the
    same within pyobjcryst. See the ParameterSet class for base attributes.

    Attributes:
    scat        --  The adapted pyobjcryst.molecule.MolAtom.
    parent      --  The ObjCrystParSet this belongs to
    element     --  Non-refinable name of the element (property).

    Managed Parameters:
    x, y, z     --  Atom position in crystal coordinates (ParameterAdapter)
    occ         --  Occupancy of the atom on its crystal location
                    (ParameterAdapter)
    Biso        --  Isotropic scattering factor (ParameterAdapter). This does
                    not exist for dummy atoms. See the 'isDummy' method.
    B11, B22, B33, B12, B21, B23, B32, B13, B31
                --  Anisotropic displacement factor for scatterer
                (ParameterAdapter or ParameterProxy). Note that the Bij and Bji
                parameters are the same.
    
    """

    def __init__(self, name, scat, parent):
        """Initialize

        name    --  The name of the scatterer
        scat    --  The Scatterer instance
        parent  --  The ObjCrystParSet this belongs to

        """
        ScattererParSet.__init__(self, name, scat, parent)
        sp = scat.GetScatteringPower()

        # Only wrap this if there is a scattering power
        if sp is not None:
            self.addParameter(ParameterAdapter("Biso", sp, attr = "Biso"))
            self.addParameter(ParameterAdapter("B11", sp, attr = "B11"))
            self.addParameter(ParameterAdapter("B22", sp, attr = "B22"))
            self.addParameter(ParameterAdapter("B33", sp, attr = "B33"))
            B12 = ParameterAdapter("B12", sp, attr = "B12")
            B21 = ParameterProxy("B21", B12)
            B13 = ParameterAdapter("B13", sp, attr = "B13")
            B31 = ParameterProxy("B31", B13)
            B23 = ParameterAdapter("B23", sp, attr = "B23")
            B32 = ParameterProxy("B32", B23)
            self.addParameter(B12)
            self.addParameter(B21)
            self.addParameter(B13)
            self.addParameter(B31)
            self.addParameter(B23)
            self.addParameter(B32)

        return

    def _getElem(self):
        """Getter for the element type."""
        sp = self.scat.GetScatteringPower()
        if sp:
            return sp.GetSymbol()
        else:
            return "dummy"

    element = property(_getElem)

    def isDummy(self):
        """Indicate whether this atom is a dummy atom."""
        return self.scat.IsDummy()

# End class MolAtomParSet

class MoleculeRestraint(object):
    """Base class for adapting pyobjcryst Molecule restraints to srfit.

    The 'penalty' method calls 'GetLogLikelihood' of the pyobjcryst restraint.
    This implements the 'penalty' method from
    diffpy.srfit.fitbase.restraint.Restraint.  The 'restrain' method is not
    needed or implemented.

    Attributes:
    res     --  The pyobjcryst Molecule restraint.
    scaled  --  A flag indicating if the restraint is scaled (multiplied) by
                the unrestrained point-average chi^2 (chi^2/numpoints) (default
                False).

    """

    def __init__(self, res, scaled = False):
        """Create a Restraint-like from a pyobjcryst Molecule restraint.

        res     --  The pyobjcryst Molecule restraint.
        scaled  --  A flag indicating if the restraint is scaled (multiplied)
                    by the unrestrained point-average chi^2 (chi^2/numpoints)
                    (default False).

        """
        self.res = res
        self.scaled = scaled
        return

    def penalty(self, w = 1.0):
        """Calculate the penalty of the restraint.

        w   --  The point-average chi^2 which is optionally used to scale the
                penalty (default 1.0).
        
        """
        penalty = self.res.GetLogLikelihood()
        if self.scaled:
            penalty *= w
        return penalty

# End class MoleculeRestraint

class BondLengthRestraint(MoleculeRestraint):
    """Restrain the distance between two atoms.

    Attributes:
    atom1   --  The first atom in the bond (MolAtomParSet)
    atom2   --  The second atom in the bond (MolAtomParSet)
    length  --  The length of the bond (Angstroms)
    sigma   --  The uncertainty of the bond length (Angstroms)
    delta   --  The width of the bond (Angstroms)
    res     --  The pyobjcryst BondLength restraint
    scaled  --  A flag indicating if the restraint is scaled (multiplied) by
                the unrestrained point-average chi^2 (chi^2/numpoints) (default
                False)

    """

    def __init__(self, atom1, atom2, length, sigma, delta, scaled = False):
        """Create a bond length restraint.

        atom1   --  First atom (MolAtomParSet) in the bond
        atom2   --  Second atom (MolAtomParSet) in the bond
        length  --  The length of the bond (Angstroms)
        sigma   --  The uncertainty of the bond length (Angstroms)
        delta   --  The width of the bond (Angstroms)
        scaled  --  A flag indicating if the restraint is scaled (multiplied)
                    by the unrestrained point-average chi^2 (chi^2/numpoints)
                    (default False)

        """
        self.atom1 = atom1
        self.atom2 = atom2

        m = self.atom1.scat.GetMolecule()
        res = m.AddBond(atom1.scat, atom2.scat, length, sigma, delta)

        MoleculeRestraint.__init__(self, res, scaled)
        return

    # Give access to the parameters of the restraint
    length = property( lambda self: self.res.GetLength0(),
                       lambda self, val: self.res.SetLength0(val))
    sigma = property( lambda self: self.res.GetLengthSigma(),
                       lambda self, val: self.res.SetLengthSigma(val))
    delta = property( lambda self: self.res.GetLengthDelta(),
                       lambda self, val: self.res.SetLengthDelta(val))

# End class BondLengthRestraint

class BondAngleRestraint(MoleculeRestraint):
    """Restrain the angle defined by three atoms.

    Attributes:
    atom1   --  The first atom in the angle (MolAtomParSet)
    atom2   --  The second atom in the angle (MolAtomParSet)
    atom3   --  The third atom in the angle (MolAtomParSet)
    angle   --  The bond angle (radians)
    sigma   --  The uncertainty of the bond angle (radians)
    delta   --  The width of the bond angle (radians)
    res     --  The pyobjcryst BondAngle restraint
    scaled  --  A flag indicating if the restraint is scaled (multiplied) by
                the unrestrained point-average chi^2 (chi^2/numpoints) (default
                False)

    """

    def __init__(self, atom1, atom2, atom3, angle, sigma, delta, scaled =
            False):
        """Create a bond angle restraint.

        atom1   --  First atom (MolAtomParSet) in the bond angle
        atom2   --  Second (central) atom (MolAtomParSet) in the bond angle
        atom3   --  Third atom (MolAtomParSet) in the bond angle
        angle   --  The bond angle (radians)
        sigma   --  The uncertainty of the bond angle (radians)
        delta   --  The width of the bond angle (radians)
        scaled  --  A flag indicating if the restraint is scaled (multiplied)
                    by the unrestrained point-average chi^2 (chi^2/numpoints)
                    (default False).

        """
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3

        m = self.atom1.scat.GetMolecule()
        res = m.AddBondAngle(atom1.scat, atom2.scat, atom3.scat, angle,
                sigma, delta)

        MoleculeRestraint.__init__(self, res, scaled)
        return

    # Give access to the parameters of the restraint
    angle = property( lambda self: self.res.GetAngle0(),
                       lambda self, val: self.res.SetAngle0(val))
    sigma = property( lambda self: self.res.GetAngleSigma(),
                       lambda self, val: self.res.SetAngleSigma(val))
    delta = property( lambda self: self.res.GetAngleDelta(),
                       lambda self, val: self.res.SetAngleDelta(val))

# End class BondAngleRestraint

class DihedralAngleRestraint(MoleculeRestraint):
    """Restrain the dihedral (torsion) angle defined by four atoms.

    Attributes:
    atom1   --  The first atom in the angle (MolAtomParSet)
    atom2   --  The second (central) atom in the angle (MolAtomParSet)
    atom3   --  The third (central) atom in the angle (MolAtomParSet)
    atom4   --  The fourth atom in the angle (MolAtomParSet)
    angle   --  The dihedral angle (radians)
    sigma   --  The uncertainty of the dihedral angle (radians)
    delta   --  The width of the dihedral angle (radians)
    res     --  The pyobjcryst DihedralAngle restraint
    scaled  --  A flag indicating if the restraint is scaled (multiplied) by
                the unrestrained point-average chi^2 (chi^2/numpoints) (default
                False)

    """

    def __init__(self, atom1, atom2, atom3, atom4, angle, sigma, delta, scaled
            = False):
        """Create a dihedral angle restraint.

        atom1   --  First atom (MolAtomParSet) in the angle
        atom2   --  Second (central) atom (MolAtomParSet) in the angle
        atom3   --  Third (central) atom (MolAtomParSet) in the angle
        atom4   --  Fourth atom in the angle (MolAtomParSet)
        angle   --  The dihedral angle (radians)
        sigma   --  The uncertainty of the dihedral angle (radians)
        delta   --  The width of the dihedral angle (radians)
        scaled  --  A flag indicating if the restraint is scaled (multiplied)
                    by the unrestrained point-average chi^2 (chi^2/numpoints)
                    (default False).

        """
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4

        m = self.atom1.scat.GetMolecule()
        res = m.AddDihedralAngle(atom1.scat, atom2.scat, atom3.scat,
                atom4.scat, angle, sigma, delta)

        MoleculeRestraint.__init__(self, res, scaled)
        return

    # Give access to the parameters of the restraint
    angle = property( lambda self: self.res.GetAngle0(),
                       lambda self, val: self.res.SetAngle0(val))
    sigma = property( lambda self: self.res.GetAngleSigma(),
                       lambda self, val: self.res.SetAngleSigma(val))
    delta = property( lambda self: self.res.GetAngleDelta(),
                       lambda self, val: self.res.SetAngleDelta(val))

# End class DihedralAngleRestraint

class StretchModeParameter(Parameter):
    """Partial Parameter class encapsulating pyobjcryst stretch modes.

    This class relies upon attributes that do not belong to it. Do not
    instantiate this class.

    Required attributes:
    matoms      --  The set of all mutated AtomParSets
    molecule    --  The MoleculeParSet the atoms belong to
    mode        --  The pyobjcryst.molecule.StretchMode used to change atomic
                    positions.

    """

    def setValue(self, val):
        """Change the value of the Parameter.

        This changes the position of the mutated atom and the absolute position
        of the molecule in such a way that the center of mass of the molecule
        does not change.

        """
        curval = self.getValue()
        val = float(val)

        if val == curval:
            return

        # The StretchMode expects the change in mutated value.
        delta = val - curval
        self.mode.Stretch(delta)

        # Let Parameter take care of the general details
        Parameter.setValue(self, val)

        return

    def addAtoms(self, atomlist):
        """Associate MolAtomParSets with the Parameter.

        This will associate additional MolAtomParSets with the Parameter. These
        will be mutated in the exact same way as the primary mutated
        MolAtomParSet. This is useful when a group of atoms should move rigidly
        in response to a change in a bond property.

        """
        if not hasattr(atomlist, "__iter__"):
            atomlist = [atomlist]
        # Record the added atoms in the Parameter
        self.matoms.update(atomlist)
        # Make sure we're observing these atoms
        for a in atomlist:
            a.x.addObserver(self._flush)
            a.y.addObserver(self._flush)
            a.z.addObserver(self._flush)

        # Record the added atoms in the StretchMode
        scatlist = [a.scat for a in atomlist]
        self.mode.AddAtoms(scatlist)
        return

    def notify(self):
        """Notify all mutated Parameters and observers.

        Some of the mutated parameters will be observing us. At the same time
        we need to observe them. Observable won't let us do both, so we notify
        the Parameters that we mutate directly.

        """
        # Notify the atoms that have moved
        for a in self.matoms:
            a.x._flush(None)
            a.y._flush(None)
            a.z._flush(None)
        # Notify the molecule position
        self.molecule.x._flush(None)
        self.molecule.y._flush(None)
        self.molecule.z._flush(None)

        # Notify observers
        Parameter.notify(self)
        return

# End class StretchModeParameter

class BondLengthParameter(StretchModeParameter):
    """Class for abstracting a bond length in a Molecule to a Parameter.

    This wraps up a pyobjcryst.molecule.StretchModeBondLength object so that
    the distance between two MolAtoms in a Molecule can be used as an
    adjustable Parameter. When a bond length is adjusted, the second MolAtom is
    moved, and the absolute position of the Molecule is altered to preserve the
    location of the center of mass within the Crystal. Thus, the x, y and z
    Parameters of the MolAtom and its parent Molecule are altered.

    This Parameter makes it possible to mutate a MolAtom multiple times in a
    single refinement step. If these mutations are not orthogonal, then this
    could lead to nonconvergence of a fit, depending on the optimizer. Consider
    mutating atom2 of a bond directly, and via a BondLengthParameter. The two
    mutations of atom2 may be determined independently by the optimizer, in
    which case the composed mutation will have an unexpected effect on the
    residual. It is best practice to either modify MolAtom positions directly,
    or thorough BondLengthParameters, BondAngleParameters and
    DihedralAngleParameters (which are mutually orthogonal).
    
    Note that by making a BondLengthParameter constant it also makes the
    underlying MolAtomParSets constant. When setting it as nonconstant, each
    MolAtomParSet is set nonconstant.  Changing the bond length changes the
    position of the second MolAtom and Molecule, even if either is set as
    constant.

    Attributes:
    atom1   --  The first MolAtomParSet in the bond
    atom2   --  The second (mutated) MolAtomParSet in the bond
    matoms  --  The set of all mutated MolAtomParSets
    molecule    --  The MoleculeParSet the MolAtomParSets belong to
    mode    --  The pyobjcryst.molecule.StretchModeBondLength for the bond

    Inherited Attributes
    name    --  A name for this Parameter.
    const   --  A flag indicating whether this is considered a constant.
    _value  --  The value of the Parameter. Modified with 'setValue'.
    value   --  Property for 'getValue' and 'setValue'.
    constraint  --  A callable that calculates the value of this Parameter. If
                this is None (None), the the Parameter is responsible for its
                own value. The callable takes no arguments.
    bounds  --  A 2-list defining the bounds on the Parameter. This can be
                used by some optimizers when the Parameter is varied.

    """

    def __init__(self, name, atom1, atom2, value = None, const = False, mode =
            None):
        """Create a BondLengthParameter.

        name    --  The name of the BondLengthParameter
        atom1   --  The first atom (MolAtomParSet) in the bond
        atom2   --  The second (mutated) atom (MolAtomParSet) in the bond
        value   --  An initial value for the bond length. If this is None
                    (default), then the current distance between the atoms will
                    be used.
        const   --  A flag indicating whether the Parameter is constant
                    (default False)
        mode    --  An extant pyobjcryst.molecule.StretchModeBondLength to use.
                    If this is None (default), then a new StretchModeBondLength
                    will be built.

        """

        # Create the mode
        self.mode = mode
        if mode is None:
            self.mode = StretchModeBondLength(atom1.scat, atom2.scat, None)
        # We only add the last atom. This is the one that will move
        self.mode.AddAtom(atom2.scat)
        self.matoms = set([atom2])

        # Observe the atom positions
        for a in [atom1, atom2]:
            a.x.addObserver(self._flush)
            a.y.addObserver(self._flush)
            a.z.addObserver(self._flush)

        self.atom1 = atom1
        self.atom2 = atom2
        self.molecule = atom1.parent

        # We do this last so the atoms are defined before we set any values.
        if value is None:
            value = GetBondLength(atom1.scat, atom2.scat)
        StretchModeParameter.__init__(self, name, value, const)
        self.setConst(const)

        return

    def setConst(self, const = True, value = None):
        """Toggle the Parameter as constant.

        This sets the underlying MolAtomParSet positions const as well.

        const   --  Flag indicating if the Parameter is constant (default
                    True).
        value   --  An optional value for the Parameter (default None). If this
                    is not None, then the Parameter will get a new value,
                    constant or otherwise.

        """
        StretchModeParameter.setConst(self, const, value)

        for a in [self.atom1, self.atom2]:
            a.x.setConst(const)
            a.y.setConst(const)
            a.z.setConst(const)
        return

    def getValue(self):
        """This calculates the value if it might have been changed.

        There is no guarantee that the MolAtomParSets underlying the bond won't
        change, so the bond length is calculated if necessary each time this is
        called.

        """
        if self._value is None:
            val = GetBondLength(self.atom1.scat, self.atom2.scat)
            Parameter.setValue(self, val)

        return self._value


# End class BondLengthParameter
                     
class BondAngleParameter(StretchModeParameter):
    """Class for abstracting a bond angle in a Molecule to a Parameter.

    This wraps up a pyobjcryst.molecule.StretchModeBondAngle object so that the
    angle defined by three MolAtoms in a Molecule can be used as an adjustable
    Parameter. When a bond angle is adjusted, the third MolAtom is moved, and
    the absolute position of the Molecule is altered to preserve the location
    of the center of mass within the crystal.

    See precautions in the BondLengthParameter class.

    Attributes
    atom1   --  The first AtomParSet in the bond angle
    atom2   --  The second (central) MolAtomParSet in the bond angle
    atom3   --  The third (mutated) MolAtomParSet in the bond angle
    matoms  --  The set of all mutated MolAtomParSets
    molecule    --  The MoleculeParSet the MolAtomParSets belong to
    mode    --  The pyobjcryst.molecule.StretchModeBondAngle for the bond angle

    Inherited Attributes
    name    --  A name for this Parameter.
    const   --  A flag indicating whether this is considered a constant.
    _value  --  The value of the Parameter. Modified with 'setValue'.
    value   --  Property for 'getValue' and 'setValue'.
    constraint  --  A callable that calculates the value of this Parameter. If
                this is None (None), the the Parameter is responsible for its
                own value. The callable takes no arguments.
    bounds  --  A 2-list defining the bounds on the Parameter. This can be
                used by some optimizers when the Parameter is varied.

    """

    def __init__(self, name, atom1, atom2, atom3, value = None, const = False,
            mode = None):
        """Create a BondAngleParameter.

        name    --  The name of the BondAngleParameter.
        atom1   --  The first atom (MolAtomParSet) in the bond angle
        atom2   --  The second (central) atom (MolAtomParSet) in the bond angle
        atom3   --  The third (mutated) atom (MolAtomParSet) in the bond angle
        value   --  An initial value for the bond length. If this is None
                    (default), then the current bond angle between the atoms
                    will be used.
        const   --  A flag indicating whether the Parameter is constant
                    (default False).
        mode    --  A pre-built mode to place in this Parameter. If this is
                    None (default), then a StretchMode will be built.

        """

        # Create the stretch mode
        self.mode = mode
        if mode is None:
            self.mode = StretchModeBondAngle(atom1.scat, atom2.scat,
                    atom3.scat, None)
        # We only add the last atom. This is the one that will move
        self.mode.AddAtom(atom3.scat)
        self.matoms = set([atom3])

        # Observe the atom positions
        for a in [atom1, atom2, atom3]:
            a.x.addObserver(self._flush)
            a.y.addObserver(self._flush)
            a.z.addObserver(self._flush)

        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.molecule = atom1.parent

        # We do this last so the atoms are defined before we set any values.
        if value is None:
            value = GetBondAngle(atom1.scat, atom2.scat, atom3.scat)
        StretchModeParameter.__init__(self, name, value, const)
        self.setConst(const)

        return

    def setConst(self, const = True, value = None):
        """Toggle the Parameter as constant.

        This sets the underlying MolAtomParSet positions const as well.

        const   --  Flag indicating if the Parameter is constant (default
                    True).
        value   --  An optional value for the Parameter (default None). If this
                    is not None, then the Parameter will get a new value,
                    constant or otherwise.

        """
        StretchModeParameter.setConst(self, const, value)
        for a in [self.atom1, self.atom2, self.atom3]:
            a.x.setConst(const)
            a.y.setConst(const)
            a.z.setConst(const)
        return

    def getValue(self):
        """This calculates the value if it might have been changed.

        There is no guarantee that the MolAtoms underlying the bond angle won't
        change, so the bond angle is calculated if necessary each time this is
        called.

        """
        if self._value is None:
            val = GetBondAngle(self.atom1.scat, self.atom2.scat,
                    self.atom3.scat)
            Parameter.setValue(self, val)

        return self._value

# End class BondAngleParameter
                     
class DihedralAngleParameter(StretchModeParameter):
    """Class for abstracting a dihedral angle in a Molecule to a Parameter.

    This wraps up a pyobjcryst.molecule.StretchModeTorsion object so that the
    angle defined by four MolAtoms ([a1-a2].[a3-a4]) in a Molecule can be used
    as an adjustable parameter. When a dihedral angle is adjusted, the fourth
    MolAtom is moved, and the absolute position of the Molecule is altered to
    preserve the location of the center of mass within the crystal.

    See precautions in the BondLengthParameter class.

    Attributes
    atom1   --  The first MolAtomParSet in the dihedral angle
    atom2   --  The second (central) MolAtomParSet in the dihedral angle
    atom3   --  The third (central) MolAtomParSet in the dihedral angle
    atom4   --  The fourth (mutated) MolAtomParSet in the dihedral angle
    matoms  --  The set of all mutated MolAtomParSets
    molecule    --  The MoleculeParSet the atoms belong to
    mode    --  The pyobjcryst.molecule.StretchModeTorsion for the dihedral
                angle

    Inherited Attributes
    name    --  A name for this Parameter.
    const   --  A flag indicating whether this is considered a constant.
    _value  --  The value of the Parameter. Modified with 'setValue'.
    value   --  Property for 'getValue' and 'setValue'.
    constraint  --  A callable that calculates the value of this Parameter. If
                this is None (None), the the Parameter is responsible for its
                own value. The callable takes no arguments.
    bounds  --  A 2-list defining the bounds on the Parameter. This can be
                used by some optimizers when the Parameter is varied.

    """

    def __init__(self, name, atom1, atom2, atom3, atom4, value = None, const =
            False, mode = None):
        """Create a DihedralAngleParameter.

        name    --  The name of the DihedralAngleParameter
        atom1   --  The first atom (MolAtomParSet) in the dihderal angle
        atom2   --  The second (central) atom (MolAtomParSet) in the dihderal
                    angle
        atom3   --  The third (central) atom (MolAtomParSet) in the dihderal
                    angle
        atom4   --  The fourth (mutated) atom (MolAtomParSet) in the dihderal
                    angle
        value   --  An initial value for the bond length. If this is None
                    (default), then the current dihedral angle between atoms
                    will be used.
        const   --  A flag indicating whether the Parameter is constant
                    (default False).
        mode    --  A pre-built mode to place in this Parameter. If this is
                    None (default), then a StretchMode will be built.

        """

        # Create the stretch mode
        self.mode = mode
        if mode is None:
            self.mode = StretchModeTorsion(atom2.scat, atom3.scat, None)
        # We only add the last atom. This is the one that will move
        self.mode.AddAtom(atom4.scat)
        self.matoms = set([atom4])

        # Observe the atom positions
        for a in [atom1, atom2, atom3, atom4]:
            a.x.addObserver(self._flush)
            a.y.addObserver(self._flush)
            a.z.addObserver(self._flush)

        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        self.molecule = atom1.parent

        # We do this last so the atoms are defined before we set any values.
        if value is None:
            value = GetDihedralAngle(atom1.scat, atom2.scat, atom3.scat,
                    atom4.scat)
        StretchModeParameter.__init__(self, name, value, const)
        self.setConst(const)

        return

    def setConst(self, const = True, value = None):
        """Toggle the Parameter as constant.

        This sets the underlying MolAtomParSet positions const as well.

        const   --  Flag indicating if the Parameter is constant (default
                    True).
        value   --  An optional value for the Parameter (default None). If this
                    is not None, then the Parameter will get a new value,
                    constant or otherwise.

        """
        StretchModeParameter.setConst(self, const, value)
        for a in [self.atom1, self.atom2, self.atom3, self.atom4]:
            a.x.setConst(const)
            a.y.setConst(const)
            a.z.setConst(const)
        return

    def getValue(self):
        """This calculates the value if it might have been changed.

        There is no guarantee that the MolAtomParSets underlying the dihedral
        angle won't change from some other Parameter, so the value is
        recalculated each time.

        """
        if self._value is None:
            val = GetDihedralAngle(self.atom1.scat, self.atom2.scat,
                    self.atom3.scat, self.atom4.scat)
            Parameter.setValue(self, val)

        return self._value

# End class DihedralAngleParameter

class ObjCrystParSet(BaseStructure):
    """A adaptor for pyobjcryst.crystal.Crystal instance.

    This class derives from diffpy.srfit.fitbase.parameterset.ParameterSet. See
    this class for base attributes.

    Attributes:
    stru        --  The adapted pyobjcryst.Crystal.
    scatterers  --  The list of aggregated ScattererParSets (either AtomParSet
                    or MoleculeParSet), provided for convenience.

    Managed Parameters:
    a, b, c, alpha, beta, gamma --  Lattice parameters (ParameterAdapter)

    Managed ParameterSets:
    <sname>     --  A ScattererParSet (either AtomParSet or MoleculeParSet),
                    where <sname> is the name of the adapted
                    pyobjcryst.atom.Atom or pyobjcryst.molecule.Molecule.

    """

    def __init__(self, cryst, name):
        """Initialize

        cryst   --  An pyobjcryst.Crystal instance.
        name    --  A name for this ParameterSet

        """
        ParameterSet.__init__(self, name)
        self.stru = cryst

        self.addParameter(ParameterAdapter("a", self.stru, attr = "a"))
        self.addParameter(ParameterAdapter("b", self.stru, attr = "b"))
        self.addParameter(ParameterAdapter("c", self.stru, attr = "c"))
        self.addParameter(ParameterAdapter("alpha", self.stru, attr =
            "alpha"))
        self.addParameter(ParameterAdapter("beta", self.stru, attr = "beta"))
        self.addParameter(ParameterAdapter("gamma", self.stru, attr =
            "gamma"))

        # Now we must loop over the scatterers and create parameter sets from
        # them.
        self.scatterers = []
        snames = []

        for j in range(self.stru.GetNbScatterer()):
            s = self.stru.GetScatt(j)
            name = s.GetName()
            if not name:
                raise AttributeError("Each Scatterer must have a name")
            if name in snames:
                raise AttributeError("MolAtom name '%s' is duplicated"%name)

            # Now create the proper object
            cname = s.GetClassName()
            if cname == "Atom":
                parset = AtomParSet(name, s, self)
            elif cname == "Molecule":
                parset = MoleculeParSet(name, s, self)
            else:
                raise AttributeError("Unrecognized scatterer '%s'"%cname)

            self.addParameterSet(parset)
            self.scatterers.append(parset)
            snames.append(name)

        # Constrain parameters to the space group
        sgname = self.getSpaceGroup()
        from diffpy.srfit.structure.sgconstraints import _constrainSpaceGroup
        _constrainSpaceGroup(self, sgname)

        return

    @classmethod
    def canAdapt(self, stru):
        """Return whether the structure can be adapted by this class."""
        from pyobjcryst.crystal import Crystal
        return isinstance(stru, Crystal)

    def getLattice(self):
        """Get the ParameterSet containing the lattice Parameters."""
        return self

    def getScatterers(self):
        """Get a list of ParameterSets that represents the scatterers.

        The site positions must be accessible from the list entries via the
        names "x", "y", and "z". The ADPs must be accessible as well, but the
        name and nature of the ADPs (U-factors, B-factors, isotropic,
        anisotropic) depends on the adapted structure.

        """
        return self.scatterers

    def getSpaceGroup(self):
        """Get the HM space group symbol for the structure."""
        sg = self.stru.GetSpaceGroup().GetName()
        extn = self.stru.GetSpaceGroup().GetExtension()
        extnstr = ":%s"%extn
        if sg.endswith(extnstr):
            sg = sg[:-len(extnstr)]
        return sg



# End class ObjCrystParSet


__id__ = "$Id$"
