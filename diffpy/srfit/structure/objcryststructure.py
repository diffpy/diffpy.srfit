#!/usr/bin/env python
"""Wrappers for adapting pyobjcryst.Crystal to a srfit ParameterSet.

This will adapt a pyobjcryst.Crystal into the ParameterSet interface. The
following classes are adapted.

pyobjcryst.Crystal  ->  ObjCrystParSet
pyobjcryst.Atom     ->  AtomParSet
pyobjcryst.Molecule ->  MoleculeParSet
pyobjcryst.MolAtom  ->  MolAtomParSet

Related to the adaptation of Molecule and MolAtom, there are adaptors for
specifying molecule restraints.
BondLengthRestraint
BondAngleRestraint
DihedralAngleRestraint

There are also Parameters for encapsulating and modifying atoms via their
relative positions. These Parameters can also act like restraints, and modify
the positions of multiple MolAtoms.
BondLengthParameter
BondAngleParameter
DihedralAngleParameter


"""


from diffpy.srfit.equation import Clicker
from diffpy.srfit.fitbase.parameter import Parameter, ParameterWrapper
from diffpy.srfit.fitbase.parameterset import ParameterSet
from diffpy.srfit.fitbase.restraint import Restraint
from diffpy.srfit.fitbase.constraint import Constraint

from pyobjcryst import GetBondLength, GetBondAngle, GetDihedralAngle
from pyobjcryst import StretchModeBondLength, StretchModeBondAngle
from pyobjcryst import StretchModeTorsion

class ScattererParSet(ParameterSet):
    """A base wrapper for an Objcryst Scatterer.

    This class derives from ParameterSet and adapts pyobjcryst.Scatterer
    derivatives (Molecule, Atom) and objects with a similar interface
    (MolAtom).

    Attributes:
    x (y, z)    --  Atom position in crystal coordinates (Parameter)
    occupancy   --  Occupancy of the atom on its crystal location (Parameter)
    parent      --  The ParameterSet this belongs to
    
    Other attributes are inherited from
    diffpy.srfit.fitbase.parameterset.ParameterSet

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

        # x, y, z, occupancy
        self.addParameter(ParameterWrapper(self.scat, "x", attr = "X"))
        self.addParameter(ParameterWrapper(self.scat, "y", attr = "Y"))
        self.addParameter(ParameterWrapper(self.scat, "z", attr = "Z"))
        self.addParameter(ParameterWrapper(self.scat, "occ", attr =
            "Occupancy"))
        return

# End class ScattererParSet

class AtomParSet(ScattererParSet):
    """A wrapper for an Objcryst Atom.

    This class derives from ParameterSet.

    Attributes:
    x (y, z)    --  Atom position in crystal coordinates (Parameter)
    occupancy   --  Occupancy of the atom on its crystal location (Parameter)
    parent      --  The ObjCrystParSet this belongs to
    biso        --  Isotropic scattering factor (Parameter).
    element     --  Non-refinable name of the element.
    
    Other attributes are inherited from diffpy.srfit.fitbase.parameterset.ParameterSet

    """

    def __init__(self, name, atom, parent):
        """Initialize

        name    --  The name of the scatterer
        scat    --  The Scatterer instance
        parent  --  The ObjCrystParSet this belongs to

        """
        ScattererParSet.__init__(self, name, atom, parent)
        sp = atom.GetScatteringPower()

        # The Biso parameter
        self.addParameter(ParameterWrapper(sp, "biso", attr = "Biso"))
        return

    def _getElem(self):
        """Getter for the element type."""
        return self.scat.GetScatteringPower().GetSymbol()

    element = property(_getElem)

# End class AtomParSet

class MoleculeParSet(ScattererParSet):
    """A wrapper for an Objcryst Molecule.

    This class derives from ParameterSet.

    Attributes:
    x (y, z)    --  Molecule position in crystal coordinates (Parameter)
    occupancy   --  Occupancy of the molecule on its crystal location
                    (Parameter)
    q0, q1, q2, q3  --  Orientational quaternion parameters (Parameter)
    parent      --  The ObjCrystParSet this belongs to
    atoms       --  A list of contained MolAtomParSets.
    
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
        self.addParameter(ParameterWrapper(self.scat, "q0", attr = "Q0"))
        self.addParameter(ParameterWrapper(self.scat, "q1", attr = "Q1"))
        self.addParameter(ParameterWrapper(self.scat, "q2", attr = "Q2"))
        self.addParameter(ParameterWrapper(self.scat, "q3", attr = "Q3"))

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
        """Wrap the stretch modes implicit to the molecule as parameters.

        This will wrap StretchModeBondLengths and StretchModeBondAngles of the
        Molecule as Parameters. Note that this requires that the MolBondAtoms
        in the Molecule came in with unique names.  Torsion angles are not
        wrapped, as there is not enough information to determine each atom in
        the angle.

        The parameters will be given the concatenated name of its constituents.

        bond lengths "bl_aname1_aname2"
        bond angles "ba_aname1_aname2_aname3"

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

        This creates a BondLengthParameter to the molecule that can be adjusted
        during the fit.

        name    --  The name of the bond length parameter
        atom1   --  The first atom (AtomParSet) in the bond
        atom2   --  The second (mutated) atom (AtomParSet) in the bond
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

        This creates a BondAngleParameter to the molecule that can be adjusted
        during the fit.

        name    --  The name of the bond length parameter
        atom1   --  The first atom (AtomParSet) in the bond angle
        atom2   --  The second (central) atom (AtomParSet) in the bond angle
        atom3   --  The third (mutated) atom (AtomParSet) in the bond angle
        value   --  An initial value for the bond length. If this is None
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

        This creates a DihedralAngleParameter to the molecule that can be
        adjusted during the fit.

        name    --  The name of the bond length parameter
        atom1   --  The first atom (AtomParSet) in the dihderal angle
        atom2   --  The second (central) atom (AtomParSet) in the dihderal
                    angle
        atom3   --  The third (central) atom (AtomParSet) in the dihderal angle
        atom4   --  The fourth (mutated) atom (AtomParSet) in the dihderal
                    angle
        value   --  An initial value for the bond length. If this is None
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
    """A wrapper for an Objcryst MolAtom.

    This class derives from ParameterSet. Note that MolAtom does not derive
    from Scatterer, but the relevant interface is the same within pyobjcryst.

    Attributes:
    x (y, z)    --  Atom position in crystal coordinates (Parameter)
    occupancy   --  Occupancy of the atom on its crystal location (Parameter)
    parent      --  The ObjCrystParSet this belongs to
    biso        --  Isotropic scattering factor (Parameter). This does not
                    exist for dummy atoms.
    element     --  Non-refinable name of the element.
    
    Other attributes are inherited from diffpy.srfit.fitbase.parameterset.ParameterSet

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
            self.addParameter(ParameterWrapper(sp, "biso", attr = "Biso"))

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
        return self.element == "dummy"

# End class MolAtomParSet

class ObjCrystParSet(ParameterSet):
    """A wrapper for ObjCryst Crystal instance.

    Attributes:
    scatterers   --  The list of aggregated ScattererParSets (either AtomParSet
                    or MoleculeParSet).
    a, b, c, alpha, beta, gamma --  Lattice parameters (Parameter)

    Other attributes are inherited from
    diffpy.srfit.fitbase.parameterset.ParameterSet
    
    """

    def __init__(self, cryst, name):
        """Initialize

        cryst   --  An pyobjcryst.Crystal instance.
        name    --  A name for this ParameterSet

        """
        ParameterSet.__init__(self, name)
        self.cryst = cryst

        self.addParameter(ParameterWrapper(self.cryst, "a", attr = "a"))
        self.addParameter(ParameterWrapper(self.cryst, "b", attr = "b"))
        self.addParameter(ParameterWrapper(self.cryst, "c", attr = "c"))
        self.addParameter(ParameterWrapper(self.cryst, "alpha", 
            attr = "alpha"))
        self.addParameter(ParameterWrapper(self.cryst, "beta", attr = "beta"))
        self.addParameter(ParameterWrapper(self.cryst, "gamma",
            attr = "gamma"))

        # Constrain the lattice before we go any further.
        sgmap = {}
        sgnum = self.cryst.GetSpaceGroup().GetSpaceGroupNumber()
        from diffpy.Structure import SpaceGroups
        sg = SpaceGroups.GetSpaceGroup(sgnum)
        system = sg.crystal_system
        if not system:
            system = "Triclinic"
        system = system.title()
        from .sgconstraints import constrainSpaceGroup
        constrainSpaceGroup(self, system)

        # Now we must loop over the scatterers and create parameter sets from
        # them.
        self.scatterers = []
        snames = []

        for j in range(self.cryst.GetNbScatterer()):
            s = self.cryst.GetScatt(j)
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


        return

# End class ObjCrystParSet

class MoleculeRestraint(object):
    """Base class for adapting pyobjcryst Molecule restraints to srfit.

    The GetLogLikelihood method of the pyobjcryst restraint is mapped to the
    penalty method. The restrain method is not present.

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
    res     --  The pyobjcryst Molecule restraint (MolAtomParSet)
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
    res     --  The pyobjcryst Molecule restraint
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
    res     --  The pyobjcryst Molecule restraint
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
    """Partial class for Parameters that encapsulate pyobjcryst.StretchModes.

    This class relies upon attributes that do not belong to it. Do not
    instantiate this class.

    Required attributes:
    matoms  --  The set of all mutated AtomParSets
    molecule    --  The molecule the atoms belong to
    mode    --  The pyobjcryst.StretchModeBondLength for the bond

    """

    def setValue(self, val):
        """Change the value of the parameter.

        This changes the position of the mutated atom and the absolute position
        of the molecule in such a way that the center of mass of the molecule
        does not change.

        """
        curval = self.getValue()
        if val == curval:
            return

        # The StretchMode expects the change in mutated value.
        delta = val - curval
        self.mode.Stretch(delta)

        self.value = val

        # Click everything that has changed
        self.click()
        # We just set the value of the parameter, so we don't have to
        # recalculate until something underneath changes.
        self.calclicker.click()
        return

    def addAtoms(self, atomlist):
        """Associate MolAtomParSets with the parameter.

        This will associate additional MolAtomParSets with the Parameter. These
        will be mutated in the exact same way as the primary mutated
        MolAtomParSet. This is useful when a group of atoms should move rigidly
        in response to a change in a bond property.

        """
        # Record the added atoms in the Parameter
        self.matoms.update(atomlist)
        for a in atomlist:
            self.clicker.addSubject(a.clicker)

        # Record the added atoms in the StretchMode
        scatlist = [a.scat for a in atomlist]
        self.mode.AddAtoms(scatlist)
        return

    def click(self):
        """Click the clickers of all the mutated parameters."""
        # Click the atoms that have moved
        for a in self.matoms:
            a.x.clicker.click()
            a.y.clicker.click()
            a.z.clicker.click()
        # Click the molecule position
        self.molecule.x.clicker.click()
        self.molecule.y.clicker.click()
        self.molecule.z.clicker.click()
        return

# End class StretchModeParameter

class BondLengthParameter(StretchModeParameter):
    """Class for abstracting a bond length in a Molecule to a Parameter.

    This wraps up a pyobjcryst.StretchModeBondLength object so that the
    distance between two atoms in a Molecule can be used as an adjustable
    parameter. When a bond length is adjusted, the second atom is moved, and
    the absolute position of the molecule is altered to preserve the location
    of the center of mass within the crystal.

    This Parameter makes it possible to mutate an atom multiple times in a
    single refinement step. If these mutations are not orthogonal, then this
    could lead to nonconvergence of a fit. Consider mutating atom2 of a bond
    directly, and via a BondLengthParameter. The two mutations of atom2 may be
    determined independently by the optimizer, in which case the composed
    mutation will have an unexpected effect on the residual. It is best
    practice to either modify atom positions directly, or thorough bond
    lengths, bond angles and dihedral angles (which are mutually orthogonal).
    
    Note that by making a bond constant it also makes the underlying atoms
    constant. When setting it as nonconstant, each atom is set nonconstant.
    Changing the bond length changes the position of the second atom and
    molecule, even if either is set as constant. 

    Attributes
    atom1   --  The first AtomParSet in the bond
    atom2   --  The second (mutated) AtomParSet in the bond
    matoms  --  The set of all mutated AtomParSets
    molecule    --  The molecule the atoms belong to
    mode    --  The pyobjcryst.StretchModeBondLength for the bond
    calclicker  --  A Clicker to record when we have calculated the latest
                    value of the parameter.

    Inherited Attributes
    name    --  A name for this Parameter. Names should be unique within a
                ParameterSet.
    clicker --  A Clicker instance for recording change in the value.
    value   --  The value of the Parameter. Modified with setValue.
    const   --  A flag indicating whether the Parameter is constant.

    """

    def __init__(self, name, atom1, atom2, value = None, const = False, mode =
            None):
        """Create a BondLengthParameter.

        name    --  The name of the bond length parameter
        atom1   --  The first atom (AtomParSet) in the bond
        atom2   --  The second (mutated) atom (AtomParSet) in the bond
        value   --  An initial value for the bond length. If this is None
                    (default), then the current distance between the atoms will
                    be used.
        const   --  A flag indicating whether the Parameter is constant
                    (default False)
        mode    --  An extant pyobjcryst.StretchModeBondLength to use. If this
                    is None (default), then a new StretchModeBondLength will be
                    built.

        """
        if value is None:
            value = GetBondLength(atom1.scat, atom2.scat)


        self.calclicker = Clicker()
        self.calclicker.click()

        self.value = value

        StretchModeParameter.__init__(self, name, value, const)

        atom1.setConst(const)
        atom2.setConst(const)

        # Create the mode
        self.mode = mode
        if mode is None:
            self.mode = StretchModeBondLength(atom1.scat, atom2.scat, None)
        # We only add the last atom. This is the one that will move
        self.mode.AddAtom(atom2.scat)

        self.clicker.addSubject(atom1.clicker)
        self.clicker.addSubject(atom2.clicker)

        self.atom1 = atom1
        self.atom2 = atom2
        self.molecule = atom1.parent
        self.matoms = set([atom2])

        return

    def setConst(self, const = True, value = None):
        """Toggle the Parameter as constant.

        This sets the underlying atoms as const as well.

        const   --  Flag indicating if the parameter is constant (default
                    True).
        value   --  An optional value for the parameter (default None). If this
                    is not None, then the parameter will get a new value,
                    constant or otherwise.

        """
        StretchModeParameter.setConst(self, const, value)
        self.atom1.setConst(const)
        self.atom2.setConst(const)
        return

    def getValue(self):
        """This calculates the value if it might have been changed.

        There is no guarantee that the atoms underlying the bond won't change,
        so the bond length is calculated if necessary each time this is called.

        """
        if self.calclicker < self.clicker:
            self.value = GetBondLength(self.atom1.scat, self.atom2.scat)
            self.calclicker.click()

        return self.value

# End class BondLengthParameter
                     
class BondAngleParameter(StretchModeParameter):
    """Class for abstracting a bond angle in a Molecule to a Parameter.

    This wraps up a pyobjcryst.StretchModeBondAngle object so that the
    angle defined by three atoms in a Molecule can be used as an adjustable
    parameter. When a bond angle is adjusted, the third atom is moved, and
    the absolute position of the molecule is altered to preserve the location
    of the center of mass within the crystal.

    This Parameter makes it possible to mutate an atom multiple times in a
    single refinement step. If these mutations are not orthogonal, then this
    could lead to nonconvergence of a fit. Consider mutating atom2 of a bond
    directly, and via a BondLengthParameter. The two mutations of atom2 may be
    determined independently by the optimizer, in which case the composed
    mutation will have an unexpected affect on the residual. It is best
    practice to either modify atom positions directly, or thorough bond
    lengths, bond angles and dihedral angles (which are orthogonal).
    
    Note that by making a bond angle constant it also makes the underlying
    atoms constant. When setting it as nonconstant, each atom is set
    nonconstant.  Changing the bond angle changes the position of the third
    atom and molecule, even if either is set as constant. 

    Attributes
    atom1   --  The first AtomParSet in the bond angle
    atom2   --  The second (central) AtomParSet in the bond angle
    atom3   --  The third (mutated) AtomParSet in the bond angle
    matoms  --  The set of all mutated AtomParSets
    molecule    --  The molecule the atoms belong to
    mode    --  The pyobjcryst.StretchModeBondAngle for the bond angle
    calclicker  --  A Clicker to record when we have calculated the latest
                    value of the parameter.

    Inherited Attributes
    name    --  A name for this Parameter. Names should be unique within a
                ParameterSet.
    clicker --  A Clicker instance for recording change in the value.
    value   --  The value of the Parameter. Modified with setValue.
    const   --  A flag indicating whether the Parameter is constant.

    """

    def __init__(self, name, atom1, atom2, atom3, value = None, const = False,
            mode = None):
        """Create a BondAngleParameter.

        name    --  The name of the bond length parameter
        atom1   --  The first atom (AtomParSet) in the bond angle
        atom2   --  The second (central) atom (AtomParSet) in the bond angle
        atom3   --  The third (mutated) atom (AtomParSet) in the bond angle
        value   --  An initial value for the bond length. If this is None
                    (default), then the current bond angle between the atoms
                    will be used.
        const   --  A flag indicating whether the Parameter is constant
                    (default False).
        mode    --  A pre-built mode to place in this Parameter. If this is
                    None (default), then a StretchMode will be built.

        """
        if value is None:
            value = GetBondAngle(atom1.scat, atom2.scat, atom3.scat)

        self.calclicker = Clicker()
        self.calclicker.click()

        StretchModeParameter.__init__(self, name, value, const)

        atom1.setConst(const)
        atom2.setConst(const)
        atom3.setConst(const)

        # Create the stretch mode
        self.mode = mode
        if mode is None:
            self.mode = StretchModeBondAngle(atom1.scat, atom2.scat,
                    atom3.scat, None)
        # We only add the last atom. This is the one that will move
        self.mode.AddAtom(atom3.scat)

        self.clicker.addSubject(atom1.clicker)
        self.clicker.addSubject(atom2.clicker)
        self.clicker.addSubject(atom3.clicker)

        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.molecule = atom1.parent
        self.matoms = set([atom3])

        return

    def setConst(self, const = True, value = None):
        """Toggle the Parameter as constant.

        This sets the underlying atoms as const as well.

        const   --  Flag indicating if the parameter is constant (default
                    True).
        value   --  An optional value for the parameter (default None). If this
                    is not None, then the parameter will get a new value,
                    constant or otherwise.

        """
        StretchModeParameter.setConst(self, const, value)
        self.atom1.setConst(const)
        self.atom2.setConst(const)
        self.atom3.setConst(const)
        return

    def getValue(self):
        """This calculates the value if it might have been changed.

        There is no guarantee that the atoms underlying the bond angle won't
        change, so the bond angle is calculated if necessary each time this is
        called.

        """
        if self.calclicker < self.clicker:
            self.value = GetBondAngle(self.atom1.scat, self.atom2.scat,
                    self.atom3.scat)
            self.calclicker.click()

        return self.value

# End class BondAngleParameter
                     
class DihedralAngleParameter(StretchModeParameter):
    """Class for abstracting a dihedral angle in a Molecule to a Parameter.

    This wraps up a pyobjcryst.StretchModeTorsion object so that the angle
    defined by four atoms ([a1-a2].[a3-a4]) in a Molecule can be used as an
    adjustable parameter. When a dihedral angle is adjusted, the fourth atom is
    moved, and the absolute position of the molecule is altered to preserve the
    location of the center of mass within the crystal.

    This Parameter makes it possible to mutate an atom multiple times in a
    single refinement step. If these mutations are not orthogonal, then this
    could lead to nonconvergence of a fit. Consider mutating atom2 of a bond
    directly, and via a BondLengthParameter. The two mutations of atom2 may be
    determined independently by the optimizer, in which case the composed
    mutation will have an unexpected affect on the residual. It is best
    practice to either modify atom positions directly, or thorough bond
    lengths, bond angles and dihedral angles (which are orthogonal).
    
    Note that by making a dihedral angle constant it also makes the underlying
    atoms constant. When setting it as nonconstant, each atom is set
    nonconstant.  Changing the bond angle changes the position of the third
    atom and molecule, even if either is set as constant. 

    Attributes
    atom1   --  The first AtomParSet in the dihedral angle
    atom2   --  The second (central) AtomParSet in the dihedral angle
    atom3   --  The third (central) AtomParSet in the dihedral angle
    atom4   --  The fourth (mutated) AtomParSet in the dihedral angle
    matoms  --  The set of all mutated AtomParSets
    molecule    --  The molecule the atoms belong to
    mode    --  The pyobjcryst.StretchModeTorsion for the dihedral angle
    calclicker  --  A Clicker to record when we have calculated the latest
                    value of the parameter.

    Inherited Attributes
    name    --  A name for this Parameter. Names should be unique within a
                ParameterSet.
    clicker --  A Clicker instance for recording change in the value.
    value   --  The value of the Parameter. Modified with setValue.
    const   --  A flag indicating whether the Parameter is constant.

    """

    def __init__(self, name, atom1, atom2, atom3, atom4, value = None, const =
            False, mode = None):
        """Create a DihedralAngleParameter.

        name    --  The name of the bond length parameter
        atom1   --  The first atom (AtomParSet) in the dihderal angle
        atom2   --  The second (central) atom (AtomParSet) in the dihderal
                    angle
        atom3   --  The third (central) atom (AtomParSet) in the dihderal angle
        atom4   --  The fourth (mutated) atom (AtomParSet) in the dihderal
                    angle
        value   --  An initial value for the bond length. If this is None
                    (default), then the current dihedral angle between atoms
                    will be used.
        const   --  A flag indicating whether the Parameter is constant
                    (default False).
        mode    --  A pre-built mode to place in this Parameter. If this is
                    None (default), then a StretchMode will be built.

        """
        if value is None:
            value = GetDihedralAngle(atom1.scat, atom2.scat, atom3.scat,
                    atom4.scat)

        self.calclicker = Clicker()
        self.calclicker.click()

        StretchModeParameter.__init__(self, name, value, const)

        atom1.setConst(const)
        atom2.setConst(const)
        atom3.setConst(const)
        atom4.setConst(const)

        # Create the stretch mode
        self.mode = mode
        if mode is None:
            self.mode = StretchModeTorsion(atom2.scat, atom3.scat, None)
        # We only add the last atom. This is the one that will move
        self.mode.AddAtom(atom4.scat)

        self.clicker.addSubject(atom1.clicker)
        self.clicker.addSubject(atom2.clicker)
        self.clicker.addSubject(atom3.clicker)
        self.clicker.addSubject(atom4.clicker)

        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        self.molecule = atom1.parent
        self.matoms = set([atom4])

        return

    def setConst(self, const = True, value = None):
        """Toggle the Parameter as constant.

        This sets the underlying atoms as const as well.

        const   --  Flag indicating if the parameter is constant (default
                    True).
        value   --  An optional value for the parameter (default None). If this
                    is not None, then the parameter will get a new value,
                    constant or otherwise.

        """
        StretchModeParameter.setConst(self, const, value)
        self.atom1.setConst(const)
        self.atom2.setConst(const)
        self.atom3.setConst(const)
        self.atom4.setConst(const)
        return

    def getValue(self):
        """This calculates the value if it might have been changed.

        There is no guarantee that the atoms underlying the dihedral angle
        won't change from some other parameter, so the value is recalculated
        each time.

        """
        if self.calclicker < self.clicker:
            self.value = GetDihedralAngle(self.atom1.scat, self.atom2.scat,
                    self.atom3.scat, self.atom4.scat)
            self.calclicker.click()

        return self.value

# End class DihedralAngleParameter

__id__ = "$Id$"
