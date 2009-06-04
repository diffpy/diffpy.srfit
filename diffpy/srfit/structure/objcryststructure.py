#!/usr/bin/env python
"""Wrappers for adapting pyobjcryst.Crystal to a srfit ParameterSet.

"""


from diffpy.srfit.equation import Clicker
from diffpy.srfit.fitbase.parameter import Parameter, ParameterProxy
from diffpy.srfit.fitbase.parameter import ParameterWrapper
from diffpy.srfit.fitbase.parameterset import ParameterSet
from diffpy.srfit.fitbase.restraint import Restraint
from diffpy.srfit.fitbase.constraint import Constraint

from pyobjcryst import GetBondLength, GetBondAngle, GetDihedralAngle
from pyobjcryst import StretchModeBondLength, StretchModeBondAngle
from pyobjcryst import StretchModeTorsion

class ScattererParSet(ParameterSet):
    """A base wrapper for an Objcryst Scatterer.

    This class derives from ParameterSet.

    Attributes:
    x (y, z)    --  Atom position in crystal coordinates (Parameter)
    occupancy   --  Occupancy of the atom on its crystal location (Parameter)
    parent      --  The ParameterSet this belongs to
    
    """

    def __init__(self, name, scat, parent):
        """Initialize

        name    --  The name of the scatterer
        scat    --  The Scatterer instance
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
    parent      --  The ParameterSet this belongs to
    biso        --  Isotropic scattering factor (Parameter).
    element     --  Non-refinable name of the element.
    
    """

    def __init__(self, name, atom, parent):
        """Initialize

        name    --  The name of the scatterer
        scat    --  The Scatterer instance
        parent  --  The Crystal this belongs to

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
    occupancy   --  Occupancy of the atom on its crystal location (Parameter)
    parent      --  The ParameterSet this belongs to
    atoms       --  A list of contained MolAtomParSets.
    
    """

    def __init__(self, name, molecule, parent):
        """Initialize

        name    --  The name of the scatterer
        molecule    --  The Molecule instance
        parent  --  The Crystal this belongs to

        """
        ScattererParSet.__init__(self, name, molecule, parent)

        self.atoms = []
        self.updateConfiguration()
        return

    def updateConfiguration(self):
        """Update the ParameterSet in response to a change in the molecule.

        This should be called whenever the number of atoms, restraints or
        constraints are changed using the molecule object itself.
        
        """
        molecule = self.scat

        self.atoms = []

        # Clear the current cache of sub-objects
        for org in self._organizers:
            self.clicker.removeSubject(org.clicker)
        self._orgdict = {}
        self._organizers = []
        self._parameters = []
        self._constraints = []
        self._restraints = set()

        # Wrap the MolAtoms within the molecule
        cdict = {}
        for a in molecule:
            name = a.GetName()
            if not name:
                sp = a.GetScatteringPower()
                name = "dummy"
                if sp:
                    name = sp.GetSymbol()

            i = cdict.get(name, 0)
            sname = name
            # Only add the number if there are atoms with the same name. We'll
            # fix the first one later.
            if i:
                sname += "_%i"%i
            cdict[name] = i+1
            atom = MolAtomParSet(sname, a, self)
            atom.molecule = self
            self.addParameterSet(atom)
            self.atoms.append(atom)

        # Fix the name of the first element in a repeated group
        for name, i in cdict.iteritems():
            if i:
                atom = getattr(self, name)
                atom.name += "_%i"%i

        # Add orientiation quaternion
        self.addParameter(ParameterWrapper(self.scat, "Q0", attr = "Q0"))
        self.addParameter(ParameterWrapper(self.scat, "Q1", attr = "Q1"))
        self.addParameter(ParameterWrapper(self.scat, "Q2", attr = "Q2"))
        self.addParameter(ParameterWrapper(self.scat, "Q3", attr = "Q3"))

        # wrap restraints
        # This turns bonds, bond angles and dihedral angles into restraints of
        # the fit with the help of the MoleculeRestraint class

        for b in self.scat.GetBondList:
            res = MoleculeRestraint(b)
            self._restraints.add(res)

        for ba in self.scat.GetBondAngleList:
            res = MoleculeRestraint(ba)
            self._restraints.add(res)

        for da in self.scat.GetDihedralAngleList:
            res = MoleculeRestraint(da)
            self._restraints.add(res)

        # FIXME - wrap contraints
        for b in self.scat.GetStretchModeBondLengthList():

        return

# End class MoleculeParSet

class MolAtomParSet(ScattererParSet):
    """A wrapper for an Objcryst MolAtom.

    This class derives from ParameterSet. Note that MolAtom does not derive
    from Scatterer, but the relevant interface is the same within pyobjcryst.

    Attributes:
    x (y, z)    --  Atom position in crystal coordinates (Parameter)
    occupancy   --  Occupancy of the atom on its crystal location (Parameter)
    parent      --  The ParameterSet this belongs to
    biso        --  Isotropic scattering factor (Parameter). This does not
                    exist for dummy atoms.
    element     --  Non-refinable name of the element.
    
    """

    def __init__(self, name, scat, parent):
        """Initialize

        name    --  The name of the scatterer
        scat    --  The Scatterer instance
        parent  --  The ParameterSet this belongs to

        """
        ScattererParSet.__init__(self, name, scat, parent)
        sp = scat.GetScatteringPower()

        # Only wrap this if there is a scattering power
        if sp:
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

# End class MolAtomParSet

class ObjCrystParSet(ParameterSet):
    """A wrapper for ObjCryst Crystal instance.

    Attributes:

    scatterers   --  The list of Scatterers (either Atoms or Molecules).
    a, b, c, alpha, beta, gamma --  Lattice parameters (Parameter)
    
    """

    def __init__(self, cryst, name):
        """Initialize

        cryst   --  An ObjCryst Crystal instance.
        name    --  A name for this
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
        from diffpy.Structure import SpaceGroups
        sgnum = self.GetSpaceGroup().GetNumber()
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

        cdict = {}
        for s in cryst:
            name = s.GetName()
            if not name:
                name = "noname"

            i = cdict.get(name, 0)
            sname = name
            # Only add the number if there are objects with the same name. We'll
            # fix the first one later.
            if i:
                sname += "_%i"%i
            cdict[name] = i+1

            # Now create the proper object
            cname = s.GetClassName()
            if cname == "Atom":
                parset = AtomParSet(sname, s, self)
            elif cname == "Molecule":
                parset = MoleculeParSet(sname, s, self)
            else:
                raise AttributeError("Unrecognized scatterer '%s'"%cname)

            self.addParameterSet(parset)
            self.atoms.append(atom)


        return

# End class ObjCrystParSet

# FIXME - flesh these out

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

        m = self.atom1.scatt.GetMolecule()
        res = m.AddBond(atom1.scatt, atom2.scatt, length, sigma, delta)

        MoleculeRestraint.__init__(self, res, scaled)
        return

# End class BondLengthRestraint

class BondAngleRestraint(MoleculeRestraint):
    """Restrain the angle defined by three atoms.

    Attributes:
    atom1   --  The first atom in the angle (MolAtomParSet)
    atom2   --  The second atom in the angle (MolAtomParSet)
    atom3   --  The third atom in the angle (MolAtomParSet)
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

        m = self.atom1.scatt.GetMolecule()
        res = m.AddBondAngle(atom1.scatt, atom2.scatt, atom3.scatt, angle,
                sigma, delta)

        MoleculeRestraint.__init__(self, res, scaled)
        return

# End class BondAngleRestraint

class DihedralAngleRestraint(MoleculeRestraint):
    """Restrain the dihedral (torsion) angle defined by four atoms.

    Attributes:
    atom1   --  The first atom in the angle (MolAtomParSet)
    atom2   --  The second (central) atom in the angle (MolAtomParSet)
    atom3   --  The third (central) atom in the angle (MolAtomParSet)
    atom4   --  The fourth atom in the angle (MolAtomParSet)
    res     --  The pyobjcryst Molecule restraint
    scaled  --  A flag indicating if the restraint is scaled (multiplied) by
                the unrestrained point-average chi^2 (chi^2/numpoints) (default
                False)
    """

    def __init__(self, atom1, atom2, atom3, angle, sigma, delta, scaled =
            False):
        """Create a dihedral angle restraint.

        atom1   --  First atom (MolAtomParSet) in the angle
        atom2   --  Second (central) atom (MolAtomParSet) in the angle
        atom3   --  Third (central) atom (MolAtomParSet) in the angle
        atom4   --  Fourth atom in the angle (MolAtomParSet)
        angle   --  The dihedral angle (radians)
        sigma   --  The uncertainty of the bond angle (radians)
        delta   --  The width of the bond angle (radians)
        scaled  --  A flag indicating if the restraint is scaled (multiplied)
                    by the unrestrained point-average chi^2 (chi^2/numpoints)
                    (default False).

        """
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4

        m = self.atom1.scatt.GetMolecule()
        res = m.AddDihedralAngle(atom1.scatt, atom2.scatt, atom3.scatt,
                atom4.scatt, angle, sigma, delta)

        MoleculeRestraint.__init__(self, res, scaled)
        return

# End class DihedralAngleRestraint

class StretchModeParameter(object):
    """Mix-in class for Parameters that encapsulate pyobjcryst.StretchModes."""

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
        self.mutated.extend(atomlist)
        scats = [a.scat for a in atomlist]
        self.mode.AddAtoms(scats)
        return

    def click():
        """Click the clickers of all the mutated parameters."""
        # Click the atoms that have moved
        for a in self.matoms:
            a.clicker.click()
        # Click the molecule position
        self.molecule.x.clicker.click()
        self.molecule.y.clicker.click()
        self.molecule.z.clicker.click()
        return

# EndClass StretchModeParameter

class BondLengthParameter(Parameter, StretchModeParameter):
    """Class for abstracting a bond length in a Molecule to a Parameter.

    This wraps up a pyobjcryst.StretchModeBondLength object so that the
    distance between two atoms in a Molecule can be used as an adjustable
    parameter. When a bond length is adjusted, the second atom is moved, and
    the absolute position of the molecule is altered to preserve the location
    of the center of mass within the crystal.
    
    Note that by making a bond constant it also makes the underlying atoms
    constant. When setting it as nonconstant, each atom is set nonconstant.
    Changing the bond length changes the position of the second atom and
    molecule, even if either is set as constant. 

    Attributes
    atom1   --  The first AtomParSet in the bond
    atom2   --  The second (mutated) AtomParSet in the bond
    matoms  --  The list of all mutated AtomParSets
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

    def __init__(self, name, atom1, atom2, value = None, const = False):
        """Create a BondLengthParameter.

        name    --  The name of the bond length parameter
        atom1   --  The first atom (AtomParSet) in the bond
        atom2   --  The second (mutated) atom (AtomParSet) in the bond
        value   --  An initial value for the bond length. If this is None
                    (default), then the current distance between the atoms will
                    be used.
        const   --  A flag indicating whether the Parameter is constant
                    (default False)

        """
        if value is None:
            value = GetBondLength(atom1.scat, atom2.scat)

        self.calclicker = Clicker()
        self.calclicker.click()

        Parameter.__init__(self, name, value, const)

        atom1.setConst(const)
        atom2.setConst(const)

        # Create the mode
        self.mode = StretchModeBondLength(atom1.scat, atom2.scat, None)
        # Add the second atom so it will be moved
        self.mode.AddAtom(atom2.scat)

        self.clicker.addSubject(atom1.clicker)
        self.clicker.addSubject(atom2.clicker)

        self.atom1 = atom1
        self.atom2 = atom2
        self.molecule = atom1.parent
        self.matoms = [atom2]

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
        Parameter.setConst(self, const, value)
        self.atom1.setConst(const)
        self.atom2.setConst(const)
        return

    def setValue(self, val):
        """Change the value of the bond length.

        This changes the position of the second atom and the absolute position
        of the molecule in such a way that the center of mass of the molecule
        does not change.

        """
        curval = self.getValue()
        if val == curval:
            return

        # The StretchModeBondLength expects the change in bond length.  This
        # will move atom2 by a distance delta, and then move the entire
        # molecule a distance 1/N * delta, where N is the number of atoms in
        # the molecule.
        delta = val - curval
        self.mode.Stretch(delta)

        self.value = val

        # Click everything that has changed
        self.click()
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
                     
class BondAngleParameter(Parameter, StretchModeParameter):
    """Class for abstracting a bond angle in a Molecule to a Parameter.

    This wraps up a pyobjcryst.StretchModeBondAngle object so that the
    angle defined by three atoms in a Molecule can be used as an adjustable
    parameter. When a bond angle is adjusted, the third atom is moved, and
    the absolute position of the molecule is altered to preserve the location
    of the center of mass within the crystal.
    
    Note that by making a bond angle constant it also makes the underlying
    atoms constant. When setting it as nonconstant, each atom is set
    nonconstant.  Changing the bond angle changes the position of the third
    atom and molecule, even if either is set as constant. 

    Attributes
    atom1   --  The first AtomParSet in the bond angle
    atom2   --  The second (central) AtomParSet in the bond angle
    atom3   --  The third (mutated) AtomParSet in the bond angle
    matoms  --  The list of all mutated AtomParSets
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

    def __init__(self, name, atom1, atom2, atom3, value = None, const = False):
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

        """
        if value is None:
            value = GetBondAngle(atom1.scat, atom2.scat, atom3.scat)

        self.calclicker = Clicker()
        self.calclicker.click()

        Parameter.__init__(self, name, value, const)

        atom1.setConst(const)
        atom2.setConst(const)
        atom3.setConst(const)

        # Create the stretch mode
        self.mode = StretchModeBondAngle(atom1.scat, atom2.scat, atom3.scat,
                None)
        # We only add the last atom. This is the one that will move
        self.mode.AddAtom(atom3.scat)

        self.clicker.addSubject(atom1.clicker)
        self.clicker.addSubject(atom2.clicker)
        self.clicker.addSubject(atom3.clicker)

        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.molecule = atom1.parent
        self.matoms = [atom3]

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
        Parameter.setConst(self, const, value)
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
                     
class DihedralAngleParameter(Parameter, StretchModeParameter):
    """Class for abstracting a dihedral angle in a Molecule to a Parameter.

    This wraps up a pyobjcryst.StretchModeTorsion object so that the angle
    defined by four atoms ([a1-a2].[a3-a4]) in a Molecule can be used as an
    adjustable parameter. When a dihedral angle is adjusted, the fourth atom is
    moved, and the absolute position of the molecule is altered to preserve the
    location of the center of mass within the crystal.
    
    Note that by making a dihedral angle constant it also makes the underlying
    atoms constant. When setting it as nonconstant, each atom is set
    nonconstant.  Changing the bond angle changes the position of the third
    atom and molecule, even if either is set as constant. 

    Attributes
    atom1   --  The first AtomParSet in the dihedral angle
    atom2   --  The second (central) AtomParSet in the dihedral angle
    atom3   --  The third (central) AtomParSet in the dihedral angle
    atom4   --  The fourth (mutated) AtomParSet in the dihedral angle
    matoms  --  The list of all mutated AtomParSets
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
            False):
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

        """
        if value is None:
            value = GetDihedralAngle(atom1.scat, atom2.scat, atom3.scat)

        self.calclicker = Clicker()
        self.calclicker.click()

        Parameter.__init__(self, name, value, const)

        atom1.setConst(const)
        atom2.setConst(const)
        atom3.setConst(const)
        atom4.setConst(const)

        # Create the stretch mode
        self.mode = StretchModeDihedralAngle(atom2.scat, atom3.scat, None)
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
        self.matoms = [atom4]

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
        Parameter.setConst(self, const, value)
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
                    self.atom3.scat)
            self.calclicker.click()

        return self.value

# End class DihedralAngleParameter

__id__ = "$Id$"
