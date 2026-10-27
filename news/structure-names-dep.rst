**Added:**

* The snake_case names of the structure adapter methods, which they have in ``diffpy.cmistructure``: ``can_adapt``, ``get_lattice``, ``get_scatterers``, ``get_space_group``, ``use_symmetry``, ``using_symmetry``, ``is_dummy``, ``has_scatterers``, ``wrap_restraints``, ``wrap_stretch_mode_parameters``, ``restrain_bond_length``, ``restrain_bond_angle``, ``restrain_dihedral_angle`` and their ``*_parameter`` forms, ``add_bond_length_parameter``, ``add_bond_angle_parameter``, ``add_dihedral_angle_parameter``, ``StretchModeParameter.add_atoms`` and ``BaseSpaceGroupParameters.add_parameter``.
* The attributes ``.angle_units``, ``.scatterer``, ``.restraint``, ``.mutated_atoms``, ``.space_group``, ``.structure_parameter_set`` and ``.index``, which they are called in ``diffpy.cmistructure``.

**Changed:**

* <news item>

**Deprecated:**

* The camelCase methods of the ``diffpy.srfit.structure`` adapters (``canAdapt``, ``getLattice``, ``getScatterers``, ``getSpaceGroup``, ``useSymmetry``, ``usingSymmetry``, ``isDummy``, ``hasScatterers``, ``wrapRestraints``, ``wrapStretchModeParameters``, the ``restrainBond*`` and ``restrainDihedral*`` methods, the ``add*Parameter`` methods of ``ObjCrystMoleculeParSet``, ``StretchModeParameter.addAtoms`` and ``BaseSpaceGroupParameters.addParameter``). Please use their snake_case names in ``diffpy.cmistructure`` instead.
* The attributes ``.angunits``, ``.scat``, ``.res``, ``.matoms``, ``.sg``, ``.strups`` and ``.idx``. Please use ``.angle_units``, ``.scatterer``, ``.restraint``, ``.mutated_atoms``, ``.space_group``, ``.structure_parameter_set`` and ``.index`` in ``diffpy.cmistructure`` instead.

**Removed:**

* <news item>

**Fixed:**

* ``ObjCrystMoleculeParSet.get_lattice`` no longer raises ``AttributeError``. It called ``newPar``, which ``ParameterSet`` doesn't have.
* ``ObjCrystMoleculeParSet.wrap_stretch_mode_parameters`` no longer raises ``AttributeError`` for a molecule with stretch modes. It called ``AddAtoms``, which ``StretchModeParameter`` doesn't have.

**Security:**

* <news item>
