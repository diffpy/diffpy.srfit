**Added:**

* ``diffpy.srfit.structure.structure_to_parameter_set`` and ``SrRealParSet.restrain_bvs``, which take the names they have in ``diffpy.cmistructure``.

**Changed:**

* The deprecation messages of ``constrainAsSpaceGroup`` and of ``setConst`` on the ObjCrystParSet bond length, bond angle and dihedral angle parameters now point to ``diffpy.cmistructure`` instead of ``diffpy.srfit.structure``.

**Deprecated:**

* ``diffpy.srfit.structure``, which moved to ``diffpy.cmistructure``. Importing it now raises a ``DeprecationWarning``.
* ``diffpy.srfit.structure.struToParameterSet``. Please use ``diffpy.cmistructure.structure_to_parameter_set`` instead.
* ``SrRealParSet.restrainBVS``. Please use ``diffpy.cmistructure.srrealparset.SrRealParSet.restrain_bvs`` instead.
* The ``.stru`` attribute of the structure parameter sets and of ``ObjCrystMoleculeParSet``. Please use ``.structure`` in ``diffpy.cmistructure`` instead.

**Removed:**

* <news item>

**Fixed:**

* <news item>

**Security:**

* <news item>
