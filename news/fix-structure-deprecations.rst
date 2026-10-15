**Added:**

* Document the ``*=``, ``%=``, ``+=``, ``|=``, and ``<<`` shorthand operators in the ``diffpy.srfit.interface`` API reference.

**Changed:**

* <news item>

**Deprecated:**

* <news item>

**Removed:**

* <news item>

**Fixed:**

* Use ``get_space_group``, ``position_formulas``, and ``u_formulas`` in place of the deprecated ``diffpy.structure`` functions ``GetSpaceGroup``, ``positionFormulas``, and ``UFormulas`` to silence deprecation warnings.
* Use ``add_constraint`` and ``add_soft_bounds`` in the ``*=`` and ``%=`` operators of ``RecipeOrganizerInterface`` so they no longer raise deprecation warnings for ``constrain`` and ``restrain``.
* Use ``pyobjcryst.crystal.create_crystal_from_cif``, ``scipy.optimize.leastsq``, and the ``particle_diameter`` argument of ``spherical_particle`` in the examples in place of deprecated or renamed APIs.

**Security:**

* <news item>
