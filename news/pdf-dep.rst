**Added:**

* ``set_qmax``, ``get_qmax``, ``set_qmin``, ``get_qmin``, ``set_scattering_type``, ``get_scattering_type``, ``set_structure`` and ``set_structure_from_parset`` to the PDF generators, and ``load_data``, ``set_calculation_range``, ``add_structure``, ``add_structure_from_parset``, ``set_qmax``, ``get_qmax``, ``set_qmin``, ``get_qmin``, ``set_scattering_type`` and ``get_scattering_type`` to ``PDFContribution``, which take the names they have in ``diffpy.cmipdf``.

**Changed:**

* The deprecation messages of the camelCase characteristic functions (``sphericalCF``, ``shellCF`` and the others) now point to ``diffpy.cmipdf.characteristicfunctions`` instead of ``diffpy.srfit.pdf.characteristicfunctions``. ``lognormalSphericalCF``, ``shellCF`` and ``shellCF2`` point to the names srfit 3.4.0 gave them, ``lognormal_spherical_particle`` and ``shell_particle``.

**Deprecated:**

* ``diffpy.srfit.pdf``, which moved to ``diffpy.cmipdf``. Importing it now raises a ``DeprecationWarning``.
* The camelCase methods of the PDF generators and ``PDFContribution`` (``setQmax``, ``getQmax``, ``setQmin``, ``getQmin``, ``setScatteringType``, ``getScatteringType``, ``setStructure``, ``setPhase``, ``loadData``, ``setCalculationRange``, ``addStructure`` and ``addPhase``). Please use their snake_case names in ``diffpy.cmipdf`` instead. ``setPhase`` and ``addPhase`` are now ``set_structure_from_parset`` and ``add_structure_from_parset``, and ``add_structure`` and ``add_structure_from_parset`` take the structure or parset before the name.

**Removed:**

* <news item>

**Fixed:**

* <news item>

**Security:**

* <news item>
