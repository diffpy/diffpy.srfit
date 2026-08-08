=============
Release notes
=============

.. current developments

3.3.0
=====

**Added:**

* Added ``add_parameter_set`` method to replace deprecated ``FitRecipe.addParameterSet``.
* Added ``add_parameter_set`` method to replace deprecated ``ParameterSet.addParameterSet``.
* Added `diffpy.structure` back to the requirements and run the test
* Added ``add_constraint`` method to ``RecipeOrganizer``.
* Added ``remove_constraint`` method to ``RecipeOrganizer``.
* Added ``add_soft_bounds`` method to ``RecipeOrganizer``.
* Added ``remove_soft_bounds`` method to ``RecipeOrganizer``.
* Added ``register_soft_bounds`` method to ``RecipeOrganizer``.
* Added ``clear_all_soft_bounds`` method to ``RecipeOrganizer``.
* Added ``get_equation_from_string`` method to ``RecipeOrganizer``.
* Add ``parse_file`` method to ``ProfileParser`` to parse a file directly with ``load_data`` from ``diffpy.utils``.
* Add ``_parse_metadata`` and ``_parse_data`` hooks to ``ProfileParser``, which subclasses override to support a new data format. ``parse_file`` is a template method that calls them.
* Add ``get_num_banks`` method to ``ProfileParser`` to replace ``getNumBanks``.
* Add ``select_bank`` method to ``ProfileParser`` to replace ``selectBank``.
* Add ``get_format`` method to ``ProfileParser`` to replace ``getFormat``.
* Add ``get_data`` method to ``ProfileParser`` to replace ``getData``.
* Add ``get_metadata`` method to ``ProfileParser`` to replace ``getMetaData``.
* Added ``get_equation`` method to ``BaseBuilder``.
* Added ``get_equation`` method to ``FitContribution``.
* Add ``plot_recipe`` method to ``FitRecipe``.
* Added ``initialize_recipe_with_results`` to ``FitRecipe``.
* Added ``iterate_over_parameters`` method.
* Added ``register_calculator`` method.
* Added ``register_function`` method.
* Added ``register_string_function`` method.
* Added ``evaluate_equation`` method.
* Added ``is_constrained`` method.
* Added ``get_constrained_parameters`` method.
* Added ``clear_all_constraints`` method.
* Added ``add_profile_generator`` to replace deprecated function.
* Added ``is_constant`` method to ``Parameter``.
* Added ``bound_range`` method to ``Parameter``.
* Added ``bound_window`` method to ``Parameter``.
* Added ``set_weight`` method to ``FitRecipe`` to replace ``setWeight``.
* Added ``get_fit_hooks`` method to ``FitRecipe`` to replace ``getFitHooks``.
* Added ``clear_fit_hooks`` method to ``FitRecipe`` to replace ``clearFitHooks``.
* Added ``pop_fit_hook`` method to ``FitRecipe`` to replace ``popFitHook``.
* Added ``push_fit_hook`` method to ``FitRecipe`` to replace ``pushFitHook``.
* Add ``spherical_particle``, ``spheroidal_particle``, ``lognormal_spherical_particle``, ``sheet_particle``, and ``shell_particle`` to ``diffpy.srfit.pdf.characteristicfunctions``, replacing ``sphericalCF``, ``spheroidalCF``, ``lognormalSphericalCF``, ``sheetCF``, and ``shellCF``.
* Add ``constrain_as_space_group`` to ``diffpy.srfit.structure.sgconstraints``, replacing ``constrainAsSpaceGroup``.
* Added ``get_residual_equation`` method to ``FitContribution``.
* Added ``set_residual_equation`` method to ``FitContribution``.
* Added ``set_value`` method to ``diffpy.srfit.fitbase.Parameter``.
* Added ``save_results`` method to ``FitResults``.
* Added ``print_results`` method to ``FitResults``.
* Added ``get_results_string`` method to ``FitResults``.
* Added ``create_new_variable`` method to ``FitRecipe``.
* Added ``delete_variable`` method to ``FitRecipe``.
* Added ``add_variable`` method to ``FitRecipe``.
* Added ``scalar_residual`` method to ``FitRecipe``.
* Added ``remove_parameter_set`` method to ``FitRecipe``.
* Added test for ``remove_parameter_set`` method in ``FitRecipe``.
* Added initialize_recipe_from_recipe to ``FitRecipe``.
* Added ``convert_bounds_to_restraints`` method to ``FitRecipe``.
* Added ``get_bounds_pairs`` method to ``FitRecipe``.
* Added ``get_bounds_array`` method to ``FitRecipe``.
* Added ``get_names`` method to ``FitRecipe`` and ``RecipeContainer``.
* Added ``get_values`` method to ``FitRecipe`` and ``RecipeContainer``.
* Added ``is_free`` method to ``FitRecipe``.
* Added ``set_equation`` method to replace ``setEquation``.
* Added ``FitResults.get_results_dictionary`` in replace of ``resultsDictionary``.
* Added ``add_contribution`` in replace of ``addContribution``.
* Support for Python 3.14
* Added ``load_parsed_data`` in replace of ``loadParsedData`` in ``Profile`` and ``SimpleRecipe``.
* Added ``set_observed_profile`` in replace of ``setObservedProfile`` in ``Profile`` and ``SimpleRecipe``.
* Added ``set_calculation_range`` in replace of ``setCalculationRange`` in ``Profile`` and ``SimpleRecipe``.
* Added ``set_calculation_points`` in replace of ``setCalculationPoints`` in ``Profile`` and ``SimpleRecipe``.

**Changed:**

* Change ``PDFParser`` to inherit ``_parse_metadata`` and ``_parse_data`` unchanged from ``ProfileParser``, since PDFgetX and PDFgetN headers are now plain ``name = value`` pairs. ``PDFParser.parse_file`` returns the same data and metadata that ``PDFParser.parseFile`` did.
* Change ``ProfileParser.parse_file`` to accept the keyword arguments of ``load_data``, such as ``usecols``, ``delimiter`` and ``comments``. Use ``usecols`` to select columns from a file with more than four of them.
* Change ``PDFContribution.loadData`` to take a file name only.
* Change the headers of the ``si-q27r60-xray.gr`` and ``ni-q27r100-neutron.gr`` test files to the modern ``diffpy.pdfgetx`` and xPDFsuite configuration formats, replacing the 2008-era PDFgetX2/PDFgetN headers. The data values are unchanged.
* Change ``spheroidalCF2`` and ``shellCF2`` in ``diffpy.srfit.pdf.characteristicfunctions`` to be deprecated aliases that convert their arguments and call ``spheroidal_particle`` and ``shell_particle``, instead of separate implementations.
* Change the current characteristic functions in ``diffpy.srfit.pdf.characteristicfunctions`` (``spherical_particle``, ``spheroidal_particle``, ``lognormal_spherical_particle``, ``sheet_particle`` and ``shell_particle``) to emit a ``RuntimeWarning`` when a non-physical shape parameter makes them return zero, since a zero return flattens the fit residual and stalls a refinement without any other sign to the user. The warning is issued once per process for each distinct problem so a refinement loop does not repeat it. The deprecated camel-case functions forward to these functions directly, so they emit the same warning and return the same result for non-physical input.
* Change ``lognormal_spherical_particle`` in ``diffpy.srfit.pdf.characteristicfunctions`` to return zero for a negative ``particle_diameter_sigma`` instead of silently returning ``spherical_particle``. A ``particle_diameter_sigma`` of zero is still the sphere limit.
* Change ``sheet_particle`` in ``diffpy.srfit.pdf.characteristicfunctions`` to return an array of zeros for a non-positive ``sheet_thickness`` when ``r`` is an array, instead of a scalar zero.
* Change ``ProfileParser.parse_file`` to set unavailable ``dx``/``dy`` uncertainties to ``None`` instead of ``0``, matching the format-specific parsers.
* Change ``Profile.set_observed_profile`` to leave ``dyobs`` as ``None`` when no uncertainties are observed, instead of silently replacing them with an array of ones. The calculated ``dy`` still defaults to 1 at every calculation point, so unweighted fits refine exactly as before.
* Change ``Profile._validate`` to accept a ``None`` ``dyobs``, since observed uncertainties are optional. ``x``, ``y``, ``dy``, ``xobs`` and ``yobs`` are still required.
* Changed `iterPars` method to match all equal-type atoms to have same ADPs

**Deprecated:**

* Deprecated ``addParameterSet`` method in ``FitRecipe``.
* Deprecate ``addParameterSet`` method in ``ParameterSet``.
* Deprecated ``constrain`` method of ``RecipeOrganizer``. Use ``add_constraint`` instead.
* Deprecated ``unconstrain`` method of ``RecipeOrganizer``. Use ``remove_constraint`` instead.
* Deprecated ``restrain`` method of ``RecipeOrganizer``. Use ``add_soft_bounds`` instead.
* Deprecated ``unrestrain`` methods of ``RecipeOrganizer``. Use ``remove_soft_bounds`` instead.
* Deprecated ``addRestraint`` method of ``RecipeOrganizer``. Use ``register_soft_bounds`` instead.
* Deprecate ``clearRestraints`` method of ``RecipeOrganizer``. Use ``clear_all_soft_bounds`` instead.
* Deprecated ``equationFromString`` method of ``RecipeOrganizer``. Use ``get_equation_from_string`` instead.
* Deprecate ``getNumBanks``, ``selectBank``, ``getFormat``, ``getData``, and ``getMetaData`` in ``ProfileParser``.
* Deprecate ``ProfileParser.parseFile``, which now delegates to ``parse_file``. Use ``parse_file`` instead.
* Deprecated ``BaseBuilder.getEquation`` method for removal in 4.0.0.
* Deprecated ``FitContributution.getEquation`` method for removal in 4.0.0.
* Deprecated ``iterPars`` method. Use ``iterate_over_parameters`` instead.
* Deprecated ``registerCalculator`` method. Use ``register_calculator`` instead.
* Deprecated ``registerFunction`` method. Use ``register_function`` instead.
* Deprecated ``registerStringFunction`` method. Use ``register_string_function`` instead.
* Deprecated ``evaluateEquation`` method. Use ``evaluate_equation`` instead.
* Deprecated ``isConstrained`` method. Use ``is_constrained`` instead.
* Deprecated ``getConstrainedPars`` method. Use ``get_constrained_parameters`` instead.
* Deprecated ``clearConstraints`` method. Use ``clear_all_constraints`` instead.
* Deprecate ``setProfile`` for removal in version 4.0.0.
* Deprecated ``addProfileGenerator`` for removal in 4.0.0.
* Deprecated ``isConst`` method of ``Parameter``. Use ``is_constant`` instead.
* Deprecated ``boundRange`` method of ``Parameter``. Use ``bound_range`` instead.
* Deprecated ``boundWindow`` method of ``Parameter``. Use ``bound_window`` instead.
* Deprecated ``setWeight`` in ``FitRecipe`` for removal in 4.0.0.
* Deprecated ``getFitHooks`` in ``FitRecipe`` for removal in 4.0.0.
* Deprecated ``clearFitHooks`` in ``FitRecipe`` for removal in 4.0.0.
* Deprecated ``popFitHook`` in ``FitRecipe`` for removal in 4.0.0.
* Deprecated ``pushFitHook`` in ``FitRecipe`` for removal in 4.0.0.
* Deprecate ``sphericalCF``, ``spheroidalCF``, ``spheroidalCF2``, ``lognormalSphericalCF``, ``sheetCF``, ``shellCF``, and ``shellCF2`` in ``diffpy.srfit.pdf.characteristicfunctions``.
* Deprecate ``constrainAsSpaceGroup`` in ``diffpy.srfit.structure.sgconstraints``.
* Deprecated ``getResidualEquation`` method of ``FitContribution`` for removal in 4.0.0.
* Deprecated ``setResidualEquation`` method of ``FitContribution`` for removal in 4.0.0.
* Deprecated ``setValue`` method in ``diffpy.srfit.fitbase.Parameter`` for removal in 4.0.0.
* Deprecated ``saveResults`` method for removal in 4.0.0.
* Deprecated ``printResults`` method for removal in 4.0.0.
* Deprecated ``formatResults`` method for removal in 4.0.0.
* Deprecated ``newVar`` method for removal in 4.0.0.
* Deprecated ``delVar`` method for removal in 4.0.0.
* Deprecated ``addVar`` method for removal in 4.0.0.
* Deprecated ``scalarResidual`` method for removal in 4.0.0.
* Deprecated ``removeParameterSet`` method for removal in 4.0.0.
* Deprecated ``boundsToRestraints`` method for removal in 4.0.0.
* Deprecated ``getBounds`` method for removal in 4.0.0.
* Deprecated ``getBounds2`` method for removal in 4.0.0.
* Deprecated ``getNames`` method for removal in 4.0.0.
* Deprecated ``getValues`` method for removal in 4.0.0.
* Deprecated ``isFree`` method for removal in 4.0.0.
* Deprecated ``setEquation`` for removal in 4.0.0.
* Deprecated ``resultsDictionary`` for removal in 4.0.0.
* Deprecated ``addContribution`` for removal in 4.0.0.
* Deprecated ``loadParsedData`` in ``Profile`` and ``SimpleRecipe`` for removal in 4.0.0.
* Deprecated ``setObservedProfile`` in ``Profile`` and ``SimpleRecipe`` for removal in 4.0.0.
* Deprecated ``setCalculationRange`` in ``Profile`` and ``SimpleRecipe`` for removal in 4.0.0.
* Deprecated ``setCalculationPoints`` in ``Profile`` and ``SimpleRecipe`` for removal in 4.0.0.

**Fixed:**

* Fixed load new parsed data with updated `Qmax` attribute
* Fix ``PDFParser.parse_file`` and ``PDFContribution.loadData`` discarding PDF metadata. They parsed PDFgetX and PDFgetN headers as generic ``name = value`` pairs, which dropped ``stype``, ``qmin``, ``qmax``, ``qdamp`` and ``qbroad``, so ``PDFGenerator`` silently fell back to its default scattering type and Q range.
* Update codebase to scikit-package 0.3.0 standards.
* Updated project to the latest scikit-package template
* Fix ``shell_particle`` in ``diffpy.srfit.pdf.characteristicfunctions`` returning a smoothly varying branch of unphysical values, some below negative one, for a negative ``thickness`` or ``radius``, which a refinement could mistake for a real solution. It now returns zero.
* Fix ``shell_particle`` in ``diffpy.srfit.pdf.characteristicfunctions`` returning one at every ``r`` for a zero ``thickness``, where the normalizing denominator vanishes. It now returns zero.
* Fix ``spheroidal_particle`` in ``diffpy.srfit.pdf.characteristicfunctions`` raising ``ZeroDivisionError`` for a zero ``equatorial_radius``. It now returns zero.
* Fix ``ProfileParser.parse_file`` producing zero-valued uncertainties for absent or invalid ``dx``/``dy`` columns, which gave non-finite residuals for any file without usable uncertainties instead of falling back to an unweighted fit.
* Fix ``Profile.dyobs`` reporting an array of ones for data that carries no uncertainties, which made an unweighted profile indistinguishable from one whose uncertainties were genuinely all 1.

**Removed:**

* Remove ``parseString`` from ``ProfileParser``, and ``PDFParser``. Parsers read from a file, so a format is now supported by overriding the ``_parse_metadata`` and ``_parse_data`` hooks.
* Remove support for passing a data string or an open file object to ``PDFContribution.loadData``. Pass a file name instead.
* Support for Python 3.11


3.2.0
=====

**Changed:**

* Temporarily removed support for SAS characteristic functions until we can migrate to the new sasview api.
* doc/ -> docs/
* requirements/test.txt -> requirements/tests.txt

**Fixed:**

* Better license designation in `pyproject.toml`


3.1.0
=====


Version 3.0.0 – 2019-03-14
--------------------------

**Added:**

* Support for Python 3.7, 3.6, 3.5 in addition to 2.7.

**Changed:**

* Always use lower-case imports from `diffpy.structure`.
* Use numeric-value sort to order variables in `PrintFitHook`.

**Deprecated:**

* Variable `__gitsha__` in the `version` module renamed to `__git_commit__`.

**Removed:**

* Optional upper and lower-bound arguments in `Parameter.setValue`.
  The bounds can be set with `Parameter.boundRange` instead.
* Unused classes `ListOperator`, `SetOperator`.

**Fixed:**

* Metadata retrieval from `PDFContribution` hierarchy.
* Refresh `PDFGenerator` when its `rgrid` is changed in-place.
* Zero division in the `nppdfsas.py` example.
* Crash in `ellipsoidsas.py` example because of bug in `Parameter.setValue`.
* Pickling of `ProfileGenerator` objects and of bound class methods.
* Invalid escape sequences in string values.
