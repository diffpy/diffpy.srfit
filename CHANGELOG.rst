=============
Release notes
=============

.. current developments

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
