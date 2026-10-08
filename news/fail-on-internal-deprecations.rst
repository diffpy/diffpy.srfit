**Added:**

* The test suite now fails if srfit's own source code calls one of srfit's
  deprecated names, so an internal caller missed during a deprecation is
  caught in CI.

**Changed:**

* <news item>

**Deprecated:**

* <news item>

**Removed:**

* <news item>

**Fixed:**

* Creating a ``FitRecipe`` no longer emits a ``DeprecationWarning`` about
  ``pushFitHook``.
* ``initializeRecipe`` now emits only its own ``DeprecationWarning`` instead of
  also warning about ``resultsDictionary``.

**Security:**

* <news item>
