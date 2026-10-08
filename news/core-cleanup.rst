**Added:**

* <news item>

**Changed:**

* <news item>

**Deprecated:**

* <news item>

**Removed:**

* Removed the undeclared ``six`` dependency and the leftover Python 2
  compatibility code that used it.

**Fixed:**

* ``SASProfile`` now updates its wrapped SAS ``DataInfo`` object when the
  observed profile is set through ``set_observed_profile`` or
  ``load_parsed_data``. Previously only the deprecated ``setObservedProfile``
  did this.
* ``FitRecipe`` tests no longer depend on the ``diffpy.srfit.pdf`` package.

**Security:**

* <news item>
