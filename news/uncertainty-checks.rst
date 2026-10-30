**Added:**

* <news item>

**Changed:**

* ``ProfileParser`` warns when it ignores a ``dy`` column whose uncertainties are zero, negative, infinite or not a number, instead of weighting every point equally without saying so.

**Deprecated:**

* <news item>

**Removed:**

* <news item>

**Fixed:**

* ``FitContribution.residual`` raises an ``SrFitError`` that counts the points whose uncertainties are zero or not a number, instead of returning an infinite or undefined residual, unless ``Profile.set_calculation_range`` excludes them.

**Security:**

* <news item>
