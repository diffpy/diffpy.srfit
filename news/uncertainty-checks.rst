**Added:**

* <news item>

**Changed:**

* ``ProfileParser`` reads the ``dy`` column as it is in the file, instead of dropping the whole column when one uncertainty is zero, negative, infinite or not a number.
* ``FitResults`` counts only the points with weight in Reduced Chi2, and says how many points of each contribution have no weight.

**Deprecated:**

* <news item>

**Removed:**

* <news item>

**Fixed:**

* A point whose uncertainty is zero, negative, infinite or not a number gets no weight in a fit, instead of making the residual infinite or undefined. When no point in the calculation range has a valid uncertainty, the profile is weighted as if it had none.

**Security:**

* <news item>
