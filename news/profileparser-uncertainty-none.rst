**Added:**

* <news item>

**Changed:**

* Change ``ProfileParser.parse_file`` to set unavailable ``dx``/``dy`` uncertainties to ``None`` instead of ``0``, matching the format-specific parsers.
* Change ``Profile.set_observed_profile`` to leave ``dyobs`` as ``None`` when no uncertainties are observed, instead of silently replacing them with an array of ones. The calculated ``dy`` still defaults to 1 at every calculation point, so unweighted fits refine exactly as before.
* Change ``Profile._validate`` to accept a ``None`` ``dyobs``, since observed uncertainties are optional. ``x``, ``y``, ``dy``, ``xobs`` and ``yobs`` are still required.

**Deprecated:**

* <news item>

**Removed:**

* <news item>

**Fixed:**

* Fix ``ProfileParser.parse_file`` producing zero-valued uncertainties for absent or invalid ``dx``/``dy`` columns, which gave non-finite residuals for any file without usable uncertainties instead of falling back to an unweighted fit.
* Fix ``Profile.dyobs`` reporting an array of ones for data that carries no uncertainties, which made an unweighted profile indistinguishable from one whose uncertainties were genuinely all 1.

**Security:**

* <news item>
