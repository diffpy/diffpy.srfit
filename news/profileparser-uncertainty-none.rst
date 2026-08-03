**Added:**

* <news item>

**Changed:**

* Change ``ProfileParser.parse_file`` to set unavailable ``dx``/``dy`` uncertainties to ``None`` instead of ``0``, matching the format-specific parsers.

**Deprecated:**

* <news item>

**Removed:**

* <news item>

**Fixed:**

* Fix ``ProfileParser.parse_file`` producing zero-valued uncertainties for absent or invalid ``dx``/``dy`` columns, which gave non-finite residuals for any file without usable uncertainties instead of falling back to an unweighted fit.

**Security:**

* <news item>
