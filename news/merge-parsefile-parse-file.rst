**Added:**

* <news item>

**Changed:**

* Change ``ProfileParser.parseFile`` to ``ProfileParser.parse_file``. Now, it dispatches to a subclass's ``parse_string`` override to perform format-specific parsing when one is provided (as ``PDFParser`` does), and otherwise reads generic column data.

**Deprecated:**

* Deprecate ``parseFile`` and ``parseString`` in ``ProfileParser``. Use ``parse_file`` and ``parse_string`` instead.

**Removed:**

* <news item>

**Fixed:**

* Fix ``PDFContribution.loadData`` discarding PDF metadata such as ``stype``, ``qmax``, ``qdamp`` and ``temperature``. It parsed with generic column data instead of dispatching to ``PDFParser``'s format-specific parsing, so ``BasePDFGenerator._process_metadata`` had nothing to configure the calculator with.

**Security:**

* <news item>
