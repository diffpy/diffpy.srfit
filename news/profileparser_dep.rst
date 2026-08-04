**Added:**

* Add ``parse_file`` method to ``ProfileParser`` to parse a file directly with ``load_data`` from ``diffpy.utils``.
* Add ``_parse_metadata`` and ``_parse_data`` hooks to ``ProfileParser``, which subclasses override to support a new data format. ``parse_file`` is a template method that calls them.
* Add ``get_num_banks`` method to ``ProfileParser`` to replace ``getNumBanks``.
* Add ``select_bank`` method to ``ProfileParser`` to replace ``selectBank``.
* Add ``get_format`` method to ``ProfileParser`` to replace ``getFormat``.
* Add ``get_data`` method to ``ProfileParser`` to replace ``getData``.
* Add ``get_metadata`` method to ``ProfileParser`` to replace ``getMetaData``.

**Changed:**

* Change ``PDFParser`` to read its data block with ``load_data`` from ``diffpy.utils``, by overriding only the ``_parse_metadata`` hook. ``PDFParser.parse_file`` returns the same data and metadata that ``PDFParser.parseFile`` did.
* Change ``ProfileParser.parse_file`` to accept the keyword arguments of ``load_data``, such as ``usecols``, ``delimiter`` and ``comments``. Use ``usecols`` to select columns from a file with more than four of them.
* Change ``PDFContribution.loadData`` to take a file name only.
* Change the headers of the ``si-q27r60-xray.gr`` and ``ni-q27r100-neutron.gr`` test files to the modern ``diffpy.pdfgetx`` and xPDFsuite configuration formats, replacing the 2008-era PDFgetX2/PDFgetN headers. The data values are unchanged.

**Deprecated:**

* Deprecate ``getNumBanks``, ``selectBank``, ``getFormat``, ``getData``, and ``getMetaData`` in ``ProfileParser``.
* Deprecate ``ProfileParser.parseFile``, which now delegates to ``parse_file``. Use ``parse_file`` instead.

**Removed:**

* Remove ``parseString`` from ``ProfileParser``, and ``PDFParser``. Parsers read from a file, so a format is now supported by overriding the ``_parse_metadata`` and ``_parse_data`` hooks.
* Remove support for passing a data string or an open file object to ``PDFContribution.loadData``. Pass a file name instead.

**Fixed:**

* Fix ``PDFParser.parse_file`` and ``PDFContribution.loadData`` discarding PDF metadata. They parsed PDFgetX and PDFgetN headers as generic ``name = value`` pairs, which dropped ``stype``, ``qmin``, ``qmax``, ``qdamp`` and ``qbroad``, so ``PDFGenerator`` silently fell back to its default scattering type and Q range.

**Security:**

* <news item>
