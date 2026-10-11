**Added:**

* Add ``add_observer``, ``remove_observer`` and ``has_observer`` to
  ``Observable``, and ``has_tags`` and ``verify_tags`` to ``TagManager``.
* Add ``is_identifier`` and ``validate_name`` to ``diffpy.srfit.util.nameutils``,
  ``convert_input_to_string`` to ``diffpy.srfit.util.inpututils``, and
  ``sort_key_for_numeric_string`` to ``diffpy.srfit.util``.

**Changed:**

* The error raised for an invalid name now says how to fix it.

**Deprecated:**

* ``Observable.addObserver``, ``removeObserver`` and ``hasObserver``,
  ``TagManager.hasTags`` and ``verifyTags``, ``isIdentifier``,
  ``validateName``, ``inputToString`` and ``sortKeyForNumericString`` are
  deprecated in favour of their snake_case replacements. The old names still
  work and warn, and will be removed in version 4.0.0.

**Removed:**

* <news item>

**Fixed:**

* <news item>

**Security:**

* <news item>
