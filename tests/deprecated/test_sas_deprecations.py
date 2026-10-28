"""Tests for the deprecation of diffpy.srfit.sas, which moved to
diffpy.cmisas.

Remove this file in version 4.0.0, together with diffpy.srfit.sas.
"""

import importlib
from types import SimpleNamespace

import pytest

import diffpy.srfit.sas
from diffpy.srfit.pdf.characteristicfunctions import SASCF


def test_sas_package_import_warns():
    # Case: diffpy.srfit.sas is imported.
    # Expected: A DeprecationWarning points to diffpy.cmisas.
    expected_msg = (
        "'diffpy.srfit.sas' is deprecated and will be removed in version "
        "4.0.0. Please use 'diffpy.cmisas' instead."
    )
    with pytest.warns(DeprecationWarning) as record:
        importlib.reload(diffpy.srfit.sas)
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg


def test_SASCF():
    # Case: SASCF is created from a sas model with no parameters.
    # Expected: A DeprecationWarning points to
    # diffpy.cmisas.characteristicfunctions.SASCF, and a SASCF is still
    # returned.
    expected_msg = (
        "'diffpy.srfit.pdf.characteristicfunctions.SASCF' is deprecated and "
        "will be removed in version 4.0.0. Please use "
        "'diffpy.cmisas.characteristicfunctions.SASCF' instead."
    )
    model = SimpleNamespace(params={}, dispersion={})
    with pytest.warns(DeprecationWarning) as record:
        characteristic_function = SASCF("sphere", model)
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg
    actual_type = type(characteristic_function)
    expected_type = SASCF
    assert actual_type == expected_type
