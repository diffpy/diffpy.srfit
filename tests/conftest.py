import importlib.resources
import logging
import sys
from functools import lru_cache

import pytest
import six

import diffpy.srfit.equation.literals as literals
from diffpy.srfit.sas.sasimport import sasimport

logger = logging.getLogger(__name__)


@lru_cache()
def has_sas():
    try:
        sasimport("sas.pr.invertor")
        sasimport("sas.models")
        return True
    except ImportError:
        return False


# diffpy.structure
@lru_cache()
def has_diffpy_structure():
    try:
        import diffpy.structure as m

        del m
        return True
    except ImportError:
        return False
        logger.warning(
            "Cannot import diffpy.structure, Structure tests skipped."
        )


@lru_cache()
def has_pyobjcryst():
    try:
        import pyobjcryst as m

        del m
        return True
    except ImportError:
        return False
        logger.warning("Cannot import pyobjcryst, pyobjcryst tests skipped.")


# diffpy.srreal


@lru_cache()
def has_diffpy_srreal():
    try:
        import diffpy.srreal.pdfcalculator as m

        del m
        return True
    except ImportError:
        return False
        logger.warning("Cannot import diffpy.srreal, PDF tests skipped.")


@pytest.fixture(scope="session")
def sas_available():
    return has_sas()


@pytest.fixture(scope="session")
def diffpy_structure_available():
    return has_diffpy_structure()


@pytest.fixture(scope="session")
def diffpy_srreal_available():
    return has_diffpy_srreal()


@pytest.fixture(scope="session")
def pyobjcryst_available():
    return has_pyobjcryst()


@pytest.fixture(scope="session")
def datafile():
    """Fixture to load a test data file from the testdata package directory."""

    def _datafile(filename):
        return importlib.resources.files("tests.testdata").joinpath(filename)

    return _datafile


@pytest.fixture(scope="session")
def make_args():
    def _makeArgs(num):
        args = []
        for i in range(num):
            j = i + 1
            args.append(literals.Argument(name="v%i" % j, value=j))
        return args

    return _makeArgs


@pytest.fixture(scope="session")
def noObserversInGlobalBuilders():
    def _noObserversInGlobalBuilders():
        """True if no observer function leaks to global builder objects.

        Ensure objects are not immortal due to a reference from static
        value.
        """
        from diffpy.srfit.equation.builder import _builders

        rv = True
        for n, b in _builders.items():
            if b.literal and b.literal._observers:
                rv = False
                break
        return rv

    return _noObserversInGlobalBuilders()


@pytest.fixture(scope="session")
def capturestdout():
    def _capturestdout(f, *args, **kwargs):
        """Capture the standard output from a call of function f."""
        savestdout = sys.stdout
        fp = six.StringIO()
        try:
            sys.stdout = fp
            f(*args, **kwargs)
        finally:
            sys.stdout = savestdout
        return fp.getvalue()

    return _capturestdout
