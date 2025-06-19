import importlib.resources
import json
import logging
from functools import lru_cache
from pathlib import Path

import pytest

logger = logging.getLogger(__name__)


@lru_cache()
def has_sas():
    try:
        __import__("sas.pr.invertor")
        __import__("sas.models")
        return True
    except ImportError:
        return False


# diffpy.structure
@lru_cache()
def has_diffpy_structure():
    _msg_nostructure = "No module named 'diffpy.structure'"
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
    _msg_nopyobjcryst = "No module named 'pyobjcryst'"
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
    _msg_nosrreal = "No module named 'diffpy.srreal'"
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


@pytest.fixture
def user_filesystem(tmp_path):
    base_dir = Path(tmp_path)
    home_dir = base_dir / "home_dir"
    home_dir.mkdir(parents=True, exist_ok=True)
    cwd_dir = base_dir / "cwd_dir"
    cwd_dir.mkdir(parents=True, exist_ok=True)

    home_config_data = {"username": "home_username", "email": "home@email.com"}
    with open(home_dir / "diffpyconfig.json", "w") as f:
        json.dump(home_config_data, f)

    yield tmp_path


@pytest.fixture
def datafile():
    """Fixture to load a test data file from the testdata package directory."""

    def _datafile(filename):
        return importlib.resources.files(
            "diffpy.srfit.tests.testdata"
        ).joinpath(filename)

    return _datafile
