import importlib.resources
import logging
import sys
from functools import lru_cache

import pytest
import six
from numpy import linspace, pi, sin

import diffpy.srfit.equation.literals as literals
from diffpy.srfit.fitbase import FitContribution, FitRecipe, Profile

logger = logging.getLogger(__name__)


@lru_cache()
def has_sas():
    try:
        import sas
        import sasmodels

        del sas
        del sasmodels
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
    """Fixture to load a test data file from the testdata package
    directory."""

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
    def _no_observers_in_global_builders():
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

    return _no_observers_in_global_builders()


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


@pytest.fixture()
def build_recipe_one_contribution():
    "helper to build a simple recipe"

    def _build_recipe():
        profile = Profile()
        x = linspace(0, pi, 10)
        y = sin(x)
        profile.set_observed_profile(x, y)
        contribution = FitContribution("c1")
        contribution.set_profile(profile)
        contribution.set_equation("amplitude*sin(wave_number*x + phase_shift)")
        recipe = FitRecipe()
        recipe.add_contribution(contribution)
        recipe.add_variable(contribution.amplitude, 1)
        recipe.add_variable(contribution.wave_number, 1)
        recipe.add_variable(contribution.phase_shift, 1)
        return recipe

    return _build_recipe


@pytest.fixture()
def build_recipe_two_contributions():
    """Helper to build a recipe with two physically related
    contributions."""
    profile1 = Profile()
    x = linspace(0, pi, 50)
    y1 = sin(x)  # amplitude=1, freq=1
    profile1.set_observed_profile(x, y1)

    contribution1 = FitContribution("c1")
    contribution1.set_profile(profile1)
    contribution1.set_equation("A*sin(k*x + c)")

    profile2 = Profile()
    y2 = 0.5 * sin(2 * x)  # amplitude=0.5, freq=2
    profile2.set_observed_profile(x, y2)

    contribution2 = FitContribution("c2")
    contribution2.set_profile(profile2)
    contribution2.set_equation("B*sin(m*x + d)")

    recipe = FitRecipe()
    recipe.add_contribution(contribution1)
    recipe.add_contribution(contribution2)

    # Add variables with reasonable initial guesses
    recipe.add_variable(contribution1.A, 0.8)
    recipe.add_variable(contribution1.k, 1.0)
    recipe.add_variable(contribution1.c, 0.1)

    recipe.add_variable(contribution2.B, 0.4)
    recipe.add_variable(contribution2.m, 2.0)
    recipe.add_variable(contribution2.d, 0.1)

    # ---- Meaningful constraints ----
    recipe.constrain(contribution2.m, "2*k")
    recipe.constrain(contribution2.d, contribution1.c)
    recipe.constrain(contribution2.B, "0.5*A")
    recipe.restrain(contribution1.A, 0.5, 1.5)
    recipe.restrain(contribution1.k, 0.8, 1.2)

    return recipe


@pytest.fixture
def temp_data_files(tmp_path):
    """
    Temporary directory containing:
    - data_with_meta.gr
    - data_without_meta.dat
    Each file contains a single line of data.
    """
    file_with_meta = tmp_path / "gr_file.gr"
    file_with_meta.write_text("1.0 2.0\n" "1.1 2.1\n" "1.2 2.2\n")

    dat_file = tmp_path / "dat_file.dat"
    dat_file.write_text("1.0 2.0\n" "1.1 2.1\n" "1.2 2.2\n")

    cgr_file = tmp_path / "cgr_file.cgr"
    cgr_file.write_text("1.0 2.0\n" "1.1 2.1\n" "1.2 2.2\n")

    results_file = tmp_path / "fit_results.res"
    results_file.write_text(
        """
Results written: Wed Feb 25 15:14:58 2026
produced by cadenmyers

Some quantities invalid due to missing profile uncertainty
Overall (Chi2 and Reduced Chi2 invalid)
------------------------------------------------------------------------------
Residual       0.00000000
Contributions  0.00000000
Restraints     0.00000000
Chi2           0.00000000
Reduced Chi2   0.00000000
Rw             0.00000000

Variables (Uncertainties invalid)
------------------------------------------------------------------------------
amplitude    1.00000000e+00 +/- 4.82804000e-01
phase_shift  -1.61291146e-18 +/- 1.00000000e+00
wave_number  1.00000000e+00 +/- 2.17496687e-01

Variable Correlations greater than 25% (Correlations invalid)
------------------------------------------------------------------------------
No correlations greater than 25%
"""
    )
    yield tmp_path
