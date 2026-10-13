"""Tests for the deprecated camelCase names in
diffpy.srfit.fitbase.fitresults.

Each test calls a deprecated name and checks that it warns. Where a new
name does the same job, the test also checks that both give the same
result on identical fresh objects.

Remove this file in version 4.0.0, together with the deprecated names.
"""

import re

import pytest

from diffpy.srfit import FitRecipe, FitResults
from diffpy.srfit.fitbase.fitresults import initializeRecipe, resultsDictionary


def expected_message(old_base, old_name, new_base, new_name):
    return re.escape(
        f"'{old_base}.{old_name}' is deprecated and will be removed in "
        f"version 4.0.0. Please use '{new_base}.{new_name}' instead."
    )


FITBASE = "diffpy.srfit.fitbase"
FIT_RESULTS = "diffpy.srfit.fitbase.FitResults"
FIT_RECIPE = "diffpy.srfit.fitbase.FitRecipe"


def _get_results_output(recipe, method_name, capsys, tmp_path):
    """Return the results text that a FitResults method gives for the
    recipe, without the save time and user that a saved file starts
    with."""
    recipe.fithooks[0].verbose = 0
    results = FitResults(recipe)
    method = getattr(results, method_name)
    if method_name in ("formatResults", "get_results_string"):
        return method(header="My Custom header")
    if method_name in ("printResults", "print_results"):
        capsys.readouterr()
        method(header="My Custom header")
        return capsys.readouterr().out
    results_file = tmp_path / f"{method_name}.res"
    method(results_file, header="My Custom header")
    return results_file.read_text().split("\n", 2)[2]


# FitResults: each deprecated method gives the same results text as its new
# name. See test_fitresults.py::test_results_output.
@pytest.mark.parametrize(
    "input_old_name, input_new_name",
    [
        # C1: formatResults returns the results text.
        # Expected: Warns and gives the same text as get_results_string.
        ("formatResults", "get_results_string"),
        # C2: printResults prints the results text.
        # Expected: Warns and prints the same text as print_results.
        ("printResults", "print_results"),
        # C3: saveResults saves the results text to a file.
        # Expected: Warns and saves the same text as save_results.
        ("saveResults", "save_results"),
    ],
)
def test_fit_results_deprecations(
    build_recipes_one_contribution,
    capsys,
    tmp_path,
    input_old_name,
    input_new_name,
):
    old_recipe, new_recipe = build_recipes_one_contribution
    with pytest.warns(
        DeprecationWarning,
        match=expected_message(
            FIT_RESULTS, input_old_name, FIT_RESULTS, input_new_name
        ),
    ):
        actual_old_result = _get_results_output(
            old_recipe, input_old_name, capsys, tmp_path
        )
    actual_new_result = _get_results_output(
        new_recipe, input_new_name, capsys, tmp_path
    )
    assert actual_old_result == actual_new_result


def test_resultsDictionary(temp_data_files):
    # C1: The results are read from a saved results file.
    # Expected: Warns and returns every "name number" pair in the file, with
    # the numbers kept as strings. This includes pairs that are not results,
    # such as "Feb 25" from the save date and "than 25" from the
    # correlations heading, and misses "Reduced Chi2", which is read as
    # "Chi2".
    with pytest.warns(
        DeprecationWarning,
        match=expected_message(
            FITBASE,
            "resultsDictionary",
            FIT_RESULTS,
            "get_results_dictionary",
        ),
    ):
        actual_results_dict = resultsDictionary(
            temp_data_files / "fit_results.res"
        )
    expected_results_dict = {
        "than": "25",
        "wave_number": "1.00000000e+00",
        "phase_shift": "-1.61291146e-18",
        "amplitude": "1.00000000e+00",
        "Rw": "0.00000000",
        "Chi2": "0.00000000",
        "Restraints": "0.00000000",
        "Contributions": "0.00000000",
        "Residual": "0.00000000",
        "Feb": "25",
    }
    assert actual_results_dict == expected_results_dict


# initializeRecipe: see
# test_fitrecipe.py::test_initialize_recipe_from_results_file.
@pytest.mark.parametrize(
    "input_kind",
    [
        # C1: The name of a results file.
        # Expected: Warns and sets A, sig and x0 to their values in the
        # file.
        "file_name",
        # C2: An open results file.
        # Expected: Warns, sets A, sig and x0 to their values in the file,
        # and leaves the file open.
        "file_object",
        # C3: The text of a results file.
        # Expected: Warns and sets A, sig and x0 to their values in the
        # text.
        "results_text",
    ],
)
def test_initializeRecipe(datafile, input_kind):
    recipe = FitRecipe("recipe")
    for name in ("A", "sig", "x0"):
        recipe.create_new_variable(name, 0)
    results_path = datafile("results.res")
    with open(results_path, "r") as results_file:
        inputs = {
            "file_name": results_path,
            "file_object": results_file,
            "results_text": results_path.read_text(),
        }
        with pytest.warns(
            DeprecationWarning,
            match=expected_message(
                FITBASE,
                "initializeRecipe",
                FIT_RECIPE,
                "initialize_recipe_with_results",
            ),
        ):
            initializeRecipe(recipe, inputs[input_kind])
        actual_file_closed = results_file.closed
    expected_file_closed = False
    assert actual_file_closed == expected_file_closed
    actual_values = [recipe.A.value, recipe.sig.value, recipe.x0.value]
    expected_values = [5.77619823e-01, -9.22758690e-01, 6.12422115e00]
    assert actual_values == pytest.approx(expected_values)
