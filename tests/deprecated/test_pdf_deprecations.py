"""Tests for the deprecation of diffpy.srfit.pdf, which moved to
diffpy.cmipdf.

Remove this file in version 4.0.0, together with diffpy.srfit.pdf.
"""

import importlib

import pytest

import diffpy.srfit.pdf
from diffpy.srfit.pdf import DebyePDFGenerator, PDFContribution, PDFGenerator
from diffpy.srfit.structure.diffpyparset import DiffpyStructureParSet
from diffpy.structure import Atom, Lattice, Structure


@pytest.fixture
def nacl_structure():
    lattice = Lattice(5.6, 5.6, 5.6, 90, 90, 90)
    atoms = [Atom("Na", [0, 0, 0]), Atom("Cl", [0.5, 0.5, 0.5])]
    structure = Structure(atoms, lattice=lattice)
    return structure


def test_pdf_package_import_warns():
    # Case: diffpy.srfit.pdf is imported.
    # Expected: A DeprecationWarning points to diffpy.cmipdf.
    expected_msg = (
        "'diffpy.srfit.pdf' is deprecated and will be removed in version "
        "4.0.0. Please use 'diffpy.cmipdf' instead."
    )
    with pytest.warns(DeprecationWarning) as record:
        importlib.reload(diffpy.srfit.pdf)
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg


@pytest.mark.parametrize(
    "old_setter, new_getter, input_value, expected_msg",
    [
        # Test the deprecated setters of BasePDFGenerator.
        # C1: setQmax sets qmax to 21.
        # Expected: Warns that set_qmax in diffpy.cmipdf replaces it, and
        # get_qmax returns 21.
        (
            "setQmax",
            "get_qmax",
            21.0,
            "'diffpy.srfit.pdf.basepdfgenerator.BasePDFGenerator.setQmax' is "
            "deprecated and will be removed in version 4.0.0. Please use "
            "'diffpy.cmipdf.basepdfgenerator.BasePDFGenerator.set_qmax' "
            "instead.",
        ),
        # C2: setQmin sets qmin to 0.5.
        # Expected: Warns that set_qmin in diffpy.cmipdf replaces it, and
        # get_qmin returns 0.5.
        (
            "setQmin",
            "get_qmin",
            0.5,
            "'diffpy.srfit.pdf.basepdfgenerator.BasePDFGenerator.setQmin' is "
            "deprecated and will be removed in version 4.0.0. Please use "
            "'diffpy.cmipdf.basepdfgenerator.BasePDFGenerator.set_qmin' "
            "instead.",
        ),
        # C3: setScatteringType sets the scattering type to neutron.
        # Expected: Warns that set_scattering_type in diffpy.cmipdf
        # replaces it, and get_scattering_type returns "N".
        (
            "setScatteringType",
            "get_scattering_type",
            "N",
            "'diffpy.srfit.pdf.basepdfgenerator.BasePDFGenerator."
            "setScatteringType' is deprecated and will be removed in version "
            "4.0.0. Please use 'diffpy.cmipdf.basepdfgenerator."
            "BasePDFGenerator.set_scattering_type' instead.",
        ),
    ],
)
def test_generator_setters(
    diffpy_srreal_available, old_setter, new_getter, input_value, expected_msg
):
    if not diffpy_srreal_available:
        pytest.skip("diffpy.srreal package not available")
    generator = PDFGenerator()
    with pytest.warns(DeprecationWarning) as record:
        getattr(generator, old_setter)(input_value)
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg
    actual_value = getattr(generator, new_getter)()
    expected_value = input_value
    assert actual_value == expected_value


@pytest.mark.parametrize(
    "old_getter, new_getter, expected_msg",
    [
        # Test the deprecated getters of BasePDFGenerator.
        # C1: getQmax is called on a new generator.
        # Expected: Warns that get_qmax in diffpy.cmipdf replaces it, and
        # returns the same value as get_qmax.
        (
            "getQmax",
            "get_qmax",
            "'diffpy.srfit.pdf.basepdfgenerator.BasePDFGenerator.getQmax' is "
            "deprecated and will be removed in version 4.0.0. Please use "
            "'diffpy.cmipdf.basepdfgenerator.BasePDFGenerator.get_qmax' "
            "instead.",
        ),
        # C2: getQmin is called on a new generator.
        # Expected: Warns that get_qmin in diffpy.cmipdf replaces it, and
        # returns the same value as get_qmin.
        (
            "getQmin",
            "get_qmin",
            "'diffpy.srfit.pdf.basepdfgenerator.BasePDFGenerator.getQmin' is "
            "deprecated and will be removed in version 4.0.0. Please use "
            "'diffpy.cmipdf.basepdfgenerator.BasePDFGenerator.get_qmin' "
            "instead.",
        ),
        # C3: getScatteringType is called on a new generator.
        # Expected: Warns that get_scattering_type in diffpy.cmipdf
        # replaces it, and returns the same value as get_scattering_type.
        (
            "getScatteringType",
            "get_scattering_type",
            "'diffpy.srfit.pdf.basepdfgenerator.BasePDFGenerator."
            "getScatteringType' is deprecated and will be removed in version "
            "4.0.0. Please use 'diffpy.cmipdf.basepdfgenerator."
            "BasePDFGenerator.get_scattering_type' instead.",
        ),
    ],
)
def test_generator_getters(
    diffpy_srreal_available, old_getter, new_getter, expected_msg
):
    if not diffpy_srreal_available:
        pytest.skip("diffpy.srreal package not available")
    generator = PDFGenerator()
    with pytest.warns(DeprecationWarning) as record:
        actual_value = getattr(generator, old_getter)()
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg
    expected_value = getattr(generator, new_getter)()
    assert actual_value == expected_value


@pytest.mark.parametrize(
    "generator_class, expected_periodic, expected_msg",
    [
        # Test the deprecated setStructure of each PDF generator.
        # C1: setStructure is called on a PDFGenerator.
        # Expected: Warns that set_structure in diffpy.cmipdf replaces it
        # and that 'stru' is now 'structure', and adapts the structure as
        # periodic.
        (
            PDFGenerator,
            True,
            "'diffpy.srfit.pdf.basepdfgenerator.BasePDFGenerator."
            "setStructure' is deprecated and will be removed in version "
            "4.0.0. Please use "
            "'diffpy.cmipdf.basepdfgenerator.BasePDFGenerator.set_structure' "
            "instead. Additionally, the signature has changed. Please pass "
            "the parameter 'stru' as 'structure'.",
        ),
        # C2: setStructure is called on a DebyePDFGenerator.
        # Expected: Warns that set_structure in diffpy.cmipdf replaces it
        # and that 'stru' is now 'structure', and adapts the structure as
        # non-periodic.
        (
            DebyePDFGenerator,
            False,
            "'diffpy.srfit.pdf.debyepdfgenerator.DebyePDFGenerator."
            "setStructure' is deprecated and will be removed in version "
            "4.0.0. Please use 'diffpy.cmipdf.debyepdfgenerator."
            "DebyePDFGenerator.set_structure' instead. Additionally, the "
            "signature has changed. Please pass the parameter 'stru' as "
            "'structure'.",
        ),
    ],
)
def test_setStructure(
    diffpy_srreal_available,
    nacl_structure,
    generator_class,
    expected_periodic,
    expected_msg,
):
    if not diffpy_srreal_available:
        pytest.skip("diffpy.srreal package not available")
    generator = generator_class()
    with pytest.warns(DeprecationWarning) as record:
        generator.setStructure(nacl_structure)
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg
    actual_structure = generator.stru
    expected_structure = nacl_structure
    assert actual_structure == expected_structure
    actual_periodic = generator.phase.using_symmetry()
    assert actual_periodic == expected_periodic


@pytest.mark.parametrize(
    "generator_class, expected_periodic, expected_msg",
    [
        # Test the deprecated setPhase of each PDF generator.
        # C1: setPhase is called on a PDFGenerator.
        # Expected: Warns that set_structure_from_parset in diffpy.cmipdf
        # replaces it, and manages the parset as periodic.
        (
            PDFGenerator,
            True,
            "'diffpy.srfit.pdf.basepdfgenerator.BasePDFGenerator.setPhase' is "
            "deprecated and will be removed in version 4.0.0. Please use "
            "'diffpy.cmipdf.basepdfgenerator.BasePDFGenerator."
            "set_structure_from_parset' instead.",
        ),
        # C2: setPhase is called on a DebyePDFGenerator.
        # Expected: Warns that set_structure_from_parset in diffpy.cmipdf
        # replaces it, and manages the parset as non-periodic.
        (
            DebyePDFGenerator,
            False,
            "'diffpy.srfit.pdf.debyepdfgenerator.DebyePDFGenerator.setPhase' "
            "is deprecated and will be removed in version 4.0.0. Please use "
            "'diffpy.cmipdf.debyepdfgenerator.DebyePDFGenerator."
            "set_structure_from_parset' instead.",
        ),
    ],
)
def test_setPhase(
    diffpy_srreal_available,
    nacl_structure,
    generator_class,
    expected_periodic,
    expected_msg,
):
    if not diffpy_srreal_available:
        pytest.skip("diffpy.srreal package not available")
    generator = generator_class()
    parameter_set = DiffpyStructureParSet("phase", nacl_structure)
    with pytest.warns(DeprecationWarning) as record:
        generator.setPhase(parameter_set)
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg
    actual_phase = generator.phase
    expected_phase = parameter_set
    assert actual_phase == expected_phase
    actual_periodic = generator.phase.using_symmetry()
    assert actual_periodic == expected_periodic


@pytest.mark.parametrize(
    "old_setter, new_getter, input_value, expected_msg",
    [
        # Test the deprecated setters of PDFContribution.
        # C1: setQmax sets qmax to 21.
        # Expected: Warns that set_qmax in diffpy.cmipdf replaces it, and
        # get_qmax returns 21.
        (
            "setQmax",
            "get_qmax",
            21.0,
            "'diffpy.srfit.pdf.pdfcontribution.PDFContribution.setQmax' is "
            "deprecated and will be removed in version 4.0.0. Please use "
            "'diffpy.cmipdf.pdfcontribution.PDFContribution.set_qmax' "
            "instead.",
        ),
        # C2: setQmin sets qmin to 0.5.
        # Expected: Warns that set_qmin in diffpy.cmipdf replaces it, and
        # get_qmin returns 0.5.
        (
            "setQmin",
            "get_qmin",
            0.5,
            "'diffpy.srfit.pdf.pdfcontribution.PDFContribution.setQmin' is "
            "deprecated and will be removed in version 4.0.0. Please use "
            "'diffpy.cmipdf.pdfcontribution.PDFContribution.set_qmin' "
            "instead.",
        ),
        # C3: setScatteringType sets the scattering type to neutron.
        # Expected: Warns that set_scattering_type in diffpy.cmipdf
        # replaces it, and get_scattering_type returns "N".
        (
            "setScatteringType",
            "get_scattering_type",
            "N",
            "'diffpy.srfit.pdf.pdfcontribution.PDFContribution."
            "setScatteringType' is deprecated and will be removed in version "
            "4.0.0. Please use 'diffpy.cmipdf.pdfcontribution."
            "PDFContribution.set_scattering_type' instead.",
        ),
    ],
)
def test_contribution_setters(
    old_setter, new_getter, input_value, expected_msg
):
    contribution = PDFContribution("pdf")
    with pytest.warns(DeprecationWarning) as record:
        getattr(contribution, old_setter)(input_value)
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg
    actual_value = getattr(contribution, new_getter)()
    expected_value = input_value
    assert actual_value == expected_value


@pytest.mark.parametrize(
    "old_getter, new_getter, expected_msg",
    [
        # Test the deprecated getters of PDFContribution.
        # C1: getQmax is called after qmax is set to 21.
        # Expected: Warns that get_qmax in diffpy.cmipdf replaces it, and
        # returns the same value as get_qmax.
        (
            "getQmax",
            "get_qmax",
            "'diffpy.srfit.pdf.pdfcontribution.PDFContribution.getQmax' is "
            "deprecated and will be removed in version 4.0.0. Please use "
            "'diffpy.cmipdf.pdfcontribution.PDFContribution.get_qmax' "
            "instead.",
        ),
        # C2: getQmin is called after qmin is set to 0.5.
        # Expected: Warns that get_qmin in diffpy.cmipdf replaces it, and
        # returns the same value as get_qmin.
        (
            "getQmin",
            "get_qmin",
            "'diffpy.srfit.pdf.pdfcontribution.PDFContribution.getQmin' is "
            "deprecated and will be removed in version 4.0.0. Please use "
            "'diffpy.cmipdf.pdfcontribution.PDFContribution.get_qmin' "
            "instead.",
        ),
        # C3: getScatteringType is called after the type is set to "N".
        # Expected: Warns that get_scattering_type in diffpy.cmipdf
        # replaces it, and returns the same value as get_scattering_type.
        (
            "getScatteringType",
            "get_scattering_type",
            "'diffpy.srfit.pdf.pdfcontribution.PDFContribution."
            "getScatteringType' is deprecated and will be removed in version "
            "4.0.0. Please use 'diffpy.cmipdf.pdfcontribution."
            "PDFContribution.get_scattering_type' instead.",
        ),
    ],
)
def test_contribution_getters(old_getter, new_getter, expected_msg):
    contribution = PDFContribution("pdf")
    contribution._meta.update({"qmax": 21.0, "qmin": 0.5, "stype": "N"})
    with pytest.warns(DeprecationWarning) as record:
        actual_value = getattr(contribution, old_getter)()
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg
    expected_value = getattr(contribution, new_getter)()
    assert actual_value == expected_value


def test_loadData(datafile):
    # Case: loadData loads a PDF data file into a PDFContribution.
    # Expected: Warns that load_data in diffpy.cmipdf replaces it, and
    # loads the same profile as load_data.
    expected_msg = (
        "'diffpy.srfit.pdf.pdfcontribution.PDFContribution.loadData' is "
        "deprecated and will be removed in version 4.0.0. Please use "
        "'diffpy.cmipdf.pdfcontribution.PDFContribution.load_data' instead."
    )
    old_contribution = PDFContribution("pdf")
    with pytest.warns(DeprecationWarning) as record:
        old_contribution.loadData(datafile("si-q27r60-xray.gr"))
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg
    new_contribution = PDFContribution("pdf")
    new_contribution.load_data(datafile("si-q27r60-xray.gr"))
    actual_metadata = old_contribution.profile.meta
    expected_metadata = new_contribution.profile.meta
    assert actual_metadata == expected_metadata
    actual_xobs = old_contribution.profile.xobs.tolist()
    expected_xobs = new_contribution.profile.xobs.tolist()
    assert actual_xobs == expected_xobs


def test_setCalculationRange(datafile):
    # Case: setCalculationRange limits a loaded profile to r = 0 to 10.
    # Expected: Warns that set_calculation_range in diffpy.cmipdf replaces
    # it, and sets the same calculation points as set_calculation_range.
    expected_msg = (
        "'diffpy.srfit.pdf.pdfcontribution.PDFContribution."
        "setCalculationRange' is deprecated and will be removed in version "
        "4.0.0. Please use 'diffpy.cmipdf.pdfcontribution.PDFContribution."
        "set_calculation_range' instead."
    )
    old_contribution = PDFContribution("pdf")
    old_contribution.load_data(datafile("si-q27r60-xray.gr"))
    with pytest.warns(DeprecationWarning) as record:
        old_contribution.setCalculationRange(0, 10)
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg
    new_contribution = PDFContribution("pdf")
    new_contribution.load_data(datafile("si-q27r60-xray.gr"))
    new_contribution.set_calculation_range(0, 10)
    actual_x = old_contribution.profile.x.tolist()
    expected_x = new_contribution.profile.x.tolist()
    assert actual_x == expected_x


def test_addStructure(diffpy_srreal_available, nacl_structure):
    # Case: addStructure adds the structure under the name "nacl".
    # Expected: Warns that add_structure in diffpy.cmipdf replaces it and
    # takes the structure before the name, and adds a generator named
    # "nacl" that adapts the structure.
    if not diffpy_srreal_available:
        pytest.skip("diffpy.srreal package not available")
    expected_msg = (
        "'diffpy.srfit.pdf.pdfcontribution.PDFContribution.addStructure' is "
        "deprecated and will be removed in version 4.0.0. Please use "
        "'diffpy.cmipdf.pdfcontribution.PDFContribution.add_structure' "
        "instead. Additionally, the signature has changed. Please pass the "
        "structure before the name, as add_structure(structure, name)."
    )
    contribution = PDFContribution("pdf")
    with pytest.warns(DeprecationWarning) as record:
        contribution.addStructure("nacl", nacl_structure)
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg
    actual_generator_names = list(contribution._generators)
    expected_generator_names = ["nacl"]
    assert actual_generator_names == expected_generator_names
    actual_structure = contribution.nacl.stru
    expected_structure = nacl_structure
    assert actual_structure == expected_structure


def test_addPhase(diffpy_srreal_available, nacl_structure):
    # Case: addPhase adds a parameter set under the name "nacl".
    # Expected: Warns that add_structure_from_parset in diffpy.cmipdf
    # replaces it and takes the parset before the name, and adds a
    # generator named "nacl" that manages the parameter set.
    if not diffpy_srreal_available:
        pytest.skip("diffpy.srreal package not available")
    expected_msg = (
        "'diffpy.srfit.pdf.pdfcontribution.PDFContribution.addPhase' is "
        "deprecated and will be removed in version 4.0.0. Please use "
        "'diffpy.cmipdf.pdfcontribution.PDFContribution."
        "add_structure_from_parset' instead. Additionally, the signature has "
        "changed. Please pass the parset before the name, as "
        "add_structure_from_parset(parset, name)."
    )
    contribution = PDFContribution("pdf")
    parameter_set = DiffpyStructureParSet("phase", nacl_structure)
    with pytest.warns(DeprecationWarning) as record:
        contribution.addPhase("nacl", parameter_set)
    actual_msg = str(record[0].message)
    assert actual_msg == expected_msg
    actual_generator_names = list(contribution._generators)
    expected_generator_names = ["nacl"]
    assert actual_generator_names == expected_generator_names
    actual_phase = contribution.nacl.phase
    expected_phase = parameter_set
    assert actual_phase == expected_phase
