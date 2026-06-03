from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main
from bijux_phylogenetics.distance import (
    compute_pairwise_genetic_distance_matrix,
    validate_distance_reference_examples,
)

FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "alignments", "metadata", "expected")


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def test_f81_alias_matches_canonical_surface_and_reports_base_frequencies() -> None:
    alias_report = compute_pairwise_genetic_distance_matrix(
        fixture("example_alignment_distance_gaps.fasta"),
        model="f81",
        gap_handling="pairwise-deletion",
    )
    canonical_report = compute_pairwise_genetic_distance_matrix(
        fixture("example_alignment_distance_gaps.fasta"),
        model="felsenstein-81",
        gap_handling="pairwise-deletion",
    )

    assert alias_report.model == "felsenstein-81"
    assert alias_report.pairs == canonical_report.pairs
    assert alias_report.model_parameters == canonical_report.model_parameters
    assert alias_report.model_parameters is not None
    assert alias_report.model_parameters.informative_base_count == 22
    assert alias_report.model_parameters.base_frequency_a == 0.363636363636364
    assert alias_report.model_parameters.base_frequency_c == 0.136363636363636
    assert alias_report.model_parameters.base_frequency_g == 0.227272727272727
    assert alias_report.model_parameters.base_frequency_t == 0.272727272727273
    assert alias_report.model_parameters.f81_limit == 0.723140495867769


def test_tn93_alias_matches_canonical_surface_and_reports_coefficients() -> None:
    alias_report = compute_pairwise_genetic_distance_matrix(
        fixture("example_alignment_duplicates.fasta"),
        model="tn93",
    )
    canonical_report = compute_pairwise_genetic_distance_matrix(
        fixture("example_alignment_duplicates.fasta"),
        model="tamura-nei-93",
    )

    assert alias_report.model == "tamura-nei-93"
    assert alias_report.pairs == canonical_report.pairs
    assert alias_report.model_parameters == canonical_report.model_parameters
    assert alias_report.model_parameters is not None
    assert alias_report.model_parameters.tn93_ag_coefficient == 0.225
    assert alias_report.model_parameters.tn93_ct_coefficient == 0.264705882352941
    assert (
        alias_report.model_parameters.tn93_transversion_coefficient == 0.254434742647059
    )


def test_tn93_invalid_composition_is_explicit_on_clean_alignment() -> None:
    report = compute_pairwise_genetic_distance_matrix(
        fixture("example_alignment_distance.fasta"),
        model="tn93",
    )
    pair = next(
        row
        for row in report.pairs
        if row.left_identifier == "A" and row.right_identifier == "B"
    )

    assert report.model_parameters is not None
    assert report.model_parameters.base_frequency_g == 0.0
    assert report.warnings == [
        "alignment-wide resolved base composition omits at least one nucleotide, so TN93 assumptions break"
    ]
    assert pair.distance is None
    assert pair.saturated is True
    assert (
        pair.saturation_reason
        == "alignment-wide resolved base composition omits at least one nucleotide class, so the TN93 correction is undefined"
    )


def test_f81_and_tn93_reference_examples_are_included() -> None:
    report = validate_distance_reference_examples()

    lookup = {observation.case: observation for observation in report.observations}
    assert lookup["dna-felsenstein-81"].passed is True
    assert lookup["dna-tamura-nei-93"].passed is True


def test_cli_alignment_distance_matrix_writes_f81_parameter_table(
    tmp_path: Path, capsys
) -> None:
    matrix_path = tmp_path / "distance.tsv"
    parameter_path = tmp_path / "parameters.tsv"

    exit_code = main(
        [
            "alignment",
            "distance-matrix",
            str(fixture("example_alignment_distance_gaps.fasta")),
            "--model",
            "f81",
            "--out",
            str(matrix_path),
            "--parameters-out",
            str(parameter_path),
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["model"] == "felsenstein-81"
    assert payload["warnings"] == []
    assert parameter_path.read_text(encoding="utf-8").splitlines()[:6] == [
        "parameter\tvalue",
        "informative_base_count\t22",
        "base_frequency_a\t0.363636363636364",
        "base_frequency_c\t0.136363636363636",
        "base_frequency_g\t0.227272727272727",
        "base_frequency_t\t0.272727272727273",
    ]
    assert matrix_path.exists()


def test_cli_alignment_distance_matrix_warns_when_tn93_assumptions_break(
    tmp_path: Path, capsys
) -> None:
    parameter_path = tmp_path / "parameters.tsv"

    exit_code = main(
        [
            "alignment",
            "distance-matrix",
            str(fixture("example_alignment_distance.fasta")),
            "--model",
            "tn93",
            "--parameters-out",
            str(parameter_path),
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["model"] == "tamura-nei-93"
    assert payload["warnings"] == [
        "alignment-wide resolved base composition omits at least one nucleotide, so TN93 assumptions break",
        "one or more pairwise distances are saturated or non-finite under the selected model",
    ]
    assert "base_frequency_g\t0" in parameter_path.read_text(encoding="utf-8")
