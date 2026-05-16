from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main
from bijux_phylogenetics.distance import compute_pairwise_genetic_distance_matrix

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


def test_k80_alias_matches_existing_kimura_two_parameter_surface() -> None:
    k80_report = compute_pairwise_genetic_distance_matrix(
        fixture("example_alignment_distance.fasta"),
        model="k80",
    )
    kimura_report = compute_pairwise_genetic_distance_matrix(
        fixture("example_alignment_distance.fasta"),
        model="kimura-2-parameter",
    )

    assert k80_report.model == "kimura-2-parameter"
    assert k80_report.identifiers == kimura_report.identifiers
    assert k80_report.pairs == kimura_report.pairs


def test_k80_alias_supports_complete_deletion_on_gapped_alignment() -> None:
    report = compute_pairwise_genetic_distance_matrix(
        fixture("example_alignment_distance_gaps.fasta"),
        model="k80",
        gap_handling="complete-deletion",
    )

    assert [
        (
            pair.left_identifier,
            pair.right_identifier,
            pair.distance,
            pair.transition_sites,
            pair.transversion_sites,
            pair.comparable_sites,
        )
        for pair in report.pairs
    ] == [
        ("A", "A", 0.0, 0.0, 0.0, 4),
        ("A", "B", 0.0, 0.0, 0.0, 4),
        ("A", "C", 0.346573590279973, 1.0, 0.0, 4),
        ("A", "D", 0.0, 0.0, 0.0, 4),
        ("B", "B", 0.0, 0.0, 0.0, 4),
        ("B", "C", 0.346573590279973, 1.0, 0.0, 4),
        ("B", "D", 0.0, 0.0, 0.0, 4),
        ("C", "C", 0.0, 0.0, 0.0, 4),
        ("C", "D", 0.346573590279973, 1.0, 0.0, 4),
        ("D", "D", 0.0, 0.0, 0.0, 4),
    ]


def test_k80_alias_tracks_ignored_ambiguity_sites_explicitly() -> None:
    report = compute_pairwise_genetic_distance_matrix(
        fixture("example_alignment_ambiguity.fasta"),
        model="k80",
        gap_handling="pairwise-deletion",
        ambiguity_policy="ignore",
    )

    pair = next(
        row
        for row in report.pairs
        if row.left_identifier == "A" and row.right_identifier == "B"
    )

    assert pair.distance == 0.0
    assert pair.transition_sites == 0.0
    assert pair.transversion_sites == 0.0
    assert pair.comparable_sites == 4
    assert pair.ambiguity_sites == 1
    assert pair.skipped_sites == 2


def test_k80_alias_distinguishes_undefined_and_infinite_saturation() -> None:
    clean_report = compute_pairwise_genetic_distance_matrix(
        fixture("example_alignment_distance.fasta"),
        model="k80",
    )
    infinite_pair = next(
        row
        for row in clean_report.pairs
        if row.left_identifier == "A" and row.right_identifier == "C"
    )
    saturated_report = compute_pairwise_genetic_distance_matrix(
        fixture("example_alignment_distance_saturated.fasta"),
        model="k80",
    )
    undefined_pair = next(
        row
        for row in saturated_report.pairs
        if row.left_identifier == "A" and row.right_identifier == "B"
    )

    assert infinite_pair.distance is None
    assert infinite_pair.saturated is True
    assert (
        infinite_pair.saturation_reason
        == "transition and transversion proportions are at the Kimura 2-parameter correction limit, so the corrected distance tends to infinity"
    )
    assert undefined_pair.distance is None
    assert undefined_pair.saturated is True
    assert (
        undefined_pair.saturation_reason
        == "transition and transversion proportions exceed the Kimura 2-parameter correction range, so the corrected distance is undefined"
    )


def test_cli_alignment_distance_matrix_accepts_k80_model_alias_and_components_out(
    tmp_path: Path, capsys
) -> None:
    matrix_path = tmp_path / "distance.tsv"
    components_path = tmp_path / "components.tsv"

    exit_code = main(
        [
            "alignment",
            "distance-matrix",
            str(fixture("example_alignment_distance.fasta")),
            "--model",
            "k80",
            "--out",
            str(matrix_path),
            "--components-out",
            str(components_path),
            "--json",
        ]
    )

    captured = capsys.readouterr()
    payload = json.loads(captured.out)

    assert exit_code == 0
    assert payload["metrics"]["model"] == "kimura-2-parameter"
    assert matrix_path.read_text(encoding="utf-8").splitlines()[:4] == [
        "left_identifier\tright_identifier\tdistance\tcomparable_sites",
        "A\tA\t0\t8",
        "A\tB\t0.14384103622589\t8",
        "A\tC\t\t8",
    ]
    assert components_path.read_text(encoding="utf-8").splitlines()[:4] == [
        "left_identifier\tright_identifier\tdistance\tcomparable_sites\tmismatch_sites\ttransition_sites\tag_transition_sites\tct_transition_sites\ttransversion_sites\tambiguity_sites\tskipped_sites\tsaturated\tsaturation_reason",
        "A\tA\t0\t8\t0\t0\t0\t0\t0\t0\t0\tfalse\t",
        "A\tB\t0.14384103622589\t8\t1\t1\t0\t1\t0\t0\t0\tfalse\t",
        "A\tC\t\t8\t4\t0\t0\t0\t4\t0\t0\ttrue\ttransition and transversion proportions are at the Kimura 2-parameter correction limit, so the corrected distance tends to infinity",
    ]
