from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.distance import (
    bootstrap_distance_trees,
    compute_pairwise_genetic_distance_matrix,
    inspect_distance_matrix_quality,
    validate_distance_reference_examples,
    write_distance_reproducibility_bundle,
)
from bijux_phylogenetics.reports.service import render_distance_report


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


def test_validate_distance_reference_examples_passes() -> None:
    report = validate_distance_reference_examples()
    assert report.all_passed is True
    assert len(report.observations) == 3


def test_compute_pairwise_genetic_distance_matrix_supports_kimura_two_parameter() -> None:
    report = compute_pairwise_genetic_distance_matrix(
        fixture("example_alignment_distance.fasta"),
        model="kimura-2-parameter",
    )
    pair = next(row for row in report.pairs if row.left_identifier == "A" and row.right_identifier == "B")
    assert pair.distance == 0.14384103622589
    assert pair.transition_sites == 1.0
    assert pair.transversion_sites == 0.0


def test_compute_pairwise_genetic_distance_matrix_supports_amino_acid_model() -> None:
    report = compute_pairwise_genetic_distance_matrix(
        fixture("example_alignment_protein.fasta"),
        model="amino-acid-p-distance",
    )
    pair = next(row for row in report.pairs if row.left_identifier == "P1" and row.right_identifier == "P2")
    assert report.inferred_alphabet == "protein"
    assert pair.distance == 0.125


def test_compute_pairwise_genetic_distance_matrix_honors_partial_match_ambiguity_policy() -> None:
    partial = compute_pairwise_genetic_distance_matrix(
        fixture("example_alignment_ambiguity.fasta"),
        model="p-distance",
        ambiguity_policy="partial-match",
    )
    ignore = compute_pairwise_genetic_distance_matrix(
        fixture("example_alignment_ambiguity.fasta"),
        model="p-distance",
        ambiguity_policy="ignore",
    )
    partial_pair = next(row for row in partial.pairs if row.left_identifier == "A" and row.right_identifier == "B")
    ignore_pair = next(row for row in ignore.pairs if row.left_identifier == "A" and row.right_identifier == "B")
    assert partial_pair.distance == 0.15
    assert partial_pair.ambiguity_sites == 1
    assert ignore_pair.distance == 0.0
    assert ignore_pair.comparable_sites == 4


def test_inspect_distance_matrix_quality_reports_saturation_risk() -> None:
    report = inspect_distance_matrix_quality(
        fixture("example_alignment_distance_saturated.fasta"),
        model="jukes-cantor",
    )
    assert report.method_assessment.decision == "risky"
    assert len(report.saturated_pairs) > 0


def test_bootstrap_distance_trees_returns_consensus_and_support() -> None:
    trees, report = bootstrap_distance_trees(
        fixture("example_alignment_distance.fasta"),
        method="neighbor-joining",
        replicates=5,
        seed=7,
    )
    assert len(trees) == 5
    assert report.tree_count == 5
    assert report.consensus_newick.endswith(";")
    assert len(report.support) > 0


def test_write_distance_reproducibility_bundle_writes_expected_files(tmp_path: Path) -> None:
    report = write_distance_reproducibility_bundle(
        tmp_path / "bundle",
        alignment_path=fixture("example_alignment_distance.fasta"),
        method="neighbor-joining",
        replicates=5,
        seed=11,
    )
    names = sorted(path.name for path in report.files)
    assert names == [
        "bootstrap-clades.tsv",
        "bootstrap-consensus.nwk",
        "bootstrap-replicates.trees",
        "bootstrap-support.tsv",
        "distance-analysis.manifest.json",
        "distance-matrix.tsv",
        "distance-tree.nwk",
    ]


def test_render_distance_report_for_alignment_embeds_quality_sections(tmp_path: Path) -> None:
    output_path = tmp_path / "distance-report.html"
    render_distance_report(
        out_path=output_path,
        alignment_path=fixture("example_alignment_distance.fasta"),
    )
    html = output_path.read_text(encoding="utf-8")
    assert "distance-quality" in html
    assert "distance-reference-validation" in html
