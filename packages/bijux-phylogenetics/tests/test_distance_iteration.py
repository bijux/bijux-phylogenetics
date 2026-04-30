from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.distance import (
    assess_distance_method_assumptions,
    assess_imported_distance_method_assumptions,
    bootstrap_distance_trees,
    compute_pairwise_genetic_distance_matrix,
    inspect_imported_distance_matrix_quality,
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
    assert len(report.observations) == 9
    assert len(report.tree_observations) == 2


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


def test_compute_pairwise_genetic_distance_matrix_supports_report_only_ambiguity_policy() -> None:
    report_only = compute_pairwise_genetic_distance_matrix(
        fixture("example_alignment_ambiguity.fasta"),
        model="p-distance",
        ambiguity_policy="report-only",
    )
    pair = next(row for row in report_only.pairs if row.left_identifier == "A" and row.right_identifier == "B")
    assert pair.distance == 0.0
    assert pair.comparable_sites == 4
    assert pair.ambiguity_sites == 1


def test_inspect_distance_matrix_quality_reports_saturation_risk() -> None:
    report = inspect_distance_matrix_quality(
        fixture("example_alignment_distance_saturated.fasta"),
        model="jukes-cantor",
    )
    assert report.method_assessment.decision == "risky"
    assert len(report.saturated_pairs) > 0


def test_assess_distance_method_assumptions_reports_ultrametric_violation() -> None:
    report = assess_distance_method_assumptions(fixture("example_alignment_distance.fasta"))
    assert report.ultrametric_compatible is False
    assert len(report.upgma_ultrametric_violations) > 0


def test_assess_imported_distance_method_assumptions_accepts_ultrametric_matrix() -> None:
    report = assess_imported_distance_method_assumptions(fixture("example_distance_matrix_ultrametric.tsv"))
    assert report.ultrametric_compatible is True
    assert report.upgma_ultrametric_violations == []


def test_inspect_imported_distance_matrix_quality_reports_structural_and_site_risks() -> None:
    report = inspect_imported_distance_matrix_quality(fixture("example_distance_matrix_nonultrametric.tsv"))
    assert report.validation.complete is True
    assert report.saturation_audit_scale == "unknown"
    assert report.low_information_pair_cutoff == 50
    assert report.low_information_pairs == []
    assert "saturation heuristics were skipped" in report.warnings[-1]


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
        "input-alignment.fasta",
    ]
    manifest = json.loads((tmp_path / "bundle" / "distance-analysis.manifest.json").read_text(encoding="utf-8"))
    assert manifest["reference_validation_passed"] is True
    assert "input-alignment.fasta" in manifest["output_checksums"]


def test_render_distance_report_for_alignment_embeds_quality_sections(tmp_path: Path) -> None:
    output_path = tmp_path / "distance-report.html"
    render_distance_report(
        out_path=output_path,
        alignment_path=fixture("example_alignment_distance.fasta"),
    )
    html = output_path.read_text(encoding="utf-8")
    assert "distance-quality" in html
    assert "distance-method-assumptions" in html
    assert "distance-reference-validation" in html
