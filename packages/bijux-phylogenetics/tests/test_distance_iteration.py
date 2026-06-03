from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.distance import (
    assess_distance_method_assumptions,
    assess_distance_method_maturity,
    assess_imported_distance_method_assumptions,
    bootstrap_distance_trees,
    build_distance_method_report,
    build_distance_tree,
    compare_distance_gap_policies,
    compare_distance_models,
    compare_distance_tree_to_reference_tree,
    compute_pairwise_genetic_distance_matrix,
    inspect_distance_matrix_quality,
    inspect_imported_distance_matrix_quality,
    list_distance_tree_method_policies,
    resolve_distance_tree_method_policy,
    summarize_distance_bootstrap_support,
    validate_distance_reference_examples,
    write_distance_reproducibility_bundle,
)
from bijux_phylogenetics.io.newick import write_newick
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
    assert len(report.observations) == 11
    assert len(report.tree_observations) == 2


def test_distance_tree_method_policies_support_bionj() -> None:
    policies = {row.method: row for row in list_distance_tree_method_policies()}
    assert set(policies) == {
        "bionj",
        "complete-linkage",
        "neighbor-joining",
        "single-linkage",
        "upgma",
        "wpgma",
    }
    assert policies["neighbor-joining"].supported is True
    assert policies["single-linkage"].supported is True
    assert policies["complete-linkage"].supported is True
    assert policies["neighbor-joining"].reference_surface == "ape::nj"
    assert policies["bionj"].supported is True
    assert policies["bionj"].reference_surface == "ape::bionj"
    assert policies["bionj"].support_scope == "owned-runtime"
    assert policies["upgma"].supported is True
    assert policies["wpgma"].supported is True


def test_resolve_distance_tree_method_policy_supports_bionj() -> None:
    policy = resolve_distance_tree_method_policy("bionj")
    assert policy.method == "bionj"
    tree, report = build_distance_tree(
        fixture("example_alignment_distance.fasta"),
        method="bionj",
    )
    assert tree.rooted is False
    assert report.method == "bionj"


def test_compute_pairwise_genetic_distance_matrix_supports_kimura_two_parameter() -> (
    None
):
    report = compute_pairwise_genetic_distance_matrix(
        fixture("example_alignment_distance.fasta"),
        model="kimura-2-parameter",
    )
    pair = next(
        row
        for row in report.pairs
        if row.left_identifier == "A" and row.right_identifier == "B"
    )
    assert pair.distance == 0.14384103622589
    assert pair.transition_sites == 1.0
    assert pair.transversion_sites == 0.0


def test_compute_pairwise_genetic_distance_matrix_supports_amino_acid_model() -> None:
    report = compute_pairwise_genetic_distance_matrix(
        fixture("example_alignment_protein.fasta"),
        model="amino-acid-p-distance",
    )
    pair = next(
        row
        for row in report.pairs
        if row.left_identifier == "P1" and row.right_identifier == "P2"
    )
    assert report.inferred_alphabet == "protein"
    assert pair.distance == 0.125


def test_compute_pairwise_genetic_distance_matrix_honors_partial_match_ambiguity_policy() -> (
    None
):
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
    partial_pair = next(
        row
        for row in partial.pairs
        if row.left_identifier == "A" and row.right_identifier == "B"
    )
    ignore_pair = next(
        row
        for row in ignore.pairs
        if row.left_identifier == "A" and row.right_identifier == "B"
    )
    assert partial_pair.distance == 0.15
    assert partial_pair.ambiguity_sites == 1
    assert ignore_pair.distance == 0.0
    assert ignore_pair.comparable_sites == 4


def test_compute_pairwise_genetic_distance_matrix_supports_report_only_ambiguity_policy() -> (
    None
):
    report_only = compute_pairwise_genetic_distance_matrix(
        fixture("example_alignment_ambiguity.fasta"),
        model="p-distance",
        ambiguity_policy="report-only",
    )
    pair = next(
        row
        for row in report_only.pairs
        if row.left_identifier == "A" and row.right_identifier == "B"
    )
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
    report = assess_distance_method_assumptions(
        fixture("example_alignment_distance.fasta")
    )
    assert report.ultrametric_compatible is False
    assert len(report.upgma_ultrametric_violations) > 0


def test_assess_imported_distance_method_assumptions_accepts_ultrametric_matrix() -> (
    None
):
    report = assess_imported_distance_method_assumptions(
        fixture("example_distance_matrix_ultrametric.tsv")
    )
    assert report.ultrametric_compatible is True
    assert report.upgma_ultrametric_violations == []


def test_inspect_imported_distance_matrix_quality_reports_structural_and_site_risks() -> (
    None
):
    report = inspect_imported_distance_matrix_quality(
        fixture("example_distance_matrix_nonultrametric.tsv")
    )
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
    assert (
        report.consensus_newick == "((A:0.0125,B:0.0125)100:0.6875,C:0.0125,D:0.0125);"
    )
    assert [row.sampled_site_indices for row in report.replicate_rows] == [
        [5, 2, 6, 0, 1, 1, 5, 0],
        [3, 0, 1, 6, 6, 1, 3, 1],
        [6, 0, 1, 3, 0, 6, 0, 3],
        [0, 2, 4, 6, 2, 1, 4, 2],
        [1, 3, 5, 1, 1, 0, 3, 7],
    ]
    assert [row.tree_newick for row in report.replicate_rows] == [
        "((A:0,B:0)Inner1:0.625,C:0,D:0)Inner2;",
        "((A:0,B:0)Inner1:0.75,C:0,D:0)Inner2;",
        "((A:0,B:0)Inner1:0.75,C:0,D:0)Inner2;",
        "((A:0,B:0)Inner1:0.625,C:0,D:0)Inner2;",
        "((A:0.0625,B:0.0625)Inner1:0.6875,C:0.0625,D:0.0625)Inner2;",
    ]
    assert [(row.clade, row.tree_count, row.frequency) for row in report.support] == [
        ("A|B", 5, 1.0)
    ]
    assert any(
        len(set(row.sampled_site_indices)) < len(row.sampled_site_indices)
        for row in report.replicate_rows
    )


def test_summarize_distance_bootstrap_support_reports_weak_clade_counts() -> None:
    report = summarize_distance_bootstrap_support(
        bootstrap_distance_trees(
            fixture("example_alignment_distance.fasta"),
            method="neighbor-joining",
            replicates=5,
            seed=7,
        )[1]
    )
    assert report.clade_count > 0
    assert report.replicates == 5


def test_compare_distance_tree_to_reference_tree_reports_matching_topology(
    tmp_path: Path,
) -> None:
    tree, _ = build_distance_tree(
        fixture("example_alignment_distance.fasta"), method="neighbor-joining"
    )
    reference_path = tmp_path / "reference.nwk"
    write_newick(reference_path, tree)
    report = compare_distance_tree_to_reference_tree(
        fixture("example_alignment_distance.fasta"),
        reference_path,
        method="neighbor-joining",
    )
    assert report.topology.topology_equal is True


def test_compare_distance_models_reports_supported_dna_models() -> None:
    report = compare_distance_models(fixture("example_alignment_distance.fasta"))
    assert [row.model for row in report.rows] == [
        "felsenstein-81",
        "jukes-cantor",
        "kimura-2-parameter",
        "p-distance",
        "tamura-nei-93",
    ]


def test_compare_distance_gap_policies_reports_changed_pairs() -> None:
    report = compare_distance_gap_policies(
        fixture("example_alignment_distance_gaps.fasta"),
        model="p-distance",
    )
    assert report.changed_pair_count > 0
    assert any(row.comparable_site_delta != 0 for row in report.rows)


def test_build_distance_method_report_combines_support_models_and_gap_sensitivity() -> (
    None
):
    report = build_distance_method_report(
        fixture("example_alignment_distance.fasta"),
        method="neighbor-joining",
        bootstrap_replicates=5,
        bootstrap_seed=3,
    )
    assert report.method_policy.method == "neighbor-joining"
    assert report.method_policy.reference_surface == "ape::nj"
    assert report.bootstrap_summary.clade_count > 0
    assert report.model_comparison.rows
    assert report.gap_policy_sensitivity.pair_count > 0


def test_assess_distance_method_maturity_validates_bundle_provenance() -> None:
    report = assess_distance_method_maturity(
        fixture("example_alignment_distance.fasta"),
        method="neighbor-joining",
        bootstrap_replicates=5,
        bootstrap_seed=3,
        validate_bundle=True,
    )
    assert report.method_policy.supported is True
    assert report.decision in {"production_candidate", "validated_with_limits"}
    assert any(
        check.name == "bundle_provenance" and check.satisfied for check in report.checks
    )


def test_write_distance_reproducibility_bundle_writes_expected_files(
    tmp_path: Path,
) -> None:
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
        "distance-summary.json",
        "distance-tree.nwk",
        "input-alignment.fasta",
    ]
    manifest = json.loads(
        (tmp_path / "bundle" / "distance-analysis.manifest.json").read_text(
            encoding="utf-8"
        )
    )
    assert manifest["reference_validation_passed"] is True
    assert "input-alignment.fasta" in manifest["output_checksums"]


def test_render_distance_report_for_alignment_embeds_quality_sections(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "distance-report.html"
    render_distance_report(
        out_path=output_path,
        alignment_path=fixture("example_alignment_distance.fasta"),
    )
    html = output_path.read_text(encoding="utf-8")
    assert "distance-quality" in html
    assert "distance-method-assumptions" in html
    assert "distance-reference-validation" in html
    assert "distance-bootstrap-summary" in html
    assert "distance-model-comparison" in html
    assert "distance-gap-policy-sensitivity" in html
    assert "distance-maturity-gate" in html
