from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.compare.influence import analyze_taxon_influence
from bijux_phylogenetics.compare.topology import (
    compare_branch_lengths,
    compare_support_values,
    compare_tree_paths,
)
from bijux_phylogenetics.engines.inference import (
    build_inference_comparison_conclusion_rows,
    build_inference_comparison_weighted_conflict_rows,
    summarize_inference_comparison_conclusions,
)


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_weighted_conflict_rows_rank_supported_conflicts_ahead_of_weak_disagreements() -> (
    None
):
    support_report = compare_support_values(
        fixture("example_tree_support_conflict_left.nwk"),
        fixture("example_tree_support_conflict_right.nwk"),
    )

    rows = build_inference_comparison_weighted_conflict_rows(support_report)

    assert [
        (row.split_id, row.severity_class, row.support_weight, row.serious_conflict)
        for row in rows
    ] == [
        ("A|B|C", "high_support_conflict", 0.92, True),
        ("C|D", "low_support_disagreement", 0.4, False),
    ]


def test_conclusion_rows_separate_stable_and_unstable_shared_clades(
    tmp_path: Path,
) -> None:
    left_path = tmp_path / "left.nwk"
    right_path = tmp_path / "right.nwk"
    left_path.write_text(
        "((A:0.1,B:0.1)95:0.2,(C:0.1,D:0.1)88:0.2);\n", encoding="utf-8"
    )
    right_path.write_text(
        "((A:0.1,B:0.1)90:0.2,(C:0.1,D:0.1)60:0.2);\n", encoding="utf-8"
    )

    support_report = compare_support_values(left_path, right_path)

    rows = build_inference_comparison_conclusion_rows(support_report)

    assert [
        (row.split_id, row.conclusion_class, row.evidence_class, row.serious_conflict)
        for row in rows
    ] == [
        ("A|B", "stable_clade", "strong_stable", False),
        ("C|D", "unstable_clade", "moderate_support_disagreement", False),
    ]


def test_inference_comparison_summary_reports_driver_taxa_where_detectable() -> None:
    left_path = fixture("example_tree_taxon_influence_left.nwk")
    right_path = fixture("example_tree_taxon_influence_right.nwk")
    topology_report = compare_tree_paths(left_path, right_path)
    branch_length_report = compare_branch_lengths(left_path, right_path)
    support_report = compare_support_values(left_path, right_path)
    conclusion_rows = build_inference_comparison_conclusion_rows(support_report)
    weighted_conflict_rows = build_inference_comparison_weighted_conflict_rows(
        support_report
    )
    taxon_influence_report = analyze_taxon_influence(left_path, right_path)

    summary = summarize_inference_comparison_conclusions(
        topology_report,
        branch_length_report,
        weighted_conflict_rows=weighted_conflict_rows,
        conclusion_rows=conclusion_rows,
        taxon_influence_report=taxon_influence_report,
    )

    assert summary.engine_specific_clade_count == 2
    assert summary.high_support_conflict_count == 1
    assert summary.top_conflict_driver_taxa == ["C"]
