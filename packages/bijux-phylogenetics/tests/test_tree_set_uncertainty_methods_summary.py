from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.trees import (
    build_tree_set_uncertainty_method_report,
    build_tree_set_uncertainty_methods_summary_text,
    write_tree_set_uncertainty_methods_summary_text,
)

FIXTURES = Path(__file__).parent / "fixtures" / "trees"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_build_tree_set_uncertainty_methods_summary_text_reports_consensus_support_and_instability() -> (
    None
):
    report = build_tree_set_uncertainty_method_report(
        fixture("example_tree_set_left.nwk")
    )

    text = build_tree_set_uncertainty_methods_summary_text(report)

    assert "Tree-Set Uncertainty Methods Summary" in text
    assert "- consensus contract: `majority-rule`" in text
    assert "- rooted topology cluster count: `2`" in text
    assert "- unstable taxon count: `4`" in text
    assert "- conflicting clade count: `0`" in text


def test_write_tree_set_uncertainty_methods_summary_text_writes_markdown(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "tree-set-uncertainty-methods-summary.md"
    report = build_tree_set_uncertainty_method_report(
        fixture("example_tree_set_left.nwk")
    )

    result = write_tree_set_uncertainty_methods_summary_text(output_path, report)

    assert result.output_path == output_path
    assert result.topology_cluster_count == 2
    assert result.unstable_taxon_count == 4
    assert result.warning_count >= 1
    assert output_path.read_text(encoding="utf-8") == result.text
