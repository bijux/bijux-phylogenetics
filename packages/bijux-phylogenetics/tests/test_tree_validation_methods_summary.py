from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.reports import write_tree_validation_methods_summary_text


def tree_fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_write_tree_validation_methods_summary_text_reports_thresholds_and_contexts(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "tree-validation-methods-summary.md"

    result = write_tree_validation_methods_summary_text(
        output_path,
        tree_path=tree_fixture("example_tree_unrooted.nwk"),
    )

    assert result.output_path == output_path
    assert result.blocked_context_count == 2
    assert result.repair_item_count == 0
    assert "Tree Validation Methods Summary" in result.text
    assert (
        "ultrametricity is evaluated with the APE-compatible tolerance" in result.text
    )
    assert "long terminal branches are flagged when they exceed" in result.text
    assert (
        "- `time tree`: blocked by time trees require ultrametric root-to-tip distances"
        in result.text
    )
    assert (
        "- `comparative methods`: blocked by comparative methods require a rooted tree"
        in result.text
    )
    assert "- no taxa were excluded or flagged for repair" in result.text
    assert output_path.read_text(encoding="utf-8") == result.text


def test_write_tree_validation_methods_summary_text_reports_repair_items(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "tree-validation-methods-summary-flags.md"

    result = write_tree_validation_methods_summary_text(
        output_path,
        tree_path=tree_fixture("example_tree_labels.nwk"),
    )

    assert result.repair_item_count >= 1
    assert "downstream-unsafe tip labels detected" in result.text
    assert "`Homo sapiens`" in result.text
    assert "`Mus musculus`" in result.text
    assert (
        "- `time tree`: blocked by time trees require ultrametric root-to-tip distances"
        in result.text
    )
