from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.compare.presentation import build_tree_comparison_report


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_build_tree_comparison_report_writes_support_conflict_section(
    tmp_path: Path,
) -> None:
    output = tmp_path / "compare.html"
    result = build_tree_comparison_report(
        fixture("example_tree_support_left.nwk"),
        fixture("example_tree_support_right.nwk"),
        out_path=output,
    )

    html = output.read_text(encoding="utf-8")

    assert result.output_path == output
    assert "Bijux Tree Comparison Report" in html
    assert "support-comparison" in html
    assert "support-conflicts" in html
