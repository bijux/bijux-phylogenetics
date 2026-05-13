from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.clades import extract_tree_clades
from bijux_phylogenetics.reports.tree_package import (
    build_tree_report_package,
    summarize_tree_support,
)


def tree_fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_build_tree_report_package_writes_html_svg_and_tsv_outputs(
    tmp_path: Path,
) -> None:
    result = build_tree_report_package(
        tree_fixture("example_tree_support_left.nwk"),
        out_dir=tmp_path / "tree-report-package",
    )

    assert result.report_path.exists()
    assert result.figure_path.exists()
    assert result.support_table_path.exists()
    assert result.clade_table_path.exists()
    assert result.branch_stats_path.exists()
    assert result.manifest_path.exists()
    assert result.figure.rendered_support_count == 2

    html = result.report_path.read_text(encoding="utf-8")
    assert "Bijux Full Tree Report" in html
    assert "<svg" in html
    assert "Support Table" in html
    assert "Clade Table" in html
    assert "Branch-Length Stats" in html

    support_lines = result.support_table_path.read_text(encoding="utf-8").splitlines()
    assert support_lines[0].startswith("node_kind\tnode\tnode_label\tdescendant_taxa")
    assert any("\tstrong\t" in line for line in support_lines[1:])

    branch_lines = result.branch_stats_path.read_text(encoding="utf-8").splitlines()
    assert branch_lines[0].startswith("branch_count\tdefined_branch_count")
    assert branch_lines[1].startswith("6\t6\t0")

    manifest = json.loads(result.manifest_path.read_text(encoding="utf-8"))
    assert manifest["report_kind"] == "tree_package"
    assert manifest["metrics"]["supported_branch_count"] == 3


def test_build_tree_report_package_reports_withheld_support_rendering(
    tmp_path: Path,
) -> None:
    result = build_tree_report_package(
        tree_fixture("example_tree_support_invalid.nwk"),
        out_dir=tmp_path / "tree-report-package",
    )

    assert result.support_audit.validated is False
    assert result.figure.rendered_support_count == 0
    assert any("withheld" in item for item in result.limitations)

    support_lines = result.support_table_path.read_text(encoding="utf-8").splitlines()
    assert len(support_lines) > 1


def test_summarize_tree_support_keeps_root_and_internal_support_rows() -> None:
    clades = extract_tree_clades(tree_fixture("example_tree_support_left.nwk"))

    rows = summarize_tree_support(clades)

    rows_by_kind = {row.node_kind: row for row in rows if row.node_kind == "root"}
    assert len(rows) == 3
    assert sorted(row.node_kind for row in rows) == ["internal", "internal", "root"]
    assert rows_by_kind["root"].support == 99.0
    internal_supports = sorted(
        row.support_fraction for row in rows if row.node_kind == "internal"
    )
    assert internal_supports == [0.88, 0.95]
