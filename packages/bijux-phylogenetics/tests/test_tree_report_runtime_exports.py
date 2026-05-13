from __future__ import annotations

import bijux_phylogenetics

from bijux_phylogenetics.reports import (
    TreeBranchStatisticsRow,
    TreeReportPackageResult,
    TreeSupportRow,
    build_tree_report_package,
    summarize_tree_branch_statistics,
    summarize_tree_support,
    write_tree_branch_statistics_table,
    write_tree_support_table,
)


def test_tree_report_package_surfaces_export_publicly() -> None:
    assert bijux_phylogenetics.TreeReportPackageResult is TreeReportPackageResult
    assert bijux_phylogenetics.TreeSupportRow is TreeSupportRow
    assert bijux_phylogenetics.TreeBranchStatisticsRow is TreeBranchStatisticsRow
    assert bijux_phylogenetics.build_tree_report_package is build_tree_report_package
    assert bijux_phylogenetics.summarize_tree_support is summarize_tree_support
    assert (
        bijux_phylogenetics.summarize_tree_branch_statistics
        is summarize_tree_branch_statistics
    )
    assert bijux_phylogenetics.write_tree_support_table is write_tree_support_table
    assert (
        bijux_phylogenetics.write_tree_branch_statistics_table
        is write_tree_branch_statistics_table
    )
