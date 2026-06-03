from __future__ import annotations

import bijux_phylogenetics.reports as reports_api
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
    assert reports_api.TreeReportPackageResult is TreeReportPackageResult
    assert reports_api.TreeSupportRow is TreeSupportRow
    assert reports_api.TreeBranchStatisticsRow is TreeBranchStatisticsRow
    assert reports_api.build_tree_report_package is build_tree_report_package
    assert reports_api.summarize_tree_support is summarize_tree_support
    assert (
        reports_api.summarize_tree_branch_statistics is summarize_tree_branch_statistics
    )
    assert reports_api.write_tree_support_table is write_tree_support_table
    assert (
        reports_api.write_tree_branch_statistics_table
        is write_tree_branch_statistics_table
    )
