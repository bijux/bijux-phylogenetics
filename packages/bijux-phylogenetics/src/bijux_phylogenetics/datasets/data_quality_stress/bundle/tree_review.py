from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows

from ..models import (
    CatarrhineDataQualityStressPanelWorkflowReport,
    DataQualityRepairAction,
)
from .shared import _format_number, _tree_warning_nodes


def _write_tree_issues_table(
    path: Path,
    report: CatarrhineDataQualityStressPanelWorkflowReport,
) -> Path:
    rows: list[dict[str, str]] = []
    if report.raw_tree_validation.zero_length_branches:
        rows.append(
            {
                "issue_code": "zero_length_branches",
                "severity": "warning",
                "affected_taxa": ",".join(report.repaired_branch_nodes),
                "affected_nodes": ",".join(report.repaired_branch_nodes),
                "raw_value": str(report.raw_tree_validation.zero_length_branches),
                "action": "apply_branch_length_floor_in_cleaned_tree",
            }
        )
    if report.raw_tree_validation.negative_branch_lengths:
        rows.append(
            {
                "issue_code": "negative_branch_lengths",
                "severity": "warning",
                "affected_taxa": "",
                "affected_nodes": ",".join(
                    _tree_warning_nodes(
                        report.raw_tree_validation,
                        warning_code="negative_branch_lengths",
                    )
                ),
                "raw_value": str(report.raw_tree_validation.negative_branch_lengths),
                "action": "apply_branch_length_floor_in_cleaned_tree",
            }
        )
    for outlier in report.raw_tree_inspection.long_branch_outliers:
        rows.append(
            {
                "issue_code": "long_branch_outlier",
                "severity": "warning",
                "affected_taxa": outlier.node
                if outlier.branch_type == "terminal"
                else "",
                "affected_nodes": outlier.node,
                "raw_value": _format_number(outlier.branch_length),
                "action": (
                    "drop_taxon_from_cleaned_tree"
                    if outlier.branch_type == "terminal"
                    else "flag_only"
                ),
            }
        )
    return write_taxon_rows(
        path,
        columns=[
            "issue_code",
            "severity",
            "affected_taxa",
            "affected_nodes",
            "raw_value",
            "action",
        ],
        rows=rows,
    )


def _write_repair_actions_table(
    path: Path,
    rows: list[DataQualityRepairAction],
) -> Path:
    return write_taxon_rows(
        path,
        columns=["action_kind", "affected_taxa", "affected_nodes", "reason", "result"],
        rows=[
            {
                "action_kind": row.action_kind,
                "affected_taxa": ",".join(row.affected_taxa),
                "affected_nodes": ",".join(row.affected_nodes),
                "reason": row.reason,
                "result": row.result,
            }
            for row in rows
        ],
    )
