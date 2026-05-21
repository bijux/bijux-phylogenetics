from __future__ import annotations

from dataclasses import asdict
from hashlib import sha256
import json
from pathlib import Path

from bijux_phylogenetics.diagnostics.validation import (
    TreeForensicReport,
    TreeInspectionReport,
    TreeValidationReport,
)
from bijux_phylogenetics.render.tree_svg import (
    SupportLabelRenderAudit,
    TreeRenderResult,
)
from bijux_phylogenetics.reports.methods import (
    TreeValidationMethodsSummaryTextResult,
)
from bijux_phylogenetics.reports.review import ReviewerAuditChecklist

from .contracts import TreeBranchStatisticsRow, TreeSupportRow


def checksum(path: Path) -> str:
    """Return the stable checksum for one report artifact or input path."""
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_machine_manifest(
    *,
    tree_path: Path,
    title: str,
    report_path: Path,
    figure_path: Path,
    methods_summary_path: Path,
    reviewer_audit_checklist_path: Path,
    support_table_path: Path,
    clade_table_path: Path,
    branch_stats_path: Path,
    validation: TreeValidationReport,
    inspection: TreeInspectionReport,
    forensic: TreeForensicReport,
    figure: TreeRenderResult,
    support_audit: SupportLabelRenderAudit,
    methods_summary: TreeValidationMethodsSummaryTextResult,
    support_rows: list[TreeSupportRow],
    branch_stats: TreeBranchStatisticsRow,
    reviewer_summary: list[str],
    limitations: list[str],
) -> dict[str, object]:
    """Build the machine manifest before the reviewer checklist is attached."""
    return {
        "report_kind": "tree_package",
        "title": title,
        "input_path": str(tree_path),
        "input_checksum": checksum(tree_path),
        "outputs": {
            "report_path": str(report_path),
            "figure_path": str(figure_path),
            "methods_summary_path": str(methods_summary_path),
            "reviewer_audit_checklist_path": str(reviewer_audit_checklist_path),
            "support_table_path": str(support_table_path),
            "clade_table_path": str(clade_table_path),
            "branch_stats_path": str(branch_stats_path),
        },
        "metrics": {
            "tip_count": inspection.tip_count,
            "node_count": inspection.node_count,
            "clade_count": inspection.clade_count,
            "supported_branch_count": sum(
                1 for row in support_rows if row.support is not None
            ),
            "long_outlier_count": branch_stats.long_outlier_count,
            "short_outlier_count": branch_stats.short_outlier_count,
            "rendered_support_count": figure.rendered_support_count,
        },
        "reviewer_summary": reviewer_summary,
        "limitations": limitations,
        "methods_summary_text": methods_summary.text,
        "validation": asdict(validation),
        "inspection": asdict(inspection),
        "forensic": asdict(forensic),
        "support_audit": asdict(support_audit),
    }


def attach_reviewer_audit_checklist(
    *,
    machine_manifest: dict[str, object],
    reviewer_audit_checklist: ReviewerAuditChecklist,
) -> dict[str, object]:
    """Attach the reviewer checklist payload to the machine manifest."""
    machine_manifest["reviewer_audit_checklist"] = asdict(reviewer_audit_checklist)
    return machine_manifest


def write_machine_manifest(path: Path, machine_manifest: dict[str, object]) -> Path:
    """Write the machine manifest artifact."""
    path.write_text(
        json.dumps(machine_manifest, default=str, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path
