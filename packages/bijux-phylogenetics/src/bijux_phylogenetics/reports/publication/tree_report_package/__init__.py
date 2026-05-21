from __future__ import annotations

from dataclasses import asdict
from hashlib import sha256
import json
from pathlib import Path

from bijux_phylogenetics.diagnostics.validation import (
    TreeForensicReport,
    TreeInspectionReport,
    TreeValidationReport,
    forensic_tree_path,
    inspect_tree_path,
    validate_tree_path,
)
from bijux_phylogenetics.evidence.provenance.method_tiers import (
    MethodTierAssessment,
    tree_report_method_tier,
)
from bijux_phylogenetics.render.tree_svg import (
    audit_support_label_rendering,
    render_tree_svg,
)
from bijux_phylogenetics.reports.methods import (
    write_tree_validation_methods_summary_text,
)
from bijux_phylogenetics.reports.review import (
    ReviewerAuditChecklist,
    write_reviewer_audit_checklist,
)
from bijux_phylogenetics.trees import (
    CladeTableRow,
    analyze_branch_length_distribution,
    extract_tree_clades,
)

from .contracts import (
    TreeBranchStatisticsRow,
    TreeReportPackageResult,
    TreeSupportRow,
)
from .artifact_outputs import (
    write_tree_branch_statistics_table,
    write_tree_support_table,
)
from .review_context import (
    build_reviewer_summary,
    summarize_tree_branch_statistics,
    summarize_tree_support,
)
from .presentation import write_tree_report_html


def _checksum(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_tree_report_package(
    tree_path: Path,
    *,
    out_dir: Path,
    title: str = "Bijux Full Tree Report",
) -> TreeReportPackageResult:
    out_dir.mkdir(parents=True, exist_ok=True)
    report_path = out_dir / "tree-report.html"
    figure_path = out_dir / "tree-image.svg"
    methods_summary_path = out_dir / "tree-validation-methods-summary.md"
    reviewer_audit_checklist_path = out_dir / "reviewer-audit-checklist.tsv"
    support_table_path = out_dir / "support-table.tsv"
    clade_table_path = out_dir / "clade-table.tsv"
    branch_stats_path = out_dir / "branch-stats.tsv"
    manifest_path = out_dir / "tree-report.manifest.json"

    validation = validate_tree_path(tree_path)
    inspection = inspect_tree_path(tree_path)
    forensic = forensic_tree_path(tree_path)
    support_audit = audit_support_label_rendering(tree_path)
    figure = render_tree_svg(
        tree_path,
        out_path=figure_path,
        layout="phylogram",
        show_support_values=support_audit.validated,
        validated_support_labels=support_audit.labels_by_node,
        support_validation_warnings=support_audit.warnings,
    )
    clades = extract_tree_clades(tree_path)
    support_rows = summarize_tree_support(clades)
    branch_lengths = analyze_branch_length_distribution(tree_path)
    branch_stats = summarize_tree_branch_statistics(branch_lengths)
    methods_summary = write_tree_validation_methods_summary_text(
        methods_summary_path,
        tree_path=tree_path,
    )
    method_tier = tree_report_method_tier()
    reviewer_summary, limitations = build_reviewer_summary(
        inspection=inspection,
        support_rows=support_rows,
        branch_stats=branch_stats,
        support_audit=support_audit,
    )

    write_tree_support_table(support_table_path, support_rows)
    from bijux_phylogenetics.trees import write_clade_table

    write_clade_table(clade_table_path, clades)
    write_tree_branch_statistics_table(branch_stats_path, branch_stats)

    machine_manifest = {
        "report_kind": "tree_package",
        "title": title,
        "input_path": str(tree_path),
        "input_checksum": _checksum(tree_path),
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
    reviewer_audit_checklist = write_reviewer_audit_checklist(
        reviewer_audit_checklist_path,
        machine_manifest,
    ).checklist
    machine_manifest["reviewer_audit_checklist"] = asdict(reviewer_audit_checklist)
    manifest_path.write_text(
        json.dumps(machine_manifest, default=str, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_tree_report_html(
        path=report_path,
        title=title,
        figure_svg=figure_path.read_text(encoding="utf-8"),
        methods_summary_text=methods_summary.text,
        reviewer_summary=reviewer_summary,
        limitations=limitations,
        support_rows=support_rows,
        clade_rows=clades.rows,
        branch_stats=branch_stats,
        method_tier=method_tier,
        reviewer_audit_checklist=reviewer_audit_checklist,
        validation=validation,
        inspection=inspection,
        forensic=forensic,
        manifest=machine_manifest,
    )
    return TreeReportPackageResult(
        output_dir=out_dir,
        report_path=report_path,
        figure_path=figure_path,
        methods_summary_path=methods_summary_path,
        reviewer_audit_checklist_path=reviewer_audit_checklist_path,
        support_table_path=support_table_path,
        clade_table_path=clade_table_path,
        branch_stats_path=branch_stats_path,
        manifest_path=manifest_path,
        validation=validation,
        inspection=inspection,
        forensic=forensic,
        figure=figure,
        support_audit=support_audit,
        clades=clades,
        branch_lengths=branch_lengths,
        support_rows=support_rows,
        branch_stats=branch_stats,
        method_tier=method_tier,
        reviewer_summary=reviewer_summary,
        limitations=limitations,
        methods_summary=methods_summary,
        reviewer_audit_checklist=reviewer_audit_checklist,
        machine_manifest=machine_manifest,
    )
