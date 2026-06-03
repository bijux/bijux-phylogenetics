from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from bijux_phylogenetics.bayesian.beast.validation import (
    detect_impossible_calibration_constraints,
    validate_fossil_calibration_table,
    validate_tip_dating_metadata,
)
from bijux_phylogenetics.evidence.provenance.method_tiers import (
    bayesian_report_method_tier,
)
from bijux_phylogenetics.render.html import write_html_report

from .contracts import CalibrationAuditReportBuildResult
from .report_policy import (
    calibration_audit_limitations,
    method_tier_section,
    method_tier_summary_metrics,
)


def render_calibration_audit_report(
    *,
    tree_path: Path,
    calibration_path: Path,
    out_path: Path,
    tip_dates_path: Path | None = None,
    alignment_path: Path | None = None,
    taxon_column: str | None = None,
    date_column: str = "date",
) -> CalibrationAuditReportBuildResult:
    """Render a deterministic HTML report for calibration and tip-date validation."""
    calibration_report = validate_fossil_calibration_table(tree_path, calibration_path)
    impossible = detect_impossible_calibration_constraints(tree_path, calibration_path)
    tip_dates = (
        validate_tip_dating_metadata(
            tree_path,
            tip_dates_path,
            alignment_path=alignment_path,
            taxon_column=taxon_column,
            date_column=date_column,
        )
        if tip_dates_path is not None
        else None
    )
    title = "Bijux Calibration Audit Report"
    method_tier = bayesian_report_method_tier("calibration-audit")
    limitations = calibration_audit_limitations(
        invalid_calibration_count=calibration_report.invalid_calibration_count,
        impossible_constraint_count=len(impossible.issues),
        invalid_tip_count=0 if tip_dates is None else tip_dates.invalid_tip_count,
    )
    sections = [
        method_tier_section(method_tier),
        (
            "fossil-calibrations",
            json.dumps(
                asdict(calibration_report), default=str, indent=2, sort_keys=True
            ),
        ),
        (
            "impossible-constraints",
            json.dumps(asdict(impossible), default=str, indent=2, sort_keys=True),
        ),
    ]
    if tip_dates is not None:
        sections.append(
            (
                "tip-dates",
                json.dumps(asdict(tip_dates), default=str, indent=2, sort_keys=True),
            )
        )
    sections.append(("limitations", json.dumps(limitations, indent=2)))
    warning_count = len(impossible.issues) + (
        0 if tip_dates is None else len(tip_dates.issues)
    )
    machine_manifest = {
        "report_kind": "calibration-audit",
        "title": title,
        "tree_path": str(tree_path),
        "calibration_path": str(calibration_path),
        "tip_dates_path": None if tip_dates_path is None else str(tip_dates_path),
        "invalid_calibration_count": calibration_report.invalid_calibration_count,
        "warning_count": warning_count,
        "limitations": limitations,
        "sections": [name for name, _ in sections],
    }
    write_html_report(
        title=title,
        sections=sections,
        out_path=out_path,
        embedded_json=machine_manifest,
        summary_metrics=method_tier_summary_metrics(method_tier),
    )
    return CalibrationAuditReportBuildResult(
        output_path=out_path,
        report_kind="calibration-audit",
        title=title,
        tree_path=tree_path,
        calibration_path=calibration_path,
        tip_dates_path=tip_dates_path,
        invalid_calibration_count=calibration_report.invalid_calibration_count,
        warning_count=warning_count,
        method_tier=method_tier,
        machine_manifest=machine_manifest,
    )
