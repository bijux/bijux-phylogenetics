from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from bijux_phylogenetics.bayesian.beast.validation import assess_time_tree_readiness
from bijux_phylogenetics.evidence.provenance.method_tiers import (
    bayesian_report_method_tier,
)
from bijux_phylogenetics.render.html import write_html_report

from .contracts import TimeTreeReadinessReportBuildResult
from .report_policy import method_tier_section, method_tier_summary_metrics


def render_time_tree_readiness_report(
    *,
    tree_path: Path,
    out_path: Path,
    calibration_path: Path | None = None,
    tip_dates_path: Path | None = None,
    alignment_path: Path | None = None,
    taxon_column: str | None = None,
    date_column: str = "date",
) -> TimeTreeReadinessReportBuildResult:
    """Render an HTML readiness decision for dated phylogenetics."""
    readiness = assess_time_tree_readiness(
        tree_path,
        calibration_path=calibration_path,
        tip_dates_path=tip_dates_path,
        alignment_path=alignment_path,
        taxon_column=taxon_column,
        date_column=date_column,
    )
    limitations: list[str] = []
    if readiness.calibration_dominance is not None:
        limitations.extend(readiness.calibration_dominance.warnings)
    if (
        readiness.tip_date_report is not None
        and readiness.tip_date_report.invalid_tip_count
    ):
        limitations.append(
            "tip-date metadata requires correction before the dated-tree workflow can be trusted"
        )
    title = "Bijux Time-Tree Readiness Report"
    method_tier = bayesian_report_method_tier("time-tree-readiness")
    sections = [
        method_tier_section(method_tier),
        (
            "readiness",
            json.dumps(asdict(readiness), default=str, indent=2, sort_keys=True),
        ),
    ]
    if readiness.calibration_report is not None:
        sections.append(
            (
                "fossil-calibrations",
                json.dumps(
                    asdict(readiness.calibration_report),
                    default=str,
                    indent=2,
                    sort_keys=True,
                ),
            )
        )
    if readiness.calibration_dominance is not None:
        sections.append(
            (
                "calibration-dominance",
                json.dumps(
                    asdict(readiness.calibration_dominance),
                    default=str,
                    indent=2,
                    sort_keys=True,
                ),
            )
        )
    if readiness.tip_date_report is not None:
        sections.append(
            (
                "tip-dates",
                json.dumps(
                    asdict(readiness.tip_date_report),
                    default=str,
                    indent=2,
                    sort_keys=True,
                ),
            )
        )
    sections.append(
        ("limitations", json.dumps(sorted(dict.fromkeys(limitations)), indent=2))
    )
    warning_count = len(readiness.blockers) + len(readiness.warnings)
    machine_manifest = {
        "report_kind": "time-tree-readiness",
        "title": title,
        "tree_path": str(tree_path),
        "decision": readiness.decision,
        "warning_count": warning_count,
        "limitations": sorted(dict.fromkeys(limitations)),
        "sections": [name for name, _ in sections],
    }
    write_html_report(
        title=title,
        sections=sections,
        out_path=out_path,
        embedded_json=machine_manifest,
        summary_metrics=method_tier_summary_metrics(method_tier),
    )
    return TimeTreeReadinessReportBuildResult(
        output_path=out_path,
        report_kind="time-tree-readiness",
        title=title,
        tree_path=tree_path,
        warning_count=warning_count,
        method_tier=method_tier,
        machine_manifest=machine_manifest,
    )
