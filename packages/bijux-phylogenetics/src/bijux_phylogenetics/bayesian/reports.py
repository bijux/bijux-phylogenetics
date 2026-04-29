from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from bijux_phylogenetics.bayesian.beast import (
    detect_impossible_calibration_constraints,
    validate_fossil_calibration_table,
    validate_tip_dating_metadata,
)
from bijux_phylogenetics.bayesian.mrbayes import (
    assess_mrbayes_convergence,
    compute_mrbayes_effective_sample_sizes,
    parse_mrbayes_parameter_traces,
    summarize_mrbayes_posterior_trees,
)
from bijux_phylogenetics.render.html import write_html_report
from bijux_phylogenetics.tree_set import compute_clade_frequency_table


@dataclass(slots=True)
class BayesianPosteriorReportBuildResult:
    output_path: Path
    report_kind: str
    title: str
    posterior_tree_path: Path
    trace_path: Path
    kept_tree_count: int
    warning_count: int
    machine_manifest: dict[str, object]


@dataclass(slots=True)
class CalibrationAuditReportBuildResult:
    output_path: Path
    report_kind: str
    title: str
    tree_path: Path
    calibration_path: Path
    tip_dates_path: Path | None
    invalid_calibration_count: int
    warning_count: int
    machine_manifest: dict[str, object]


def render_bayesian_posterior_report(
    *,
    posterior_tree_path: Path,
    trace_path: Path,
    out_path: Path,
    burnin_fraction: float = 0.25,
    ess_threshold: float = 200.0,
    mean_shift_threshold: float = 0.5,
) -> BayesianPosteriorReportBuildResult:
    """Render a deterministic HTML report for one posterior tree and trace pair."""
    _, posterior_summary = summarize_mrbayes_posterior_trees(
        posterior_tree_path,
        burnin_fraction=burnin_fraction,
    )
    trace_report = parse_mrbayes_parameter_traces(trace_path)
    ess_report = compute_mrbayes_effective_sample_sizes(trace_path)
    convergence = assess_mrbayes_convergence(
        trace_path,
        ess_threshold=ess_threshold,
        mean_shift_threshold=mean_shift_threshold,
    )
    clade_frequencies = compute_clade_frequency_table(posterior_summary.filtered_tree_set_path)
    title = "Bijux Bayesian Posterior Report"
    sections = [
        ("posterior-summary", json.dumps(asdict(posterior_summary), default=str, indent=2, sort_keys=True)),
        ("trace-summary", json.dumps(asdict(trace_report), default=str, indent=2, sort_keys=True)),
        ("effective-sample-sizes", json.dumps(asdict(ess_report), default=str, indent=2, sort_keys=True)),
        ("convergence", json.dumps(asdict(convergence), default=str, indent=2, sort_keys=True)),
        ("clade-frequencies", json.dumps(asdict(clade_frequencies), default=str, indent=2, sort_keys=True)),
    ]
    machine_manifest = {
        "report_kind": "bayesian-posterior",
        "title": title,
        "posterior_tree_path": str(posterior_tree_path),
        "trace_path": str(trace_path),
        "filtered_tree_set_path": str(posterior_summary.filtered_tree_set_path),
        "kept_tree_count": posterior_summary.kept_tree_count,
        "warning_count": len(convergence.warnings),
        "sections": [name for name, _ in sections],
    }
    write_html_report(title=title, sections=sections, out_path=out_path, embedded_json=machine_manifest)
    return BayesianPosteriorReportBuildResult(
        output_path=out_path,
        report_kind="bayesian-posterior",
        title=title,
        posterior_tree_path=posterior_tree_path,
        trace_path=trace_path,
        kept_tree_count=posterior_summary.kept_tree_count,
        warning_count=len(convergence.warnings),
        machine_manifest=machine_manifest,
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
    sections = [
        ("fossil-calibrations", json.dumps(asdict(calibration_report), default=str, indent=2, sort_keys=True)),
        ("impossible-constraints", json.dumps(asdict(impossible), default=str, indent=2, sort_keys=True)),
    ]
    if tip_dates is not None:
        sections.append(("tip-dates", json.dumps(asdict(tip_dates), default=str, indent=2, sort_keys=True)))
    warning_count = len(impossible.issues) + (0 if tip_dates is None else len(tip_dates.issues))
    machine_manifest = {
        "report_kind": "calibration-audit",
        "title": title,
        "tree_path": str(tree_path),
        "calibration_path": str(calibration_path),
        "tip_dates_path": None if tip_dates_path is None else str(tip_dates_path),
        "invalid_calibration_count": calibration_report.invalid_calibration_count,
        "warning_count": warning_count,
        "sections": [name for name, _ in sections],
    }
    write_html_report(title=title, sections=sections, out_path=out_path, embedded_json=machine_manifest)
    return CalibrationAuditReportBuildResult(
        output_path=out_path,
        report_kind="calibration-audit",
        title=title,
        tree_path=tree_path,
        calibration_path=calibration_path,
        tip_dates_path=tip_dates_path,
        invalid_calibration_count=calibration_report.invalid_calibration_count,
        warning_count=warning_count,
        machine_manifest=machine_manifest,
    )
