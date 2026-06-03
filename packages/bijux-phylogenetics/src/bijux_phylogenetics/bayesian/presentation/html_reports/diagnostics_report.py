from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path
import tempfile

from bijux_phylogenetics.bayesian.beast.logs import (
    assess_beast_burnin_sensitivity,
    assess_beast_chain_mixing,
    validate_beast_posterior_log,
)
from bijux_phylogenetics.bayesian.beast.validation import (
    detect_impossible_calibration_constraints,
    validate_fossil_calibration_table,
    validate_tip_dating_metadata,
)
from bijux_phylogenetics.bayesian.beast.xml_analysis import summarize_beast_analysis_xml
from bijux_phylogenetics.evidence.provenance.method_tiers import (
    bayesian_report_method_tier,
)
from bijux_phylogenetics.render.html import write_html_report

from ..posterior_uncertainty import (
    summarize_beast_workflow_evidence,
    write_bayesian_limitations_text,
    write_bayesian_methods_summary_text,
    write_supplementary_bayesian_diagnostics_table,
)
from .contracts import BayesianDiagnosticsReportBuildResult
from .report_policy import method_tier_section, method_tier_summary_metrics


def render_bayesian_diagnostics_report(
    *,
    posterior_tree_path: Path,
    primary_log_path: Path,
    out_path: Path,
    additional_log_paths: list[Path] | None = None,
    analysis_xml_path: Path | None = None,
    tree_path: Path | None = None,
    calibration_path: Path | None = None,
    tip_dates_path: Path | None = None,
    alignment_path: Path | None = None,
    taxon_column: str | None = None,
    date_column: str = "date",
    burnin_fractions: tuple[float, ...] = (0.1, 0.25, 0.5),
    required_columns: tuple[str, ...] = ("posterior", "likelihood"),
    ess_threshold: float = 200.0,
    mean_shift_threshold: float = 0.5,
    cross_chain_mean_shift_threshold: float = 0.75,
) -> BayesianDiagnosticsReportBuildResult:
    """Render a reviewer-facing Bayesian diagnostics report for posterior logs and calibrations."""
    validation = validate_beast_posterior_log(
        primary_log_path, required_columns=required_columns
    )
    burnin = assess_beast_burnin_sensitivity(
        posterior_tree_path,
        log_path=primary_log_path,
        burnin_fractions=burnin_fractions,
    )
    log_paths = [primary_log_path, *(additional_log_paths or [])]
    resolved_analysis_xml_path = _resolve_beast_analysis_xml_path(
        analysis_xml_path=analysis_xml_path,
        primary_log_path=primary_log_path,
    )
    analysis_summary = (
        None
        if resolved_analysis_xml_path is None
        else summarize_beast_analysis_xml(resolved_analysis_xml_path)
    )
    workflow_evidence = summarize_beast_workflow_evidence(
        posterior_tree_path=posterior_tree_path,
        primary_log_path=primary_log_path,
        analysis_xml_path=resolved_analysis_xml_path,
    )
    mixing = assess_beast_chain_mixing(
        log_paths,
        ess_threshold=ess_threshold,
        mean_shift_threshold=mean_shift_threshold,
        cross_chain_mean_shift_threshold=cross_chain_mean_shift_threshold,
    )
    calibration_report = None
    impossible_report = None
    tip_date_report = None
    if tree_path is not None and calibration_path is not None:
        calibration_report = validate_fossil_calibration_table(
            tree_path, calibration_path
        )
        impossible_report = detect_impossible_calibration_constraints(
            tree_path, calibration_path
        )
    if tree_path is not None and tip_dates_path is not None:
        tip_date_report = validate_tip_dating_metadata(
            tree_path,
            tip_dates_path,
            alignment_path=alignment_path,
            taxon_column=taxon_column,
            date_column=date_column,
        )
    supplementary_table_path = Path(
        tempfile.mkstemp(prefix="bijux-bayesian-diagnostics-", suffix=".tsv")[1]
    )
    methods_text_path = Path(
        tempfile.mkstemp(prefix="bijux-bayesian-methods-", suffix=".md")[1]
    )
    supplementary_table = write_supplementary_bayesian_diagnostics_table(
        supplementary_table_path,
        posterior_tree_path=posterior_tree_path,
        primary_log_path=primary_log_path,
        additional_log_paths=additional_log_paths,
        burnin_fractions=burnin_fractions,
        required_columns=required_columns,
        ess_threshold=ess_threshold,
        mean_shift_threshold=mean_shift_threshold,
        cross_chain_mean_shift_threshold=cross_chain_mean_shift_threshold,
    )
    methods_summary = write_bayesian_methods_summary_text(
        methods_text_path,
        posterior_tree_path=posterior_tree_path,
        primary_log_path=primary_log_path,
        additional_log_paths=additional_log_paths,
        analysis_xml_path=resolved_analysis_xml_path,
        burnin_fractions=burnin_fractions,
        ess_threshold=ess_threshold,
        mean_shift_threshold=mean_shift_threshold,
        cross_chain_mean_shift_threshold=cross_chain_mean_shift_threshold,
        calibration_path=calibration_path,
        tip_dates_path=tip_dates_path,
    )
    limitations_text_path = Path(
        tempfile.mkstemp(prefix="bijux-bayesian-limitations-", suffix=".md")[1]
    )
    limitations_summary = write_bayesian_limitations_text(
        limitations_text_path,
        posterior_tree_path=posterior_tree_path,
        primary_log_path=primary_log_path,
        additional_log_paths=additional_log_paths,
        tree_path=tree_path,
        calibration_path=calibration_path,
        tip_dates_path=tip_dates_path,
        alignment_path=alignment_path,
        burnin_fractions=burnin_fractions,
        ess_threshold=ess_threshold,
        mean_shift_threshold=mean_shift_threshold,
        cross_chain_mean_shift_threshold=cross_chain_mean_shift_threshold,
    )
    title = "Bijux Bayesian Diagnostics Report"
    method_tier = bayesian_report_method_tier("bayesian-diagnostics")
    sections = [
        method_tier_section(method_tier),
        *(
            [
                (
                    "workflow-evidence",
                    json.dumps(
                        asdict(workflow_evidence), default=str, indent=2, sort_keys=True
                    ),
                ),
                (
                    "analysis-assumptions",
                    json.dumps(
                        asdict(analysis_summary), default=str, indent=2, sort_keys=True
                    ),
                ),
            ]
            if analysis_summary is not None
            else []
        ),
        *(
            []
            if analysis_summary is not None
            else [
                (
                    "workflow-evidence",
                    json.dumps(
                        asdict(workflow_evidence), default=str, indent=2, sort_keys=True
                    ),
                )
            ]
        ),
        (
            "posterior-log-validation",
            json.dumps(asdict(validation), default=str, indent=2, sort_keys=True),
        ),
        (
            "burnin-sensitivity",
            json.dumps(asdict(burnin), default=str, indent=2, sort_keys=True),
        ),
        (
            "chain-mixing",
            json.dumps(asdict(mixing), default=str, indent=2, sort_keys=True),
        ),
        (
            "supplementary-diagnostics-table",
            supplementary_table.output_path.read_text(encoding="utf-8"),
        ),
        ("methods-summary-text", methods_summary.text),
        ("limitations-text", limitations_summary.text),
    ]
    warning_count = len(validation.issues) + len(burnin.warnings) + len(mixing.issues)
    if analysis_summary is not None:
        warning_count += len(analysis_summary.issues)
    if calibration_report is not None:
        sections.append(
            (
                "fossil-calibrations",
                json.dumps(
                    asdict(calibration_report), default=str, indent=2, sort_keys=True
                ),
            )
        )
        warning_count += calibration_report.invalid_calibration_count
    if impossible_report is not None:
        sections.append(
            (
                "impossible-constraints",
                json.dumps(
                    asdict(impossible_report), default=str, indent=2, sort_keys=True
                ),
            )
        )
        warning_count += len(impossible_report.issues)
    if tip_date_report is not None:
        sections.append(
            (
                "tip-dates",
                json.dumps(
                    asdict(tip_date_report), default=str, indent=2, sort_keys=True
                ),
            )
        )
        warning_count += len(tip_date_report.issues)
    machine_manifest = {
        "report_kind": "bayesian-diagnostics",
        "title": title,
        "posterior_tree_path": str(posterior_tree_path),
        "primary_log_path": str(primary_log_path),
        "analysis_xml_path": None
        if resolved_analysis_xml_path is None
        else str(resolved_analysis_xml_path),
        "chain_count": len(log_paths),
        "warning_count": warning_count,
        "limitations": [
            line.removeprefix("- ").strip()
            for line in limitations_summary.text.splitlines()
            if line.startswith("- ")
        ],
        "sections": [name for name, _ in sections],
    }
    write_html_report(
        title=title,
        sections=sections,
        out_path=out_path,
        embedded_json=machine_manifest,
        summary_metrics=method_tier_summary_metrics(method_tier),
    )
    supplementary_table_path.unlink(missing_ok=True)
    methods_text_path.unlink(missing_ok=True)
    limitations_text_path.unlink(missing_ok=True)
    return BayesianDiagnosticsReportBuildResult(
        output_path=out_path,
        report_kind="bayesian-diagnostics",
        title=title,
        posterior_tree_path=posterior_tree_path,
        primary_log_path=primary_log_path,
        chain_count=len(log_paths),
        warning_count=warning_count,
        method_tier=method_tier,
        machine_manifest=machine_manifest,
    )


def _resolve_beast_analysis_xml_path(
    *,
    analysis_xml_path: Path | None,
    primary_log_path: Path,
) -> Path | None:
    if analysis_xml_path is not None:
        return analysis_xml_path
    direct_candidate = primary_log_path.with_suffix(".xml")
    if direct_candidate.exists():
        return direct_candidate
    seeded_candidate = primary_log_path.with_name(
        _strip_beast_output_suffix(primary_log_path.name) + ".xml"
    )
    if seeded_candidate.exists():
        return seeded_candidate
    return None


def _strip_beast_output_suffix(file_name: str) -> str:
    for suffix in (".log", ".trees"):
        if not file_name.endswith(suffix):
            continue
        stem = file_name[: -len(suffix)]
        if stem.endswith(".$(seed)"):
            return stem[: -len(".$(seed)")]
        parts = stem.rsplit(".", 1)
        if len(parts) == 2 and parts[1].isdigit():
            return parts[0]
        return stem
    return file_name
