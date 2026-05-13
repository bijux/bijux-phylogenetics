from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path
import tempfile

from bijux_phylogenetics.bayesian.beast import (
    assess_beast_burnin_sensitivity,
    assess_beast_chain_mixing,
    assess_time_tree_readiness,
    detect_impossible_calibration_constraints,
    summarize_beast_analysis_xml,
    validate_beast_posterior_log,
    validate_fossil_calibration_table,
    validate_tip_dating_metadata,
)
from bijux_phylogenetics.bayesian.comparison import (
    compare_independent_bayesian_runs,
    compare_ml_tree_to_bayesian_posterior,
)
from bijux_phylogenetics.bayesian.mrbayes import (
    assess_mrbayes_convergence,
    compute_mrbayes_effective_sample_sizes,
    parse_mrbayes_parameter_traces,
    summarize_mrbayes_posterior_trees,
)
from bijux_phylogenetics.bayesian.uncertainty import (
    write_bayesian_limitations_text,
    write_bayesian_methods_summary_text,
    write_supplementary_bayesian_diagnostics_table,
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


@dataclass(slots=True)
class BayesianRunComparisonReportBuildResult:
    output_path: Path
    report_kind: str
    title: str
    left_tree_set_path: Path
    right_tree_set_path: Path
    left_trace_path: Path
    right_trace_path: Path
    trace_kind: str
    warning_count: int
    machine_manifest: dict[str, object]


@dataclass(slots=True)
class BayesianDiagnosticsReportBuildResult:
    output_path: Path
    report_kind: str
    title: str
    posterior_tree_path: Path
    primary_log_path: Path
    chain_count: int
    warning_count: int
    machine_manifest: dict[str, object]


@dataclass(slots=True)
class BayesianMlComparisonReportBuildResult:
    output_path: Path
    report_kind: str
    title: str
    ml_tree_path: Path
    posterior_tree_path: Path
    warning_count: int
    machine_manifest: dict[str, object]


@dataclass(slots=True)
class TimeTreeReadinessReportBuildResult:
    output_path: Path
    report_kind: str
    title: str
    tree_path: Path
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
    clade_frequencies = compute_clade_frequency_table(
        posterior_summary.filtered_tree_set_path
    )
    title = "Bijux Bayesian Posterior Report"
    sections = [
        (
            "posterior-summary",
            json.dumps(
                asdict(posterior_summary), default=str, indent=2, sort_keys=True
            ),
        ),
        (
            "trace-summary",
            json.dumps(asdict(trace_report), default=str, indent=2, sort_keys=True),
        ),
        (
            "effective-sample-sizes",
            json.dumps(asdict(ess_report), default=str, indent=2, sort_keys=True),
        ),
        (
            "convergence",
            json.dumps(asdict(convergence), default=str, indent=2, sort_keys=True),
        ),
        (
            "clade-frequencies",
            json.dumps(
                asdict(clade_frequencies), default=str, indent=2, sort_keys=True
            ),
        ),
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
    write_html_report(
        title=title,
        sections=sections,
        out_path=out_path,
        embedded_json=machine_manifest,
    )
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
        "sections": [name for name, _ in sections],
    }
    write_html_report(
        title=title,
        sections=sections,
        out_path=out_path,
        embedded_json=machine_manifest,
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
        machine_manifest=machine_manifest,
    )


def render_bayesian_run_comparison_report(
    *,
    left_tree_set_path: Path,
    right_tree_set_path: Path,
    left_trace_path: Path,
    right_trace_path: Path,
    out_path: Path,
    trace_kind: str = "mrbayes",
    burnin_fraction: float = 0.25,
    ess_threshold: float = 200.0,
    mean_shift_threshold: float = 0.5,
) -> BayesianRunComparisonReportBuildResult:
    """Render a deterministic HTML report comparing two Bayesian runs."""
    comparison = compare_independent_bayesian_runs(
        left_tree_set_path,
        right_tree_set_path,
        left_trace_path=left_trace_path,
        right_trace_path=right_trace_path,
        trace_kind=trace_kind,
        burnin_fraction=burnin_fraction,
        ess_threshold=ess_threshold,
        mean_shift_threshold=mean_shift_threshold,
    )
    title = "Bijux Bayesian Run Comparison Report"
    sections = [
        (
            "run-comparison",
            json.dumps(asdict(comparison), default=str, indent=2, sort_keys=True),
        ),
        (
            "tree-comparison",
            json.dumps(
                asdict(comparison.tree_comparison),
                default=str,
                indent=2,
                sort_keys=True,
            ),
        ),
        (
            "left-convergence",
            json.dumps(
                asdict(comparison.left_convergence),
                default=str,
                indent=2,
                sort_keys=True,
            ),
        ),
        (
            "right-convergence",
            json.dumps(
                asdict(comparison.right_convergence),
                default=str,
                indent=2,
                sort_keys=True,
            ),
        ),
        (
            "parameter-differences",
            json.dumps(
                [asdict(row) for row in comparison.parameter_differences],
                indent=2,
                sort_keys=True,
            ),
        ),
    ]
    machine_manifest = {
        "report_kind": "bayesian-run-comparison",
        "title": title,
        "left_tree_set_path": str(left_tree_set_path),
        "right_tree_set_path": str(right_tree_set_path),
        "left_trace_path": str(left_trace_path),
        "right_trace_path": str(right_trace_path),
        "trace_kind": trace_kind,
        "warning_count": len(comparison.warnings),
        "sections": [name for name, _ in sections],
    }
    write_html_report(
        title=title,
        sections=sections,
        out_path=out_path,
        embedded_json=machine_manifest,
    )
    return BayesianRunComparisonReportBuildResult(
        output_path=out_path,
        report_kind="bayesian-run-comparison",
        title=title,
        left_tree_set_path=left_tree_set_path,
        right_tree_set_path=right_tree_set_path,
        left_trace_path=left_trace_path,
        right_trace_path=right_trace_path,
        trace_kind=trace_kind,
        warning_count=len(comparison.warnings),
        machine_manifest=machine_manifest,
    )


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
    sections = [
        *(
            [
                (
                    "analysis-assumptions",
                    json.dumps(
                        asdict(analysis_summary), default=str, indent=2, sort_keys=True
                    ),
                )
            ]
            if analysis_summary is not None
            else []
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
        "sections": [name for name, _ in sections],
    }
    write_html_report(
        title=title,
        sections=sections,
        out_path=out_path,
        embedded_json=machine_manifest,
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


def render_ml_vs_bayesian_tree_report(
    *,
    ml_tree_path: Path,
    posterior_tree_path: Path,
    out_path: Path,
    burnin_fraction: float = 0.25,
) -> BayesianMlComparisonReportBuildResult:
    """Render a reviewer-facing comparison between one ML tree and one Bayesian MCC summary."""
    comparison = compare_ml_tree_to_bayesian_posterior(
        ml_tree_path,
        posterior_tree_path,
        burnin_fraction=burnin_fraction,
    )
    title = "Bijux ML Versus Bayesian Tree Report"
    sections = [
        (
            "ml-versus-bayesian-summary",
            json.dumps(asdict(comparison), default=str, indent=2, sort_keys=True),
        ),
        (
            "topology-comparison",
            json.dumps(
                asdict(comparison.topology), default=str, indent=2, sort_keys=True
            ),
        ),
        (
            "branch-length-comparison",
            json.dumps(
                asdict(comparison.branch_lengths), default=str, indent=2, sort_keys=True
            ),
        ),
    ]
    machine_manifest = {
        "report_kind": "ml-vs-bayesian-tree",
        "title": title,
        "ml_tree_path": str(ml_tree_path),
        "posterior_tree_path": str(posterior_tree_path),
        "warning_count": len(comparison.warnings),
        "sections": [name for name, _ in sections],
    }
    write_html_report(
        title=title,
        sections=sections,
        out_path=out_path,
        embedded_json=machine_manifest,
    )
    return BayesianMlComparisonReportBuildResult(
        output_path=out_path,
        report_kind="ml-vs-bayesian-tree",
        title=title,
        ml_tree_path=ml_tree_path,
        posterior_tree_path=posterior_tree_path,
        warning_count=len(comparison.warnings),
        machine_manifest=machine_manifest,
    )


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
    sections = [
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
        "sections": [name for name, _ in sections],
    }
    write_html_report(
        title=title,
        sections=sections,
        out_path=out_path,
        embedded_json=machine_manifest,
    )
    return TimeTreeReadinessReportBuildResult(
        output_path=out_path,
        report_kind="time-tree-readiness",
        title=title,
        tree_path=tree_path,
        warning_count=warning_count,
        machine_manifest=machine_manifest,
    )
