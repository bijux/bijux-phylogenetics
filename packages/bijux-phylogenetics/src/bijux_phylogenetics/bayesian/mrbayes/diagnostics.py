from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.bayesian.burnin import (
    DEFAULT_BURNIN_FRACTIONS,
    normalize_burnin_fractions,
    summarize_burnin_clade_shifts,
    summarize_burnin_parameter_shifts,
)
from bijux_phylogenetics.bayesian.diagnostics import (
    TraceConvergenceReport,
    summarize_trace_convergence,
    summarize_trace_parameters,
)
from bijux_phylogenetics.core.metadata import write_taxon_rows
from bijux_phylogenetics.trees import compute_clade_frequency_table

from .models import (
    EffectiveSampleSize,
    MrBayesBurninSensitivityReport,
    MrBayesBurninSensitivitySlice,
    MrBayesConvergenceReport,
    MrBayesESSReport,
    MrBayesParameterDiagnosticsReport,
    MrBayesParameterSummary,
    MrBayesTraceReport,
    MrBayesTraceRow,
)
from .posterior_trees import summarize_mrbayes_posterior_trees
from .tabular import parse_mrbayes_parameter_traces


def compute_mrbayes_effective_sample_sizes(path: Path) -> MrBayesESSReport:
    """Compute per-parameter effective sample sizes from a MrBayes trace file."""
    report = parse_mrbayes_parameter_traces(path)
    convergence = summarize_trace_convergence(
        path=path,
        rows=[row.values for row in report.rows],
        columns=report.columns,
    )
    effective_sample_sizes = [
        EffectiveSampleSize(
            parameter=summary.parameter,
            sample_count=summary.sample_count,
            effective_sample_size=summary.effective_sample_size,
        )
        for summary in convergence.series
    ]
    return MrBayesESSReport(
        path=path,
        sample_count=report.row_count,
        effective_sample_sizes=effective_sample_sizes,
    )


def summarize_mrbayes_parameter_diagnostics(
    path: Path,
    *,
    burnin_fraction: float = 0.0,
) -> MrBayesParameterDiagnosticsReport:
    """Summarize burn-in-aware posterior parameter diagnostics from MrBayes traces."""
    report = parse_mrbayes_parameter_traces(path)
    burnin_row_count, kept_rows = _split_mrbayes_trace_rows(
        report, burnin_fraction=burnin_fraction
    )
    diagnostics = summarize_trace_parameters(
        path=path,
        rows=[row.values for row in kept_rows],
        columns=report.columns,
    )
    return MrBayesParameterDiagnosticsReport(
        path=path,
        burnin_fraction=burnin_fraction,
        burnin_row_count=burnin_row_count,
        kept_row_count=len(kept_rows),
        first_kept_generation=kept_rows[0].generation,
        last_kept_generation=kept_rows[-1].generation,
        parameter_summaries=[
            MrBayesParameterSummary(
                parameter=summary.parameter,
                sample_count=summary.sample_count,
                effective_sample_size=summary.effective_sample_size,
                mean=summary.mean,
                median=summary.median,
                standard_deviation=summary.standard_deviation,
                minimum=summary.minimum,
                maximum=summary.maximum,
                hpd_95_lower=summary.hpd_95_lower,
                hpd_95_upper=summary.hpd_95_upper,
                first_half_mean=summary.first_half_mean,
                second_half_mean=summary.second_half_mean,
                standardized_mean_shift=summary.standardized_mean_shift,
            )
            for summary in diagnostics.series
        ],
    )


def write_mrbayes_parameter_summary_table(
    path: Path,
    report: MrBayesParameterDiagnosticsReport,
) -> Path:
    """Write a reviewer-facing TSV summary of MrBayes posterior parameter diagnostics."""
    return write_taxon_rows(
        path,
        columns=[
            "parameter",
            "sample_count",
            "effective_sample_size",
            "mean",
            "median",
            "standard_deviation",
            "minimum",
            "maximum",
            "hpd_95_lower",
            "hpd_95_upper",
            "first_half_mean",
            "second_half_mean",
            "standardized_mean_shift",
            "burnin_fraction",
            "burnin_row_count",
            "kept_row_count",
            "first_kept_generation",
            "last_kept_generation",
        ],
        rows=[
            {
                "parameter": summary.parameter,
                "sample_count": str(summary.sample_count),
                "effective_sample_size": format(summary.effective_sample_size, ".15g"),
                "mean": format(summary.mean, ".15g"),
                "median": format(summary.median, ".15g"),
                "standard_deviation": format(summary.standard_deviation, ".15g"),
                "minimum": format(summary.minimum, ".15g"),
                "maximum": format(summary.maximum, ".15g"),
                "hpd_95_lower": format(summary.hpd_95_lower, ".15g"),
                "hpd_95_upper": format(summary.hpd_95_upper, ".15g"),
                "first_half_mean": format(summary.first_half_mean, ".15g"),
                "second_half_mean": format(summary.second_half_mean, ".15g"),
                "standardized_mean_shift": format(
                    summary.standardized_mean_shift, ".15g"
                ),
                "burnin_fraction": format(report.burnin_fraction, ".15g"),
                "burnin_row_count": str(report.burnin_row_count),
                "kept_row_count": str(report.kept_row_count),
                "first_kept_generation": str(report.first_kept_generation),
                "last_kept_generation": str(report.last_kept_generation),
            }
            for summary in report.parameter_summaries
        ],
    )


def assess_mrbayes_convergence(
    path: Path,
    *,
    ess_threshold: float = 200.0,
    mean_shift_threshold: float = 0.5,
) -> MrBayesConvergenceReport:
    """Flag low-ESS or unstable MrBayes trace parameters."""
    report = parse_mrbayes_parameter_traces(path)
    convergence = summarize_trace_convergence(
        path=path,
        rows=[row.values for row in report.rows],
        columns=report.columns,
        ess_threshold=ess_threshold,
        mean_shift_threshold=mean_shift_threshold,
    )
    return _build_mrbayes_convergence_report(convergence)


def _build_mrbayes_convergence_report(
    convergence: TraceConvergenceReport,
) -> MrBayesConvergenceReport:
    return MrBayesConvergenceReport(
        path=convergence.path,
        sample_count=convergence.sample_count,
        converged=convergence.converged,
        ess_threshold=convergence.ess_threshold,
        mean_shift_threshold=convergence.mean_shift_threshold,
        warnings=[
            {
                "parameter": warning.parameter,
                "code": warning.code,
                "message": warning.message,
                "observed_value": warning.observed_value,
                "threshold": warning.threshold,
            }
            for warning in convergence.warnings
        ],
        parameter_summaries=[
            {
                "parameter": summary.parameter,
                "sample_count": summary.sample_count,
                "effective_sample_size": summary.effective_sample_size,
                "mean": summary.mean,
                "median": summary.median,
                "standard_deviation": summary.standard_deviation,
                "minimum": summary.minimum,
                "maximum": summary.maximum,
                "hpd_95_lower": summary.hpd_95_lower,
                "hpd_95_upper": summary.hpd_95_upper,
                "first_half_mean": summary.first_half_mean,
                "second_half_mean": summary.second_half_mean,
                "standardized_mean_shift": summary.standardized_mean_shift,
            }
            for summary in convergence.series
        ],
    )


def _split_mrbayes_trace_rows(
    report: MrBayesTraceReport,
    *,
    burnin_fraction: float,
) -> tuple[int, list[MrBayesTraceRow]]:
    if not 0.0 <= burnin_fraction < 1.0:
        raise ValueError(
            f"burnin_fraction must be between 0 and 1, got {burnin_fraction}"
        )
    burnin_row_count = int(report.row_count * burnin_fraction)
    kept_rows = report.rows[burnin_row_count:]
    if not kept_rows:
        raise ValueError(
            "burnin_fraction discards every MrBayes trace row; reduce the burn-in"
        )
    return burnin_row_count, kept_rows


def assess_mrbayes_burnin_sensitivity(
    posterior_tree_path: Path,
    *,
    trace_path: Path | None = None,
    burnin_fractions: tuple[float, ...] = DEFAULT_BURNIN_FRACTIONS,
) -> MrBayesBurninSensitivityReport:
    """Compare MrBayes posterior summaries across multiple burn-in fractions."""
    ordered_fractions = normalize_burnin_fractions(burnin_fractions)
    slices: list[MrBayesBurninSensitivitySlice] = []
    previous_consensus: str | None = None
    changed_consensus_count = 0
    parameter_summaries_by_fraction: dict[float, list[MrBayesParameterSummary]] = {}
    clade_frequencies_by_fraction: dict[float, list[object]] = {}
    for fraction in ordered_fractions:
        _, posterior_summary = summarize_mrbayes_posterior_trees(
            posterior_tree_path,
            burnin_fraction=fraction,
        )
        clade_report = compute_clade_frequency_table(
            posterior_summary.filtered_tree_set_path
        )
        kept_row_count = None
        first_kept_generation = None
        last_kept_generation = None
        lnl_mean = None
        tree_length_mean = None
        if trace_path is not None:
            trace_summary = summarize_mrbayes_parameter_diagnostics(
                trace_path,
                burnin_fraction=fraction,
            )
            parameter_summaries_by_fraction[fraction] = (
                trace_summary.parameter_summaries
            )
            kept_row_count = trace_summary.kept_row_count
            first_kept_generation = trace_summary.first_kept_generation
            last_kept_generation = trace_summary.last_kept_generation
            lnl_mean = _mean_mrbayes_parameter(trace_summary, "LnL")
            tree_length_mean = _mean_mrbayes_parameter(trace_summary, "TL")
            if tree_length_mean is None:
                tree_length_mean = _mean_mrbayes_parameter(trace_summary, "TL{all}")
        clade_frequencies_by_fraction[fraction] = list(clade_report.clade_frequencies)
        slices.append(
            MrBayesBurninSensitivitySlice(
                burnin_fraction=fraction,
                burnin_tree_count=posterior_summary.burnin_tree_count,
                kept_tree_count=posterior_summary.kept_tree_count,
                rooted_topology_count=posterior_summary.rooted_topology_count,
                clade_frequency_count=posterior_summary.clade_frequency_count,
                consensus_newick=posterior_summary.consensus_newick,
                kept_row_count=kept_row_count,
                first_kept_generation=first_kept_generation,
                last_kept_generation=last_kept_generation,
                lnl_mean=lnl_mean,
                tree_length_mean=tree_length_mean,
            )
        )
        if (
            previous_consensus is not None
            and previous_consensus != posterior_summary.consensus_newick
        ):
            changed_consensus_count += 1
        previous_consensus = posterior_summary.consensus_newick
    parameter_shifts = summarize_burnin_parameter_shifts(
        parameter_summaries_by_fraction
    )
    clade_shifts = summarize_burnin_clade_shifts(clade_frequencies_by_fraction)
    warnings: list[str] = []
    if changed_consensus_count:
        warnings.append(
            "majority-rule consensus topology changes across tested burn-in fractions"
        )
    if any(shift.unstable for shift in parameter_shifts):
        warnings.append(
            "one or more posterior parameter 95% HPD intervals do not overlap across tested burn-in fractions"
        )
    if any(shift.unstable for shift in clade_shifts):
        warnings.append(
            "one or more posterior clade probabilities cross the majority-rule threshold across tested burn-in fractions"
        )
    return MrBayesBurninSensitivityReport(
        posterior_tree_path=posterior_tree_path,
        trace_path=trace_path,
        slices=slices,
        changed_consensus_count=changed_consensus_count,
        parameter_shifts=parameter_shifts,
        clade_shifts=clade_shifts,
        unstable_parameter_count=sum(1 for shift in parameter_shifts if shift.unstable),
        unstable_clade_count=sum(1 for shift in clade_shifts if shift.unstable),
        warnings=warnings,
    )


def write_mrbayes_burnin_sensitivity_slice_table(
    path: Path,
    report: MrBayesBurninSensitivityReport,
) -> Path:
    """Write one row per tested MrBayes burn-in fraction."""
    return write_taxon_rows(
        path,
        columns=[
            "burnin_fraction",
            "burnin_tree_count",
            "kept_tree_count",
            "rooted_topology_count",
            "clade_frequency_count",
            "kept_row_count",
            "first_kept_generation",
            "last_kept_generation",
            "lnl_mean",
            "tree_length_mean",
            "consensus_newick",
        ],
        rows=[
            {
                "burnin_fraction": format(row.burnin_fraction, ".15g"),
                "burnin_tree_count": str(row.burnin_tree_count),
                "kept_tree_count": str(row.kept_tree_count),
                "rooted_topology_count": str(row.rooted_topology_count),
                "clade_frequency_count": str(row.clade_frequency_count),
                "kept_row_count": ""
                if row.kept_row_count is None
                else str(row.kept_row_count),
                "first_kept_generation": ""
                if row.first_kept_generation is None
                else str(row.first_kept_generation),
                "last_kept_generation": ""
                if row.last_kept_generation is None
                else str(row.last_kept_generation),
                "lnl_mean": ""
                if row.lnl_mean is None
                else format(row.lnl_mean, ".15g"),
                "tree_length_mean": ""
                if row.tree_length_mean is None
                else format(row.tree_length_mean, ".15g"),
                "consensus_newick": row.consensus_newick,
            }
            for row in report.slices
        ],
    )


def _mean_mrbayes_parameter(
    report: MrBayesParameterDiagnosticsReport,
    parameter: str,
) -> float | None:
    for summary in report.parameter_summaries:
        if summary.parameter == parameter:
            return summary.mean
    return None
