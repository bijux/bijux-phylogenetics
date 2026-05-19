from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.ancestral.common import write_ancestral_rows

from .models import (
    DiscreteAncestralExclusion,
    DiscreteAncestralReport,
    DiscreteAncestralSummary,
)


def summarize_discrete_ancestral_report(
    report: DiscreteAncestralReport,
) -> DiscreteAncestralSummary:
    """Summarize the main review facts for one discrete ancestral report."""
    internal_estimates = [
        estimate for estimate in report.estimates if not estimate.is_tip
    ]
    if not internal_estimates:
        raise ValueError(
            "discrete ancestral summary requires at least one internal-node estimate"
        )
    ambiguous_internal_node_count = sum(
        1 for estimate in internal_estimates if estimate.ambiguous
    )
    root_estimate = max(
        internal_estimates,
        key=lambda estimate: (
            len(estimate.descendant_taxa),
            estimate.node,
        ),
    )
    return DiscreteAncestralSummary(
        trait=report.trait,
        taxon_column=report.taxon_column,
        model=report.model,
        state_ordering=report.state_ordering,
        root_prior_mode=report.root_prior_mode,
        fixed_root_state=report.fixed_root_state,
        analyzed_taxon_count=report.taxon_count,
        excluded_taxon_count=len(report.dropped_missing_taxa),
        internal_node_count=len(internal_estimates),
        ambiguous_internal_node_count=ambiguous_internal_node_count,
        unstable_node_count=len(report.unstable_nodes),
        weak_support_node_count=len(report.weak_support_nodes),
        observed_state_count=len(report.observed_states),
        sparse_state_count=len(report.sparse_states),
        minimal_change_count=report.minimal_change_count,
        parsimonious_root_state_count=report.parsimonious_root_state_count,
        root_node=root_estimate.node,
        root_most_likely_state=root_estimate.most_likely_state,
        root_confidence=root_estimate.confidence,
        phytools_rerooting_method_comparable=(
            report.rerooting_method_compatibility.comparable
        ),
        log_likelihood=report.log_likelihood,
        parameter_count=report.parameter_count,
        aic=report.aic,
        optimizer_converged=(
            None
            if report.optimizer_diagnostics is None
            else report.optimizer_diagnostics.converged
        ),
        optimizer_iteration_count=(
            None
            if report.optimizer_diagnostics is None
            else report.optimizer_diagnostics.iteration_count
        ),
        optimizer_function_evaluation_count=(
            None
            if report.optimizer_diagnostics is None
            else report.optimizer_diagnostics.function_evaluation_count
        ),
        overparameterized=report.overparameterized,
        baseline_model=(
            None
            if report.baseline_comparison is None
            else report.baseline_comparison.baseline_model
        ),
        baseline_delta_aic=(
            None
            if report.baseline_comparison is None
            else report.baseline_comparison.delta_aic
        ),
        preferred_model_by_aic=(
            None
            if report.baseline_comparison is None
            else report.baseline_comparison.preferred_model_by_aic
        ),
        warning_count=len(report.warnings),
    )


def discrete_ancestral_exclusions(
    report: DiscreteAncestralReport,
) -> list[DiscreteAncestralExclusion]:
    """Return one explicit exclusion row per dropped tip taxon."""
    return [
        DiscreteAncestralExclusion(
            taxon=taxon,
            reason="missing_discrete_trait_state",
        )
        for taxon in report.dropped_missing_taxa
    ]


def write_discrete_ancestral_summary_table(
    path: Path,
    report: DiscreteAncestralReport,
) -> Path:
    """Write one summary ledger for a discrete ancestral reconstruction."""
    summary = summarize_discrete_ancestral_report(report)
    return write_ancestral_rows(
        path,
        columns=[
            "trait",
            "taxon_column",
            "model",
            "state_ordering",
            "root_prior_mode",
            "fixed_root_state",
            "analyzed_taxon_count",
            "excluded_taxon_count",
            "internal_node_count",
            "ambiguous_internal_node_count",
            "unstable_node_count",
            "weak_support_node_count",
            "observed_state_count",
            "sparse_state_count",
            "minimal_change_count",
            "parsimonious_root_state_count",
            "root_node",
            "root_most_likely_state",
            "root_confidence",
            "phytools_rerooting_method_comparable",
            "log_likelihood",
            "parameter_count",
            "aic",
            "optimizer_converged",
            "optimizer_iteration_count",
            "optimizer_function_evaluation_count",
            "overparameterized",
            "baseline_model",
            "baseline_delta_aic",
            "preferred_model_by_aic",
            "warning_count",
        ],
        rows=[
            {
                "trait": summary.trait,
                "taxon_column": summary.taxon_column,
                "model": summary.model,
                "state_ordering": summary.state_ordering,
                "root_prior_mode": summary.root_prior_mode or "",
                "fixed_root_state": summary.fixed_root_state or "",
                "analyzed_taxon_count": str(summary.analyzed_taxon_count),
                "excluded_taxon_count": str(summary.excluded_taxon_count),
                "internal_node_count": str(summary.internal_node_count),
                "ambiguous_internal_node_count": str(
                    summary.ambiguous_internal_node_count
                ),
                "unstable_node_count": str(summary.unstable_node_count),
                "weak_support_node_count": str(summary.weak_support_node_count),
                "observed_state_count": str(summary.observed_state_count),
                "sparse_state_count": str(summary.sparse_state_count),
                "minimal_change_count": _format_optional_int(
                    summary.minimal_change_count
                ),
                "parsimonious_root_state_count": _format_optional_int(
                    summary.parsimonious_root_state_count
                ),
                "root_node": summary.root_node,
                "root_most_likely_state": summary.root_most_likely_state,
                "root_confidence": str(summary.root_confidence),
                "phytools_rerooting_method_comparable": str(
                    summary.phytools_rerooting_method_comparable
                ).lower(),
                "log_likelihood": _format_optional_float(summary.log_likelihood),
                "parameter_count": _format_optional_int(summary.parameter_count),
                "aic": _format_optional_float(summary.aic),
                "optimizer_converged": _format_optional_bool(
                    summary.optimizer_converged
                ),
                "optimizer_iteration_count": _format_optional_int(
                    summary.optimizer_iteration_count
                ),
                "optimizer_function_evaluation_count": _format_optional_int(
                    summary.optimizer_function_evaluation_count
                ),
                "overparameterized": str(summary.overparameterized).lower(),
                "baseline_model": summary.baseline_model or "",
                "baseline_delta_aic": _format_optional_float(
                    summary.baseline_delta_aic
                ),
                "preferred_model_by_aic": summary.preferred_model_by_aic or "",
                "warning_count": str(summary.warning_count),
            }
        ],
    )


def write_discrete_ancestral_probability_table(
    path: Path,
    report: DiscreteAncestralReport,
) -> Path:
    """Write one internal-node marginal-probability ledger for a discrete reconstruction."""
    return write_ancestral_rows(
        path,
        columns=[
            "node",
            "node_name",
            "descendant_taxa",
            "most_likely_state",
            "state_set",
            "state_probabilities",
            "confidence",
            "ambiguous",
            "unstable",
            "interpretation",
        ],
        rows=[
            {
                "node": estimate.node,
                "node_name": estimate.node_name or "",
                "descendant_taxa": ",".join(estimate.descendant_taxa),
                "most_likely_state": estimate.most_likely_state,
                "state_set": ",".join(estimate.state_set),
                "state_probabilities": json.dumps(
                    estimate.state_probabilities,
                    sort_keys=True,
                ),
                "confidence": str(estimate.confidence),
                "ambiguous": str(estimate.ambiguous).lower(),
                "unstable": str(estimate.unstable).lower(),
                "interpretation": estimate.interpretation,
            }
            for estimate in report.estimates
            if not estimate.is_tip
        ],
    )


def write_discrete_ancestral_exclusion_table(
    path: Path,
    report: DiscreteAncestralReport,
) -> Path:
    """Write one explicit excluded-tip ledger for a discrete reconstruction."""
    exclusions = discrete_ancestral_exclusions(report)
    return write_ancestral_rows(
        path,
        columns=["taxon", "reason"],
        rows=[
            {
                "taxon": row.taxon,
                "reason": row.reason,
            }
            for row in exclusions
        ],
    )


def write_discrete_ancestral_transition_table(
    path: Path,
    report: DiscreteAncestralReport,
) -> Path:
    """Write one fitted transition-rate ledger for a discrete likelihood reconstruction."""
    return write_ancestral_rows(
        path,
        columns=[
            "source_state",
            "target_state",
            "transition_allowed",
            "step_distance",
            "rate",
        ],
        rows=[
            {
                "source_state": row.source_state,
                "target_state": row.target_state,
                "transition_allowed": str(row.transition_allowed).lower(),
                "step_distance": str(row.step_distance),
                "rate": str(row.rate),
            }
            for row in report.transition_rate_rows
        ],
    )


def write_discrete_ancestral_fit_table(
    path: Path,
    report: DiscreteAncestralReport,
) -> Path:
    """Write one optimizer and model-fit ledger for a discrete likelihood reconstruction."""
    diagnostics = report.optimizer_diagnostics
    baseline = report.baseline_comparison
    return write_ancestral_rows(
        path,
        columns=[
            "model",
            "taxon_count",
            "state_count",
            "parameter_count",
            "log_likelihood",
            "aic",
            "overparameterized",
            "optimizer_name",
            "optimizer_converged",
            "optimizer_iteration_count",
            "optimizer_function_evaluation_count",
            "simplex_shrink_count",
            "initial_candidate_count",
            "best_initial_scale",
            "hit_lower_parameter_bound",
            "hit_upper_parameter_bound",
            "baseline_model",
            "baseline_parameter_count",
            "baseline_log_likelihood",
            "baseline_aic",
            "delta_log_likelihood",
            "delta_aic",
            "preferred_model_by_aic",
        ],
        rows=[
            {
                "model": report.model,
                "taxon_count": str(report.taxon_count),
                "state_count": str(len(report.observed_states)),
                "parameter_count": _format_optional_int(report.parameter_count),
                "log_likelihood": _format_optional_float(report.log_likelihood),
                "aic": _format_optional_float(report.aic),
                "overparameterized": str(report.overparameterized).lower(),
                "optimizer_name": ""
                if diagnostics is None
                else diagnostics.optimizer_name,
                "optimizer_converged": _format_optional_bool(
                    None if diagnostics is None else diagnostics.converged
                ),
                "optimizer_iteration_count": _format_optional_int(
                    None if diagnostics is None else diagnostics.iteration_count
                ),
                "optimizer_function_evaluation_count": _format_optional_int(
                    None
                    if diagnostics is None
                    else diagnostics.function_evaluation_count
                ),
                "simplex_shrink_count": _format_optional_int(
                    None if diagnostics is None else diagnostics.simplex_shrink_count
                ),
                "initial_candidate_count": _format_optional_int(
                    None if diagnostics is None else diagnostics.initial_candidate_count
                ),
                "best_initial_scale": _format_optional_float(
                    None if diagnostics is None else diagnostics.best_initial_scale
                ),
                "hit_lower_parameter_bound": _format_optional_bool(
                    None
                    if diagnostics is None
                    else diagnostics.hit_lower_parameter_bound
                ),
                "hit_upper_parameter_bound": _format_optional_bool(
                    None
                    if diagnostics is None
                    else diagnostics.hit_upper_parameter_bound
                ),
                "baseline_model": "" if baseline is None else baseline.baseline_model,
                "baseline_parameter_count": _format_optional_int(
                    None if baseline is None else baseline.baseline_parameter_count
                ),
                "baseline_log_likelihood": _format_optional_float(
                    None if baseline is None else baseline.baseline_log_likelihood
                ),
                "baseline_aic": _format_optional_float(
                    None if baseline is None else baseline.baseline_aic
                ),
                "delta_log_likelihood": _format_optional_float(
                    None if baseline is None else baseline.delta_log_likelihood
                ),
                "delta_aic": _format_optional_float(
                    None if baseline is None else baseline.delta_aic
                ),
                "preferred_model_by_aic": (
                    "" if baseline is None else baseline.preferred_model_by_aic
                ),
            }
        ],
    )


def _format_optional_int(value: int | None) -> str:
    if value is None:
        return ""
    return str(value)


def _format_optional_float(value: float | None) -> str:
    if value is None:
        return ""
    return str(value)


def _format_optional_bool(value: bool | None) -> str:
    if value is None:
        return ""
    return str(value).lower()
