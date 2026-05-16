from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path

from bijux_phylogenetics.ancestral.common import (
    AncestralDiscreteDataset,
    load_discrete_dataset,
    write_ancestral_rows,
)
from bijux_phylogenetics.ancestral.discrete import (
    DiscreteModelBaselineComparison,
    DiscreteOptimizerDiagnostics,
    DiscreteTransitionRateRow,
    _build_transition_rate_rows,
    _detect_discrete_overparameterization,
    _fit_discrete_mk_model,
    _parameter_count,
    _resolve_allowed_transition_pairs,
    _resolve_discrete_model_name,
    _resolve_state_order,
    _tree_log_likelihood,
)
from bijux_phylogenetics.comparative.common import tip_root_depths
from bijux_phylogenetics.core.ultrametric import summarize_ultrametric_tip_depths


@dataclass(slots=True)
class DiscreteMkInputAudit:
    """Owned input-policy audit for one discrete Mk model fit."""

    tree_path: Path
    traits_path: Path
    trait: str
    taxon_count: int
    taxa: list[str]
    observed_states: list[str]
    state_counts: dict[str, int]
    sparse_states: list[str]
    tree_is_ultrametric: bool
    minimum_root_to_tip_depth: float
    maximum_root_to_tip_depth: float
    ultrametric_policy: str
    missing_value_policy: str
    pruned_missing_value_taxa: list[str]
    warnings: list[str]


@dataclass(slots=True)
class DiscreteMkFitReport:
    """Discrete Mk trait-evolution fit over one rooted tree."""

    tree_path: Path
    traits_path: Path
    taxon_column: str
    trait: str
    model: str
    state_ordering: str
    state_order: list[str]
    taxon_count: int
    input_audit: DiscreteMkInputAudit
    log_likelihood: float
    parameter_count: int
    aic: float
    aicc: float
    transition_rate_rows: list[DiscreteTransitionRateRow]
    allowed_transition_pairs: list[tuple[str, str]]
    optimizer_diagnostics: DiscreteOptimizerDiagnostics
    overparameterized: bool
    baseline_comparison: DiscreteModelBaselineComparison | None


def _normalize_transition_rate_rows(
    rows: list[DiscreteTransitionRateRow],
    *,
    state_ordering: str,
) -> list[DiscreteTransitionRateRow]:
    if state_ordering != "unordered":
        return rows
    return [
        DiscreteTransitionRateRow(
            source_state=row.source_state,
            target_state=row.target_state,
            transition_allowed=row.transition_allowed,
            step_distance=(1 if row.transition_allowed else row.step_distance),
            rate=row.rate,
        )
        for row in rows
    ]


def fit_discrete_mk_model(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    model: str = "equal-rates",
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
    allowed_transition_pairs: list[tuple[str, str]] | None = None,
) -> DiscreteMkFitReport:
    """Fit one Mk discrete-trait model on a rooted tree."""
    dataset = load_discrete_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    return fit_discrete_mk_model_from_dataset(
        dataset,
        model=model,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
        allowed_transition_pairs=allowed_transition_pairs,
    )


def fit_discrete_mk_model_from_dataset(
    dataset: AncestralDiscreteDataset,
    *,
    model: str = "equal-rates",
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
    allowed_transition_pairs: list[tuple[str, str]] | None = None,
) -> DiscreteMkFitReport:
    """Fit one Mk discrete-trait model from a native discrete dataset."""
    resolved_model = _resolve_discrete_model_name(model)
    state_order = _resolve_state_order(
        dataset.observed_states,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
    )
    resolved_allowed_transition_pairs = _resolve_allowed_transition_pairs(
        state_order,
        model=resolved_model,
        state_ordering=state_ordering,
        allowed_transition_pairs=allowed_transition_pairs,
    )
    rate_matrix, root_prior, optimizer_diagnostics = _fit_discrete_mk_model(
        dataset.tree,
        dataset.states_by_taxon,
        state_order=state_order,
        model=resolved_model,
        state_ordering=state_ordering,
        allowed_transition_pairs=resolved_allowed_transition_pairs,
    )
    log_likelihood = _tree_log_likelihood(
        dataset.tree,
        dataset.states_by_taxon,
        state_order=state_order,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
    )
    parameter_count = _parameter_count(
        len(state_order),
        model=resolved_model,
        state_ordering=state_ordering,
        allowed_transition_pairs=resolved_allowed_transition_pairs,
    )
    aic = (2.0 * parameter_count) - (2.0 * log_likelihood)
    aicc = _fit_aicc(
        aic, sample_size=len(dataset.taxa), parameter_count=parameter_count
    )
    baseline_comparison: DiscreteModelBaselineComparison | None = None
    if (
        resolved_model != "equal-rates"
        and allowed_transition_pairs is None
        and state_ordering == "unordered"
    ):
        baseline_fit = fit_discrete_mk_model_from_dataset(
            dataset,
            model="equal-rates",
            state_ordering=state_ordering,
            ordered_states=state_order,
            allowed_transition_pairs=None,
        )
        baseline_comparison = DiscreteModelBaselineComparison(
            baseline_model="equal-rates",
            baseline_log_likelihood=baseline_fit.log_likelihood,
            baseline_parameter_count=baseline_fit.parameter_count,
            baseline_aic=baseline_fit.aic,
            delta_log_likelihood=log_likelihood - baseline_fit.log_likelihood,
            delta_aic=aic - baseline_fit.aic,
            preferred_model_by_aic=(
                resolved_model if aic <= baseline_fit.aic else "equal-rates"
            ),
        )
    overparameterized = _detect_discrete_overparameterization(
        taxon_count=len(dataset.taxa),
        parameter_count=parameter_count,
    )
    warnings = list(dataset.warnings)
    if overparameterized:
        warnings.append(
            "the discrete Mk likelihood fit is likely overparameterized relative to the analyzed taxon count"
        )
    if not optimizer_diagnostics.converged:
        warnings.append(
            "the discrete Mk optimizer did not converge and should be interpreted cautiously"
        )
    if (
        optimizer_diagnostics.hit_lower_parameter_bound
        or optimizer_diagnostics.hit_upper_parameter_bound
    ):
        warnings.append(
            "one or more discrete Mk rate parameters hit an optimizer bound and should be interpreted as weakly identified"
        )
    if (
        baseline_comparison is not None
        and baseline_comparison.preferred_model_by_aic == "equal-rates"
    ):
        warnings.append(
            "the equal-rates baseline remains preferred by AIC over the requested discrete Mk model"
        )
    input_audit = _build_discrete_mk_input_audit(dataset, warnings=warnings)
    return DiscreteMkFitReport(
        tree_path=dataset.tree_path,
        traits_path=dataset.traits_path,
        taxon_column=dataset.taxon_column,
        trait=dataset.trait,
        model=resolved_model,
        state_ordering=state_ordering,
        state_order=state_order,
        taxon_count=len(dataset.taxa),
        input_audit=input_audit,
        log_likelihood=log_likelihood,
        parameter_count=parameter_count,
        aic=aic,
        aicc=aicc,
        transition_rate_rows=_normalize_transition_rate_rows(
            _build_transition_rate_rows(
                state_order=state_order,
                state_ordering=state_ordering,
                rate_matrix=rate_matrix,
                allowed_transition_pairs=resolved_allowed_transition_pairs,
            ),
            state_ordering=state_ordering,
        ),
        allowed_transition_pairs=[
            (state_order[left_index], state_order[right_index])
            for left_index, right_index in sorted(resolved_allowed_transition_pairs)
        ],
        optimizer_diagnostics=optimizer_diagnostics,
        overparameterized=overparameterized,
        baseline_comparison=baseline_comparison,
    )


def write_discrete_mk_summary_table(path: Path, report: DiscreteMkFitReport) -> Path:
    """Write one flat summary ledger for a discrete Mk fit."""
    baseline = report.baseline_comparison
    diagnostics = report.optimizer_diagnostics
    return write_ancestral_rows(
        path,
        columns=[
            "trait",
            "taxon_column",
            "model",
            "state_ordering",
            "analyzed_taxon_count",
            "excluded_taxon_count",
            "observed_state_count",
            "sparse_state_count",
            "log_likelihood",
            "parameter_count",
            "aic",
            "aicc",
            "optimizer_name",
            "optimizer_converged",
            "optimizer_iteration_count",
            "optimizer_function_evaluation_count",
            "optimizer_simplex_shrink_count",
            "optimizer_initial_candidate_count",
            "optimizer_best_initial_scale",
            "optimizer_hit_lower_parameter_bound",
            "optimizer_hit_upper_parameter_bound",
            "overparameterized",
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
                "trait": report.trait,
                "taxon_column": report.taxon_column,
                "model": report.model,
                "state_ordering": report.state_ordering,
                "analyzed_taxon_count": str(report.taxon_count),
                "excluded_taxon_count": str(
                    len(report.input_audit.pruned_missing_value_taxa)
                ),
                "observed_state_count": str(len(report.input_audit.observed_states)),
                "sparse_state_count": str(len(report.input_audit.sparse_states)),
                "log_likelihood": format(report.log_likelihood, ".15g"),
                "parameter_count": str(report.parameter_count),
                "aic": format(report.aic, ".15g"),
                "aicc": "inf"
                if math.isinf(report.aicc)
                else format(report.aicc, ".15g"),
                "optimizer_name": diagnostics.optimizer_name,
                "optimizer_converged": str(diagnostics.converged).lower(),
                "optimizer_iteration_count": str(diagnostics.iteration_count),
                "optimizer_function_evaluation_count": str(
                    diagnostics.function_evaluation_count
                ),
                "optimizer_simplex_shrink_count": str(diagnostics.simplex_shrink_count),
                "optimizer_initial_candidate_count": str(
                    diagnostics.initial_candidate_count
                ),
                "optimizer_best_initial_scale": format(
                    diagnostics.best_initial_scale,
                    ".15g",
                ),
                "optimizer_hit_lower_parameter_bound": str(
                    diagnostics.hit_lower_parameter_bound
                ).lower(),
                "optimizer_hit_upper_parameter_bound": str(
                    diagnostics.hit_upper_parameter_bound
                ).lower(),
                "overparameterized": str(report.overparameterized).lower(),
                "baseline_model": "" if baseline is None else baseline.baseline_model,
                "baseline_parameter_count": ""
                if baseline is None
                else str(baseline.baseline_parameter_count),
                "baseline_log_likelihood": ""
                if baseline is None
                else format(baseline.baseline_log_likelihood, ".15g"),
                "baseline_aic": ""
                if baseline is None
                else format(baseline.baseline_aic, ".15g"),
                "delta_log_likelihood": ""
                if baseline is None
                else format(baseline.delta_log_likelihood, ".15g"),
                "delta_aic": ""
                if baseline is None
                else format(baseline.delta_aic, ".15g"),
                "preferred_model_by_aic": (
                    "" if baseline is None else baseline.preferred_model_by_aic
                ),
            }
        ],
    )


def write_discrete_mk_rate_table(path: Path, report: DiscreteMkFitReport) -> Path:
    """Write one directed rate-matrix ledger for a discrete Mk fit."""
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
                "rate": format(row.rate, ".15g"),
            }
            for row in report.transition_rate_rows
        ],
    )


def _build_discrete_mk_input_audit(
    dataset: AncestralDiscreteDataset,
    *,
    warnings: list[str],
) -> DiscreteMkInputAudit:
    ultrametric_summary = summarize_ultrametric_tip_depths(
        tip_root_depths(dataset.tree, dataset.taxa),
        tolerance=1e-12,
    )
    return DiscreteMkInputAudit(
        tree_path=dataset.tree_path,
        traits_path=dataset.traits_path,
        trait=dataset.trait,
        taxon_count=len(dataset.taxa),
        taxa=list(dataset.taxa),
        observed_states=list(dataset.observed_states),
        state_counts=dict(dataset.state_counts),
        sparse_states=list(dataset.sparse_states),
        tree_is_ultrametric=ultrametric_summary.ultrametric,
        minimum_root_to_tip_depth=ultrametric_summary.minimum_tip_depth,
        maximum_root_to_tip_depth=ultrametric_summary.maximum_tip_depth,
        ultrametric_policy="accept-rooted-trees-and-report-ultrametricity",
        missing_value_policy="prune-overlapping-missing-values",
        pruned_missing_value_taxa=list(dataset.dropped_missing_taxa),
        warnings=warnings,
    )


def _fit_aicc(aic: float, *, sample_size: int, parameter_count: int) -> float:
    denominator = sample_size - parameter_count - 1
    if denominator <= 0:
        return math.inf
    return aic + ((2.0 * parameter_count * (parameter_count + 1.0)) / denominator)
