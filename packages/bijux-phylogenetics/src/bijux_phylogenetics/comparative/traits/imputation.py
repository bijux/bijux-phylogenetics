from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path

from bijux_phylogenetics.comparative._math import (
    dot,
    invert_matrix,
    log_determinant,
    matrix_vector_multiply,
    quadratic_form,
    stable_covariance,
)
from bijux_phylogenetics.comparative.common import (
    ComparativeReadinessReport,
    build_brownian_covariance_matrix,
    load_comparative_dataset,
    summarize_numeric_trait_readiness,
)
from bijux_phylogenetics.comparative.continuous.brownian import (
    BrownianTraitEvolutionSummaryReport,
    summarize_brownian_trait_evolution,
)
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.phylo.pruning import prune_tree_to_requested_taxa

_Z_95 = 1.959963984540054


@dataclass(slots=True)
class TraitImputationExclusion:
    """One taxon excluded from Brownian trait imputation."""

    taxon: str
    reason: str


@dataclass(slots=True)
class TraitImputationRow:
    """One missing trait value imputed from observed taxa under Brownian motion."""

    taxon: str
    missing_reason: str
    observed_support_taxon_count: int
    predicted_value: float
    conditional_variance: float
    conditional_standard_error: float
    lower_95_confidence_interval: float
    upper_95_confidence_interval: float


@dataclass(slots=True)
class TraitImputationHoldoutRow:
    """One leave-one-observed-out Brownian validation row."""

    taxon: str
    observed_value: float
    predicted_value: float
    residual: float
    absolute_error: float
    conditional_variance: float
    conditional_standard_error: float
    lower_95_confidence_interval: float
    upper_95_confidence_interval: float
    covered_by_95_confidence_interval: bool
    observed_support_taxon_count: int
    rank: int


@dataclass(slots=True)
class TraitImputationSummaryReport:
    """Reviewer-facing Brownian trait imputation report with holdout evidence."""

    tree_path: Path
    traits_path: Path
    taxon_column: str
    trait: str
    model: str
    tree_taxon_count: int
    observed_taxa: list[str]
    observed_taxon_count: int
    imputation_rows: list[TraitImputationRow]
    excluded_taxa: list[TraitImputationExclusion]
    root_state: float
    sigma_squared: float
    log_likelihood: float
    aic: float
    aicc: float
    holdout_validation_status: str
    holdout_rows: list[TraitImputationHoldoutRow]
    holdout_mean_absolute_error: float | None
    holdout_root_mean_squared_error: float | None
    holdout_interval_coverage: float | None
    assumptions: list[str]
    warnings: list[str]
    readiness: ComparativeReadinessReport
    brownian_fit: BrownianTraitEvolutionSummaryReport


def summarize_trait_imputation(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
) -> TraitImputationSummaryReport:
    """Predict missing continuous traits from observed taxa under Brownian motion."""
    readiness = summarize_numeric_trait_readiness(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    observed_dataset = load_comparative_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        minimum_taxa=3,
        require_rooted=True,
        require_binary=False,
    )
    brownian_fit = summarize_brownian_trait_evolution(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    imputation_candidates = _imputation_candidates(readiness, observed_dataset.taxa)
    imputation_rows = _build_imputation_rows(
        tree_path,
        observed_taxa=observed_dataset.taxa,
        observed_values=observed_dataset.trait_values,
        target_rows=imputation_candidates,
        root_state=brownian_fit.root_state,
        sigma_squared=brownian_fit.sigma_squared,
    )
    holdout_rows, holdout_validation_status = _build_holdout_rows(
        tree_path,
        observed_taxa=observed_dataset.taxa,
        observed_values=observed_dataset.trait_values,
    )
    holdout_mean_absolute_error = (
        sum(row.absolute_error for row in holdout_rows) / len(holdout_rows)
        if holdout_rows
        else None
    )
    holdout_root_mean_squared_error = (
        math.sqrt(
            sum(row.residual * row.residual for row in holdout_rows) / len(holdout_rows)
        )
        if holdout_rows
        else None
    )
    holdout_interval_coverage = (
        sum(1 for row in holdout_rows if row.covered_by_95_confidence_interval)
        / len(holdout_rows)
        if holdout_rows
        else None
    )
    warnings = _warnings(
        readiness,
        imputation_rows=imputation_rows,
        holdout_validation_status=holdout_validation_status,
    )
    return TraitImputationSummaryReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=observed_dataset.taxon_column,
        trait=trait,
        model="brownian",
        tree_taxon_count=readiness.tree_taxa,
        observed_taxa=list(observed_dataset.taxa),
        observed_taxon_count=len(observed_dataset.taxa),
        imputation_rows=imputation_rows,
        excluded_taxa=_excluded_taxa(readiness),
        root_state=brownian_fit.root_state,
        sigma_squared=brownian_fit.sigma_squared,
        log_likelihood=brownian_fit.log_likelihood,
        aic=brownian_fit.aic,
        aicc=brownian_fit.aicc,
        holdout_validation_status=holdout_validation_status,
        holdout_rows=holdout_rows,
        holdout_mean_absolute_error=holdout_mean_absolute_error,
        holdout_root_mean_squared_error=holdout_root_mean_squared_error,
        holdout_interval_coverage=holdout_interval_coverage,
        assumptions=[
            *brownian_fit.assumptions,
            "Missing values are imputed from the Brownian conditional distribution given all observed taxa retained for the trait",
            "Holdout validation removes one observed taxon at a time, refits the Brownian mean and diffusion rate on the remaining observed taxa, and predicts the held-out tip",
        ],
        warnings=warnings,
        readiness=readiness,
        brownian_fit=brownian_fit,
    )


def write_trait_imputation_summary_table(
    path: Path,
    report: TraitImputationSummaryReport,
) -> Path:
    """Write one summary ledger for Brownian trait imputation."""
    return write_taxon_rows(
        path,
        columns=[
            "trait",
            "taxon_column",
            "model",
            "tree_taxon_count",
            "observed_taxon_count",
            "imputed_taxon_count",
            "excluded_taxon_count",
            "root_state",
            "sigma_squared",
            "log_likelihood",
            "aic",
            "aicc",
            "holdout_validation_status",
            "holdout_count",
            "holdout_mean_absolute_error",
            "holdout_root_mean_squared_error",
            "holdout_interval_coverage",
        ],
        rows=[
            {
                "trait": report.trait,
                "taxon_column": report.taxon_column,
                "model": report.model,
                "tree_taxon_count": str(report.tree_taxon_count),
                "observed_taxon_count": str(report.observed_taxon_count),
                "imputed_taxon_count": str(len(report.imputation_rows)),
                "excluded_taxon_count": str(len(report.excluded_taxa)),
                "root_state": format(report.root_state, ".15g"),
                "sigma_squared": format(report.sigma_squared, ".15g"),
                "log_likelihood": format(report.log_likelihood, ".15g"),
                "aic": format(report.aic, ".15g"),
                "aicc": format(report.aicc, ".15g"),
                "holdout_validation_status": report.holdout_validation_status,
                "holdout_count": str(len(report.holdout_rows)),
                "holdout_mean_absolute_error": _format_optional(
                    report.holdout_mean_absolute_error
                ),
                "holdout_root_mean_squared_error": _format_optional(
                    report.holdout_root_mean_squared_error
                ),
                "holdout_interval_coverage": _format_optional(
                    report.holdout_interval_coverage
                ),
            }
        ],
    )


def write_trait_imputation_table(
    path: Path,
    report: TraitImputationSummaryReport,
) -> Path:
    """Write one imputed-value ledger for Brownian trait imputation."""
    return write_taxon_rows(
        path,
        columns=[
            "taxon",
            "missing_reason",
            "observed_support_taxon_count",
            "predicted_value",
            "conditional_variance",
            "conditional_standard_error",
            "lower_95_confidence_interval",
            "upper_95_confidence_interval",
        ],
        rows=[
            {
                "taxon": row.taxon,
                "missing_reason": row.missing_reason,
                "observed_support_taxon_count": str(row.observed_support_taxon_count),
                "predicted_value": format(row.predicted_value, ".15g"),
                "conditional_variance": format(row.conditional_variance, ".15g"),
                "conditional_standard_error": format(
                    row.conditional_standard_error, ".15g"
                ),
                "lower_95_confidence_interval": format(
                    row.lower_95_confidence_interval, ".15g"
                ),
                "upper_95_confidence_interval": format(
                    row.upper_95_confidence_interval, ".15g"
                ),
            }
            for row in report.imputation_rows
        ],
    )


def write_trait_imputation_holdout_table(
    path: Path,
    report: TraitImputationSummaryReport,
) -> Path:
    """Write one leave-one-observed-out validation ledger for Brownian imputation."""
    return write_taxon_rows(
        path,
        columns=[
            "taxon",
            "observed_value",
            "predicted_value",
            "residual",
            "absolute_error",
            "conditional_variance",
            "conditional_standard_error",
            "lower_95_confidence_interval",
            "upper_95_confidence_interval",
            "covered_by_95_confidence_interval",
            "observed_support_taxon_count",
            "rank",
        ],
        rows=[
            {
                "taxon": row.taxon,
                "observed_value": format(row.observed_value, ".15g"),
                "predicted_value": format(row.predicted_value, ".15g"),
                "residual": format(row.residual, ".15g"),
                "absolute_error": format(row.absolute_error, ".15g"),
                "conditional_variance": format(row.conditional_variance, ".15g"),
                "conditional_standard_error": format(
                    row.conditional_standard_error, ".15g"
                ),
                "lower_95_confidence_interval": format(
                    row.lower_95_confidence_interval, ".15g"
                ),
                "upper_95_confidence_interval": format(
                    row.upper_95_confidence_interval, ".15g"
                ),
                "covered_by_95_confidence_interval": str(
                    row.covered_by_95_confidence_interval
                ).lower(),
                "observed_support_taxon_count": str(row.observed_support_taxon_count),
                "rank": str(row.rank),
            }
            for row in report.holdout_rows
        ],
    )


def write_trait_imputation_exclusion_table(
    path: Path,
    report: TraitImputationSummaryReport,
) -> Path:
    """Write one excluded-taxa ledger for Brownian trait imputation."""
    return write_taxon_rows(
        path,
        columns=["taxon", "reason"],
        rows=[
            {
                "taxon": row.taxon,
                "reason": row.reason,
            }
            for row in report.excluded_taxa
        ],
    )


@dataclass(slots=True)
class _ImputationTarget:
    taxon: str
    missing_reason: str


def _imputation_candidates(
    readiness: ComparativeReadinessReport,
    observed_taxa: list[str],
) -> list[_ImputationTarget]:
    observed_set = set(observed_taxa)
    rows: list[_ImputationTarget] = []
    rows.extend(
        _ImputationTarget(taxon=taxon, missing_reason="missing_from_trait_table")
        for taxon in readiness.missing_from_traits
        if taxon not in observed_set
    )
    rows.extend(
        _ImputationTarget(taxon=taxon, missing_reason="missing_trait_value")
        for taxon in readiness.pruned_missing_value_taxa
        if taxon not in observed_set
    )
    return rows


def _build_imputation_rows(
    tree_path: Path,
    *,
    observed_taxa: list[str],
    observed_values: list[float],
    target_rows: list[_ImputationTarget],
    root_state: float,
    sigma_squared: float,
) -> list[TraitImputationRow]:
    if not target_rows:
        return []
    target_lookup = {row.taxon: row for row in target_rows}
    requested_taxa = [*observed_taxa, *[row.taxon for row in target_rows]]
    pruned_tree, _ = prune_tree_to_requested_taxa(tree_path, requested_taxa)
    model_taxa = pruned_tree.tip_names
    covariance = stable_covariance(
        [
            [value * sigma_squared for value in row]
            for row in build_brownian_covariance_matrix(pruned_tree, model_taxa)
        ]
    )
    observed_indices = [model_taxa.index(taxon) for taxon in observed_taxa]
    observed_covariance = stable_covariance(
        [
            [covariance[row][column] for column in observed_indices]
            for row in observed_indices
        ]
    )
    inverse_observed_covariance = invert_matrix(observed_covariance)
    centered_observed = [value - root_state for value in observed_values]
    weighted_observed = matrix_vector_multiply(
        inverse_observed_covariance,
        centered_observed,
    )
    rows: list[TraitImputationRow] = []
    for taxon in model_taxa:
        target = target_lookup.get(taxon)
        if target is None:
            continue
        target_index = model_taxa.index(taxon)
        cross_covariance = [
            covariance[target_index][index] for index in observed_indices
        ]
        predicted_value = root_state + dot(cross_covariance, weighted_observed)
        weighted_cross = matrix_vector_multiply(
            inverse_observed_covariance,
            cross_covariance,
        )
        conditional_variance = max(
            covariance[target_index][target_index]
            - dot(cross_covariance, weighted_cross),
            0.0,
        )
        conditional_standard_error = math.sqrt(max(conditional_variance, 1e-12))
        rows.append(
            TraitImputationRow(
                taxon=taxon,
                missing_reason=target.missing_reason,
                observed_support_taxon_count=len(observed_taxa),
                predicted_value=predicted_value,
                conditional_variance=conditional_variance,
                conditional_standard_error=conditional_standard_error,
                lower_95_confidence_interval=predicted_value
                - (_Z_95 * conditional_standard_error),
                upper_95_confidence_interval=predicted_value
                + (_Z_95 * conditional_standard_error),
            )
        )
    return rows


def _build_holdout_rows(
    tree_path: Path,
    *,
    observed_taxa: list[str],
    observed_values: list[float],
) -> tuple[list[TraitImputationHoldoutRow], str]:
    if len(observed_taxa) < 4:
        return [], "insufficient_observed_taxa"
    observed_lookup = dict(zip(observed_taxa, observed_values, strict=True))
    rows: list[TraitImputationHoldoutRow] = []
    for taxon in observed_taxa:
        support_taxa = [candidate for candidate in observed_taxa if candidate != taxon]
        support_values = [observed_lookup[candidate] for candidate in support_taxa]
        pruned_tree, _ = prune_tree_to_requested_taxa(tree_path, [*support_taxa, taxon])
        support_covariance = stable_covariance(
            build_brownian_covariance_matrix(pruned_tree, support_taxa)
        )
        fit = _fit_brownian_parameters(
            support_values,
            support_covariance,
        )
        prediction = _build_imputation_rows(
            tree_path,
            observed_taxa=support_taxa,
            observed_values=support_values,
            target_rows=[_ImputationTarget(taxon=taxon, missing_reason="holdout")],
            root_state=fit.root_state,
            sigma_squared=fit.sigma_squared,
        )[0]
        observed_value = observed_lookup[taxon]
        residual = observed_value - prediction.predicted_value
        rows.append(
            TraitImputationHoldoutRow(
                taxon=taxon,
                observed_value=observed_value,
                predicted_value=prediction.predicted_value,
                residual=residual,
                absolute_error=abs(residual),
                conditional_variance=prediction.conditional_variance,
                conditional_standard_error=prediction.conditional_standard_error,
                lower_95_confidence_interval=prediction.lower_95_confidence_interval,
                upper_95_confidence_interval=prediction.upper_95_confidence_interval,
                covered_by_95_confidence_interval=(
                    prediction.lower_95_confidence_interval
                    <= observed_value
                    <= prediction.upper_95_confidence_interval
                ),
                observed_support_taxon_count=len(support_taxa),
                rank=0,
            )
        )
    rows.sort(key=lambda row: (-row.absolute_error, row.taxon))
    for rank, row in enumerate(rows, start=1):
        row.rank = rank
    return rows, "performed"


@dataclass(slots=True)
class _BrownianParameterFit:
    root_state: float
    sigma_squared: float
    log_likelihood: float


def _fit_brownian_parameters(
    values: list[float],
    covariance: list[list[float]],
) -> _BrownianParameterFit:
    inverse_covariance = invert_matrix(covariance)
    ones = [1.0] * len(values)
    denominator = quadratic_form(ones, inverse_covariance)
    root_state = (
        sum(
            inverse_covariance[row_index][column_index] * values[column_index]
            for row_index in range(len(values))
            for column_index in range(len(values))
        )
        / denominator
    )
    residuals = [value - root_state for value in values]
    sigma_squared = quadratic_form(residuals, inverse_covariance) / len(values)
    log_likelihood = -0.5 * (
        len(values) * math.log(2.0 * math.pi * sigma_squared)
        + log_determinant(covariance)
        + len(values)
    )
    return _BrownianParameterFit(
        root_state=root_state,
        sigma_squared=sigma_squared,
        log_likelihood=log_likelihood,
    )


def _excluded_taxa(
    readiness: ComparativeReadinessReport,
) -> list[TraitImputationExclusion]:
    rows: list[TraitImputationExclusion] = []
    rows.extend(
        TraitImputationExclusion(taxon=taxon, reason="non_numeric_trait_value")
        for taxon in readiness.pruned_non_numeric_taxa
    )
    rows.extend(
        TraitImputationExclusion(taxon=taxon, reason="absent_from_tree")
        for taxon in readiness.extra_trait_taxa
    )
    return rows


def _warnings(
    readiness: ComparativeReadinessReport,
    *,
    imputation_rows: list[TraitImputationRow],
    holdout_validation_status: str,
) -> list[str]:
    warnings: list[str] = []
    if readiness.missing_from_traits:
        warnings.append(
            "trait table is missing one or more tree taxa and those taxa were imputed from Brownian phylogenetic context"
        )
    if readiness.pruned_missing_value_taxa:
        warnings.append(
            "one or more overlapping taxa have missing trait values and were imputed from Brownian phylogenetic context"
        )
    if readiness.pruned_non_numeric_taxa:
        warnings.append(
            "one or more overlapping taxa have non-numeric trait values and were excluded from Brownian imputation"
        )
    if readiness.extra_trait_taxa:
        warnings.append("trait table contains taxa absent from the tree")
    if not imputation_rows:
        warnings.append("no tree taxa required Brownian trait imputation")
    if holdout_validation_status != "performed":
        warnings.append(
            "holdout validation requires at least four observed taxa so every reduced Brownian fit retains at least three taxa"
        )
    return warnings


def _format_optional(value: float | None) -> str:
    if value is None:
        return ""
    return format(value, ".15g")
