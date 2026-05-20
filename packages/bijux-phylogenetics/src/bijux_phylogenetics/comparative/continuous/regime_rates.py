from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path

from bijux_phylogenetics.comparative._math import (
    invert_matrix,
    log_determinant,
    matrix_multiply,
    quadratic_form,
    stable_covariance,
    transpose,
)
from bijux_phylogenetics.comparative.common import (
    ComparativeDataset,
    ComparativeReadinessReport,
    load_comparative_dataset,
    summarize_numeric_trait_readiness,
)
from bijux_phylogenetics.comparative.models import (
    ComparativeModelComparisonRow,
    ComparativeParameterInterval,
    ComparativeResidualOutlier,
    ComparativeResidualSummary,
    _comparison_row,
    _estimate_lambda_for_values,
    fit_brownian_motion_model,
)
from bijux_phylogenetics.comparative.traits.regime_mapping import (
    build_branch_identity_lookup,
    resolve_branch_regime_id_column,
)
from bijux_phylogenetics.datasets.study_inputs import load_taxon_table, write_taxon_rows
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import ComparativeMethodError
from bijux_phylogenetics.io.trees import load_tree

_PROFILE_CONFIDENCE_DELTA = 1.920729410347062
_Z_95 = 1.959963984540054


@dataclass(slots=True)
class BrownianRegimeExclusion:
    """One taxon excluded before multi-rate Brownian fitting."""

    taxon: str
    reason: str


@dataclass(slots=True)
class BrownianRegimeBranchRow:
    """One non-root branch assigned to a regime for multi-rate Brownian fitting."""

    branch_id: str
    regime: str
    branch_length: float
    descendant_taxa: list[str]
    analyzed_descendant_taxa: list[str]
    contributes_to_analysis: bool


@dataclass(slots=True)
class BrownianRegimeRateRow:
    """One regime-specific Brownian rate estimate."""

    regime: str
    branch_count: int
    contributing_branch_count: int
    total_branch_length: float
    contributing_branch_length: float
    sigma_squared: float
    lower_95: float
    upper_95: float
    interval_method: str


@dataclass(slots=True)
class BrownianRegimeProfileRow:
    """One conditional likelihood-profile row for one regime rate."""

    regime: str
    sigma_squared: float
    log_likelihood: float
    delta_log_likelihood: float
    in_support_interval: bool
    selected: bool


@dataclass(slots=True)
class BrownianRegimeIdentifiabilityWarning:
    """Warning that one regime-specific Brownian rate is weakly identified."""

    regime: str
    kind: str
    message: str


@dataclass(slots=True)
class BrownianRegimeFitSummaryReport:
    """Reviewer-facing multi-rate Brownian fit driven by a branch regime map."""

    tree_path: Path
    traits_path: Path
    regime_map_path: Path
    taxon_column: str
    branch_id_column: str
    regime_column: str
    trait: str
    tree_taxon_count: int
    analyzed_taxa: list[str]
    analyzed_taxon_count: int
    excluded_taxa: list[BrownianRegimeExclusion]
    branch_rows: list[BrownianRegimeBranchRow]
    regime_rows: list[BrownianRegimeRateRow]
    profile_rows: list[BrownianRegimeProfileRow]
    root_state: float
    root_state_interval: ComparativeParameterInterval
    log_likelihood: float
    aic: float
    aicc: float
    comparison_rows: list[ComparativeModelComparisonRow]
    better_model: str
    likelihood_ratio_statistic: float
    likelihood_ratio_degrees_of_freedom: int
    likelihood_ratio_p_value: float
    likelihood_ratio_p_value_method: str
    identifiability_warnings: list[BrownianRegimeIdentifiabilityWarning]
    residual_diagnostics: ComparativeResidualSummary
    assumptions: list[str]
    warnings: list[str]
    readiness: ComparativeReadinessReport


def summarize_brownian_regime_rates(
    tree_path: Path,
    traits_path: Path,
    regime_map_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    branch_id_column: str | None = None,
    regime_column: str = "regime",
) -> BrownianRegimeFitSummaryReport:
    """Summarize a multi-rate Brownian model from a user-provided branch regime map."""
    readiness = summarize_numeric_trait_readiness(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    dataset = load_comparative_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        minimum_taxa=3,
        require_rooted=True,
        require_binary=False,
    )
    source_tree = load_tree(tree_path)
    branch_rows, resolved_branch_id_column = _load_branch_regime_rows(
        regime_map_path,
        tree=source_tree,
        analyzed_taxa=dataset.taxa,
        branch_id_column=branch_id_column,
        regime_column=regime_column,
    )
    regime_matrices = _build_regime_covariance_components(dataset.taxa, branch_rows)
    if len(regime_matrices) < 2:
        raise ComparativeMethodError(
            "multi-rate Brownian fitting requires at least two contributing regimes"
        )
    baseline = fit_brownian_motion_model(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    fit = _fit_multirate_brownian_model(
        dataset,
        regime_matrices=regime_matrices,
        baseline_rate=baseline.rate,
    )
    comparison_rows = _build_comparison_rows(
        taxon_count=dataset.taxon_count
        if hasattr(dataset, "taxon_count")
        else len(dataset.taxa),
        baseline_log_likelihood=baseline.log_likelihood,
        multirate_log_likelihood=fit.log_likelihood,
        regime_count=len(regime_matrices),
    )
    better_model = next(row.model for row in comparison_rows if row.selected)
    likelihood_ratio_degrees_of_freedom = len(regime_matrices) - 1
    likelihood_ratio_statistic = max(
        0.0,
        -2.0 * (baseline.log_likelihood - fit.log_likelihood),
    )
    likelihood_ratio_p_value = _chi_square_survival(
        likelihood_ratio_statistic,
        likelihood_ratio_degrees_of_freedom,
    )
    profile_rows = _build_profile_rows(
        dataset,
        regime_matrices=regime_matrices,
        fitted_rates=fit.regime_rates,
        baseline_rate=baseline.rate,
    )
    regime_rows = _build_regime_rows(branch_rows, fit.regime_rates, profile_rows)
    identifiability_warnings = _build_identifiability_warnings(
        regime_rows,
        profile_rows,
        better_model=better_model,
    )
    residual_diagnostics = _build_multirate_residual_diagnostics(
        dataset,
        fit.covariance,
        fit.residuals,
    )
    warnings = list(
        dict.fromkeys(
            [
                *readiness.warnings,
                *residual_diagnostics.warnings,
                *[warning.message for warning in identifiability_warnings],
            ]
        )
    )
    return BrownianRegimeFitSummaryReport(
        tree_path=tree_path,
        traits_path=traits_path,
        regime_map_path=regime_map_path,
        taxon_column=dataset.taxon_column,
        branch_id_column=resolved_branch_id_column,
        regime_column=regime_column,
        trait=trait,
        tree_taxon_count=readiness.tree_taxa,
        analyzed_taxa=list(dataset.taxa),
        analyzed_taxon_count=len(dataset.taxa),
        excluded_taxa=_build_excluded_taxa(readiness),
        branch_rows=branch_rows,
        regime_rows=regime_rows,
        profile_rows=profile_rows,
        root_state=fit.root_state,
        root_state_interval=fit.root_state_interval,
        log_likelihood=fit.log_likelihood,
        aic=next(row.aic for row in comparison_rows if row.model == "brownian-regimes"),
        aicc=next(
            row.aicc for row in comparison_rows if row.model == "brownian-regimes"
        ),
        comparison_rows=comparison_rows,
        better_model=better_model,
        likelihood_ratio_statistic=likelihood_ratio_statistic,
        likelihood_ratio_degrees_of_freedom=likelihood_ratio_degrees_of_freedom,
        likelihood_ratio_p_value=likelihood_ratio_p_value,
        likelihood_ratio_p_value_method="wilson-hilferty-approximation",
        identifiability_warnings=identifiability_warnings,
        residual_diagnostics=residual_diagnostics,
        assumptions=[
            "Every non-root branch must be assigned to one user-provided regime.",
            "Each regime contributes its own Brownian sigma-squared rate to shared-path covariance.",
            "Regime-specific uncertainty is reported from conditional likelihood profiles with other regimes fixed at their best-supported values.",
        ],
        warnings=warnings,
        readiness=readiness,
    )


@dataclass(slots=True)
class _MultirateFit:
    regime_rates: dict[str, float]
    root_state: float
    root_state_interval: ComparativeParameterInterval
    log_likelihood: float
    fitted_values: list[float]
    residuals: list[float]
    covariance: list[list[float]]


def _fit_multirate_brownian_model(
    dataset: ComparativeDataset,
    *,
    regime_matrices: dict[str, list[list[float]]],
    baseline_rate: float,
) -> _MultirateFit:
    lower = max(baseline_rate * 0.02, 1e-6)
    upper = max(baseline_rate * 50.0, lower * 10.0)
    regimes = sorted(regime_matrices)
    bounds = dict.fromkeys(regimes, (lower, upper))
    regime_rates = dict.fromkeys(regimes, baseline_rate)
    for _ in range(6):
        for regime in regimes:
            search = _logspace(bounds[regime][0], bounds[regime][1], 41)
            best_index = 0
            best_log_likelihood = -math.inf
            best_rate = search[0]
            for index, candidate in enumerate(search):
                trial_rates = dict(regime_rates)
                trial_rates[regime] = candidate
                fit = _fit_multirate_covariance(dataset, regime_matrices, trial_rates)
                if fit.log_likelihood > best_log_likelihood:
                    best_log_likelihood = fit.log_likelihood
                    best_index = index
                    best_rate = candidate
            regime_rates[regime] = best_rate
            bounds[regime] = (
                search[max(0, best_index - 3)],
                search[min(len(search) - 1, best_index + 3)],
            )
    fit = _fit_multirate_covariance(dataset, regime_matrices, regime_rates)
    return _MultirateFit(
        regime_rates=regime_rates,
        root_state=fit.root_state,
        root_state_interval=_root_state_interval(fit.root_state, fit.covariance),
        log_likelihood=fit.log_likelihood,
        fitted_values=fit.fitted_values,
        residuals=fit.residuals,
        covariance=fit.covariance,
    )


@dataclass(slots=True)
class _CovarianceFit:
    root_state: float
    log_likelihood: float
    fitted_values: list[float]
    residuals: list[float]
    covariance: list[list[float]]


def _fit_multirate_covariance(
    dataset: ComparativeDataset,
    regime_matrices: dict[str, list[list[float]]],
    regime_rates: dict[str, float],
) -> _CovarianceFit:
    covariance = stable_covariance(
        _combine_regime_covariance(regime_matrices, regime_rates)
    )
    inverse_covariance = invert_matrix(covariance)
    ones = [1.0] * len(dataset.trait_values)
    denom = quadratic_form(ones, inverse_covariance)
    root_state = (
        sum(
            sum(
                inverse_covariance[row_index][column_index]
                * dataset.trait_values[column_index]
                for column_index in range(len(dataset.trait_values))
            )
            for row_index in range(len(dataset.trait_values))
        )
        / denom
    )
    fitted_values = [root_state] * len(dataset.trait_values)
    residuals = [value - root_state for value in dataset.trait_values]
    log_likelihood = -0.5 * (
        len(dataset.trait_values) * math.log(2.0 * math.pi)
        + log_determinant(covariance)
        + quadratic_form(residuals, inverse_covariance)
    )
    return _CovarianceFit(
        root_state=root_state,
        log_likelihood=log_likelihood,
        fitted_values=fitted_values,
        residuals=residuals,
        covariance=covariance,
    )


def _combine_regime_covariance(
    regime_matrices: dict[str, list[list[float]]],
    regime_rates: dict[str, float],
) -> list[list[float]]:
    size = len(next(iter(regime_matrices.values())))
    combined = [[0.0] * size for _ in range(size)]
    for regime, matrix in regime_matrices.items():
        rate = regime_rates[regime]
        for row_index in range(size):
            for column_index in range(size):
                combined[row_index][column_index] += (
                    rate * matrix[row_index][column_index]
                )
    return combined


def _root_state_interval(
    root_state: float,
    covariance: list[list[float]],
) -> ComparativeParameterInterval:
    inverse_covariance = invert_matrix(covariance)
    ones = [1.0] * len(covariance)
    denom = quadratic_form(ones, inverse_covariance)
    root_se = math.sqrt(max(1.0 / denom, 1e-12))
    return ComparativeParameterInterval(
        name="root_state",
        estimate=root_state,
        lower_95=root_state - (_Z_95 * root_se),
        upper_95=root_state + (_Z_95 * root_se),
        method="wald",
    )


def _build_comparison_rows(
    *,
    taxon_count: int,
    baseline_log_likelihood: float,
    multirate_log_likelihood: float,
    regime_count: int,
) -> list[ComparativeModelComparisonRow]:
    rows = [
        _comparison_row("brownian", 2, baseline_log_likelihood, taxon_count),
        _comparison_row(
            "brownian-regimes",
            regime_count + 1,
            multirate_log_likelihood,
            taxon_count,
        ),
    ]
    best_aicc = min(row.aicc for row in rows)
    for row in rows:
        row.selected = math.isclose(row.aicc, best_aicc, rel_tol=0.0, abs_tol=1e-12)
    return rows


def _build_profile_rows(
    dataset: ComparativeDataset,
    *,
    regime_matrices: dict[str, list[list[float]]],
    fitted_rates: dict[str, float],
    baseline_rate: float,
) -> list[BrownianRegimeProfileRow]:
    lower = max(baseline_rate * 0.02, 1e-6)
    upper = max(baseline_rate * 50.0, lower * 10.0)
    rows: list[BrownianRegimeProfileRow] = []
    for regime in sorted(regime_matrices):
        candidates = sorted(
            {
                *_logspace(lower, upper, 81),
                fitted_rates[regime],
            }
        )
        profile = []
        best_log_likelihood = -math.inf
        for candidate in candidates:
            trial_rates = dict(fitted_rates)
            trial_rates[regime] = candidate
            fit = _fit_multirate_covariance(dataset, regime_matrices, trial_rates)
            profile.append((candidate, fit.log_likelihood))
            best_log_likelihood = max(best_log_likelihood, fit.log_likelihood)
        for candidate, log_likelihood in profile:
            rows.append(
                BrownianRegimeProfileRow(
                    regime=regime,
                    sigma_squared=candidate,
                    log_likelihood=log_likelihood,
                    delta_log_likelihood=best_log_likelihood - log_likelihood,
                    in_support_interval=(
                        (best_log_likelihood - log_likelihood)
                        <= _PROFILE_CONFIDENCE_DELTA
                    ),
                    selected=math.isclose(
                        candidate,
                        fitted_rates[regime],
                        rel_tol=0.0,
                        abs_tol=1e-12,
                    ),
                )
            )
    return rows


def _build_regime_rows(
    branch_rows: list[BrownianRegimeBranchRow],
    regime_rates: dict[str, float],
    profile_rows: list[BrownianRegimeProfileRow],
) -> list[BrownianRegimeRateRow]:
    by_regime: dict[str, list[BrownianRegimeBranchRow]] = {}
    for row in branch_rows:
        by_regime.setdefault(row.regime, []).append(row)
    profile_by_regime: dict[str, list[BrownianRegimeProfileRow]] = {}
    for row in profile_rows:
        profile_by_regime.setdefault(row.regime, []).append(row)
    rows: list[BrownianRegimeRateRow] = []
    for regime in sorted(regime_rates):
        branches = by_regime[regime]
        support = [
            row.sigma_squared
            for row in profile_by_regime.get(regime, [])
            if row.in_support_interval
        ]
        rows.append(
            BrownianRegimeRateRow(
                regime=regime,
                branch_count=len(branches),
                contributing_branch_count=sum(
                    1 for branch in branches if branch.contributes_to_analysis
                ),
                total_branch_length=sum(branch.branch_length for branch in branches),
                contributing_branch_length=sum(
                    branch.branch_length
                    for branch in branches
                    if branch.contributes_to_analysis
                ),
                sigma_squared=regime_rates[regime],
                lower_95=min(support) if support else regime_rates[regime],
                upper_95=max(support) if support else regime_rates[regime],
                interval_method="conditional_profile_likelihood_95",
            )
        )
    return rows


def _build_identifiability_warnings(
    regime_rows: list[BrownianRegimeRateRow],
    profile_rows: list[BrownianRegimeProfileRow],
    *,
    better_model: str,
) -> list[BrownianRegimeIdentifiabilityWarning]:
    warnings: list[BrownianRegimeIdentifiabilityWarning] = []
    profile_by_regime: dict[str, list[BrownianRegimeProfileRow]] = {}
    for row in profile_rows:
        profile_by_regime.setdefault(row.regime, []).append(row)
    for row in regime_rows:
        profile = profile_by_regime.get(row.regime, [])
        if not profile:
            continue
        supported_count = sum(
            1 for candidate in profile if candidate.in_support_interval
        )
        if supported_count >= math.ceil(len(profile) / 2.0):
            warnings.append(
                BrownianRegimeIdentifiabilityWarning(
                    regime=row.regime,
                    kind="flat_profile",
                    message=(
                        f"regime '{row.regime}' has a broad conditional likelihood "
                        "profile, so its Brownian rate is weakly identified"
                    ),
                )
            )
        if math.isclose(
            row.lower_95,
            min(candidate.sigma_squared for candidate in profile),
            rel_tol=0.0,
            abs_tol=1e-12,
        ) or math.isclose(
            row.upper_95,
            max(candidate.sigma_squared for candidate in profile),
            rel_tol=0.0,
            abs_tol=1e-12,
        ):
            warnings.append(
                BrownianRegimeIdentifiabilityWarning(
                    regime=row.regime,
                    kind="boundary_rate",
                    message=(
                        f"regime '{row.regime}' reaches a profile search boundary, "
                        "so its Brownian rate is weakly identified"
                    ),
                )
            )
        if row.contributing_branch_length <= 0.5:
            warnings.append(
                BrownianRegimeIdentifiabilityWarning(
                    regime=row.regime,
                    kind="low_branch_exposure",
                    message=(
                        f"regime '{row.regime}' contributes very little analyzed "
                        "branch length, so its Brownian rate is weakly identified"
                    ),
                )
            )
    if better_model != "brownian-regimes":
        warnings.append(
            BrownianRegimeIdentifiabilityWarning(
                regime="*",
                kind="comparison_not_preferred",
                message=(
                    "single-rate Brownian remains preferred over the multi-rate "
                    "regime model by AICc on the analyzed taxon set"
                ),
            )
        )
    return warnings


def _build_multirate_residual_diagnostics(
    dataset: ComparativeDataset,
    covariance: list[list[float]],
    residuals: list[float],
) -> ComparativeResidualSummary:
    residual_mean = sum(residuals) / len(residuals)
    residual_variance = sum((value - residual_mean) ** 2 for value in residuals) / max(
        1,
        len(residuals) - 1,
    )
    residual_skewness = 0.0
    if residual_variance > 0.0:
        residual_sd = math.sqrt(residual_variance)
        residual_skewness = sum(
            ((value - residual_mean) / residual_sd) ** 3 for value in residuals
        ) / len(residuals)
    inverse_covariance = invert_matrix(covariance)
    hat = _hat_matrix(len(dataset.taxa), covariance, inverse_covariance)
    outliers: list[ComparativeResidualOutlier] = []
    max_abs_standardized = 0.0
    for index, taxon in enumerate(dataset.taxa):
        leverage = min(max(hat[index][index], 0.0), 0.999999)
        denominator = math.sqrt(max(covariance[index][index] * (1.0 - leverage), 1e-12))
        standardized = residuals[index] / denominator
        max_abs_standardized = max(max_abs_standardized, abs(standardized))
        if abs(standardized) >= 2.0:
            outliers.append(
                ComparativeResidualOutlier(
                    taxon=taxon,
                    residual=residuals[index],
                    standardized_residual=standardized,
                    note="absolute standardized residual exceeds 2.0",
                )
            )
    residual_lambda = _estimate_lambda_for_values(dataset, residuals)
    warnings: list[str] = []
    if abs(residual_skewness) > 1.0:
        warnings.append("residual distribution is noticeably skewed")
    if residual_lambda > 0.5:
        warnings.append("residuals retain moderate phylogenetic structure")
    if outliers:
        warnings.append("one or more taxa have unusually large residuals")
    return ComparativeResidualSummary(
        residual_mean=residual_mean,
        residual_variance=residual_variance,
        residual_skewness=residual_skewness,
        max_abs_standardized_residual=max_abs_standardized,
        phylogenetic_residual_lambda=residual_lambda,
        outlier_taxa=outliers,
        warnings=warnings,
    )


def _hat_matrix(
    sample_size: int,
    covariance: list[list[float]],
    inverse_covariance: list[list[float]],
) -> list[list[float]]:
    design = [[1.0] for _ in range(sample_size)]
    xt_vinv = matrix_multiply(transpose(design), inverse_covariance)
    xt_vinv_x_inverse = invert_matrix(matrix_multiply(xt_vinv, design))
    return matrix_multiply(design, matrix_multiply(xt_vinv_x_inverse, xt_vinv))


def _load_branch_regime_rows(
    regime_map_path: Path,
    *,
    tree: PhyloTree,
    analyzed_taxa: list[str],
    branch_id_column: str | None,
    regime_column: str,
) -> tuple[list[BrownianRegimeBranchRow], str]:
    resolved_branch_id_column = resolve_branch_regime_id_column(
        regime_map_path,
        requested=branch_id_column,
    )
    table = load_taxon_table(
        regime_map_path,
        taxon_column=resolved_branch_id_column,
    )
    if regime_column not in table.columns:
        raise ComparativeMethodError(
            f"regime map does not contain column '{regime_column}'"
        )
    branch_lookup = build_branch_identity_lookup(tree, analyzed_taxa=analyzed_taxa)
    mapped_branch_ids = {row[table.taxon_column] for row in table.rows}
    expected_branch_ids = set(branch_lookup)
    missing = sorted(expected_branch_ids - mapped_branch_ids)
    if missing:
        raise ComparativeMethodError(
            "regime map is missing one or more non-root branches: " + ", ".join(missing)
        )
    extra = sorted(mapped_branch_ids - expected_branch_ids)
    if extra:
        raise ComparativeMethodError(
            "regime map contains branches absent from the tree: " + ", ".join(extra)
        )
    rows: list[BrownianRegimeBranchRow] = []
    for row in table.rows:
        branch = branch_lookup[row[table.taxon_column]]
        regime = row[regime_column]
        if not regime:
            raise ComparativeMethodError(
                f"regime map branch '{branch.branch_id}' has an empty '{regime_column}' value"
            )
        rows.append(
            BrownianRegimeBranchRow(
                branch_id=branch.branch_id,
                regime=regime,
                branch_length=branch.branch_length,
                descendant_taxa=branch.descendant_taxa,
                analyzed_descendant_taxa=branch.analyzed_descendant_taxa,
                contributes_to_analysis=branch.contributes_to_analysis,
            )
        )
    contributing_regimes = {
        row.regime
        for row in rows
        if row.contributes_to_analysis and row.branch_length > 0.0
    }
    if len(contributing_regimes) < 2:
        raise ComparativeMethodError(
            "regime map must expose at least two regimes that contribute analyzed branch length"
        )
    return sorted(rows, key=lambda row: row.branch_id), resolved_branch_id_column


def _build_regime_covariance_components(
    taxa: list[str],
    branch_rows: list[BrownianRegimeBranchRow],
) -> dict[str, list[list[float]]]:
    index = {taxon: position for position, taxon in enumerate(taxa)}
    components: dict[str, list[list[float]]] = {}
    for branch in branch_rows:
        if not branch.contributes_to_analysis or branch.branch_length <= 0.0:
            continue
        matrix = components.setdefault(
            branch.regime,
            [[0.0] * len(taxa) for _ in range(len(taxa))],
        )
        for left_taxon in branch.analyzed_descendant_taxa:
            for right_taxon in branch.analyzed_descendant_taxa:
                matrix[index[left_taxon]][index[right_taxon]] += branch.branch_length
    return components


def _build_excluded_taxa(
    readiness: ComparativeReadinessReport,
) -> list[BrownianRegimeExclusion]:
    rows: list[BrownianRegimeExclusion] = []
    rows.extend(
        BrownianRegimeExclusion(taxon=taxon, reason="missing_from_trait_table")
        for taxon in readiness.missing_from_traits
    )
    rows.extend(
        BrownianRegimeExclusion(taxon=taxon, reason="missing_trait_value")
        for taxon in readiness.pruned_missing_value_taxa
    )
    rows.extend(
        BrownianRegimeExclusion(taxon=taxon, reason="non_numeric_trait_value")
        for taxon in readiness.pruned_non_numeric_taxa
    )
    rows.extend(
        BrownianRegimeExclusion(taxon=taxon, reason="absent_from_tree")
        for taxon in readiness.extra_trait_taxa
    )
    return rows


def write_brownian_regime_summary_table(
    path: Path,
    report: BrownianRegimeFitSummaryReport,
) -> Path:
    """Write one overall summary ledger for a multi-rate Brownian fit."""
    return write_taxon_rows(
        path,
        columns=[
            "trait",
            "taxon_column",
            "branch_id_column",
            "regime_column",
            "tree_taxon_count",
            "analyzed_taxon_count",
            "excluded_taxon_count",
            "regime_count",
            "root_state",
            "root_state_lower_95",
            "root_state_upper_95",
            "log_likelihood",
            "aic",
            "aicc",
            "better_model",
            "likelihood_ratio_statistic",
            "likelihood_ratio_degrees_of_freedom",
            "likelihood_ratio_p_value",
            "identifiability_warning_count",
            "residual_variance",
            "max_abs_standardized_residual",
            "phylogenetic_residual_lambda",
        ],
        rows=[
            {
                "trait": report.trait,
                "taxon_column": report.taxon_column,
                "branch_id_column": report.branch_id_column,
                "regime_column": report.regime_column,
                "tree_taxon_count": report.tree_taxon_count,
                "analyzed_taxon_count": report.analyzed_taxon_count,
                "excluded_taxon_count": len(report.excluded_taxa),
                "regime_count": len(report.regime_rows),
                "root_state": format(report.root_state, ".15g"),
                "root_state_lower_95": format(
                    report.root_state_interval.lower_95, ".15g"
                ),
                "root_state_upper_95": format(
                    report.root_state_interval.upper_95, ".15g"
                ),
                "log_likelihood": format(report.log_likelihood, ".15g"),
                "aic": format(report.aic, ".15g"),
                "aicc": format(report.aicc, ".15g"),
                "better_model": report.better_model,
                "likelihood_ratio_statistic": format(
                    report.likelihood_ratio_statistic,
                    ".15g",
                ),
                "likelihood_ratio_degrees_of_freedom": report.likelihood_ratio_degrees_of_freedom,
                "likelihood_ratio_p_value": format(
                    report.likelihood_ratio_p_value,
                    ".15g",
                ),
                "identifiability_warning_count": len(report.identifiability_warnings),
                "residual_variance": format(
                    report.residual_diagnostics.residual_variance,
                    ".15g",
                ),
                "max_abs_standardized_residual": format(
                    report.residual_diagnostics.max_abs_standardized_residual,
                    ".15g",
                ),
                "phylogenetic_residual_lambda": format(
                    report.residual_diagnostics.phylogenetic_residual_lambda,
                    ".15g",
                ),
            }
        ],
    )


def write_brownian_regime_comparison_table(
    path: Path,
    report: BrownianRegimeFitSummaryReport,
) -> Path:
    """Write one single-rate versus multi-rate Brownian comparison ledger."""
    best_aicc = min(row.aicc for row in report.comparison_rows)
    rows = [
        {
            "row_kind": "model_fit",
            "model": row.model,
            "comparison_id": "",
            "parameter_count": row.parameter_count,
            "log_likelihood": format(row.log_likelihood, ".15g"),
            "aic": format(row.aic, ".15g"),
            "aicc": format(row.aicc, ".15g"),
            "delta_aicc": format(row.aicc - best_aicc, ".15g"),
            "selected": str(row.selected).lower(),
            "left_model": "",
            "right_model": "",
            "statistic": "",
            "degrees_of_freedom": "",
            "p_value": "",
            "p_value_method": "",
        }
        for row in report.comparison_rows
    ]
    rows.append(
        {
            "row_kind": "likelihood_ratio_test",
            "model": "",
            "comparison_id": "brownian-vs-brownian-regimes",
            "parameter_count": "",
            "log_likelihood": "",
            "aic": "",
            "aicc": "",
            "delta_aicc": "",
            "selected": "",
            "left_model": "brownian",
            "right_model": "brownian-regimes",
            "statistic": format(report.likelihood_ratio_statistic, ".15g"),
            "degrees_of_freedom": report.likelihood_ratio_degrees_of_freedom,
            "p_value": format(report.likelihood_ratio_p_value, ".15g"),
            "p_value_method": report.likelihood_ratio_p_value_method,
        }
    )
    return write_taxon_rows(
        path,
        columns=[
            "row_kind",
            "model",
            "comparison_id",
            "parameter_count",
            "log_likelihood",
            "aic",
            "aicc",
            "delta_aicc",
            "selected",
            "left_model",
            "right_model",
            "statistic",
            "degrees_of_freedom",
            "p_value",
            "p_value_method",
        ],
        rows=rows,
    )


def write_brownian_regime_rate_table(
    path: Path,
    report: BrownianRegimeFitSummaryReport,
) -> Path:
    """Write one per-regime rate ledger for a multi-rate Brownian fit."""
    return write_taxon_rows(
        path,
        columns=[
            "regime",
            "branch_count",
            "contributing_branch_count",
            "total_branch_length",
            "contributing_branch_length",
            "sigma_squared",
            "sigma_squared_lower_95",
            "sigma_squared_upper_95",
            "interval_method",
        ],
        rows=[
            {
                "regime": row.regime,
                "branch_count": row.branch_count,
                "contributing_branch_count": row.contributing_branch_count,
                "total_branch_length": format(row.total_branch_length, ".15g"),
                "contributing_branch_length": format(
                    row.contributing_branch_length,
                    ".15g",
                ),
                "sigma_squared": format(row.sigma_squared, ".15g"),
                "sigma_squared_lower_95": format(row.lower_95, ".15g"),
                "sigma_squared_upper_95": format(row.upper_95, ".15g"),
                "interval_method": row.interval_method,
            }
            for row in report.regime_rows
        ],
    )


def write_brownian_regime_profile_table(
    path: Path,
    report: BrownianRegimeFitSummaryReport,
) -> Path:
    """Write one conditional rate-profile ledger for a multi-rate Brownian fit."""
    return write_taxon_rows(
        path,
        columns=[
            "regime",
            "sigma_squared",
            "log_likelihood",
            "delta_log_likelihood",
            "in_support_interval",
            "selected",
        ],
        rows=[
            {
                "regime": row.regime,
                "sigma_squared": format(row.sigma_squared, ".15g"),
                "log_likelihood": format(row.log_likelihood, ".15g"),
                "delta_log_likelihood": format(row.delta_log_likelihood, ".15g"),
                "in_support_interval": str(row.in_support_interval).lower(),
                "selected": str(row.selected).lower(),
            }
            for row in report.profile_rows
        ],
    )


def write_brownian_regime_branch_table(
    path: Path,
    report: BrownianRegimeFitSummaryReport,
) -> Path:
    """Write one normalized branch-regime assignment ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "branch_id",
            "regime",
            "branch_length",
            "descendant_taxa",
            "analyzed_descendant_taxa",
            "contributes_to_analysis",
        ],
        rows=[
            {
                "branch_id": row.branch_id,
                "regime": row.regime,
                "branch_length": format(row.branch_length, ".15g"),
                "descendant_taxa": ",".join(row.descendant_taxa),
                "analyzed_descendant_taxa": ",".join(row.analyzed_descendant_taxa),
                "contributes_to_analysis": str(row.contributes_to_analysis).lower(),
            }
            for row in report.branch_rows
        ],
    )


def write_brownian_regime_exclusion_table(
    path: Path,
    report: BrownianRegimeFitSummaryReport,
) -> Path:
    """Write one excluded-taxon ledger for a multi-rate Brownian fit."""
    return write_taxon_rows(
        path,
        columns=["taxon", "reason"],
        rows=[
            {"taxon": row.taxon, "reason": row.reason} for row in report.excluded_taxa
        ],
    )


def _logspace(lower: float, upper: float, count: int) -> list[float]:
    if count < 2:
        return [lower]
    low_log = math.log(lower)
    high_log = math.log(upper)
    return [
        math.exp(low_log + ((high_log - low_log) * (index / float(count - 1))))
        for index in range(count)
    ]


def _chi_square_survival(statistic: float, degrees_of_freedom: int) -> float:
    if statistic <= 0.0:
        return 1.0
    if degrees_of_freedom <= 0:
        return 1.0
    if degrees_of_freedom == 1:
        return math.erfc(math.sqrt(statistic / 2.0))
    z_score = (
        ((statistic / degrees_of_freedom) ** (1.0 / 3.0))
        - (1.0 - (2.0 / (9.0 * degrees_of_freedom)))
    ) / math.sqrt(2.0 / (9.0 * degrees_of_freedom))
    return 0.5 * math.erfc(z_score / math.sqrt(2.0))
