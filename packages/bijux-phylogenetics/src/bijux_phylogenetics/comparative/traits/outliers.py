from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path

from bijux_phylogenetics.comparative._math import (
    dot,
    invert_matrix,
    matrix_vector_multiply,
    stable_covariance,
)
from bijux_phylogenetics.comparative.common import (
    ComparativeReadinessReport,
    build_brownian_covariance_matrix,
    build_ou_covariance_matrix,
    descendant_taxa,
    load_comparative_dataset,
    node_signature,
    summarize_numeric_trait_readiness,
)
from bijux_phylogenetics.comparative.continuous.model_fitting import (
    compare_brownian_and_ou_models,
    fit_brownian_motion_model,
    fit_ornstein_uhlenbeck_model,
)
from bijux_phylogenetics.comparative.model_selection import (
    ComparativeModelComparisonRow,
)
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.phylo.pruning import prune_tree_to_requested_taxa
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

_OUTLIER_THRESHOLD = 2.0


@dataclass(slots=True)
class TraitOutlierExclusion:
    """One taxon excluded before phylogenetic outlier scoring."""

    taxon: str
    reason: str


@dataclass(slots=True)
class TraitOutlierTaxonRow:
    """One taxon ranked by leave-one-taxon-out phylogenetic residual size."""

    taxon: str
    observed_value: float
    conditional_expected_value: float
    residual: float
    residual_direction: str
    conditional_variance: float
    conditional_standard_error: float
    standardized_residual: float
    abs_standardized_residual: float
    context_clade_id: str
    context_node_label: str | None
    context_taxon_count: int
    context_taxa: list[str]
    context_mean: float | None
    sibling_context_id: str | None
    sibling_taxon_count: int
    sibling_taxa: list[str]
    sibling_mean: float | None
    context_mean_shift: float | None
    outlier: bool
    rank: int


@dataclass(slots=True)
class TraitOutlierSummaryReport:
    """Reviewer-facing trait outlier review under the selected phylogenetic model."""

    tree_path: Path
    traits_path: Path
    taxon_column: str
    trait: str
    tree_taxon_count: int
    analyzed_taxa: list[str]
    analyzed_taxon_count: int
    excluded_taxa: list[TraitOutlierExclusion]
    selected_model: str
    selected_mean_parameter: str
    selected_mean_value: float
    selected_alpha: float | None
    selected_sigma_squared: float
    model_comparison_rows: list[ComparativeModelComparisonRow]
    outlier_threshold: float
    taxon_rows: list[TraitOutlierTaxonRow]
    outlier_taxa: list[str]
    top_outlier_taxon: str | None
    top_abs_standardized_residual: float | None
    assumptions: list[str]
    warnings: list[str]
    readiness: ComparativeReadinessReport


def summarize_trait_outliers(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
) -> TraitOutlierSummaryReport:
    """Rank continuous-trait outliers by conditional phylogenetic residual size."""
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
    model_comparison = compare_brownian_and_ou_models(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    (
        selected_model,
        selected_mean_parameter,
        selected_mean_value,
        selected_alpha,
        selected_sigma_squared,
        assumptions,
        warnings,
    ) = _selected_model_parameters(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        selected_model=model_comparison.better_model,
        readiness=readiness,
    )
    pruned_tree, _ = prune_tree_to_requested_taxa(tree_path, dataset.taxa)
    covariance = _scaled_selected_covariance(
        pruned_tree,
        dataset.taxa,
        model=selected_model,
        sigma_squared=selected_sigma_squared,
        alpha=selected_alpha,
    )
    values_by_taxon = dict(zip(dataset.taxa, dataset.trait_values, strict=True))
    context_by_taxon = _build_context_by_taxon(pruned_tree, values_by_taxon)
    taxon_rows = _build_taxon_rows(
        dataset.taxa,
        dataset.trait_values,
        covariance,
        selected_mean_value,
        context_by_taxon,
    )
    outlier_taxa = [row.taxon for row in taxon_rows if row.outlier]
    if outlier_taxa:
        warnings = [
            *warnings,
            "one or more taxa have unusually large conditional phylogenetic residuals",
        ]
    return TraitOutlierSummaryReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=dataset.taxon_column,
        trait=trait,
        tree_taxon_count=readiness.tree_taxa,
        analyzed_taxa=list(dataset.taxa),
        analyzed_taxon_count=len(dataset.taxa),
        excluded_taxa=_build_excluded_taxa(readiness),
        selected_model=selected_model,
        selected_mean_parameter=selected_mean_parameter,
        selected_mean_value=selected_mean_value,
        selected_alpha=selected_alpha,
        selected_sigma_squared=selected_sigma_squared,
        model_comparison_rows=list(model_comparison.rows),
        outlier_threshold=_OUTLIER_THRESHOLD,
        taxon_rows=taxon_rows,
        outlier_taxa=outlier_taxa,
        top_outlier_taxon=taxon_rows[0].taxon if taxon_rows else None,
        top_abs_standardized_residual=(
            taxon_rows[0].abs_standardized_residual if taxon_rows else None
        ),
        assumptions=assumptions,
        warnings=list(dict.fromkeys(warnings)),
        readiness=readiness,
    )


def write_trait_outlier_summary_table(
    path: Path,
    report: TraitOutlierSummaryReport,
) -> Path:
    """Write one summary ledger for phylogenetic trait outlier review."""
    aicc_by_model = {row.model: row.aicc for row in report.model_comparison_rows}
    return write_taxon_rows(
        path,
        columns=[
            "trait",
            "taxon_column",
            "tree_taxon_count",
            "analyzed_taxon_count",
            "excluded_taxon_count",
            "selected_model",
            "selected_mean_parameter",
            "selected_mean_value",
            "selected_alpha",
            "selected_sigma_squared",
            "brownian_aicc",
            "ou_aicc",
            "outlier_threshold",
            "outlier_count",
            "top_outlier_taxon",
            "top_abs_standardized_residual",
        ],
        rows=[
            {
                "trait": report.trait,
                "taxon_column": report.taxon_column,
                "tree_taxon_count": str(report.tree_taxon_count),
                "analyzed_taxon_count": str(report.analyzed_taxon_count),
                "excluded_taxon_count": str(len(report.excluded_taxa)),
                "selected_model": report.selected_model,
                "selected_mean_parameter": report.selected_mean_parameter,
                "selected_mean_value": format(report.selected_mean_value, ".15g"),
                "selected_alpha": _format_optional(report.selected_alpha),
                "selected_sigma_squared": format(report.selected_sigma_squared, ".15g"),
                "brownian_aicc": _format_optional(aicc_by_model.get("brownian")),
                "ou_aicc": _format_optional(aicc_by_model.get("ou")),
                "outlier_threshold": format(report.outlier_threshold, ".15g"),
                "outlier_count": str(len(report.outlier_taxa)),
                "top_outlier_taxon": report.top_outlier_taxon or "",
                "top_abs_standardized_residual": _format_optional(
                    report.top_abs_standardized_residual
                ),
            }
        ],
    )


def write_trait_outlier_taxon_table(
    path: Path,
    report: TraitOutlierSummaryReport,
) -> Path:
    """Write one ranked taxon ledger for phylogenetic trait outlier review."""
    return write_taxon_rows(
        path,
        columns=[
            "taxon",
            "observed_value",
            "conditional_expected_value",
            "residual",
            "residual_direction",
            "conditional_variance",
            "conditional_standard_error",
            "standardized_residual",
            "abs_standardized_residual",
            "context_clade_id",
            "context_node_label",
            "context_taxon_count",
            "context_taxa",
            "context_mean",
            "sibling_context_id",
            "sibling_taxon_count",
            "sibling_taxa",
            "sibling_mean",
            "context_mean_shift",
            "outlier",
            "rank",
        ],
        rows=[
            {
                "taxon": row.taxon,
                "observed_value": format(row.observed_value, ".15g"),
                "conditional_expected_value": format(
                    row.conditional_expected_value, ".15g"
                ),
                "residual": format(row.residual, ".15g"),
                "residual_direction": row.residual_direction,
                "conditional_variance": format(row.conditional_variance, ".15g"),
                "conditional_standard_error": format(
                    row.conditional_standard_error, ".15g"
                ),
                "standardized_residual": format(row.standardized_residual, ".15g"),
                "abs_standardized_residual": format(
                    row.abs_standardized_residual, ".15g"
                ),
                "context_clade_id": row.context_clade_id,
                "context_node_label": row.context_node_label or "",
                "context_taxon_count": str(row.context_taxon_count),
                "context_taxa": ",".join(row.context_taxa),
                "context_mean": _format_optional(row.context_mean),
                "sibling_context_id": row.sibling_context_id or "",
                "sibling_taxon_count": str(row.sibling_taxon_count),
                "sibling_taxa": ",".join(row.sibling_taxa),
                "sibling_mean": _format_optional(row.sibling_mean),
                "context_mean_shift": _format_optional(row.context_mean_shift),
                "outlier": str(row.outlier).lower(),
                "rank": str(row.rank),
            }
            for row in report.taxon_rows
        ],
    )


def write_trait_outlier_exclusion_table(
    path: Path,
    report: TraitOutlierSummaryReport,
) -> Path:
    """Write one excluded-taxa ledger for phylogenetic trait outlier review."""
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


def _selected_model_parameters(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None,
    selected_model: str,
    readiness: ComparativeReadinessReport,
) -> tuple[str, str, float, float | None, float, list[str], list[str]]:
    if selected_model == "brownian":
        fit = fit_brownian_motion_model(
            tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
        )
        warnings = [*readiness.warnings, *fit.residual_diagnostics.warnings]
        assumptions = [
            *fit.assumptions,
            "Taxon-level outlier scores condition each tip on all other retained tips under the selected covariance model",
        ]
        return (
            "brownian",
            "root_state",
            fit.root_state,
            None,
            fit.rate,
            assumptions,
            warnings,
        )
    if selected_model == "ou":
        fit = fit_ornstein_uhlenbeck_model(
            tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
        )
        warnings = [
            *readiness.warnings,
            *fit.residual_diagnostics.warnings,
            *[warning.message for warning in fit.identifiability_warnings],
        ]
        assumptions = [
            *fit.assumptions,
            "Taxon-level outlier scores condition each tip on all other retained tips under the selected covariance model",
        ]
        return (
            "ou",
            "theta",
            fit.theta,
            fit.alpha,
            fit.sigma_squared,
            assumptions,
            warnings,
        )
    raise ValueError(f"unsupported selected outlier model: {selected_model}")


def _scaled_selected_covariance(
    tree: PhyloTree,
    taxa: list[str],
    *,
    model: str,
    sigma_squared: float,
    alpha: float | None,
) -> list[list[float]]:
    if model == "brownian":
        base = build_brownian_covariance_matrix(tree, taxa)
    elif model == "ou":
        if alpha is None:
            raise ValueError("OU trait outlier scoring requires an alpha value")
        base = build_ou_covariance_matrix(tree, taxa, alpha=alpha)
    else:
        raise ValueError(f"unsupported covariance model: {model}")
    return stable_covariance([[value * sigma_squared for value in row] for row in base])


def _build_taxon_rows(
    taxa: list[str],
    values: list[float],
    covariance: list[list[float]],
    mean_value: float,
    context_by_taxon: dict[str, _TaxonContext],
) -> list[TraitOutlierTaxonRow]:
    rows: list[TraitOutlierTaxonRow] = []
    for index, taxon in enumerate(taxa):
        expected_value, conditional_variance = _conditional_prediction(
            values,
            covariance,
            mean_value,
            index=index,
        )
        residual = values[index] - expected_value
        conditional_standard_error = math.sqrt(max(conditional_variance, 1e-12))
        standardized_residual = residual / conditional_standard_error
        context = context_by_taxon[taxon]
        rows.append(
            TraitOutlierTaxonRow(
                taxon=taxon,
                observed_value=values[index],
                conditional_expected_value=expected_value,
                residual=residual,
                residual_direction=(
                    "higher_than_expected" if residual >= 0.0 else "lower_than_expected"
                ),
                conditional_variance=conditional_variance,
                conditional_standard_error=conditional_standard_error,
                standardized_residual=standardized_residual,
                abs_standardized_residual=abs(standardized_residual),
                context_clade_id=context.context_clade_id,
                context_node_label=context.context_node_label,
                context_taxon_count=len(context.context_taxa),
                context_taxa=context.context_taxa,
                context_mean=context.context_mean,
                sibling_context_id=context.sibling_context_id,
                sibling_taxon_count=len(context.sibling_taxa),
                sibling_taxa=context.sibling_taxa,
                sibling_mean=context.sibling_mean,
                context_mean_shift=context.context_mean_shift,
                outlier=abs(standardized_residual) >= _OUTLIER_THRESHOLD,
                rank=0,
            )
        )
    rows.sort(
        key=lambda row: (
            -row.abs_standardized_residual,
            -abs(row.residual),
            row.taxon,
        )
    )
    for rank, row in enumerate(rows, start=1):
        row.rank = rank
    return rows


def _conditional_prediction(
    values: list[float],
    covariance: list[list[float]],
    mean_value: float,
    *,
    index: int,
) -> tuple[float, float]:
    other_indices = [
        candidate for candidate in range(len(values)) if candidate != index
    ]
    cross_covariance = [covariance[index][candidate] for candidate in other_indices]
    other_covariance = stable_covariance(
        [[covariance[row][column] for column in other_indices] for row in other_indices]
    )
    centered_other = [values[candidate] - mean_value for candidate in other_indices]
    inverse_other = invert_matrix(other_covariance)
    weighted_other = matrix_vector_multiply(inverse_other, centered_other)
    expected_value = mean_value + dot(cross_covariance, weighted_other)
    weighted_cross = matrix_vector_multiply(inverse_other, cross_covariance)
    conditional_variance = max(
        covariance[index][index] - dot(cross_covariance, weighted_cross),
        0.0,
    )
    return expected_value, conditional_variance


@dataclass(slots=True)
class _TaxonContext:
    context_clade_id: str
    context_node_label: str | None
    context_taxa: list[str]
    context_mean: float | None
    sibling_context_id: str | None
    sibling_taxa: list[str]
    sibling_mean: float | None
    context_mean_shift: float | None


def _build_context_by_taxon(
    tree: PhyloTree,
    values_by_taxon: dict[str, float],
) -> dict[str, _TaxonContext]:
    contexts: dict[str, _TaxonContext] = {}
    for leaf in tree.iter_leaves():
        if leaf.name is None or leaf.name not in values_by_taxon:
            continue
        parent = leaf.parent
        if parent is None:
            contexts[leaf.name] = _TaxonContext(
                context_clade_id=leaf.name,
                context_node_label=None,
                context_taxa=[leaf.name],
                context_mean=values_by_taxon[leaf.name],
                sibling_context_id=None,
                sibling_taxa=[],
                sibling_mean=None,
                context_mean_shift=None,
            )
            continue
        context_taxa = [
            taxon for taxon in descendant_taxa(parent) if taxon in values_by_taxon
        ]
        sibling_taxa = sorted(
            taxon
            for child in parent.children
            if child is not leaf
            for taxon in descendant_taxa(child)
            if taxon in values_by_taxon
        )
        context_mean = _mean_for_taxa(context_taxa, values_by_taxon)
        sibling_mean = _mean_for_taxa(sibling_taxa, values_by_taxon)
        contexts[leaf.name] = _TaxonContext(
            context_clade_id=node_signature(parent),
            context_node_label=parent.name,
            context_taxa=context_taxa,
            context_mean=context_mean,
            sibling_context_id=_taxa_signature(sibling_taxa),
            sibling_taxa=sibling_taxa,
            sibling_mean=sibling_mean,
            context_mean_shift=(
                None
                if context_mean is None or sibling_mean is None
                else context_mean - sibling_mean
            ),
        )
    return contexts


def _mean_for_taxa(
    taxa: list[str],
    values_by_taxon: dict[str, float],
) -> float | None:
    if not taxa:
        return None
    return sum(values_by_taxon[taxon] for taxon in taxa) / len(taxa)


def _taxa_signature(taxa: list[str]) -> str | None:
    if not taxa:
        return None
    return "|".join(sorted(taxa))


def _build_excluded_taxa(
    readiness: ComparativeReadinessReport,
) -> list[TraitOutlierExclusion]:
    rows: list[TraitOutlierExclusion] = []
    rows.extend(
        TraitOutlierExclusion(taxon=taxon, reason="missing_from_trait_table")
        for taxon in readiness.missing_from_traits
    )
    rows.extend(
        TraitOutlierExclusion(taxon=taxon, reason="missing_trait_value")
        for taxon in readiness.pruned_missing_value_taxa
    )
    rows.extend(
        TraitOutlierExclusion(taxon=taxon, reason="non_numeric_trait_value")
        for taxon in readiness.pruned_non_numeric_taxa
    )
    rows.extend(
        TraitOutlierExclusion(taxon=taxon, reason="absent_from_tree")
        for taxon in readiness.extra_trait_taxa
    )
    return rows


def _format_optional(value: float | None) -> str:
    if value is None:
        return ""
    return format(value, ".15g")
