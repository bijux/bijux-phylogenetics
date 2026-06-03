from __future__ import annotations

from pathlib import Path
import tempfile

from bijux_phylogenetics.comparative.pgls import PGLSResult
from bijux_phylogenetics.comparative.reporting.contracts import (
    ComparativeCoefficientDeltaRow,
    ComparativePruningComparisonReport,
    ComparativeTreeComparisonReport,
)
from bijux_phylogenetics.comparative.reporting.snapshot import (
    build_comparative_model_snapshot,
)
from bijux_phylogenetics.datasets.study_inputs import load_taxon_table, write_taxon_rows
from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.pruning import prune_tree_to_requested_taxa


def compare_comparative_results_across_trees(
    left_tree_path: Path,
    right_tree_path: Path,
    traits_path: Path,
    *,
    response: str | None = None,
    predictors: list[str] | None = None,
    formula: str | None = None,
    taxon_column: str | None = None,
    lambda_value: float | str = "estimate",
) -> ComparativeTreeComparisonReport:
    """Compare comparative fits across two alternative trees."""
    left = build_comparative_model_snapshot(
        left_tree_path,
        traits_path,
        response=response,
        predictors=predictors,
        formula=formula,
        taxon_column=taxon_column,
        lambda_value=lambda_value,
    )
    right = build_comparative_model_snapshot(
        right_tree_path,
        traits_path,
        response=response,
        predictors=predictors,
        formula=formula,
        taxon_column=taxon_column,
        lambda_value=lambda_value,
    )
    coefficient_deltas = _coefficient_delta_rows(left.pgls_model, right.pgls_model)
    sign_changed_terms = [row.term for row in coefficient_deltas if row.sign_changed]
    warnings = _comparative_sensitivity_warnings(
        left_selected_model=left.model_comparison.better_model,
        right_selected_model=right.model_comparison.better_model,
        sign_changed_terms=sign_changed_terms,
        delta_blombergs_k=right.signal_k.k - left.signal_k.k,
        delta_pagels_lambda=right.signal_lambda.lambda_value
        - left.signal_lambda.lambda_value,
        change_context="alternative-tree comparison",
    )
    return ComparativeTreeComparisonReport(
        left_tree_path=left_tree_path,
        right_tree_path=right_tree_path,
        response=left.response,
        predictors=left.formula.predictors,
        left_selected_model=left.model_comparison.better_model,
        right_selected_model=right.model_comparison.better_model,
        delta_blombergs_k=right.signal_k.k - left.signal_k.k,
        delta_pagels_lambda=right.signal_lambda.lambda_value
        - left.signal_lambda.lambda_value,
        delta_brownian_rate=right.brownian.rate - left.brownian.rate,
        delta_ou_alpha=right.ou.alpha - left.ou.alpha,
        coefficient_deltas=coefficient_deltas,
        sign_changed_terms=sign_changed_terms,
        conclusion_changed=bool(
            sign_changed_terms
            or left.model_comparison.better_model != right.model_comparison.better_model
        ),
        warnings=warnings,
    )


def compare_comparative_results_across_pruning(
    tree_path: Path,
    traits_path: Path,
    *,
    response: str | None = None,
    predictors: list[str] | None = None,
    formula: str | None = None,
    drop_taxa: list[str] | None = None,
    keep_taxa: list[str] | None = None,
    taxon_column: str | None = None,
    lambda_value: float | str = "estimate",
) -> ComparativePruningComparisonReport:
    """Compare comparative fits before and after explicit tree/trait pruning."""
    if bool(drop_taxa) == bool(keep_taxa):
        raise ValueError(
            "provide exactly one of drop_taxa or keep_taxa for pruning comparison"
        )
    baseline = build_comparative_model_snapshot(
        tree_path,
        traits_path,
        response=response,
        predictors=predictors,
        formula=formula,
        taxon_column=taxon_column,
        lambda_value=lambda_value,
    )
    tree = load_tree(tree_path)
    kept_taxa = sorted(set(keep_taxa or tree.tip_names) - set(drop_taxa or []))
    pruned_tree, pruning_report = prune_tree_to_requested_taxa(tree_path, kept_taxa)
    pruned_tree_path = Path(
        tempfile.mkstemp(prefix="bijux-comparative-pruned-", suffix=".nwk")[1]
    )
    pruned_traits_path = Path(
        tempfile.mkstemp(
            prefix="bijux-comparative-pruned-", suffix=Path(traits_path).suffix
        )[1]
    )
    try:
        write_newick(pruned_tree_path, pruned_tree)
        _write_pruned_trait_table(
            traits_path,
            pruned_traits_path,
            kept_taxa=pruning_report.kept_taxa,
            taxon_column=taxon_column,
        )
        pruned = build_comparative_model_snapshot(
            pruned_tree_path,
            pruned_traits_path,
            response=response,
            predictors=predictors,
            formula=formula,
            taxon_column=taxon_column,
            lambda_value=lambda_value,
        )
    finally:
        pruned_tree_path.unlink(missing_ok=True)
        pruned_traits_path.unlink(missing_ok=True)
    coefficient_deltas = _coefficient_delta_rows(baseline.pgls_model, pruned.pgls_model)
    sign_changed_terms = [row.term for row in coefficient_deltas if row.sign_changed]
    warnings = _comparative_sensitivity_warnings(
        left_selected_model=baseline.model_comparison.better_model,
        right_selected_model=pruned.model_comparison.better_model,
        sign_changed_terms=sign_changed_terms,
        delta_blombergs_k=pruned.signal_k.k - baseline.signal_k.k,
        delta_pagels_lambda=pruned.signal_lambda.lambda_value
        - baseline.signal_lambda.lambda_value,
        change_context="taxon-pruning comparison",
    )
    return ComparativePruningComparisonReport(
        tree_path=tree_path,
        response=baseline.response,
        predictors=baseline.formula.predictors,
        baseline_taxa=baseline.pgls_model.taxa,
        pruned_taxa=pruned.pgls_model.taxa,
        dropped_taxa=sorted(
            set(baseline.pgls_model.taxa) - set(pruned.pgls_model.taxa)
        ),
        delta_blombergs_k=pruned.signal_k.k - baseline.signal_k.k,
        delta_pagels_lambda=pruned.signal_lambda.lambda_value
        - baseline.signal_lambda.lambda_value,
        coefficient_deltas=coefficient_deltas,
        baseline_selected_model=baseline.model_comparison.better_model,
        pruned_selected_model=pruned.model_comparison.better_model,
        sign_changed_terms=sign_changed_terms,
        conclusion_changed=bool(
            sign_changed_terms
            or baseline.model_comparison.better_model
            != pruned.model_comparison.better_model
        ),
        warnings=warnings,
    )


def _coefficient_delta_rows(
    left: PGLSResult, right: PGLSResult
) -> list[ComparativeCoefficientDeltaRow]:
    left_coefficients = {row.name: row.estimate for row in left.coefficients}
    right_coefficients = {row.name: row.estimate for row in right.coefficients}
    rows: list[ComparativeCoefficientDeltaRow] = []
    for term in sorted(set(left_coefficients) | set(right_coefficients)):
        if term == "intercept":
            continue
        left_estimate = left_coefficients.get(term)
        right_estimate = right_coefficients.get(term)
        delta = None
        sign_changed = False
        if left_estimate is not None and right_estimate is not None:
            delta = right_estimate - left_estimate
            sign_changed = (left_estimate < 0.0 < right_estimate) or (
                right_estimate < 0.0 < left_estimate
            )
        rows.append(
            ComparativeCoefficientDeltaRow(
                term=term,
                left_estimate=left_estimate,
                right_estimate=right_estimate,
                delta=delta,
                sign_changed=sign_changed,
            )
        )
    return rows


def _comparative_sensitivity_warnings(
    *,
    left_selected_model: str,
    right_selected_model: str,
    sign_changed_terms: list[str],
    delta_blombergs_k: float,
    delta_pagels_lambda: float,
    change_context: str,
) -> list[str]:
    warnings: list[str] = []
    if left_selected_model != right_selected_model:
        warnings.append(
            f"{change_context} changes the selected comparative process model from {left_selected_model} to {right_selected_model}"
        )
    if sign_changed_terms:
        warnings.append(
            f"{change_context} flips the sign of one or more coefficient estimates: {', '.join(sign_changed_terms)}"
        )
    if abs(delta_blombergs_k) >= 0.2:
        warnings.append(f"{change_context} materially changes Blomberg's K")
    if abs(delta_pagels_lambda) >= 0.2:
        warnings.append(f"{change_context} materially changes Pagel's lambda")
    return warnings


def _write_pruned_trait_table(
    source_path: Path,
    out_path: Path,
    *,
    kept_taxa: list[str],
    taxon_column: str | None,
) -> Path:
    table = load_taxon_table(source_path, taxon_column=taxon_column)
    kept = set(kept_taxa)
    rows = [row for row in table.rows if row[table.taxon_column] in kept]
    return write_taxon_rows(out_path, columns=table.columns, rows=rows)
