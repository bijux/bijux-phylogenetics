from __future__ import annotations

from pathlib import Path
import tempfile

from bijux_phylogenetics.comparative.pgls import run_pgls
from bijux_phylogenetics.datasets.study_inputs import load_taxon_table, write_taxon_rows
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.pruning import prune_tree_to_requested_taxa
from bijux_phylogenetics.runtime.errors import ComparativeMethodError

from .contracts import (
    MULTIVARIATE_MISSING_VALUE_POLICY,
    MULTIVARIATE_NUMERICAL_TOLERANCE,
    MultivariateComparativeRegressionReport,
    MultivariateTaxonExclusion,
)
from .input_curation import (
    build_response_formula,
    build_shared_taxon_exclusions,
    inspect_response_model,
    raise_for_input_blockers,
    shared_analysis_taxa,
)
from .report_rows import build_response_coefficient_rows, build_response_model_rows
from .residual_analysis import (
    build_multivariate_warnings,
    build_residual_association_rows,
    build_residual_correlation_rows,
    build_residual_covariance_diagnostics,
    build_residual_covariance_rows,
)


def run_multivariate_comparative_regression(
    tree_path: Path,
    traits_path: Path,
    *,
    responses: list[str],
    predictors: list[str],
    taxon_column: str | None = None,
    lambda_value: float | str = "estimate",
) -> MultivariateComparativeRegressionReport:
    """Fit one shared-taxon comparative workflow across multiple response traits."""
    if len(responses) < 2:
        raise ComparativeMethodError(
            "multivariate comparative regression requires at least two response traits"
        )
    if not predictors:
        raise ComparativeMethodError(
            "multivariate comparative regression requires at least one predictor"
        )

    table = load_taxon_table(traits_path, taxon_column=taxon_column)
    tree = load_tree(tree_path)
    rows_by_taxon = {row[table.taxon_column]: row for row in table.rows}
    tree_taxa = set(tree.tip_names)
    table_taxa = set(table.taxa)
    excluded_taxa = build_initial_taxon_exclusions(
        responses=responses,
        predictors=predictors,
        table_taxa=table_taxa,
        tree_taxa=tree_taxa,
    )

    overlap_taxa = sorted(tree_taxa & table_taxa)
    if len(overlap_taxa) < 2:
        raise ComparativeMethodError(
            "multivariate comparative regression requires at least two taxa shared between the tree and trait table"
        )

    overlap_tree, _ = prune_tree_to_requested_taxa(tree_path, overlap_taxa)
    overlap_rows = [rows_by_taxon[taxon] for taxon in overlap_taxa]
    with tempfile.TemporaryDirectory(
        prefix="bijux-phylogenetics-multivariate-"
    ) as tmp_dir:
        overlap_tree_path = Path(tmp_dir) / "multivariate-overlap-tree.nwk"
        overlap_table_path = Path(tmp_dir) / "multivariate-overlap-traits.tsv"
        overlap_tree_path.write_text(
            dumps_newick(overlap_tree) + "\n", encoding="utf-8"
        )
        write_taxon_rows(
            overlap_table_path,
            columns=table.columns,
            rows=overlap_rows,
        )

        overlap_reports = [
            inspect_response_model(
                overlap_tree_path,
                overlap_table_path,
                response=response,
                predictors=predictors,
                taxon_column=table.taxon_column,
            )
            for response in responses
        ]
        raise_for_input_blockers(overlap_reports)

        analysis_taxa = shared_analysis_taxa(overlap_reports)
        excluded_taxa.extend(
            build_shared_taxon_exclusions(
                overlap_taxa=overlap_taxa,
                responses=responses,
                overlap_reports=overlap_reports,
            )
        )
        if len(analysis_taxa) < 2:
            raise ComparativeMethodError(
                "multivariate comparative regression does not retain enough shared complete-case taxa"
            )

        reduced_tree, _ = prune_tree_to_requested_taxa(overlap_tree_path, analysis_taxa)
        reduced_rows = [rows_by_taxon[taxon] for taxon in analysis_taxa]
        reduced_tree_path = Path(tmp_dir) / "multivariate-tree.nwk"
        reduced_table_path = Path(tmp_dir) / "multivariate-traits.tsv"
        reduced_tree_path.write_text(
            dumps_newick(reduced_tree) + "\n", encoding="utf-8"
        )
        write_taxon_rows(
            reduced_table_path,
            columns=table.columns,
            rows=reduced_rows,
        )

        final_input_reports = [
            inspect_response_model(
                reduced_tree_path,
                reduced_table_path,
                response=response,
                predictors=predictors,
                taxon_column=table.taxon_column,
            )
            for response in responses
        ]
        raise_for_input_blockers(final_input_reports)

        response_models = [
            run_pgls(
                reduced_tree_path,
                reduced_table_path,
                formula=build_response_formula(response, predictors),
                taxon_column=table.taxon_column,
                lambda_value=lambda_value,
            )
            for response in responses
        ]

    for model in response_models:
        model.tree_path = tree_path
        model.traits_path = traits_path

    response_model_rows = build_response_model_rows(response_models)
    coefficient_rows = build_response_coefficient_rows(response_models)
    covariance_rows = build_residual_covariance_rows(response_models)
    correlation_rows = build_residual_correlation_rows(covariance_rows)
    association_rows = build_residual_association_rows(response_models)
    covariance_diagnostics = build_residual_covariance_diagnostics(covariance_rows)
    warnings = build_multivariate_warnings(
        final_input_reports=final_input_reports,
        response_models=response_models,
        covariance_diagnostics=covariance_diagnostics,
    )
    return MultivariateComparativeRegressionReport(
        tree_path=tree_path,
        traits_path=traits_path,
        responses=list(responses),
        predictors=list(predictors),
        taxon_column=table.taxon_column,
        missing_value_policy=MULTIVARIATE_MISSING_VALUE_POLICY,
        numerical_tolerance=MULTIVARIATE_NUMERICAL_TOLERANCE,
        analysis_taxa=list(response_models[0].taxa),
        excluded_taxa=excluded_taxa,
        response_models=response_models,
        response_model_rows=response_model_rows,
        coefficient_rows=coefficient_rows,
        covariance_rows=covariance_rows,
        correlation_rows=correlation_rows,
        association_rows=association_rows,
        covariance_diagnostics=covariance_diagnostics,
        warnings=warnings,
    )


def build_initial_taxon_exclusions(
    *,
    responses: list[str],
    predictors: list[str],
    table_taxa: set[str],
    tree_taxa: set[str],
) -> list[MultivariateTaxonExclusion]:
    excluded_taxa: list[MultivariateTaxonExclusion] = []
    required_terms = [*responses, *predictors]
    for taxon in sorted(tree_taxa - table_taxa):
        excluded_taxa.append(
            MultivariateTaxonExclusion(
                taxon=taxon,
                reason="missing_from_trait_table",
                missing_columns=list(required_terms),
                blocking_responses=list(responses),
                details="taxon is present in the tree but absent from the trait table",
            )
        )
    for taxon in sorted(table_taxa - tree_taxa):
        excluded_taxa.append(
            MultivariateTaxonExclusion(
                taxon=taxon,
                reason="missing_from_tree",
                missing_columns=[],
                blocking_responses=[],
                details="taxon is present in the trait table but absent from the tree",
            )
        )
    return excluded_taxa


__all__ = [
    "run_multivariate_comparative_regression",
]
