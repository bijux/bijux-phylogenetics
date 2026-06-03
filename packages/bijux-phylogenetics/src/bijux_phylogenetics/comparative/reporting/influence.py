from __future__ import annotations

from bijux_phylogenetics.comparative.reporting.contracts import (
    ComparativeInfluenceReport,
    ComparativePredictorInfluenceRow,
    ComparativeTaxonInfluenceRow,
)
from bijux_phylogenetics.comparative.reporting.snapshot import (
    build_comparative_model_snapshot,
)


def build_trait_influence_report(
    tree_path,
    traits_path,
    *,
    response: str | None = None,
    predictors: list[str] | None = None,
    formula: str | None = None,
    taxon_column: str | None = None,
    lambda_value: float | str = "estimate",
) -> ComparativeInfluenceReport:
    """Identify predictor terms and taxa driving one comparative result."""
    snapshot = build_comparative_model_snapshot(
        tree_path,
        traits_path,
        response=response,
        predictors=predictors,
        formula=formula,
        taxon_column=taxon_column,
        lambda_value=lambda_value,
    )
    predictor_rows = sorted(
        [
            ComparativePredictorInfluenceRow(
                term=coefficient.name,
                estimate=coefficient.estimate,
                test_statistic=coefficient.test_statistic,
                p_value=coefficient.p_value,
                significant=coefficient.p_value <= 0.05,
            )
            for coefficient in snapshot.pgls_model.coefficients
            if coefficient.name != "intercept"
        ],
        key=lambda row: (abs(row.test_statistic), abs(row.estimate), row.term),
        reverse=True,
    )
    leverage_by_taxon = {
        row.taxon: row for row in snapshot.pgls_model.diagnostics.leverage_rows
    }
    sensitivity_by_taxon = {row.dropped_taxon: row for row in snapshot.sensitivity.rows}
    taxon_rows = [
        ComparativeTaxonInfluenceRow(
            taxon=taxon,
            leverage=leverage_by_taxon[taxon].leverage,
            standardized_residual=leverage_by_taxon[taxon].standardized_residual,
            sensitivity_delta_log_likelihood=sensitivity_by_taxon[
                taxon
            ].delta_log_likelihood,
            sensitivity_delta_primary_parameter=sensitivity_by_taxon[
                taxon
            ].delta_primary_parameter,
            influence_score=(
                abs(leverage_by_taxon[taxon].standardized_residual)
                + leverage_by_taxon[taxon].leverage
                + abs(sensitivity_by_taxon[taxon].delta_log_likelihood)
            ),
        )
        for taxon in snapshot.pgls_model.taxa
    ]
    taxon_rows.sort(key=lambda row: (row.influence_score, row.taxon), reverse=True)
    top_predictor_terms = [row.term for row in predictor_rows[:3]]
    top_taxa = [row.taxon for row in taxon_rows[:3]]
    warnings: list[str] = []
    if any(row.significant for row in predictor_rows):
        warnings.append(
            "one or more predictor terms show nominal coefficient-level significance"
        )
    if any(abs(row.standardized_residual) >= 2.0 for row in taxon_rows):
        warnings.append(
            "one or more taxa contribute large standardized residuals to the fitted model"
        )
    if any(
        row.leverage
        >= (2.0 * len(snapshot.pgls_model.encoded_columns))
        / len(snapshot.pgls_model.taxa)
        for row in taxon_rows
    ):
        warnings.append(
            "one or more taxa have high leverage relative to model complexity"
        )
    return ComparativeInfluenceReport(
        tree_path=tree_path,
        traits_path=traits_path,
        response=snapshot.response,
        selected_model=snapshot.model_comparison.better_model,
        predictor_rows=predictor_rows,
        taxon_rows=taxon_rows,
        top_predictor_terms=top_predictor_terms,
        top_taxa=top_taxa,
        warnings=warnings,
    )
