from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path

from bijux_phylogenetics.comparative.pgls import (
    ComparativeFormulaSpecification,
    inspect_pgls_inputs,
    run_pgls,
)
from bijux_phylogenetics.comparative.regression import (
    summarize_phylogenetic_logistic,
)
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.phylo.pruning import prune_tree_to_requested_taxa
from bijux_phylogenetics.phylo.topology.tree import TreeNode
from bijux_phylogenetics.runtime.errors import ComparativeMethodError


@dataclass(slots=True)
class ComparativeResidualTaxonRow:
    """One analyzed taxon with fitted comparative residual diagnostics."""

    taxon: str
    observed_value: float
    fitted_value: float
    residual: float
    standardized_residual: float


@dataclass(slots=True)
class ComparativeResidualCladeRow:
    """Residual aggregation across one internal comparative clade."""

    clade_id: str
    node_label: str | None
    taxon_count: int
    taxa: list[str]
    mean_residual: float
    mean_abs_residual: float
    mean_standardized_residual: float
    mean_abs_standardized_residual: float
    max_abs_standardized_residual: float
    residual_sum_of_squares: float
    residual_sum_of_squares_share: float
    positive_residual_taxa: int
    negative_residual_taxa: int
    influence_score: float
    residual_heavy: bool
    rank: int


@dataclass(slots=True)
class ComparativeCladeResidualReport:
    """Clade-aware residual aggregation for one comparative model fit."""

    tree_path: Path
    traits_path: Path
    response: str
    formula: ComparativeFormulaSpecification
    model_family: str
    standardized_residual_method: str
    taxon_rows: list[ComparativeResidualTaxonRow]
    clade_rows: list[ComparativeResidualCladeRow]
    residual_heavy_clades: list[str]
    warnings: list[str]


def analyze_comparative_residual_clades(
    tree_path: Path,
    traits_path: Path,
    *,
    response: str | None = None,
    predictors: list[str] | None = None,
    formula: str | None = None,
    taxon_column: str | None = None,
    lambda_value: float | str = "estimate",
) -> ComparativeCladeResidualReport:
    """Aggregate fitted comparative residuals across internal clades."""
    input_report = inspect_pgls_inputs(
        tree_path,
        traits_path,
        response=response,
        predictors=predictors,
        formula=formula,
        taxon_column=taxon_column,
    )
    if not input_report.ready:
        raise ComparativeMethodError("; ".join(input_report.blockers))

    response_values = [row.response_value for row in input_report.model_matrix.rows]
    model_family = _shared_response_family(response_values)
    if model_family == "logistic":
        if lambda_value == "estimate":
            raise ComparativeMethodError(
                "clade-aware residual analysis requires a numeric lambda value for binary-response comparative models"
            )
        taxon_rows = _logistic_taxon_rows(
            summarize_phylogenetic_logistic(
                tree_path,
                traits_path,
                response=response,
                predictors=predictors,
                formula=formula,
                taxon_column=taxon_column,
                lambda_value=float(lambda_value),
            )
        )
        standardized_residual_method = "pearson-binomial-residual"
    else:
        taxon_rows = _pgls_taxon_rows(
            run_pgls(
                tree_path,
                traits_path,
                response=response,
                predictors=predictors,
                formula=formula,
                taxon_column=taxon_column,
                lambda_value=lambda_value,
            )
        )
        standardized_residual_method = "leveraged-gls-residual"

    if len(taxon_rows) < 3:
        raise ComparativeMethodError(
            "clade-aware residual analysis requires at least three analyzed taxa"
        )
    pruned_tree, _ = prune_tree_to_requested_taxa(
        tree_path,
        [row.taxon for row in taxon_rows],
    )
    clade_rows = _build_clade_rows(pruned_tree.root, taxon_rows=taxon_rows)
    _rank_clade_rows(clade_rows)
    residual_heavy_clades = [row.clade_id for row in clade_rows if row.residual_heavy]
    warnings: list[str] = []
    if residual_heavy_clades:
        warnings.append(
            "one or more internal clades concentrate elevated residual burden"
        )
    return ComparativeCladeResidualReport(
        tree_path=tree_path,
        traits_path=traits_path,
        response=input_report.formula.response,
        formula=input_report.formula,
        model_family=model_family,
        standardized_residual_method=standardized_residual_method,
        taxon_rows=taxon_rows,
        clade_rows=clade_rows,
        residual_heavy_clades=residual_heavy_clades,
        warnings=warnings,
    )


def write_comparative_residual_taxon_table(
    path: Path,
    report: ComparativeCladeResidualReport,
) -> Path:
    """Write one fitted residual ledger across analyzed taxa."""
    return write_taxon_rows(
        path,
        columns=[
            "taxon",
            "observed_value",
            "fitted_value",
            "residual",
            "standardized_residual",
        ],
        rows=[
            {
                "taxon": row.taxon,
                "observed_value": format(row.observed_value, ".15g"),
                "fitted_value": format(row.fitted_value, ".15g"),
                "residual": format(row.residual, ".15g"),
                "standardized_residual": format(row.standardized_residual, ".15g"),
            }
            for row in report.taxon_rows
        ],
    )


def write_comparative_residual_clade_table(
    path: Path,
    report: ComparativeCladeResidualReport,
) -> Path:
    """Write one clade-level residual aggregation ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "clade_id",
            "node_label",
            "taxon_count",
            "taxa",
            "mean_residual",
            "mean_abs_residual",
            "mean_standardized_residual",
            "mean_abs_standardized_residual",
            "max_abs_standardized_residual",
            "residual_sum_of_squares",
            "residual_sum_of_squares_share",
            "positive_residual_taxa",
            "negative_residual_taxa",
            "influence_score",
            "residual_heavy",
            "rank",
        ],
        rows=[
            {
                "clade_id": row.clade_id,
                "node_label": row.node_label or "",
                "taxon_count": str(row.taxon_count),
                "taxa": ",".join(row.taxa),
                "mean_residual": format(row.mean_residual, ".15g"),
                "mean_abs_residual": format(row.mean_abs_residual, ".15g"),
                "mean_standardized_residual": format(
                    row.mean_standardized_residual, ".15g"
                ),
                "mean_abs_standardized_residual": format(
                    row.mean_abs_standardized_residual, ".15g"
                ),
                "max_abs_standardized_residual": format(
                    row.max_abs_standardized_residual, ".15g"
                ),
                "residual_sum_of_squares": format(row.residual_sum_of_squares, ".15g"),
                "residual_sum_of_squares_share": format(
                    row.residual_sum_of_squares_share, ".15g"
                ),
                "positive_residual_taxa": str(row.positive_residual_taxa),
                "negative_residual_taxa": str(row.negative_residual_taxa),
                "influence_score": format(row.influence_score, ".15g"),
                "residual_heavy": str(row.residual_heavy).lower(),
                "rank": str(row.rank),
            }
            for row in report.clade_rows
        ],
    )


def _shared_response_family(response_values: list[float]) -> str:
    if all(
        math.isclose(value, round(value), abs_tol=1e-12) for value in response_values
    ) and {int(round(value)) for value in response_values} <= {0, 1}:
        return "logistic"
    return "pgls"


def _pgls_taxon_rows(report: object) -> list[ComparativeResidualTaxonRow]:
    fitted_by_taxon = {
        row.taxon: row for row in report.diagnostics.fitted_observed_rows
    }
    standardized_by_taxon = {
        row.taxon: row.standardized_residual for row in report.diagnostics.leverage_rows
    }
    return [
        ComparativeResidualTaxonRow(
            taxon=taxon,
            observed_value=fitted_by_taxon[taxon].observed,
            fitted_value=fitted_by_taxon[taxon].fitted,
            residual=fitted_by_taxon[taxon].residual,
            standardized_residual=standardized_by_taxon[taxon],
        )
        for taxon in report.taxa
    ]


def _logistic_taxon_rows(report: object) -> list[ComparativeResidualTaxonRow]:
    rows: list[ComparativeResidualTaxonRow] = []
    for fitted_row in report.fitted_rows:
        variance = max(
            fitted_row.fitted_probability * (1.0 - fitted_row.fitted_probability),
            1e-12,
        )
        rows.append(
            ComparativeResidualTaxonRow(
                taxon=fitted_row.taxon,
                observed_value=float(fitted_row.observed_response),
                fitted_value=fitted_row.fitted_probability,
                residual=fitted_row.residual,
                standardized_residual=fitted_row.residual / math.sqrt(variance),
            )
        )
    return rows


def _build_clade_rows(
    root: TreeNode,
    *,
    taxon_rows: list[ComparativeResidualTaxonRow],
) -> list[ComparativeResidualCladeRow]:
    row_by_taxon = {row.taxon: row for row in taxon_rows}
    total_residual_sum_of_squares = sum(row.residual**2 for row in taxon_rows)
    clade_rows: list[ComparativeResidualCladeRow] = []

    def visit(node: TreeNode, *, is_root: bool) -> list[str]:
        if node.is_leaf():
            return [node.name] if node.name is not None else []
        taxa: list[str] = []
        for child in node.children:
            taxa.extend(visit(child, is_root=False))
        ordered_taxa = sorted(taxa)
        if not is_root:
            member_rows = [row_by_taxon[taxon] for taxon in ordered_taxa]
            mean_residual = sum(row.residual for row in member_rows) / len(member_rows)
            mean_abs_residual = sum(abs(row.residual) for row in member_rows) / len(
                member_rows
            )
            mean_standardized_residual = sum(
                row.standardized_residual for row in member_rows
            ) / len(member_rows)
            mean_abs_standardized_residual = sum(
                abs(row.standardized_residual) for row in member_rows
            ) / len(member_rows)
            max_abs_standardized_residual = max(
                abs(row.standardized_residual) for row in member_rows
            )
            residual_sum_of_squares = sum(row.residual**2 for row in member_rows)
            residual_sum_of_squares_share = (
                residual_sum_of_squares / total_residual_sum_of_squares
                if total_residual_sum_of_squares > 0.0
                else 0.0
            )
            clade_rows.append(
                ComparativeResidualCladeRow(
                    clade_id="|".join(ordered_taxa),
                    node_label=node.name,
                    taxon_count=len(member_rows),
                    taxa=ordered_taxa,
                    mean_residual=mean_residual,
                    mean_abs_residual=mean_abs_residual,
                    mean_standardized_residual=mean_standardized_residual,
                    mean_abs_standardized_residual=mean_abs_standardized_residual,
                    max_abs_standardized_residual=max_abs_standardized_residual,
                    residual_sum_of_squares=residual_sum_of_squares,
                    residual_sum_of_squares_share=residual_sum_of_squares_share,
                    positive_residual_taxa=sum(
                        1 for row in member_rows if row.residual > 0.0
                    ),
                    negative_residual_taxa=sum(
                        1 for row in member_rows if row.residual < 0.0
                    ),
                    influence_score=(
                        residual_sum_of_squares_share * mean_abs_standardized_residual
                    ),
                    residual_heavy=(
                        mean_abs_standardized_residual >= 1.5
                        or (
                            max_abs_standardized_residual >= 2.0
                            and residual_sum_of_squares_share >= 0.2
                        )
                    ),
                    rank=0,
                )
            )
        return ordered_taxa

    visit(root, is_root=True)
    return clade_rows


def _rank_clade_rows(rows: list[ComparativeResidualCladeRow]) -> None:
    ranked_rows = sorted(
        rows,
        key=lambda row: (
            -row.influence_score,
            -row.residual_sum_of_squares_share,
            -row.mean_abs_standardized_residual,
            -row.taxon_count,
            row.clade_id,
        ),
    )
    for rank, row in enumerate(ranked_rows, start=1):
        row.rank = rank
