from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Never

from bijux_phylogenetics.comparative._math import invert_matrix, log_determinant
from bijux_phylogenetics.comparative.common import (
    build_brownian_covariance_matrix,
    load_comparative_dataset,
    tip_root_depths,
)
from bijux_phylogenetics.comparative.pgls.design import inspect_pgls_inputs
from bijux_phylogenetics.comparative.pgls.fitting import run_pgls
from bijux_phylogenetics.comparative.pgls.models import (
    ComparativeFormulaSpecification,
    PGLSResult,
)
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.phylo.branch_lengths.ultrametric import (
    summarize_ultrametric_tip_depths,
)
from bijux_phylogenetics.runtime.errors import ComparativeMethodError


@dataclass(slots=True)
class BrownianCovarianceRow:
    """One pairwise shared-path observation in a Brownian covariance matrix."""

    left_taxon: str
    right_taxon: str
    is_diagonal: bool
    shared_path_length: float
    left_root_depth: float
    right_root_depth: float


@dataclass(slots=True)
class BrownianCovariancePGLSReport:
    """Reviewer-facing Brownian-covariance audit plus the fitted PGLS model."""

    tree_path: Path
    traits_path: Path
    response: str
    formula: ComparativeFormulaSpecification
    taxon_count: int
    tree_is_ultrametric: bool
    minimum_root_to_tip_depth: float
    maximum_root_to_tip_depth: float
    minimum_branch_length: float
    maximum_branch_length: float
    raw_log_determinant: float
    positive_definite_before_stabilization: bool
    rows: list[BrownianCovarianceRow]
    model: PGLSResult


def summarize_brownian_covariance_pgls(
    tree_path: Path,
    traits_path: Path,
    *,
    response: str | None = None,
    predictors: list[str] | None = None,
    formula: str | None = None,
    taxon_column: str | None = None,
) -> BrownianCovariancePGLSReport:
    """Fit one PGLS model under fixed Brownian covariance and audit the raw matrix."""
    try:
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

        dataset = load_comparative_dataset(
            tree_path,
            traits_path,
            trait=input_report.formula_audit.response_column,
            taxon_column=taxon_column,
            minimum_taxa=len(input_report.encoded_columns) + 1,
            require_rooted=True,
            require_binary=False,
        )
        taxa = list(input_report.analysis_taxa)
        raw_covariance = build_brownian_covariance_matrix(dataset.tree, taxa)
        root_depths = tip_root_depths(dataset.tree, taxa)
        ultrametric_summary = summarize_ultrametric_tip_depths(
            root_depths,
            tolerance=1e-12,
        )
        minimum_branch_length, maximum_branch_length = _branch_length_range(
            dataset.tree_path
        )
        raw_log_determinant = _validate_raw_brownian_covariance(
            tree_path=tree_path,
            taxa=taxa,
            covariance_matrix=raw_covariance,
            minimum_branch_length=minimum_branch_length,
        )
        rows = _build_covariance_rows(
            taxa=taxa,
            covariance_matrix=raw_covariance,
            root_depths=root_depths,
        )
        return BrownianCovariancePGLSReport(
            tree_path=tree_path,
            traits_path=traits_path,
            response=input_report.response,
            formula=input_report.formula,
            taxon_count=len(taxa),
            tree_is_ultrametric=ultrametric_summary.ultrametric,
            minimum_root_to_tip_depth=ultrametric_summary.minimum_tip_depth,
            maximum_root_to_tip_depth=ultrametric_summary.maximum_tip_depth,
            minimum_branch_length=minimum_branch_length,
            maximum_branch_length=maximum_branch_length,
            raw_log_determinant=raw_log_determinant,
            positive_definite_before_stabilization=True,
            rows=rows,
            model=run_pgls(
                tree_path,
                traits_path,
                response=response,
                predictors=predictors,
                formula=formula,
                taxon_column=taxon_column,
                lambda_value=1.0,
            ),
        )
    except ComparativeMethodError as error:
        _reraise_brownian_input_error(tree_path, error)


def _reraise_brownian_input_error(
    tree_path: Path,
    error: ComparativeMethodError,
) -> Never:
    details = error.details or {}
    failure_reason = details.get("failure_reason")
    evidence = details.get("evidence", {})
    if failure_reason == "comparative_negative_branch_lengths":
        raise ComparativeMethodError(
            "Brownian covariance is invalid: tree contains negative branch lengths",
            details={
                "failure_reason": "brownian_covariance_negative_branch_lengths",
                "scientific_explanation": (
                    "Brownian shared-path covariance is invalid on a tree with negative branch lengths because shared evolutionary distance cannot be negative."
                ),
                "likely_causes": [
                    "the tree file contains one or more negative branch lengths",
                ],
                "actionable_fixes": [
                    "repair or re-estimate the tree so every non-root branch length is non-negative",
                    "inspect the tree for scaling or export errors that introduced negative lengths",
                ],
                "evidence": {
                    "tree_path": str(tree_path),
                    "minimum_branch_length": evidence.get("minimum_branch_length"),
                },
            },
        ) from error
    if failure_reason == "comparative_branch_lengths_incomplete":
        raise ComparativeMethodError(
            "Brownian covariance requires complete branch lengths",
            details={
                "failure_reason": "brownian_covariance_branch_lengths_incomplete",
                "scientific_explanation": (
                    "Brownian covariance needs complete numeric branch lengths because shared evolutionary path length defines the covariance."
                ),
                "likely_causes": [
                    "the tree was exported without complete branch lengths",
                    "one or more branches have blank or missing lengths",
                ],
                "actionable_fixes": [
                    "rerun tree inference or export with branch lengths preserved",
                    "inspect the tree file for missing branch-length fields",
                ],
                "evidence": {"tree_path": str(tree_path)},
            },
        ) from error
    raise error


def write_brownian_covariance_table(
    path: Path, report: BrownianCovariancePGLSReport
) -> Path:
    """Write the pairwise Brownian covariance ledger as TSV or CSV."""
    return write_taxon_rows(
        path,
        columns=[
            "left_taxon",
            "right_taxon",
            "is_diagonal",
            "shared_path_length",
            "left_root_depth",
            "right_root_depth",
            "tree_is_ultrametric",
            "minimum_root_to_tip_depth",
            "maximum_root_to_tip_depth",
            "minimum_branch_length",
            "maximum_branch_length",
            "raw_log_determinant",
            "positive_definite_before_stabilization",
        ],
        rows=[
            {
                "left_taxon": row.left_taxon,
                "right_taxon": row.right_taxon,
                "is_diagonal": "true" if row.is_diagonal else "false",
                "shared_path_length": format(row.shared_path_length, ".15g"),
                "left_root_depth": format(row.left_root_depth, ".15g"),
                "right_root_depth": format(row.right_root_depth, ".15g"),
                "tree_is_ultrametric": (
                    "true" if report.tree_is_ultrametric else "false"
                ),
                "minimum_root_to_tip_depth": format(
                    report.minimum_root_to_tip_depth, ".15g"
                ),
                "maximum_root_to_tip_depth": format(
                    report.maximum_root_to_tip_depth, ".15g"
                ),
                "minimum_branch_length": format(report.minimum_branch_length, ".15g"),
                "maximum_branch_length": format(report.maximum_branch_length, ".15g"),
                "raw_log_determinant": format(report.raw_log_determinant, ".15g"),
                "positive_definite_before_stabilization": (
                    "true" if report.positive_definite_before_stabilization else "false"
                ),
            }
            for row in report.rows
        ],
    )


def _branch_length_range(tree_path: Path) -> tuple[float, float]:
    from bijux_phylogenetics.io.trees import load_tree

    tree = load_tree(tree_path)
    branch_lengths = [
        node.branch_length
        for node in tree.iter_nodes()
        if node is not tree.root and node.branch_length is not None
    ]
    if not branch_lengths:
        raise ComparativeMethodError(
            "Brownian covariance requires complete branch lengths",
            details={
                "failure_reason": "brownian_covariance_branch_lengths_incomplete",
                "scientific_explanation": (
                    "Brownian covariance needs complete numeric branch lengths because shared evolutionary path length defines the covariance."
                ),
                "likely_causes": [
                    "the tree was exported without complete branch lengths",
                    "one or more branches have blank or missing lengths",
                ],
                "actionable_fixes": [
                    "rerun tree inference or export with branch lengths preserved",
                    "inspect the tree file for missing branch-length fields",
                ],
                "evidence": {"tree_path": str(tree_path)},
            },
        )
    return min(branch_lengths), max(branch_lengths)


def _validate_raw_brownian_covariance(
    *,
    tree_path: Path,
    taxa: list[str],
    covariance_matrix: list[list[float]],
    minimum_branch_length: float,
) -> float:
    if minimum_branch_length < 0.0:
        raise ComparativeMethodError(
            "Brownian covariance is invalid: tree contains negative branch lengths",
            details={
                "failure_reason": "brownian_covariance_negative_branch_lengths",
                "scientific_explanation": (
                    "Brownian shared-path covariance is invalid on a tree with negative branch lengths because shared evolutionary distance cannot be negative."
                ),
                "likely_causes": [
                    "the tree file contains one or more negative branch lengths",
                ],
                "actionable_fixes": [
                    "repair or re-estimate the tree so every non-root branch length is non-negative",
                    "inspect the tree for scaling or export errors that introduced negative lengths",
                ],
                "evidence": {
                    "tree_path": str(tree_path),
                    "minimum_branch_length": minimum_branch_length,
                },
            },
        )
    diagonal = [covariance_matrix[index][index] for index in range(len(taxa))]
    non_positive_taxa = [
        taxon for taxon, value in zip(taxa, diagonal, strict=True) if value <= 0.0
    ]
    if non_positive_taxa:
        raise ComparativeMethodError(
            "Brownian covariance is invalid: non-positive root-to-tip path length for "
            + ", ".join(non_positive_taxa)
        )
    try:
        invert_matrix(covariance_matrix)
        return log_determinant(covariance_matrix)
    except ValueError as error:
        raise ComparativeMethodError(
            "Brownian covariance is invalid before stabilization: "
            f"{tree_path.name}: {error}"
        ) from error


def _build_covariance_rows(
    *,
    taxa: list[str],
    covariance_matrix: list[list[float]],
    root_depths: dict[str, float],
) -> list[BrownianCovarianceRow]:
    rows: list[BrownianCovarianceRow] = []
    for row_index, left_taxon in enumerate(taxa):
        for column_index, right_taxon in enumerate(taxa):
            rows.append(
                BrownianCovarianceRow(
                    left_taxon=left_taxon,
                    right_taxon=right_taxon,
                    is_diagonal=row_index == column_index,
                    shared_path_length=covariance_matrix[row_index][column_index],
                    left_root_depth=root_depths[left_taxon],
                    right_root_depth=root_depths[right_taxon],
                )
            )
    return rows
