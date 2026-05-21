from __future__ import annotations

import csv
from dataclasses import dataclass
import math
from pathlib import Path

from bijux_phylogenetics.comparative._math import (
    invert_matrix,
    log_determinant,
    matrix_rank,
    symmetric_matrix_condition_number,
)
from bijux_phylogenetics.comparative.common import (
    build_brownian_covariance_matrix,
    tip_root_depths,
)
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.branch_lengths.ultrametric import (
    summarize_ultrametric_tip_depths,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import ComparativeMethodError

BROWNIAN_COVARIANCE_CONDITION_THRESHOLD = 1e12


@dataclass(slots=True)
class BrownianCovarianceRow:
    """One ordered shared-ancestry covariance observation."""

    left_taxon: str
    right_taxon: str
    is_diagonal: bool
    shared_ancestry_covariance: float
    left_root_depth: float
    right_root_depth: float


@dataclass(slots=True)
class BrownianCovarianceReport:
    """One explicit Brownian shared-ancestry covariance matrix report."""

    tree_path: Path
    taxa: list[str]
    tree_is_rooted: bool
    tree_is_ultrametric: bool
    minimum_root_to_tip_depth: float
    maximum_root_to_tip_depth: float
    minimum_branch_length: float
    maximum_branch_length: float
    matrix_dimension: int
    matrix_rank: int
    singular: bool
    near_singular: bool
    positive_definite: bool
    condition_number: float
    raw_log_determinant: float | None
    rows: list[BrownianCovarianceRow]


def summarize_brownian_covariance(
    tree_path: Path,
    *,
    taxa: list[str] | None = None,
) -> BrownianCovarianceReport:
    """Summarize one Brownian shared-ancestry covariance matrix for an explicit tip order."""
    tree = load_tree(tree_path)
    return summarize_brownian_covariance_from_tree(
        tree,
        tree_path=tree_path,
        taxa=taxa,
    )


def summarize_brownian_covariance_from_tree(
    tree: PhyloTree,
    *,
    tree_path: Path | None = None,
    taxa: list[str] | None = None,
) -> BrownianCovarianceReport:
    """Summarize one Brownian shared-ancestry covariance matrix from a native tree."""
    ordered_taxa = _resolve_taxa(tree, taxa)
    effective_tree_path = Path("<in-memory-tree>") if tree_path is None else tree_path
    minimum_branch_length, maximum_branch_length = _branch_length_range(
        tree, effective_tree_path
    )
    root_depths = tip_root_depths(tree, ordered_taxa)
    ultrametric_summary = summarize_ultrametric_tip_depths(
        root_depths,
        tolerance=1e-12,
    )
    covariance_matrix = build_brownian_covariance_matrix(tree, ordered_taxa)
    covariance_rank = matrix_rank(covariance_matrix, tolerance=1e-12)
    singular = covariance_rank < len(ordered_taxa)
    positive_definite, raw_log_determinant = _matrix_positive_definite_diagnostics(
        covariance_matrix
    )
    condition_number = math.inf
    if not singular:
        condition_number = symmetric_matrix_condition_number(covariance_matrix)
    near_singular = (
        singular or condition_number >= BROWNIAN_COVARIANCE_CONDITION_THRESHOLD
    )
    return BrownianCovarianceReport(
        tree_path=effective_tree_path,
        taxa=ordered_taxa,
        tree_is_rooted=_tree_is_rooted(tree),
        tree_is_ultrametric=ultrametric_summary.ultrametric,
        minimum_root_to_tip_depth=ultrametric_summary.minimum_tip_depth,
        maximum_root_to_tip_depth=ultrametric_summary.maximum_tip_depth,
        minimum_branch_length=minimum_branch_length,
        maximum_branch_length=maximum_branch_length,
        matrix_dimension=len(ordered_taxa),
        matrix_rank=covariance_rank,
        singular=singular,
        near_singular=near_singular,
        positive_definite=positive_definite,
        condition_number=condition_number,
        raw_log_determinant=raw_log_determinant,
        rows=_build_covariance_rows(
            taxa=ordered_taxa,
            covariance_matrix=covariance_matrix,
            root_depths=root_depths,
        ),
    )


def write_brownian_covariance_long_table(
    path: Path, report: BrownianCovarianceReport
) -> Path:
    """Write one long-form Brownian shared-ancestry covariance ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "left_taxon",
            "right_taxon",
            "is_diagonal",
            "shared_ancestry_covariance",
            "left_root_depth",
            "right_root_depth",
            "tree_is_rooted",
            "tree_is_ultrametric",
            "matrix_dimension",
            "matrix_rank",
            "singular",
            "near_singular",
            "positive_definite",
            "condition_number",
            "minimum_root_to_tip_depth",
            "maximum_root_to_tip_depth",
            "minimum_branch_length",
            "maximum_branch_length",
            "raw_log_determinant",
        ],
        rows=[
            {
                "left_taxon": row.left_taxon,
                "right_taxon": row.right_taxon,
                "is_diagonal": row.is_diagonal,
                "shared_ancestry_covariance": format(
                    row.shared_ancestry_covariance, ".15g"
                ),
                "left_root_depth": format(row.left_root_depth, ".15g"),
                "right_root_depth": format(row.right_root_depth, ".15g"),
                "tree_is_rooted": report.tree_is_rooted,
                "tree_is_ultrametric": report.tree_is_ultrametric,
                "matrix_dimension": report.matrix_dimension,
                "matrix_rank": report.matrix_rank,
                "singular": report.singular,
                "near_singular": report.near_singular,
                "positive_definite": report.positive_definite,
                "condition_number": (
                    "inf"
                    if math.isinf(report.condition_number)
                    else format(report.condition_number, ".15g")
                ),
                "minimum_root_to_tip_depth": format(
                    report.minimum_root_to_tip_depth, ".15g"
                ),
                "maximum_root_to_tip_depth": format(
                    report.maximum_root_to_tip_depth, ".15g"
                ),
                "minimum_branch_length": format(report.minimum_branch_length, ".15g"),
                "maximum_branch_length": format(report.maximum_branch_length, ".15g"),
                "raw_log_determinant": (
                    ""
                    if report.raw_log_determinant is None
                    else format(report.raw_log_determinant, ".15g")
                ),
            }
            for row in report.rows
        ],
    )


def write_brownian_covariance_matrix_table(
    path: Path, report: BrownianCovarianceReport
) -> Path:
    """Write one wide Brownian covariance matrix with explicit tip order."""
    delimiter = "," if path.suffix == ".csv" else "\t"
    matrix_lookup = {
        (row.left_taxon, row.right_taxon): row.shared_ancestry_covariance
        for row in report.rows
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter=delimiter)
        writer.writerow(["taxon", *report.taxa])
        for left_taxon in report.taxa:
            writer.writerow(
                [
                    left_taxon,
                    *[
                        format(matrix_lookup[(left_taxon, right_taxon)], ".15g")
                        for right_taxon in report.taxa
                    ],
                ]
            )
    return path


def _resolve_taxa(tree: PhyloTree, taxa: list[str] | None) -> list[str]:
    ordered_taxa = list(tree.tip_names if taxa is None else taxa)
    duplicate_taxa = sorted(
        {taxon for taxon in ordered_taxa if ordered_taxa.count(taxon) > 1}
    )
    if duplicate_taxa:
        raise ComparativeMethodError(
            "Brownian covariance requires each requested taxon at most once",
            details={
                "failure_reason": "brownian_covariance_taxa_duplicated",
                "scientific_explanation": (
                    "A Brownian covariance matrix must use one unique row and one unique column per retained tip."
                ),
                "likely_causes": [
                    "the requested taxon order repeats one or more names",
                ],
                "actionable_fixes": [
                    "deduplicate the requested taxon order before building the covariance matrix",
                ],
                "evidence": {
                    "duplicate_taxa": duplicate_taxa,
                },
            },
        )
    tree_taxa = set(tree.tip_names)
    missing_taxa = [taxon for taxon in ordered_taxa if taxon not in tree_taxa]
    if missing_taxa:
        raise ComparativeMethodError(
            "Brownian covariance requires every requested taxon to be present in the tree",
            details={
                "failure_reason": "brownian_covariance_taxa_missing",
                "scientific_explanation": (
                    "The shared-ancestry covariance matrix can only be computed for taxa that are present in the phylogeny."
                ),
                "likely_causes": [
                    "the requested taxon order includes names absent from the tree",
                ],
                "actionable_fixes": [
                    "align the requested taxon order to the tree tip labels before rerunning the covariance summary",
                ],
                "evidence": {
                    "missing_taxa": missing_taxa,
                    "available_taxa": list(tree.tip_names),
                },
            },
        )
    return ordered_taxa


def _tree_is_rooted(tree: PhyloTree) -> bool:
    return tree.rooted is True or len(tree.root.children) == 2


def _branch_length_range(tree: PhyloTree, tree_path: Path) -> tuple[float, float]:
    branch_lengths = [
        node.branch_length
        for node in tree.iter_nodes()
        if node is not tree.root and node.branch_length is not None
    ]
    branch_count = sum(1 for node in tree.iter_nodes() if node is not tree.root)
    if len(branch_lengths) != branch_count:
        raise ComparativeMethodError(
            "Brownian covariance requires complete branch lengths",
            details={
                "failure_reason": "brownian_covariance_branch_lengths_incomplete",
                "scientific_explanation": (
                    "Brownian shared-ancestry covariance needs complete numeric branch lengths because shared evolutionary path length defines every covariance entry."
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
    minimum_branch_length = min(branch_lengths)
    if minimum_branch_length < 0.0:
        raise ComparativeMethodError(
            "Brownian covariance is invalid: tree contains negative branch lengths",
            details={
                "failure_reason": "brownian_covariance_negative_branch_lengths",
                "scientific_explanation": (
                    "Brownian shared-ancestry covariance is invalid on a tree with negative branch lengths because shared evolutionary distance cannot be negative."
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
    return minimum_branch_length, max(branch_lengths)


def _matrix_positive_definite_diagnostics(
    covariance_matrix: list[list[float]],
) -> tuple[bool, float | None]:
    try:
        invert_matrix(covariance_matrix)
        return True, log_determinant(covariance_matrix)
    except ValueError:
        return False, None


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
                    shared_ancestry_covariance=covariance_matrix[row_index][
                        column_index
                    ],
                    left_root_depth=root_depths[left_taxon],
                    right_root_depth=root_depths[right_taxon],
                )
            )
    return rows
