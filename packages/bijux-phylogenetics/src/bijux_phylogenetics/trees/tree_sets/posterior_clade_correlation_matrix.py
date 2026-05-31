from __future__ import annotations

import csv
import math
from pathlib import Path

from bijux_phylogenetics.phylo.topology.clades import informative_rooted_clades

from .contracts import PosteriorCladeCorrelationReport, PosteriorCladeCorrelationRow
from .inventory import _analyze_tree_set, _require_exact_taxa
from .topology import _clades_conflict, _format_clade


def _compatibility_relation(
    left: frozenset[str],
    right: frozenset[str],
) -> tuple[str, str]:
    if left == right:
        return ("identical", "same-clade")
    if _clades_conflict(left, right):
        return ("conflict", "overlap-without-containment")
    if left <= right or right <= left:
        return ("compatible", "nested")
    return ("compatible", "disjoint")


def _binary_correlation(
    *,
    left_count: int,
    right_count: int,
    joint_count: int,
    total_tree_count: int,
) -> float | None:
    if total_tree_count <= 0:
        raise ValueError("total_tree_count must be positive")
    if left_count == right_count == joint_count:
        return 1.0
    left_frequency = left_count / total_tree_count
    right_frequency = right_count / total_tree_count
    joint_frequency = joint_count / total_tree_count
    denominator = math.sqrt(
        left_frequency
        * (1.0 - left_frequency)
        * right_frequency
        * (1.0 - right_frequency)
    )
    if denominator <= 0.0:
        return None
    return round(
        (joint_frequency - (left_frequency * right_frequency)) / denominator,
        15,
    )


def compute_posterior_clade_correlation_matrix(
    path: Path,
) -> PosteriorCladeCorrelationReport:
    """Compute clade co-occurrence and binary correlation across posterior trees."""
    analysis = _analyze_tree_set(path)
    exact_taxa = _require_exact_taxa(analysis)
    exact_taxa_set = set(exact_taxa)
    clade_counts = analysis.clade_counts or {}
    tree_count = len(analysis.trees)
    tree_clade_sets = [
        informative_rooted_clades(tree, exact_taxa_set) for tree in analysis.trees
    ]
    ordered_clades = sorted(
        clade_counts,
        key=lambda clade: (-clade_counts[clade], _format_clade(clade)),
    )
    rows: list[PosteriorCladeCorrelationRow] = []
    formatted_clades = [_format_clade(clade) for clade in ordered_clades]

    for left_index, left_clade in enumerate(ordered_clades):
        left_count = clade_counts[left_clade]
        left_frequency = round(left_count / tree_count, 15)
        for right_clade in ordered_clades[left_index:]:
            right_count = clade_counts[right_clade]
            right_frequency = round(right_count / tree_count, 15)
            joint_count = sum(
                1
                for tree_clades in tree_clade_sets
                if left_clade in tree_clades and right_clade in tree_clades
            )
            joint_frequency = round(joint_count / tree_count, 15)
            expected_joint_frequency = round(left_frequency * right_frequency, 15)
            relation, reason = _compatibility_relation(left_clade, right_clade)
            rows.append(
                PosteriorCladeCorrelationRow(
                    left_clade=_format_clade(left_clade),
                    right_clade=_format_clade(right_clade),
                    compatibility_relation=relation,
                    compatibility_reason=reason,
                    left_tree_count=left_count,
                    right_tree_count=right_count,
                    left_frequency=left_frequency,
                    right_frequency=right_frequency,
                    cooccurrence_tree_count=joint_count,
                    cooccurrence_frequency=joint_frequency,
                    expected_cooccurrence_frequency=expected_joint_frequency,
                    binary_correlation=_binary_correlation(
                        left_count=left_count,
                        right_count=right_count,
                        joint_count=joint_count,
                        total_tree_count=tree_count,
                    ),
                )
            )

    return PosteriorCladeCorrelationReport(
        path=analysis.path,
        tree_count=tree_count,
        processing=analysis.processing,
        shared_taxa=exact_taxa,
        clade_count=len(formatted_clades),
        pair_count=len(rows),
        clade_order=formatted_clades,
        rows=rows,
    )


def write_posterior_clade_correlation_pair_table(
    path: Path,
    report: PosteriorCladeCorrelationReport,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "left_clade",
                "right_clade",
                "compatibility_relation",
                "compatibility_reason",
                "left_tree_count",
                "right_tree_count",
                "left_frequency",
                "right_frequency",
                "cooccurrence_tree_count",
                "cooccurrence_frequency",
                "expected_cooccurrence_frequency",
                "binary_correlation",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.rows:
            writer.writerow(
                {
                    "left_clade": row.left_clade,
                    "right_clade": row.right_clade,
                    "compatibility_relation": row.compatibility_relation,
                    "compatibility_reason": row.compatibility_reason,
                    "left_tree_count": row.left_tree_count,
                    "right_tree_count": row.right_tree_count,
                    "left_frequency": format(row.left_frequency, ".15g"),
                    "right_frequency": format(row.right_frequency, ".15g"),
                    "cooccurrence_tree_count": row.cooccurrence_tree_count,
                    "cooccurrence_frequency": format(
                        row.cooccurrence_frequency,
                        ".15g",
                    ),
                    "expected_cooccurrence_frequency": format(
                        row.expected_cooccurrence_frequency,
                        ".15g",
                    ),
                    "binary_correlation": (
                        ""
                        if row.binary_correlation is None
                        else format(row.binary_correlation, ".15g")
                    ),
                }
            )
    return path


def write_posterior_clade_correlation_matrix_table(
    path: Path,
    report: PosteriorCladeCorrelationReport,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    row_lookup: dict[tuple[str, str], float | None] = {}
    for row in report.rows:
        row_lookup[(row.left_clade, row.right_clade)] = row.binary_correlation
        row_lookup[(row.right_clade, row.left_clade)] = row.binary_correlation
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["clade", *report.clade_order])
        for left_clade in report.clade_order:
            matrix_row = [left_clade]
            for right_clade in report.clade_order:
                value = row_lookup[(left_clade, right_clade)]
                matrix_row.append("" if value is None else format(value, ".15g"))
            writer.writerow(matrix_row)
    return path


def write_posterior_clade_correlation_artifacts(
    out_dir: Path,
    report: PosteriorCladeCorrelationReport,
) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    matrix_path = write_posterior_clade_correlation_matrix_table(
        out_dir / "posterior-clade-correlation-matrix.tsv",
        report,
    )
    pair_path = write_posterior_clade_correlation_pair_table(
        out_dir / "posterior-clade-correlation-pairs.tsv",
        report,
    )
    return {
        "posterior_clade_correlation_matrix_path": matrix_path,
        "posterior_clade_correlation_pair_path": pair_path,
    }
