from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.topology import summarize_tree_tip_distances
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

from .fixed_topology_policy import validate_fixed_topology_distance_input
from .imported import (
    _distance_lookup_from_imported,
    load_imported_distance_matrix,
    validate_imported_distance_matrix,
)
from .models import PatristicResidualDiagnosticsReport, PatristicResidualRow

_RESIDUAL_ROUND_DIGITS = 12


def compute_patristic_residual_diagnostics(
    tree: PhyloTree,
    identifiers: list[str],
    distance_lookup: dict[tuple[str, str], float],
    *,
    matrix_path: Path | None = None,
    tree_path: Path | None = None,
) -> PatristicResidualDiagnosticsReport:
    """Compare one observed distance matrix against one tree's patristic distances."""
    validate_fixed_topology_distance_input(tree, identifiers, distance_lookup)
    tip_distance_report = summarize_tree_tip_distances(tree)
    fitted_lookup = {
        (row.left_identifier, row.right_identifier): row.distance
        for row in tip_distance_report.pairs
    }
    rows: list[PatristicResidualRow] = []
    for left_index, left_identifier in enumerate(identifiers):
        for right_index in range(left_index + 1, len(identifiers)):
            right_identifier = identifiers[right_index]
            observed_distance = float(
                distance_lookup[(left_identifier, right_identifier)]
            )
            fitted_distance = float(fitted_lookup[(left_identifier, right_identifier)])
            residual = round(
                observed_distance - fitted_distance,
                _RESIDUAL_ROUND_DIGITS,
            )
            rows.append(
                PatristicResidualRow(
                    left_identifier=left_identifier,
                    right_identifier=right_identifier,
                    observed_distance=round(
                        observed_distance,
                        _RESIDUAL_ROUND_DIGITS,
                    ),
                    fitted_distance=round(
                        fitted_distance,
                        _RESIDUAL_ROUND_DIGITS,
                    ),
                    residual=residual,
                    absolute_residual=round(abs(residual), _RESIDUAL_ROUND_DIGITS),
                    rank=0,
                )
            )
    ranked_rows = sorted(
        rows,
        key=lambda row: (
            -row.absolute_residual,
            row.left_identifier,
            row.right_identifier,
        ),
    )
    ranked_rows = [
        PatristicResidualRow(
            left_identifier=row.left_identifier,
            right_identifier=row.right_identifier,
            observed_distance=row.observed_distance,
            fitted_distance=row.fitted_distance,
            residual=row.residual,
            absolute_residual=row.absolute_residual,
            rank=index,
        )
        for index, row in enumerate(ranked_rows, start=1)
    ]
    residual_sum_squares = round(
        sum(row.residual * row.residual for row in ranked_rows),
        _RESIDUAL_ROUND_DIGITS,
    )
    max_absolute_residual = round(
        max((row.absolute_residual for row in ranked_rows), default=0.0),
        _RESIDUAL_ROUND_DIGITS,
    )
    return PatristicResidualDiagnosticsReport(
        matrix_path=matrix_path,
        tree_path=tree_path,
        taxa=list(identifiers),
        pair_count=len(ranked_rows),
        residual_sum_squares=residual_sum_squares,
        max_absolute_residual=max_absolute_residual,
        rows=ranked_rows,
    )


def compute_patristic_residual_diagnostics_from_imported_distance_matrix(
    matrix_path: Path,
    tree_path: Path,
) -> PatristicResidualDiagnosticsReport:
    """Compare one imported distance matrix against one on-disk tree."""
    entries = load_imported_distance_matrix(matrix_path)
    validation = validate_imported_distance_matrix(matrix_path)
    distance_lookup, _missing_distance_policy_report = _distance_lookup_from_imported(
        validation,
        entries,
    )
    tree = load_tree(tree_path)
    return compute_patristic_residual_diagnostics(
        tree,
        validation.identifiers,
        distance_lookup,
        matrix_path=matrix_path,
        tree_path=tree_path,
    )


def write_patristic_residual_table(
    path: Path,
    report: PatristicResidualDiagnosticsReport,
) -> Path:
    """Write one ranked long-form patristic residual table."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "left_identifier",
        "right_identifier",
        "observed_distance",
        "fitted_distance",
        "residual",
        "absolute_residual",
        "rank",
    ]
    lines = ["\t".join(columns)]
    for row in report.rows:
        lines.append(
            "\t".join(
                [
                    row.left_identifier,
                    row.right_identifier,
                    format(row.observed_distance, ".12g"),
                    format(row.fitted_distance, ".12g"),
                    format(row.residual, ".12g"),
                    format(row.absolute_residual, ".12g"),
                    str(row.rank),
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_patristic_residual_run_json(
    path: Path,
    report: PatristicResidualDiagnosticsReport,
) -> Path:
    """Write the complete patristic residual diagnostics payload as JSON."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "matrix_path": None if report.matrix_path is None else str(report.matrix_path),
        "tree_path": None if report.tree_path is None else str(report.tree_path),
        "taxa": report.taxa,
        "pair_count": report.pair_count,
        "residual_sum_squares": report.residual_sum_squares,
        "max_absolute_residual": report.max_absolute_residual,
        "rows": [asdict(row) for row in report.rows],
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return path


def write_patristic_residual_artifacts(
    out_dir: Path,
    report: PatristicResidualDiagnosticsReport,
) -> dict[str, Path]:
    """Write governed patristic residual diagnostics artifacts."""
    out_dir.mkdir(parents=True, exist_ok=True)
    residual_table_path = write_patristic_residual_table(
        out_dir / "distance_residuals.tsv",
        report,
    )
    run_json_path = write_patristic_residual_run_json(out_dir / "run.json", report)
    return {
        "distance_residuals": residual_table_path,
        "run_json": run_json_path,
    }
