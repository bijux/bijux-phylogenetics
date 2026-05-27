from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from bijux_phylogenetics.phylo.pruning import prune_tree_to_requested_taxa
from bijux_phylogenetics.phylo.topology.clades import robinson_foulds_metrics
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

from .imported import (
    load_imported_distance_matrix,
    validate_imported_distance_matrix,
)
from .missing_distance_policy import apply_missing_distance_policy
from .models import (
    DistanceTaxonInfluenceReport,
    DistanceTaxonInfluenceRow,
    ImportedDistanceEntry,
    MissingDistancePolicy,
)
from .patristic_residuals import compute_patristic_residual_diagnostics
from .shared import _build_distance_tree_from_lookup, _pair_key

_ROUND_DIGITS = 12


def _defined_distance_lookup_from_entries(
    entries: list[ImportedDistanceEntry],
) -> dict[tuple[str, str], float]:
    defined_lookup: dict[tuple[str, str], float] = {}
    for entry in entries:
        if entry.left_identifier == entry.right_identifier:
            continue
        defined_lookup.setdefault(
            _pair_key(entry.left_identifier, entry.right_identifier),
            float(entry.distance),
        )
    return defined_lookup


def _raw_missing_pair_counts(
    identifiers: list[str],
    defined_lookup: dict[tuple[str, str], float],
) -> dict[str, int]:
    missing_pair_counts = {identifier: 0 for identifier in identifiers}
    for left_index, left_identifier in enumerate(identifiers):
        for right_index in range(left_index + 1, len(identifiers)):
            right_identifier = identifiers[right_index]
            if _pair_key(left_identifier, right_identifier) in defined_lookup:
                continue
            missing_pair_counts[left_identifier] += 1
            missing_pair_counts[right_identifier] += 1
    return missing_pair_counts


def _subset_defined_lookup(
    identifiers: list[str],
    defined_lookup: dict[tuple[str, str], float],
) -> dict[tuple[str, str], float]:
    retained = set(identifiers)
    return {
        pair: distance
        for pair, distance in defined_lookup.items()
        if pair[0] in retained and pair[1] in retained
    }


def _prune_reference_tree(
    reference_tree_path: Path,
    identifiers: list[str],
) -> PhyloTree:
    pruned_tree, pruning_report = prune_tree_to_requested_taxa(
        reference_tree_path,
        identifiers,
    )
    if pruning_report.absent_requested_taxa:
        missing = ", ".join(pruning_report.absent_requested_taxa)
        raise ValueError(
            "reference tree is missing one or more matrix taxa required for taxon influence analysis: "
            + missing
        )
    return pruned_tree


def _rooted_rf_metrics(
    tree: PhyloTree,
    reference_tree: PhyloTree,
) -> tuple[int, float]:
    shared_taxa = set(tree.tip_names) & set(reference_tree.tip_names)
    metrics = robinson_foulds_metrics(
        tree,
        reference_tree,
        shared_taxa,
        rf_mode="rooted",
    )
    return metrics.distance, round(metrics.normalized_distance, _ROUND_DIGITS)


def _resolve_tree_and_residuals(
    identifiers: list[str],
    defined_lookup: dict[tuple[str, str], float],
    *,
    method: str,
    missing_distance_policy: MissingDistancePolicy,
) -> tuple[PhyloTree, float]:
    resolved_lookup, _policy_report = apply_missing_distance_policy(
        identifiers,
        defined_lookup,
        policy=missing_distance_policy,
    )
    tree = _build_distance_tree_from_lookup(
        identifiers,
        resolved_lookup,
        method=method,
    )
    residuals = compute_patristic_residual_diagnostics(
        tree,
        identifiers,
        resolved_lookup,
    )
    return tree, residuals.residual_sum_squares


def _rank_taxon_influence_rows(
    rows: list[DistanceTaxonInfluenceRow],
) -> list[DistanceTaxonInfluenceRow]:
    ranked = sorted(
        rows,
        key=lambda row: (
            -row.residual_sum_squares_improvement,
            -row.rooted_robinson_foulds_improvement,
            -row.rooted_normalized_robinson_foulds_improvement,
            row.taxon,
        ),
    )
    return [
        DistanceTaxonInfluenceRow(
            taxon=row.taxon,
            retained_taxa=row.retained_taxa,
            raw_missing_pair_count=row.raw_missing_pair_count,
            baseline_residual_sum_squares=row.baseline_residual_sum_squares,
            leave_one_out_residual_sum_squares=row.leave_one_out_residual_sum_squares,
            residual_sum_squares_improvement=row.residual_sum_squares_improvement,
            baseline_rooted_robinson_foulds_distance=(
                row.baseline_rooted_robinson_foulds_distance
            ),
            leave_one_out_rooted_robinson_foulds_distance=(
                row.leave_one_out_rooted_robinson_foulds_distance
            ),
            rooted_robinson_foulds_improvement=row.rooted_robinson_foulds_improvement,
            baseline_rooted_normalized_robinson_foulds=(
                row.baseline_rooted_normalized_robinson_foulds
            ),
            leave_one_out_rooted_normalized_robinson_foulds=(
                row.leave_one_out_rooted_normalized_robinson_foulds
            ),
            rooted_normalized_robinson_foulds_improvement=(
                row.rooted_normalized_robinson_foulds_improvement
            ),
            topology_improved=row.topology_improved,
            residual_improved=row.residual_improved,
            influence_rank=index,
        )
        for index, row in enumerate(ranked, start=1)
    ]


def analyze_distance_taxon_influence_from_imported_distance_matrix(
    matrix_path: Path,
    reference_tree_path: Path,
    *,
    method: str,
    missing_distance_policy: MissingDistancePolicy = "reject",
) -> DistanceTaxonInfluenceReport:
    """Rank taxon influence by leave-one-out distance-tree fit and RF improvement."""
    validation = validate_imported_distance_matrix(matrix_path)
    if len(validation.identifiers) < 4:
        raise ValueError(
            "distance taxon influence analysis requires at least four taxa"
        )
    entries = load_imported_distance_matrix(matrix_path)
    defined_lookup = _defined_distance_lookup_from_entries(entries)
    raw_missing_pair_counts = _raw_missing_pair_counts(
        validation.identifiers,
        defined_lookup,
    )
    baseline_reference_tree = _prune_reference_tree(
        reference_tree_path,
        validation.identifiers,
    )
    baseline_tree, baseline_rss = _resolve_tree_and_residuals(
        validation.identifiers,
        defined_lookup,
        method=method,
        missing_distance_policy=missing_distance_policy,
    )
    (
        baseline_rooted_rf_distance,
        baseline_rooted_normalized_rf,
    ) = _rooted_rf_metrics(baseline_tree, baseline_reference_tree)

    rows: list[DistanceTaxonInfluenceRow] = []
    for taxon in validation.identifiers:
        retained_taxa = [
            identifier for identifier in validation.identifiers if identifier != taxon
        ]
        reduced_lookup = _subset_defined_lookup(retained_taxa, defined_lookup)
        reduced_reference_tree = _prune_reference_tree(reference_tree_path, retained_taxa)
        leave_one_out_tree, leave_one_out_rss = _resolve_tree_and_residuals(
            retained_taxa,
            reduced_lookup,
            method=method,
            missing_distance_policy=missing_distance_policy,
        )
        (
            leave_one_out_rooted_rf_distance,
            leave_one_out_rooted_normalized_rf,
        ) = _rooted_rf_metrics(leave_one_out_tree, reduced_reference_tree)
        rss_improvement = round(baseline_rss - leave_one_out_rss, _ROUND_DIGITS)
        rooted_rf_improvement = (
            baseline_rooted_rf_distance - leave_one_out_rooted_rf_distance
        )
        normalized_rf_improvement = round(
            baseline_rooted_normalized_rf - leave_one_out_rooted_normalized_rf,
            _ROUND_DIGITS,
        )
        rows.append(
            DistanceTaxonInfluenceRow(
                taxon=taxon,
                retained_taxa=retained_taxa,
                raw_missing_pair_count=raw_missing_pair_counts[taxon],
                baseline_residual_sum_squares=baseline_rss,
                leave_one_out_residual_sum_squares=leave_one_out_rss,
                residual_sum_squares_improvement=rss_improvement,
                baseline_rooted_robinson_foulds_distance=baseline_rooted_rf_distance,
                leave_one_out_rooted_robinson_foulds_distance=(
                    leave_one_out_rooted_rf_distance
                ),
                rooted_robinson_foulds_improvement=rooted_rf_improvement,
                baseline_rooted_normalized_robinson_foulds=baseline_rooted_normalized_rf,
                leave_one_out_rooted_normalized_robinson_foulds=(
                    leave_one_out_rooted_normalized_rf
                ),
                rooted_normalized_robinson_foulds_improvement=(
                    normalized_rf_improvement
                ),
                topology_improved=(
                    rooted_rf_improvement > 0 or normalized_rf_improvement > 0.0
                ),
                residual_improved=rss_improvement > 0.0,
                influence_rank=0,
            )
        )

    return DistanceTaxonInfluenceReport(
        source_path=matrix_path,
        source_kind="imported-distance-matrix",
        reference_tree_path=reference_tree_path,
        method=method,
        missing_distance_policy=missing_distance_policy,
        taxa=list(validation.identifiers),
        baseline_residual_sum_squares=baseline_rss,
        baseline_rooted_robinson_foulds_distance=baseline_rooted_rf_distance,
        baseline_rooted_normalized_robinson_foulds=baseline_rooted_normalized_rf,
        rows=_rank_taxon_influence_rows(rows),
    )


def write_distance_taxon_influence_table(
    path: Path,
    report: DistanceTaxonInfluenceReport,
) -> Path:
    """Write one ranked leave-one-out taxon influence table."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "influence_rank",
        "taxon",
        "retained_taxa",
        "raw_missing_pair_count",
        "baseline_residual_sum_squares",
        "leave_one_out_residual_sum_squares",
        "residual_sum_squares_improvement",
        "baseline_rooted_robinson_foulds_distance",
        "leave_one_out_rooted_robinson_foulds_distance",
        "rooted_robinson_foulds_improvement",
        "baseline_rooted_normalized_robinson_foulds",
        "leave_one_out_rooted_normalized_robinson_foulds",
        "rooted_normalized_robinson_foulds_improvement",
        "topology_improved",
        "residual_improved",
    ]
    lines = ["\t".join(columns)]
    for row in report.rows:
        lines.append(
            "\t".join(
                [
                    str(row.influence_rank),
                    row.taxon,
                    "|".join(row.retained_taxa),
                    str(row.raw_missing_pair_count),
                    format(row.baseline_residual_sum_squares, ".12g"),
                    format(row.leave_one_out_residual_sum_squares, ".12g"),
                    format(row.residual_sum_squares_improvement, ".12g"),
                    str(row.baseline_rooted_robinson_foulds_distance),
                    str(row.leave_one_out_rooted_robinson_foulds_distance),
                    str(row.rooted_robinson_foulds_improvement),
                    format(row.baseline_rooted_normalized_robinson_foulds, ".12g"),
                    format(row.leave_one_out_rooted_normalized_robinson_foulds, ".12g"),
                    format(row.rooted_normalized_robinson_foulds_improvement, ".12g"),
                    "true" if row.topology_improved else "false",
                    "true" if row.residual_improved else "false",
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_distance_taxon_influence_run_json(
    path: Path,
    report: DistanceTaxonInfluenceReport,
) -> Path:
    """Write the complete distance taxon influence payload as JSON."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "source_path": str(report.source_path),
        "source_kind": report.source_kind,
        "reference_tree_path": str(report.reference_tree_path),
        "method": report.method,
        "missing_distance_policy": report.missing_distance_policy,
        "taxa": report.taxa,
        "baseline_residual_sum_squares": report.baseline_residual_sum_squares,
        "baseline_rooted_robinson_foulds_distance": (
            report.baseline_rooted_robinson_foulds_distance
        ),
        "baseline_rooted_normalized_robinson_foulds": (
            report.baseline_rooted_normalized_robinson_foulds
        ),
        "rows": [asdict(row) for row in report.rows],
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def write_distance_taxon_influence_artifacts(
    out_dir: Path,
    report: DistanceTaxonInfluenceReport,
) -> dict[str, Path]:
    """Write governed distance taxon influence artifacts."""
    out_dir.mkdir(parents=True, exist_ok=True)
    influence_table_path = write_distance_taxon_influence_table(
        out_dir / "taxon_influence.tsv",
        report,
    )
    run_json_path = write_distance_taxon_influence_run_json(
        out_dir / "run.json",
        report,
    )
    return {
        "taxon_influence": influence_table_path,
        "run_json": run_json_path,
    }
