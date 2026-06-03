from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from .imported import (
    load_imported_distance_matrix,
    validate_imported_distance_matrix,
)
from .models import (
    DistanceTaxonInfluenceReport,
    DistanceTaxonInfluenceRow,
    MissingDistancePolicy,
)
from .shared import _pair_key
from .taxon_leave_one_out import (
    ROUND_DIGITS,
    defined_distance_lookup_from_entries,
    prune_tree_to_identifiers_or_raise,
    resolve_tree_and_residuals,
    rooted_rf_metrics,
    subset_defined_lookup,
)


def _raw_missing_pair_counts(
    identifiers: list[str],
    defined_lookup: dict[tuple[str, str], float],
) -> dict[str, int]:
    missing_pair_counts = dict.fromkeys(identifiers, 0)
    for left_index, left_identifier in enumerate(identifiers):
        for right_index in range(left_index + 1, len(identifiers)):
            right_identifier = identifiers[right_index]
            if _pair_key(left_identifier, right_identifier) in defined_lookup:
                continue
            missing_pair_counts[left_identifier] += 1
            missing_pair_counts[right_identifier] += 1
    return missing_pair_counts


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
    defined_lookup = defined_distance_lookup_from_entries(entries)
    raw_missing_pair_counts = _raw_missing_pair_counts(
        validation.identifiers,
        defined_lookup,
    )
    baseline_reference_tree = prune_tree_to_identifiers_or_raise(
        reference_tree_path,
        validation.identifiers,
        missing_taxa_message_prefix=(
            "reference tree is missing one or more matrix taxa required for taxon influence analysis: "
        ),
    )
    baseline_tree, baseline_rss = resolve_tree_and_residuals(
        validation.identifiers,
        defined_lookup,
        method=method,
        missing_distance_policy=missing_distance_policy,
    )
    (
        baseline_rooted_rf_distance,
        baseline_rooted_normalized_rf,
    ) = rooted_rf_metrics(baseline_tree, baseline_reference_tree)

    rows: list[DistanceTaxonInfluenceRow] = []
    for taxon in validation.identifiers:
        retained_taxa = [
            identifier for identifier in validation.identifiers if identifier != taxon
        ]
        reduced_lookup = subset_defined_lookup(retained_taxa, defined_lookup)
        reduced_reference_tree = prune_tree_to_identifiers_or_raise(
            reference_tree_path,
            retained_taxa,
            missing_taxa_message_prefix=(
                "reference tree is missing one or more matrix taxa required for taxon influence analysis: "
            ),
        )
        leave_one_out_tree, leave_one_out_rss = resolve_tree_and_residuals(
            retained_taxa,
            reduced_lookup,
            method=method,
            missing_distance_policy=missing_distance_policy,
        )
        (
            leave_one_out_rooted_rf_distance,
            leave_one_out_rooted_normalized_rf,
        ) = rooted_rf_metrics(leave_one_out_tree, reduced_reference_tree)
        rss_improvement = round(baseline_rss - leave_one_out_rss, ROUND_DIGITS)
        rooted_rf_improvement = (
            baseline_rooted_rf_distance - leave_one_out_rooted_rf_distance
        )
        normalized_rf_improvement = round(
            baseline_rooted_normalized_rf - leave_one_out_rooted_normalized_rf,
            ROUND_DIGITS,
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
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
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
