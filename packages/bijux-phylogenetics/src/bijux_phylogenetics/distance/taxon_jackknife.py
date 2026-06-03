from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.phylo.pruning import prune_tree_object_to_requested_taxa

from .imported import load_imported_distance_matrix, validate_imported_distance_matrix
from .models import (
    DistanceTaxonJackknifeReport,
    DistanceTaxonJackknifeRow,
    MissingDistancePolicy,
)
from .taxon_leave_one_out import (
    ROUND_DIGITS,
    compute_tree_residual_sum_squares,
    defined_distance_lookup_from_entries,
    resolve_distance_lookup,
    resolve_tree_and_residuals,
    rooted_rf_clade_differences,
    subset_defined_lookup,
)


def analyze_distance_taxon_jackknife_from_imported_distance_matrix(
    matrix_path: Path,
    *,
    method: str,
    missing_distance_policy: MissingDistancePolicy = "reject",
) -> DistanceTaxonJackknifeReport:
    """Rebuild one distance tree after removing each taxon and compare it to the pruned baseline tree."""
    validation = validate_imported_distance_matrix(matrix_path)
    if len(validation.identifiers) < 4:
        raise ValueError("distance taxon jackknife requires at least four taxa")
    entries = load_imported_distance_matrix(matrix_path)
    defined_lookup = defined_distance_lookup_from_entries(entries)
    baseline_tree, baseline_residual_sum_squares = resolve_tree_and_residuals(
        validation.identifiers,
        defined_lookup,
        method=method,
        missing_distance_policy=missing_distance_policy,
    )

    rows: list[DistanceTaxonJackknifeRow] = []
    for removed_taxon in validation.identifiers:
        retained_taxa = [
            identifier
            for identifier in validation.identifiers
            if identifier != removed_taxon
        ]
        reduced_lookup = subset_defined_lookup(retained_taxa, defined_lookup)
        resolved_reduced_lookup = resolve_distance_lookup(
            retained_taxa,
            reduced_lookup,
            missing_distance_policy=missing_distance_policy,
        )
        pruned_baseline_tree = prune_tree_object_to_requested_taxa(
            baseline_tree,
            retained_taxa,
        )
        pruned_baseline_residual_sum_squares = compute_tree_residual_sum_squares(
            pruned_baseline_tree,
            retained_taxa,
            resolved_reduced_lookup,
        )
        rebuilt_tree, rebuilt_residual_sum_squares = resolve_tree_and_residuals(
            retained_taxa,
            reduced_lookup,
            method=method,
            missing_distance_policy=missing_distance_policy,
        )
        (
            reference_only_clades,
            rebuilt_only_clades,
            affected_clades,
            rooted_robinson_foulds_distance,
            rooted_normalized_robinson_foulds,
        ) = rooted_rf_clade_differences(
            rebuilt_tree,
            pruned_baseline_tree,
        )
        rows.append(
            DistanceTaxonJackknifeRow(
                removed_taxon=removed_taxon,
                retained_taxa=retained_taxa,
                pruned_baseline_tree_newick=dumps_newick(pruned_baseline_tree),
                rebuilt_tree_newick=dumps_newick(rebuilt_tree),
                pruned_baseline_residual_sum_squares=pruned_baseline_residual_sum_squares,
                rebuilt_residual_sum_squares=rebuilt_residual_sum_squares,
                residual_sum_squares_change=round(
                    pruned_baseline_residual_sum_squares - rebuilt_residual_sum_squares,
                    ROUND_DIGITS,
                ),
                rooted_robinson_foulds_distance=rooted_robinson_foulds_distance,
                rooted_normalized_robinson_foulds=rooted_normalized_robinson_foulds,
                reference_only_clades=reference_only_clades,
                rebuilt_only_clades=rebuilt_only_clades,
                affected_clades=affected_clades,
                topology_changed=rooted_robinson_foulds_distance > 0,
            )
        )

    return DistanceTaxonJackknifeReport(
        source_path=matrix_path,
        source_kind="imported-distance-matrix",
        method=method,
        missing_distance_policy=missing_distance_policy,
        taxa=list(validation.identifiers),
        baseline_tree_newick=dumps_newick(baseline_tree),
        baseline_residual_sum_squares=baseline_residual_sum_squares,
        rows=rows,
    )


def write_distance_taxon_jackknife_table(
    path: Path,
    report: DistanceTaxonJackknifeReport,
) -> Path:
    """Write one leave-one-taxon-out jackknife ledger."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "removed_taxon",
        "retained_taxa",
        "rooted_robinson_foulds_distance",
        "rooted_normalized_robinson_foulds",
        "pruned_baseline_residual_sum_squares",
        "rebuilt_residual_sum_squares",
        "residual_sum_squares_change",
        "affected_clades",
        "topology_changed",
        "pruned_baseline_tree_newick",
        "rebuilt_tree_newick",
    ]
    lines = ["\t".join(columns)]
    for row in report.rows:
        lines.append(
            "\t".join(
                [
                    row.removed_taxon,
                    "|".join(row.retained_taxa),
                    str(row.rooted_robinson_foulds_distance),
                    format(row.rooted_normalized_robinson_foulds, ".12g"),
                    format(row.pruned_baseline_residual_sum_squares, ".12g"),
                    format(row.rebuilt_residual_sum_squares, ".12g"),
                    format(row.residual_sum_squares_change, ".12g"),
                    ";".join(row.affected_clades),
                    "true" if row.topology_changed else "false",
                    row.pruned_baseline_tree_newick,
                    row.rebuilt_tree_newick,
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_distance_taxon_jackknife_run_json(
    path: Path,
    report: DistanceTaxonJackknifeReport,
) -> Path:
    """Write the complete taxon jackknife payload as JSON."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "source_path": str(report.source_path),
        "source_kind": report.source_kind,
        "method": report.method,
        "missing_distance_policy": report.missing_distance_policy,
        "taxa": report.taxa,
        "baseline_tree_newick": report.baseline_tree_newick,
        "baseline_residual_sum_squares": report.baseline_residual_sum_squares,
        "rows": [asdict(row) for row in report.rows],
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return path


def write_distance_taxon_jackknife_artifacts(
    out_dir: Path,
    report: DistanceTaxonJackknifeReport,
) -> dict[str, Path]:
    """Write governed taxon jackknife artifacts."""
    out_dir.mkdir(parents=True, exist_ok=True)
    baseline_tree_path = out_dir / "baseline_tree.nwk"
    baseline_tree_path.write_text(report.baseline_tree_newick + "\n", encoding="utf-8")
    table_path = write_distance_taxon_jackknife_table(
        out_dir / "taxon_jackknife.tsv",
        report,
    )
    run_json_path = write_distance_taxon_jackknife_run_json(
        out_dir / "run.json",
        report,
    )
    return {
        "baseline_tree": baseline_tree_path,
        "taxon_jackknife": table_path,
        "run_json": run_json_path,
    }
