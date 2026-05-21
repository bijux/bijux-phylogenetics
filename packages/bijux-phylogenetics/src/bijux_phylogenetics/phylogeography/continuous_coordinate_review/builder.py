from __future__ import annotations

from pathlib import Path
import tempfile

from bijux_phylogenetics.ancestral.continuous import (
    reconstruct_continuous_ancestral_states,
)
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.phylo.pruning import prune_tree_to_requested_taxa
from bijux_phylogenetics.runtime.errors import AncestralReconstructionError

from .contracts import PhylogeographicCoordinateReport
from .inputs import prepare_coordinate_dataset, write_filtered_coordinate_table
from .movement_analysis import (
    build_branch_rows,
    build_estimate_rows,
    build_summary,
    review_branch_outliers,
)


def summarize_continuous_phylogeography(
    tree_path: Path,
    table_path: Path,
    *,
    latitude_column: str,
    longitude_column: str,
    taxon_column: str | None = None,
    model: str = "brownian",
    alpha: float = 1.0,
) -> PhylogeographicCoordinateReport:
    prepared = prepare_coordinate_dataset(
        tree_path,
        table_path,
        latitude_column=latitude_column,
        longitude_column=longitude_column,
        taxon_column=taxon_column,
    )
    if len(prepared.included_taxa) < 2:
        raise AncestralReconstructionError(
            "continuous phylogeography requires at least two taxa with usable coordinates"
        )

    with tempfile.TemporaryDirectory(prefix="bijux-phylogeography-") as temp_dir:
        temp_root = Path(temp_dir)
        pruned_tree, _ = prune_tree_to_requested_taxa(tree_path, prepared.included_taxa)
        pruned_tree_path = temp_root / "analysis-tree.nwk"
        pruned_tree_path.write_text(f"{dumps_newick(pruned_tree)}\n", encoding="utf-8")
        filtered_table_path = temp_root / "coordinates.tsv"
        write_filtered_coordinate_table(filtered_table_path, prepared)
        latitude_report = reconstruct_continuous_ancestral_states(
            pruned_tree_path,
            filtered_table_path,
            trait=latitude_column,
            taxon_column=prepared.taxon_column,
            model=model,
            alpha=alpha,
        )
        longitude_report = reconstruct_continuous_ancestral_states(
            pruned_tree_path,
            filtered_table_path,
            trait=longitude_column,
            taxon_column=prepared.taxon_column,
            model=model,
            alpha=alpha,
        )

    warnings = list(
        dict.fromkeys([*latitude_report.warnings, *longitude_report.warnings])
    )
    estimate_rows = build_estimate_rows(latitude_report, longitude_report)
    branch_rows = build_branch_rows(latitude_report, longitude_report)
    outlier_rows, reviewed_branch_rows = review_branch_outliers(branch_rows)
    summary = build_summary(
        report=latitude_report,
        prepared=prepared,
        latitude_column=latitude_column,
        longitude_column=longitude_column,
        estimate_rows=estimate_rows,
        branch_rows=reviewed_branch_rows,
        outlier_rows=outlier_rows,
        warnings=warnings,
    )
    return PhylogeographicCoordinateReport(
        tree_path=tree_path,
        table_path=table_path,
        taxon_column=prepared.taxon_column,
        latitude_column=latitude_column,
        longitude_column=longitude_column,
        model=model,
        alpha=alpha,
        summary=summary,
        estimate_rows=estimate_rows,
        branch_rows=reviewed_branch_rows,
        outlier_rows=outlier_rows,
        exclusion_rows=prepared.exclusion_rows,
        warnings=warnings,
    )
