from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

from bijux_phylogenetics.io.newick import dumps_newick

from ...artifacts import preview_report_rows, section
from .instability_sections import build_instability_sections
from .shared import artifact_link, preview_payload, truncate_dataclass_rows
from .topology_sections import build_topology_sections


def build_tree_uncertainty_sections(
    *,
    budget,
    summary,
    consensus_tree,
    consensus,
    clade_frequencies,
    diversity,
    clusters,
    unstable_taxa,
    unstable_clades,
    clade_conflicts,
    conclusion_summary,
    storage_risk,
    thinning_sensitivity,
    consensus_sensitivity,
    benchmark,
    benchmark_tree_count,
    benchmark_taxon_count,
    multimodality,
    maturity,
    scaled_report_note: dict[str, object],
    limitations: list[str],
    methods_summary_path: Path,
    artifact_paths: dict[str, Path],
    out_path: Path,
) -> tuple[list[tuple[str, object]], list[str]]:
    truncated_sections: list[str] = []
    preview_limit = 5
    clade_frequency_rows, clade_frequency_truncated = truncate_dataclass_rows(
        rows=clade_frequencies.clade_frequencies,
        limit=budget.max_report_table_rows,
        section_name="clade-frequencies",
        truncated_sections=truncated_sections,
    )
    rf_rows, rf_truncated = truncate_dataclass_rows(
        rows=diversity.rf_distribution,
        limit=budget.max_report_table_rows,
        section_name="rf-distance-distribution",
        truncated_sections=truncated_sections,
    )
    cluster_rows, cluster_truncated = truncate_dataclass_rows(
        rows=clusters.clusters,
        limit=budget.max_report_table_rows,
        section_name="topology-clusters",
        truncated_sections=truncated_sections,
    )
    unstable_taxa_rows, unstable_taxa_truncated = truncate_dataclass_rows(
        rows=unstable_taxa.taxa,
        limit=budget.max_report_table_rows,
        section_name="unstable-taxa",
        truncated_sections=truncated_sections,
    )
    unstable_clade_rows, unstable_clade_truncated = truncate_dataclass_rows(
        rows=unstable_clades.clades,
        limit=budget.max_report_table_rows,
        section_name="unstable-clades",
        truncated_sections=truncated_sections,
    )
    conflict_rows, conflict_truncated = truncate_dataclass_rows(
        rows=[] if clade_conflicts is None else clade_conflicts.conflicts,
        limit=budget.max_report_table_rows,
        section_name="clade-credibility-conflicts",
        truncated_sections=truncated_sections,
    )
    robust_rows, robust_truncated = truncate_dataclass_rows(
        rows=[] if conclusion_summary is None else conclusion_summary.robust_clades,
        limit=budget.max_report_table_rows,
        section_name="uncertainty-aware-conclusions.robust",
        truncated_sections=truncated_sections,
    )
    uncertain_rows, uncertain_truncated = truncate_dataclass_rows(
        rows=[] if conclusion_summary is None else conclusion_summary.uncertain_clades,
        limit=budget.max_report_table_rows,
        section_name="uncertainty-aware-conclusions.uncertain",
        truncated_sections=truncated_sections,
    )
    conflicting_rows, conflicting_truncated = truncate_dataclass_rows(
        rows=[]
        if conclusion_summary is None
        else conclusion_summary.conflicting_clades,
        limit=budget.max_report_table_rows,
        section_name="uncertainty-aware-conclusions.conflicting",
        truncated_sections=truncated_sections,
    )
    thinning_rows, thinning_truncated = truncate_dataclass_rows(
        rows=[] if thinning_sensitivity is None else thinning_sensitivity.rows,
        limit=budget.max_report_table_rows,
        section_name="thinning-sensitivity",
        truncated_sections=truncated_sections,
    )
    consensus_rows, consensus_truncated = truncate_dataclass_rows(
        rows=[] if consensus_sensitivity is None else consensus_sensitivity.rows,
        limit=budget.max_report_table_rows,
        section_name="consensus-threshold-sensitivity",
        truncated_sections=truncated_sections,
    )
    benchmark_rows, benchmark_truncated = truncate_dataclass_rows(
        rows=[] if benchmark is None else benchmark.rows,
        limit=budget.max_report_table_rows,
        section_name="tree-set-benchmark",
        truncated_sections=truncated_sections,
    )
    sections = [
        section("methods-summary-text", methods_summary_path.read_text(encoding="utf-8")),
        section("limitations", limitations),
        section("tree-set-summary", asdict(summary)),
        section(
            "consensus-tree",
            {"newick": dumps_newick(consensus_tree), "report": asdict(consensus)},
        ),
        *build_topology_sections(
            clade_frequencies=clade_frequencies,
            diversity=diversity,
            clusters=clusters,
            multimodality=multimodality,
            scaled_report_note=scaled_report_note,
            artifact_paths=artifact_paths,
            out_path=out_path,
            preview_limit=preview_limit,
            clade_frequency_rows=clade_frequency_rows,
            clade_frequency_truncated=clade_frequency_truncated,
            rf_rows=rf_rows,
            rf_truncated=rf_truncated,
            cluster_rows=cluster_rows,
            cluster_truncated=cluster_truncated,
        ),
        *build_instability_sections(
            unstable_taxa=unstable_taxa,
            unstable_clades=unstable_clades,
            clade_conflicts=clade_conflicts,
            conclusion_summary=conclusion_summary,
            scaled_report_note=scaled_report_note,
            artifact_paths=artifact_paths,
            out_path=out_path,
            preview_limit=preview_limit,
            unstable_taxa_rows=unstable_taxa_rows,
            unstable_taxa_truncated=unstable_taxa_truncated,
            unstable_clade_rows=unstable_clade_rows,
            unstable_clade_truncated=unstable_clade_truncated,
            conflict_rows=conflict_rows,
            conflict_truncated=conflict_truncated,
            robust_rows=robust_rows,
            robust_truncated=robust_truncated,
            uncertain_rows=uncertain_rows,
            uncertain_truncated=uncertain_truncated,
            conflicting_rows=conflicting_rows,
            conflicting_truncated=conflicting_truncated,
        ),
        section(
            "storage-risk",
            {
                **asdict(storage_risk),
                "artifact_path": (
                    artifact_paths["storage_risk"].relative_to(out_path.parent).as_posix()
                ),
            },
        ),
        section(
            "thinning-sensitivity",
            {
                **(
                    {
                        "path": str(thinning_sensitivity.path),
                        "original_tree_count": thinning_sensitivity.original_tree_count,
                        "original_rooted_topology_count": thinning_sensitivity.original_rooted_topology_count,
                        "original_dominant_topology_frequency": thinning_sensitivity.original_dominant_topology_frequency,
                        "warning_count": len(thinning_sensitivity.warnings),
                        "warnings": thinning_sensitivity.warnings,
                        "row_count": len(thinning_sensitivity.rows),
                        "truncated_row_count": thinning_truncated,
                        "preview_row_count": min(len(thinning_rows), preview_limit),
                        "preview_rows": preview_report_rows(
                            thinning_rows, limit=preview_limit
                        ),
                    }
                    if thinning_sensitivity is not None
                    else scaled_report_note
                ),
                "artifact_path": (
                    artifact_paths["thinning_sensitivity"]
                    .relative_to(out_path.parent)
                    .as_posix()
                ),
            },
        ),
        section(
            "consensus-threshold-sensitivity",
            {
                **(
                    {
                        "path": str(consensus_sensitivity.path),
                        "tree_count": consensus_sensitivity.tree_count,
                        "warning_count": len(consensus_sensitivity.warnings),
                        "warnings": consensus_sensitivity.warnings,
                        "row_count": len(consensus_sensitivity.rows),
                        "truncated_row_count": consensus_truncated,
                        "preview_row_count": min(len(consensus_rows), preview_limit),
                        "preview_rows": preview_report_rows(
                            consensus_rows, limit=preview_limit
                        ),
                    }
                    if consensus_sensitivity is not None
                    else scaled_report_note
                ),
                "artifact_path": (
                    artifact_paths["consensus_threshold_sensitivity"]
                    .relative_to(out_path.parent)
                    .as_posix()
                ),
            },
        ),
        section(
            "tree-set-benchmark",
            {
                **(
                    {
                        "tree_counts": benchmark.tree_counts,
                        "taxon_counts": benchmark.taxon_counts,
                        "sampled_tree_count": benchmark_tree_count,
                        "sampled_taxon_count": benchmark_taxon_count,
                        "benchmark_capped": benchmark_tree_count != summary.tree_count
                        or benchmark_taxon_count != max(len(summary.shared_taxa), 2),
                        "row_count": len(benchmark.rows),
                        "truncated_row_count": benchmark_truncated,
                        "preview_row_count": min(len(benchmark_rows), preview_limit),
                        "preview_rows": preview_report_rows(
                            benchmark_rows, limit=preview_limit
                        ),
                    }
                    if benchmark is not None
                    else scaled_report_note
                ),
                "artifact_path": (
                    artifact_paths["tree_set_benchmark"]
                    .relative_to(out_path.parent)
                    .as_posix()
                ),
            },
        ),
        section(
            "maturity-gate",
            (
                {
                    **asdict(maturity),
                    "artifact_path": (
                        artifact_paths["maturity_gate"]
                        .relative_to(out_path.parent)
                        .as_posix()
                    ),
                }
                if maturity is not None
                else {
                    **scaled_report_note,
                    "artifact_path": (
                        artifact_paths["maturity_gate"]
                        .relative_to(out_path.parent)
                        .as_posix()
                    ),
                }
            ),
        ),
    ]
    return sections, truncated_sections
