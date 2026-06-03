from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

from bijux_phylogenetics.io.newick import dumps_newick

from ....artifacts import section
from .instability_sections import build_instability_sections
from .sensitivity_sections import build_sensitivity_sections
from .shared import truncate_dataclass_rows
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
    """Build reviewer-facing sections for the tree uncertainty report."""
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
        section(
            "methods-summary-text", methods_summary_path.read_text(encoding="utf-8")
        ),
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
        *build_sensitivity_sections(
            storage_risk=storage_risk,
            thinning_sensitivity=thinning_sensitivity,
            consensus_sensitivity=consensus_sensitivity,
            benchmark=benchmark,
            benchmark_tree_count=benchmark_tree_count,
            benchmark_taxon_count=benchmark_taxon_count,
            maturity=maturity,
            summary=summary,
            scaled_report_note=scaled_report_note,
            artifact_paths=artifact_paths,
            out_path=out_path,
            preview_limit=preview_limit,
            thinning_rows=thinning_rows,
            thinning_truncated=thinning_truncated,
            consensus_rows=consensus_rows,
            consensus_truncated=consensus_truncated,
            benchmark_rows=benchmark_rows,
            benchmark_truncated=benchmark_truncated,
        ),
    ]
    return sections, truncated_sections
