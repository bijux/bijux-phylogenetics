from __future__ import annotations

import csv
from pathlib import Path

from bijux_phylogenetics.io.newick import loads_newick, write_newick

from ..tree_sets.clade_support import (
    _build_clade_frequency_report,
    write_clade_frequency_table,
)
from ..tree_sets.inventory import _analyze_tree_set
from ..tree_sets.quartet_concordance import (
    _build_gene_tree_quartet_concordance_report,
)
from .instability import (
    summarize_clade_credibility_conflicts,
    write_clade_credibility_conflict_table,
)
from .models import (
    GeneTreeConflictArtifactReport,
    GeneTreeConflictQuartetSummary,
    GeneTreeConflictReferenceTree,
    GeneTreeConflictSummaryReport,
)
from .rogue_taxa import detect_rogue_taxa, write_rogue_taxon_table
from .topology_diversity import _build_topology_cluster_report


def summarize_gene_tree_conflicts(
    path: Path,
    *,
    credibility_threshold: float = 0.5,
    rogue_consensus_threshold: float = 0.5,
) -> GeneTreeConflictSummaryReport:
    """Bundle clade, quartet, rogue-taxon, and conflict evidence for one gene-tree set."""
    analysis = _analyze_tree_set(path)
    topology_clusters = _build_topology_cluster_report(analysis)
    reference_cluster = topology_clusters.clusters[0]
    reference_tree = analysis.rooted_representatives[
        reference_cluster.rooted_topology_id
    ][2]
    clade_frequencies = _build_clade_frequency_report(analysis)
    quartet_report = _build_gene_tree_quartet_concordance_report(
        species_tree_path=path,
        species_tree=reference_tree,
        analysis=analysis,
    )
    rogue_taxa = detect_rogue_taxa(
        path,
        consensus_threshold=rogue_consensus_threshold,
    )
    clade_conflicts = summarize_clade_credibility_conflicts(
        path,
        credibility_threshold=credibility_threshold,
    )
    return GeneTreeConflictSummaryReport(
        path=path,
        tree_count=len(analysis.trees),
        processing=analysis.processing,
        shared_taxa=analysis.shared_taxa,
        credibility_threshold=credibility_threshold,
        rogue_consensus_threshold=rogue_consensus_threshold,
        reference_tree=GeneTreeConflictReferenceTree(
            selection_method="dominant-rooted-topology-representative",
            rooted_topology_id=reference_cluster.rooted_topology_id,
            frequency=reference_cluster.frequency,
            newick=reference_cluster.representative_newick,
        ),
        clade_frequencies=clade_frequencies,
        quartet_concordance=GeneTreeConflictQuartetSummary(
            branch_count=quartet_report.branch_count,
            total_quartet_count=quartet_report.total_quartet_count,
            concordant_quartet_count=quartet_report.concordant_quartet_count,
            discordant_first_quartet_count=(
                quartet_report.discordant_first_quartet_count
            ),
            discordant_second_quartet_count=(
                quartet_report.discordant_second_quartet_count
            ),
            uninformative_quartet_count=quartet_report.uninformative_quartet_count,
            informative_quartet_count=quartet_report.informative_quartet_count,
            rows=quartet_report.rows,
        ),
        rogue_taxa=rogue_taxa,
        clade_conflicts=clade_conflicts,
    )


def write_gene_tree_conflict_summary_table(
    path: Path,
    report: GeneTreeConflictSummaryReport,
) -> Path:
    """Write one summary row for a gene-tree conflict review bundle."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "tree_count",
                "runtime_seconds",
                "peak_memory_bytes",
                "skipped_malformed_tree_count",
                "shared_taxon_count",
                "reference_rooted_topology_id",
                "reference_tree_frequency",
                "clade_count",
                "quartet_branch_count",
                "total_quartet_count",
                "conflict_count",
                "conflicting_clade_count",
                "rogue_taxon_count",
                "top_ranked_rogue_taxon",
                "credibility_threshold",
                "rogue_consensus_threshold",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        conflicting_clades = {
            row.left_clade for row in report.clade_conflicts.conflicts
        } | {row.right_clade for row in report.clade_conflicts.conflicts}
        writer.writerow(
            {
                "tree_count": report.tree_count,
                "runtime_seconds": format(report.processing.runtime_seconds, ".15g"),
                "peak_memory_bytes": report.processing.peak_memory_bytes,
                "skipped_malformed_tree_count": (
                    report.processing.skipped_malformed_tree_count
                ),
                "shared_taxon_count": len(report.shared_taxa),
                "reference_rooted_topology_id": (
                    report.reference_tree.rooted_topology_id
                ),
                "reference_tree_frequency": format(
                    report.reference_tree.frequency,
                    ".15g",
                ),
                "clade_count": len(report.clade_frequencies.clade_frequencies),
                "quartet_branch_count": report.quartet_concordance.branch_count,
                "total_quartet_count": report.quartet_concordance.total_quartet_count,
                "conflict_count": report.clade_conflicts.conflict_count,
                "conflicting_clade_count": len(conflicting_clades),
                "rogue_taxon_count": len(report.rogue_taxa.rows),
                "top_ranked_rogue_taxon": report.rogue_taxa.rows[0].taxon,
                "credibility_threshold": format(
                    report.credibility_threshold,
                    ".15g",
                ),
                "rogue_consensus_threshold": format(
                    report.rogue_consensus_threshold,
                    ".15g",
                ),
            }
        )
    return path


def write_gene_tree_conflict_quartet_table(
    path: Path,
    report: GeneTreeConflictSummaryReport,
) -> Path:
    """Write branch-level quartet concordance for the selected reference topology."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "branch_id",
                "left_taxa",
                "right_taxa",
                "quartet_count_per_tree",
                "concordant_quartet_count",
                "discordant_first_quartet_count",
                "discordant_second_quartet_count",
                "uninformative_quartet_count",
                "informative_quartet_count",
                "concordance_factor",
                "concordant_frequency",
                "discordant_first_frequency",
                "discordant_second_frequency",
                "uninformative_frequency",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.quartet_concordance.rows:
            writer.writerow(
                {
                    "branch_id": row.branch_id,
                    "left_taxa": "|".join(row.left_taxa),
                    "right_taxa": "|".join(row.right_taxa),
                    "quartet_count_per_tree": row.quartet_count_per_tree,
                    "concordant_quartet_count": row.concordant_quartet_count,
                    "discordant_first_quartet_count": (
                        row.discordant_first_quartet_count
                    ),
                    "discordant_second_quartet_count": (
                        row.discordant_second_quartet_count
                    ),
                    "uninformative_quartet_count": row.uninformative_quartet_count,
                    "informative_quartet_count": row.informative_quartet_count,
                    "concordance_factor": (
                        ""
                        if row.concordance_factor is None
                        else format(row.concordance_factor, ".15g")
                    ),
                    "concordant_frequency": format(
                        row.concordant_frequency,
                        ".15g",
                    ),
                    "discordant_first_frequency": format(
                        row.discordant_first_frequency,
                        ".15g",
                    ),
                    "discordant_second_frequency": format(
                        row.discordant_second_frequency,
                        ".15g",
                    ),
                    "uninformative_frequency": format(
                        row.uninformative_frequency,
                        ".15g",
                    ),
                }
            )
    return path


def write_gene_tree_conflict_artifacts(
    tree_set_path: Path,
    *,
    out_dir: Path,
    prefix: str = "gene-tree-conflicts",
    credibility_threshold: float = 0.5,
    rogue_consensus_threshold: float = 0.5,
) -> GeneTreeConflictArtifactReport:
    """Write a governed artifact bundle for one gene-tree conflict review."""
    report = summarize_gene_tree_conflicts(
        tree_set_path,
        credibility_threshold=credibility_threshold,
        rogue_consensus_threshold=rogue_consensus_threshold,
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    base_path = out_dir / prefix
    output_paths = {
        "summary_table": write_gene_tree_conflict_summary_table(
            base_path.with_suffix(".summary.tsv"),
            report,
        ),
        "reference_tree": write_newick(
            base_path.with_suffix(".reference-tree.nwk"),
            loads_newick(report.reference_tree.newick),
        ),
        "clade_frequencies": write_clade_frequency_table(
            base_path.with_suffix(".clade-frequencies.tsv"),
            report.clade_frequencies,
        ),
        "quartet_concordance": write_gene_tree_conflict_quartet_table(
            base_path.with_suffix(".quartet-concordance.tsv"),
            report,
        ),
        "rogue_taxa": write_rogue_taxon_table(
            base_path.with_suffix(".rogue-taxa.tsv"),
            report.rogue_taxa,
        ),
        "clade_conflicts": write_clade_credibility_conflict_table(
            base_path.with_suffix(".clade-conflicts.tsv"),
            report.clade_conflicts,
        ),
    }
    return GeneTreeConflictArtifactReport(
        input_path=tree_set_path,
        out_dir=out_dir,
        prefix=prefix,
        summary_report=report,
        output_paths=output_paths,
    )
