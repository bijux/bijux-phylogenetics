from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ..tree_sets import (
    CladeFrequencyReport,
    ConsensusTreeReport,
    TreeSetReport,
    compute_clade_frequency_table,
    compute_consensus_tree,
    load_tree_set,
)
from .instability import (
    CladeCredibilityConflictReport,
    UncertaintyAwareConclusionSummaryReport,
    UnstableTaxaReport,
    detect_unstable_taxa,
    summarize_clade_credibility_conflicts,
    summarize_uncertainty_aware_conclusions,
)
from .models import (
    PosteriorTopologyMultimodalityReport,
    TreeTopologyClusterReport,
)
from .topology_diversity import (
    cluster_trees_by_topology,
    detect_posterior_topology_multimodality,
)


@dataclass(slots=True)
class TreeSetUncertaintyMethodReport:
    tree_set_path: Path
    summary: TreeSetReport
    consensus: ConsensusTreeReport
    clade_frequencies: CladeFrequencyReport
    topology_clusters: TreeTopologyClusterReport
    unstable_taxa: UnstableTaxaReport
    multimodality: PosteriorTopologyMultimodalityReport
    clade_conflicts: CladeCredibilityConflictReport
    conclusion_summary: UncertaintyAwareConclusionSummaryReport


@dataclass(slots=True)
class TreeSetUncertaintyMethodsSummaryTextResult:
    output_path: Path
    title: str
    warning_count: int
    warnings: tuple[str, ...]
    topology_cluster_count: int
    unstable_taxon_count: int
    text: str
    report: TreeSetUncertaintyMethodReport


def _bullet_list(values: list[str]) -> str:
    if not values:
        return "none"
    return ", ".join(f"`{value}`" for value in values)


def _deduplicate_text(values: list[str]) -> list[str]:
    seen: set[str] = set()
    ordered: list[str] = []
    for value in values:
        normalized = value.strip()
        if not normalized or normalized in seen:
            continue
        seen.add(normalized)
        ordered.append(normalized)
    return ordered


def build_tree_set_uncertainty_method_report(
    tree_set_path: Path,
) -> TreeSetUncertaintyMethodReport:
    """Build an integrated reviewer-facing evidence model for one tree-set uncertainty analysis."""
    summary = load_tree_set(tree_set_path)
    _consensus_tree, consensus = compute_consensus_tree(tree_set_path)
    clade_frequencies = compute_clade_frequency_table(tree_set_path)
    topology_clusters = cluster_trees_by_topology(tree_set_path)
    unstable_taxa = detect_unstable_taxa(tree_set_path)
    multimodality = detect_posterior_topology_multimodality(tree_set_path)
    clade_conflicts = summarize_clade_credibility_conflicts(tree_set_path)
    conclusion_summary = summarize_uncertainty_aware_conclusions(tree_set_path)
    return TreeSetUncertaintyMethodReport(
        tree_set_path=tree_set_path,
        summary=summary,
        consensus=consensus,
        clade_frequencies=clade_frequencies,
        topology_clusters=topology_clusters,
        unstable_taxa=unstable_taxa,
        multimodality=multimodality,
        clade_conflicts=clade_conflicts,
        conclusion_summary=conclusion_summary,
    )


def _tree_set_uncertainty_methods_summary_warnings(
    report: TreeSetUncertaintyMethodReport,
) -> list[str]:
    warnings = list(report.multimodality.warnings)
    if report.unstable_taxa.taxa:
        warnings.append(
            f"{len(report.unstable_taxa.taxa)} taxa change placement across the tree set"
        )
    if report.clade_conflicts.conflict_count:
        warnings.append(
            f"{report.clade_conflicts.conflict_count} high-credibility clade conflicts remain unresolved"
        )
    if report.multimodality.multimodal:
        warnings.append(
            f"tree-set topology is multimodal across {report.multimodality.mode_count} modes"
        )
    return _deduplicate_text(warnings)


def build_tree_set_uncertainty_methods_summary_text(
    report: TreeSetUncertaintyMethodReport,
) -> str:
    """Build reviewer-facing Markdown methods text for one tree-set uncertainty analysis."""
    warnings = _tree_set_uncertainty_methods_summary_warnings(report)
    summary = report.summary
    consensus = report.consensus
    multimodality = report.multimodality
    return (
        "# Tree-Set Uncertainty Methods Summary\n\n"
        f"This tree-set uncertainty analysis reviewed Newick tree-set input `{report.tree_set_path.name}` through one generic workflow that applies to posterior samples, bootstrap replicates, or other tree collections when the trees share interpretable taxa. The current input retained `{summary.tree_count}` trees with `{len(summary.shared_taxa)}` shared taxa and `{summary.rooted_topology_count}` rooted topology signatures.\n\n"
        "## Tree-Set Input And Consensus Construction\n\n"
        f"- tree-set format: `{summary.source_format}`\n"
        f"- retained tree count: `{summary.tree_count}`\n"
        f"- skipped malformed tree count: `{summary.processing.skipped_malformed_tree_count}`\n"
        f"- shared taxon count: `{len(summary.shared_taxa)}`\n"
        f"- taxa represented anywhere in the tree set: `{len(summary.taxa_union)}`\n"
        f"- consensus contract: `{consensus.consensus_method}`\n"
        f"- consensus threshold: `{consensus.consensus_threshold}`\n"
        f"- included consensus clade count: `{consensus.included_clade_count}`\n"
        f"- consensus tree newick: `{consensus.consensus_newick}`\n\n"
        "## Support And Topology Dispersion\n\n"
        f"- clade support rows retained: `{len(report.clade_frequencies.clade_frequencies)}`\n"
        f"- dominant rooted topology frequency: `{format(multimodality.dominant_mode_frequency, '.15g')}`\n"
        f"- rooted topology cluster count: `{len(report.topology_clusters.clusters)}`\n"
        f"- multimodal topology distribution: `{'yes' if multimodality.multimodal else 'no'}`\n"
        f"- topology mode count: `{multimodality.mode_count}`\n"
        f"- multimodality warnings: {_bullet_list(multimodality.warnings)}\n\n"
        "## Instability And Conflict Review\n\n"
        f"- unstable taxon count: `{len(report.unstable_taxa.taxa)}`\n"
        f"- high-credibility clade conflict count: `{report.clade_conflicts.conflict_count}`\n"
        f"- robust clade count: `{report.conclusion_summary.robust_clade_count}`\n"
        f"- uncertain clade count: `{report.conclusion_summary.uncertain_clade_count}`\n"
        f"- conflicting clade count: `{report.conclusion_summary.conflicting_clade_count}`\n\n"
        "## Assumptions And Interpretation Boundaries\n\n"
        "- this workflow treats the input as a generic tree set and does not infer posterior versus bootstrap provenance unless that provenance is supplied elsewhere in the study record\n"
        "- consensus and clade support are computed on the shared taxon set, so trees with non-identical taxon coverage are summarized through their overlap rather than through taxon-imputed expansions\n"
        "- instability scores are descriptive placement summaries and do not on their own justify biological process claims without reference to the underlying inference model\n"
        "- topology-cluster and conflict summaries show structural disagreement across trees but do not resolve which mode should be treated as biologically preferred\n\n"
        "## Reviewer Warnings\n\n"
        f"- combined warning count: `{len(warnings)}`\n"
        f"- warning details: {_bullet_list(warnings)}\n"
    )


def write_tree_set_uncertainty_methods_summary_text(
    path: Path,
    report: TreeSetUncertaintyMethodReport,
) -> TreeSetUncertaintyMethodsSummaryTextResult:
    """Write reviewer-facing Markdown methods text for one tree-set uncertainty analysis."""
    text = build_tree_set_uncertainty_methods_summary_text(report)
    warnings = _tree_set_uncertainty_methods_summary_warnings(report)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return TreeSetUncertaintyMethodsSummaryTextResult(
        output_path=path,
        title="Tree-Set Uncertainty Methods Summary",
        warning_count=len(warnings),
        warnings=tuple(warnings),
        topology_cluster_count=len(report.topology_clusters.clusters),
        unstable_taxon_count=len(report.unstable_taxa.taxa),
        text=text,
        report=report,
    )
