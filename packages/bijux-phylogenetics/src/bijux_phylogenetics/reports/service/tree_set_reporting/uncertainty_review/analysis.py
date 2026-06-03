from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from bijux_phylogenetics.trees import (
    assess_tree_set_maturity,
    assess_tree_set_storage_risk,
    assess_tree_set_thinning_sensitivity,
    benchmark_tree_set_uncertainty,
    cluster_trees_by_topology,
    compare_consensus_thresholds,
    compute_clade_frequency_table,
    compute_consensus_tree,
    detect_posterior_topology_multimodality,
    detect_unstable_clades,
    detect_unstable_taxa,
    enforce_tree_set_tree_budget,
    load_tree_set,
    summarize_clade_credibility_conflicts,
    summarize_posterior_topology_diversity,
    summarize_uncertainty_aware_conclusions,
    write_tree_set_uncertainty_methods_summary_text,
)
from bijux_phylogenetics.trees.uncertainty import (
    build_tree_set_uncertainty_method_report,
)


@dataclass(slots=True)
class TreeUncertaintyReviewAnalysis:
    title: str
    artifact_root: Path
    summary: Any
    scaled_report_mode: bool
    scaled_report_note: dict[str, object]
    methods_summary_result: Any
    consensus_tree: Any
    consensus: Any
    clade_frequencies: Any
    diversity: Any
    clusters: Any
    unstable_taxa: Any
    unstable_clades: Any
    clade_conflicts: Any | None
    conclusion_summary: Any | None
    thinning_sensitivity: Any | None
    consensus_sensitivity: Any | None
    benchmark: Any | None
    benchmark_tree_count: int | None
    benchmark_taxon_count: int | None
    multimodality: Any | None
    storage_risk: Any
    maturity: Any | None
    limitations: list[str]


def build_tree_uncertainty_review_analysis(
    *,
    tree_set_path: Path,
    out_path: Path,
    budget,
) -> TreeUncertaintyReviewAnalysis:
    """Resolve the analysis state that drives the tree uncertainty review."""
    summary = load_tree_set(tree_set_path)
    scaled_report_mode = summary.tree_count >= 1000
    scaled_report_note = {
        "status": "summary-only",
        "reason": (
            "supplemental sensitivity analyses were replaced with linked note artifacts "
            "because the report input exceeds the large tree-set scaling threshold"
        ),
        "tree_count": summary.tree_count,
    }
    enforce_tree_set_tree_budget(
        tree_count=summary.tree_count,
        budget=budget,
        workflow_name="tree uncertainty report",
        source_path=tree_set_path,
    )
    methods_report = build_tree_set_uncertainty_method_report(tree_set_path)
    consensus_tree, consensus = compute_consensus_tree(tree_set_path)
    clade_frequencies = compute_clade_frequency_table(tree_set_path)
    clusters = cluster_trees_by_topology(tree_set_path)
    diversity = summarize_posterior_topology_diversity(tree_set_path)
    unstable_taxa = detect_unstable_taxa(tree_set_path)
    unstable_clades = detect_unstable_clades(tree_set_path)
    storage_risk = assess_tree_set_storage_risk(tree_set_path)
    if scaled_report_mode:
        multimodality = None
        clade_conflicts = None
        conclusion_summary = None
        thinning_sensitivity = None
        consensus_sensitivity = None
        maturity = None
        benchmark = None
        benchmark_tree_count = None
        benchmark_taxon_count = None
    else:
        multimodality = detect_posterior_topology_multimodality(tree_set_path)
        clade_conflicts = summarize_clade_credibility_conflicts(tree_set_path)
        conclusion_summary = summarize_uncertainty_aware_conclusions(tree_set_path)
        thinning_sensitivity = assess_tree_set_thinning_sensitivity(tree_set_path)
        consensus_sensitivity = compare_consensus_thresholds(tree_set_path)
        maturity = assess_tree_set_maturity(tree_set_path)
        benchmark_tree_count = min(summary.tree_count, 128)
        benchmark_taxon_count = min(max(len(summary.shared_taxa), 2), 64)
        benchmark = benchmark_tree_set_uncertainty(
            tree_counts=[benchmark_tree_count],
            taxon_counts=[benchmark_taxon_count],
        )
    title = "Bijux Tree Uncertainty Report"
    artifact_root = out_path.parent / f"{out_path.stem}.artifacts"
    methods_summary_result = write_tree_set_uncertainty_methods_summary_text(
        artifact_root / "tree-set-uncertainty-methods-summary.md",
        methods_report,
    )
    limitations = sorted(
        {
            "consensus support and topology summaries describe the supplied tree set and should not be treated as direct proof of one true history",
            "alternative rooted modes, unstable taxa, and conflict-prone clades must remain part of interpretation instead of being collapsed into the consensus tree alone",
            *methods_summary_result.warnings,
            *([str(scaled_report_note["reason"])] if scaled_report_mode else []),
        }
    )
    return TreeUncertaintyReviewAnalysis(
        title=title,
        artifact_root=artifact_root,
        summary=summary,
        scaled_report_mode=scaled_report_mode,
        scaled_report_note=scaled_report_note,
        methods_summary_result=methods_summary_result,
        consensus_tree=consensus_tree,
        consensus=consensus,
        clade_frequencies=clade_frequencies,
        diversity=diversity,
        clusters=clusters,
        unstable_taxa=unstable_taxa,
        unstable_clades=unstable_clades,
        clade_conflicts=clade_conflicts,
        conclusion_summary=conclusion_summary,
        thinning_sensitivity=thinning_sensitivity,
        consensus_sensitivity=consensus_sensitivity,
        benchmark=benchmark,
        benchmark_tree_count=benchmark_tree_count,
        benchmark_taxon_count=benchmark_taxon_count,
        multimodality=multimodality,
        storage_risk=storage_risk,
        maturity=maturity,
        limitations=limitations,
    )
