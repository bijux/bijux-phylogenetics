from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

from bijux_phylogenetics.render.html import write_html_report
from bijux_phylogenetics.trees import (
    cluster_trees_by_topology,
    compare_posterior_topological_diversity,
    compare_posterior_tree_sets,
    detect_posterior_topology_multimodality,
    detect_unstable_clades,
    detect_unstable_taxa,
    load_tree_set,
    summarize_clade_credibility_conflicts,
    summarize_uncertainty_aware_conclusions,
)

from ..artifacts import section
from ..models import TreeSetComparisonReportBuildResult


def render_tree_set_comparison_report(
    *,
    left_tree_set_path: Path,
    right_tree_set_path: Path,
    out_path: Path,
) -> TreeSetComparisonReportBuildResult:
    """Render an HTML comparison report for two tree sets."""
    comparison = compare_posterior_tree_sets(left_tree_set_path, right_tree_set_path)
    left_summary = load_tree_set(left_tree_set_path)
    right_summary = load_tree_set(right_tree_set_path)
    left_clusters = cluster_trees_by_topology(left_tree_set_path)
    right_clusters = cluster_trees_by_topology(right_tree_set_path)
    diversity = compare_posterior_topological_diversity(
        left_tree_set_path, right_tree_set_path
    )
    left_multimodality = detect_posterior_topology_multimodality(left_tree_set_path)
    right_multimodality = detect_posterior_topology_multimodality(right_tree_set_path)
    left_unstable_taxa = detect_unstable_taxa(left_tree_set_path)
    right_unstable_taxa = detect_unstable_taxa(right_tree_set_path)
    left_unstable_clades = detect_unstable_clades(left_tree_set_path)
    right_unstable_clades = detect_unstable_clades(right_tree_set_path)
    left_conflicts = summarize_clade_credibility_conflicts(left_tree_set_path)
    right_conflicts = summarize_clade_credibility_conflicts(right_tree_set_path)
    left_conclusions = summarize_uncertainty_aware_conclusions(left_tree_set_path)
    right_conclusions = summarize_uncertainty_aware_conclusions(right_tree_set_path)
    sections = [
        section("tree-set-comparison", asdict(comparison)),
        section("topological-diversity-comparison", asdict(diversity)),
        section("left-tree-set-summary", asdict(left_summary)),
        section("right-tree-set-summary", asdict(right_summary)),
        section("left-topology-clusters", asdict(left_clusters)),
        section("right-topology-clusters", asdict(right_clusters)),
        section("left-topology-multimodality", asdict(left_multimodality)),
        section("right-topology-multimodality", asdict(right_multimodality)),
        section("left-unstable-taxa", asdict(left_unstable_taxa)),
        section("right-unstable-taxa", asdict(right_unstable_taxa)),
        section("left-unstable-clades", asdict(left_unstable_clades)),
        section("right-unstable-clades", asdict(right_unstable_clades)),
        section("left-clade-credibility-conflicts", asdict(left_conflicts)),
        section("right-clade-credibility-conflicts", asdict(right_conflicts)),
        section("left-uncertainty-aware-conclusions", asdict(left_conclusions)),
        section("right-uncertainty-aware-conclusions", asdict(right_conclusions)),
    ]
    title = "Bijux Tree-Set Comparison Report"
    machine_manifest = {
        "report_kind": "tree-set-comparison",
        "title": title,
        "left_path": str(left_tree_set_path),
        "right_path": str(right_tree_set_path),
        "shared_rooted_topology_count": comparison.shared_rooted_topology_count,
        "sections": [name for name, _ in sections],
    }
    write_html_report(
        title=title,
        sections=sections,
        out_path=out_path,
        embedded_json=machine_manifest,
    )
    return TreeSetComparisonReportBuildResult(
        output_path=out_path,
        report_kind="tree-set-comparison",
        title=title,
        left_path=left_tree_set_path,
        right_path=right_tree_set_path,
        shared_rooted_topology_count=comparison.shared_rooted_topology_count,
        machine_manifest=machine_manifest,
    )
