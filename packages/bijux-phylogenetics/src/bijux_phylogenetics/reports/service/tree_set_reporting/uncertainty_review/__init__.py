from __future__ import annotations

from pathlib import Path
from time import perf_counter
import tracemalloc

from bijux_phylogenetics.trees import (
    TreeSetProcessingSummary,
    assess_tree_set_maturity,
    assess_tree_set_storage_risk,
    assess_tree_set_thinning_sensitivity,
    benchmark_tree_set_uncertainty,
    build_tree_set_budget_report,
    build_tree_set_workflow_budget,
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

from ...models import TreeUncertaintyReportBuildResult
from .artifact_outputs import (
    finalize_tree_uncertainty_outputs,
    write_tree_uncertainty_artifacts,
)
from .machine_manifest import build_tree_uncertainty_manifest
from .sections import build_tree_uncertainty_sections


def render_tree_uncertainty_report(
    *,
    tree_set_path: Path,
    out_path: Path,
    max_tree_count: int | None = None,
    max_report_table_rows: int | None = None,
    memory_warning_threshold_bytes: int | None = None,
) -> TreeUncertaintyReportBuildResult:
    """Build a deterministic HTML report for consensus and uncertainty across a tree set."""
    budget = build_tree_set_workflow_budget(
        max_tree_count=max_tree_count,
        max_report_table_rows=max_report_table_rows,
        memory_warning_threshold_bytes=memory_warning_threshold_bytes,
    )
    started = perf_counter()
    started_tracing = tracemalloc.is_tracing()
    if not started_tracing:
        tracemalloc.start()
    try:
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
        _current, peak = tracemalloc.get_traced_memory()
        processing = TreeSetProcessingSummary(
            runtime_seconds=round(perf_counter() - started, 6),
            peak_memory_bytes=peak,
            skipped_malformed_tree_count=summary.processing.skipped_malformed_tree_count,
        )
        artifact_paths = write_tree_uncertainty_artifacts(
            artifact_root=artifact_root,
            summary=summary,
            methods_summary_result=methods_summary_result,
            consensus_tree=consensus_tree,
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
            multimodality=multimodality,
            storage_risk=storage_risk,
            maturity=maturity,
            scaled_report_note=scaled_report_note,
        )
        sections, truncated_sections = build_tree_uncertainty_sections(
            budget=budget,
            summary=summary,
            consensus_tree=consensus_tree,
            consensus=consensus,
            clade_frequencies=clade_frequencies,
            diversity=diversity,
            clusters=clusters,
            unstable_taxa=unstable_taxa,
            unstable_clades=unstable_clades,
            clade_conflicts=clade_conflicts,
            conclusion_summary=conclusion_summary,
            storage_risk=storage_risk,
            thinning_sensitivity=thinning_sensitivity,
            consensus_sensitivity=consensus_sensitivity,
            benchmark=benchmark,
            benchmark_tree_count=benchmark_tree_count,
            benchmark_taxon_count=benchmark_taxon_count,
            multimodality=multimodality,
            maturity=maturity,
            scaled_report_note=scaled_report_note,
            limitations=limitations,
            methods_summary_path=artifact_paths["methods_summary"],
            artifact_paths=artifact_paths,
            out_path=out_path,
        )
        core_sections = sections[:11]
        supplemental_sections = sections[11:]
        budget_report = build_tree_set_budget_report(
            budget=budget,
            peak_memory_bytes=processing.peak_memory_bytes,
            truncated_section_names=truncated_sections,
        )
        machine_manifest, artifact_links, artifact_manifest_path = (
            build_tree_uncertainty_manifest(
                title=title,
                tree_set_path=tree_set_path,
                out_path=out_path,
                artifact_root=artifact_root,
                summary=summary,
                processing=processing,
                budget_report=budget_report,
                scaled_report_mode=scaled_report_mode,
                methods_summary_result=methods_summary_result,
                limitations=limitations,
                artifact_paths=artifact_paths,
                core_sections=core_sections,
                supplemental_sections=supplemental_sections,
            )
        )
        summary_metrics = [
            ("tree count", summary.tree_count),
            ("rooted topologies", summary.rooted_topology_count),
            ("report mode", "scaled-summary" if scaled_report_mode else "full-review"),
            ("runtime seconds", processing.runtime_seconds),
            ("peak memory bytes", processing.peak_memory_bytes),
            ("linked artifacts", len(artifact_paths)),
        ]
        machine_manifest, html_size_bytes, linked_artifact_bytes = (
            finalize_tree_uncertainty_outputs(
                title=title,
                out_path=out_path,
                sections=sections,
                summary_metrics=summary_metrics,
                artifact_links=artifact_links,
                artifact_manifest_path=artifact_manifest_path,
                artifact_paths=artifact_paths,
                machine_manifest=machine_manifest,
            )
        )
        total_output_bytes = int(machine_manifest["total_output_bytes"])
        return TreeUncertaintyReportBuildResult(
            output_path=out_path,
            artifact_root=artifact_root,
            artifact_manifest_path=artifact_manifest_path,
            report_kind="tree-uncertainty",
            title=title,
            source_path=tree_set_path,
            tree_count=summary.tree_count,
            rooted_topology_count=summary.rooted_topology_count,
            processing=processing,
            budget_report=budget_report,
            methods_summary_path=artifact_paths["methods_summary"],
            methods_summary_warning_count=methods_summary_result.warning_count,
            limitations=limitations,
            linked_artifact_count=len(artifact_paths) + 1,
            html_size_bytes=html_size_bytes,
            linked_artifact_bytes=linked_artifact_bytes,
            total_output_bytes=total_output_bytes,
            machine_manifest=machine_manifest,
        )
    finally:
        if not started_tracing:
            tracemalloc.stop()
