from __future__ import annotations

from pathlib import Path
from time import perf_counter
import tracemalloc

from bijux_phylogenetics.trees import (
    TreeSetProcessingSummary,
    build_tree_set_budget_report,
    build_tree_set_workflow_budget,
)

from ...models import TreeUncertaintyReportBuildResult
from .analysis import build_tree_uncertainty_review_analysis
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
        analysis = build_tree_uncertainty_review_analysis(
            tree_set_path=tree_set_path,
            out_path=out_path,
            budget=budget,
        )
        _current, peak = tracemalloc.get_traced_memory()
        processing = TreeSetProcessingSummary(
            runtime_seconds=round(perf_counter() - started, 6),
            peak_memory_bytes=peak,
            skipped_malformed_tree_count=analysis.summary.processing.skipped_malformed_tree_count,
        )
        artifact_paths = write_tree_uncertainty_artifacts(
            artifact_root=analysis.artifact_root,
            summary=analysis.summary,
            methods_summary_result=analysis.methods_summary_result,
            consensus_tree=analysis.consensus_tree,
            clade_frequencies=analysis.clade_frequencies,
            diversity=analysis.diversity,
            clusters=analysis.clusters,
            unstable_taxa=analysis.unstable_taxa,
            unstable_clades=analysis.unstable_clades,
            clade_conflicts=analysis.clade_conflicts,
            conclusion_summary=analysis.conclusion_summary,
            thinning_sensitivity=analysis.thinning_sensitivity,
            consensus_sensitivity=analysis.consensus_sensitivity,
            benchmark=analysis.benchmark,
            multimodality=analysis.multimodality,
            storage_risk=analysis.storage_risk,
            maturity=analysis.maturity,
            scaled_report_note=analysis.scaled_report_note,
        )
        sections, truncated_sections = build_tree_uncertainty_sections(
            budget=budget,
            summary=analysis.summary,
            consensus_tree=analysis.consensus_tree,
            consensus=analysis.consensus,
            clade_frequencies=analysis.clade_frequencies,
            diversity=analysis.diversity,
            clusters=analysis.clusters,
            unstable_taxa=analysis.unstable_taxa,
            unstable_clades=analysis.unstable_clades,
            clade_conflicts=analysis.clade_conflicts,
            conclusion_summary=analysis.conclusion_summary,
            storage_risk=analysis.storage_risk,
            thinning_sensitivity=analysis.thinning_sensitivity,
            consensus_sensitivity=analysis.consensus_sensitivity,
            benchmark=analysis.benchmark,
            benchmark_tree_count=analysis.benchmark_tree_count,
            benchmark_taxon_count=analysis.benchmark_taxon_count,
            multimodality=analysis.multimodality,
            maturity=analysis.maturity,
            scaled_report_note=analysis.scaled_report_note,
            limitations=analysis.limitations,
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
                title=analysis.title,
                tree_set_path=tree_set_path,
                out_path=out_path,
                artifact_root=analysis.artifact_root,
                summary=analysis.summary,
                processing=processing,
                budget_report=budget_report,
                scaled_report_mode=analysis.scaled_report_mode,
                methods_summary_result=analysis.methods_summary_result,
                limitations=analysis.limitations,
                artifact_paths=artifact_paths,
                core_sections=core_sections,
                supplemental_sections=supplemental_sections,
            )
        )
        summary_metrics = [
            ("tree count", analysis.summary.tree_count),
            ("rooted topologies", analysis.summary.rooted_topology_count),
            (
                "report mode",
                "scaled-summary" if analysis.scaled_report_mode else "full-review",
            ),
            ("runtime seconds", processing.runtime_seconds),
            ("peak memory bytes", processing.peak_memory_bytes),
            ("linked artifacts", len(artifact_paths)),
        ]
        machine_manifest, html_size_bytes, linked_artifact_bytes = (
            finalize_tree_uncertainty_outputs(
                title=analysis.title,
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
            artifact_root=analysis.artifact_root,
            artifact_manifest_path=artifact_manifest_path,
            report_kind="tree-uncertainty",
            title=analysis.title,
            source_path=tree_set_path,
            tree_count=analysis.summary.tree_count,
            rooted_topology_count=analysis.summary.rooted_topology_count,
            processing=processing,
            budget_report=budget_report,
            methods_summary_path=artifact_paths["methods_summary"],
            methods_summary_warning_count=analysis.methods_summary_result.warning_count,
            limitations=analysis.limitations,
            linked_artifact_count=len(artifact_paths) + 1,
            html_size_bytes=html_size_bytes,
            linked_artifact_bytes=linked_artifact_bytes,
            total_output_bytes=total_output_bytes,
            machine_manifest=machine_manifest,
        )
    finally:
        if not started_tracing:
            tracemalloc.stop()
