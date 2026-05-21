from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
from time import perf_counter
import tracemalloc

from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.render.html import write_html_report
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
    write_clade_frequency_table,
    write_consensus_tree,
    write_topology_cluster_table,
    write_tree_distance_distribution_table,
    write_tree_set_uncertainty_methods_summary_text,
    write_unstable_clade_table,
)
from bijux_phylogenetics.trees.uncertainty import (
    build_tree_set_uncertainty_method_report,
)

from ..artifacts import (
    preview_report_rows,
    section,
    truncate_report_rows,
    write_json_artifact,
    write_tabular_artifact,
)
from ..ledger import sha256
from ..models import TreeUncertaintyReportBuildResult
from .comparison_report import render_tree_set_comparison_report


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
        truncated_sections: list[str] = []
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
        preview_limit = 5
        clade_frequency_rows, clade_frequency_truncated = truncate_report_rows(
            [asdict(row) for row in clade_frequencies.clade_frequencies],
            limit=budget.max_report_table_rows,
            section_name="clade-frequencies",
            truncated_sections=truncated_sections,
        )
        rf_rows, rf_truncated = truncate_report_rows(
            [asdict(row) for row in diversity.rf_distribution],
            limit=budget.max_report_table_rows,
            section_name="rf-distance-distribution",
            truncated_sections=truncated_sections,
        )
        cluster_rows, cluster_truncated = truncate_report_rows(
            [asdict(row) for row in clusters.clusters],
            limit=budget.max_report_table_rows,
            section_name="topology-clusters",
            truncated_sections=truncated_sections,
        )
        unstable_taxa_rows, unstable_taxa_truncated = truncate_report_rows(
            [asdict(row) for row in unstable_taxa.taxa],
            limit=budget.max_report_table_rows,
            section_name="unstable-taxa",
            truncated_sections=truncated_sections,
        )
        unstable_clade_rows, unstable_clade_truncated = truncate_report_rows(
            [asdict(row) for row in unstable_clades.clades],
            limit=budget.max_report_table_rows,
            section_name="unstable-clades",
            truncated_sections=truncated_sections,
        )
        conflict_rows, conflict_truncated = truncate_report_rows(
            []
            if clade_conflicts is None
            else [asdict(row) for row in clade_conflicts.conflicts],
            limit=budget.max_report_table_rows,
            section_name="clade-credibility-conflicts",
            truncated_sections=truncated_sections,
        )
        robust_rows, robust_truncated = truncate_report_rows(
            []
            if conclusion_summary is None
            else [asdict(row) for row in conclusion_summary.robust_clades],
            limit=budget.max_report_table_rows,
            section_name="uncertainty-aware-conclusions.robust",
            truncated_sections=truncated_sections,
        )
        uncertain_rows, uncertain_truncated = truncate_report_rows(
            []
            if conclusion_summary is None
            else [asdict(row) for row in conclusion_summary.uncertain_clades],
            limit=budget.max_report_table_rows,
            section_name="uncertainty-aware-conclusions.uncertain",
            truncated_sections=truncated_sections,
        )
        conflicting_rows, conflicting_truncated = truncate_report_rows(
            []
            if conclusion_summary is None
            else [asdict(row) for row in conclusion_summary.conflicting_clades],
            limit=budget.max_report_table_rows,
            section_name="uncertainty-aware-conclusions.conflicting",
            truncated_sections=truncated_sections,
        )
        thinning_rows, thinning_truncated = truncate_report_rows(
            []
            if thinning_sensitivity is None
            else [asdict(row) for row in thinning_sensitivity.rows],
            limit=budget.max_report_table_rows,
            section_name="thinning-sensitivity",
            truncated_sections=truncated_sections,
        )
        consensus_rows, consensus_truncated = truncate_report_rows(
            []
            if consensus_sensitivity is None
            else [asdict(row) for row in consensus_sensitivity.rows],
            limit=budget.max_report_table_rows,
            section_name="consensus-threshold-sensitivity",
            truncated_sections=truncated_sections,
        )
        benchmark_rows, benchmark_truncated = truncate_report_rows(
            [] if benchmark is None else [asdict(row) for row in benchmark.rows],
            limit=budget.max_report_table_rows,
            section_name="tree-set-benchmark",
            truncated_sections=truncated_sections,
        )
        artifact_paths = {
            "tree_set_summary": write_json_artifact(
                artifact_root / "tree-set-summary.json", asdict(summary)
            ),
            "methods_summary": methods_summary_result.output_path,
            "consensus_tree": write_consensus_tree(
                artifact_root / "consensus-tree.nwk", consensus_tree
            ),
            "clade_frequencies": write_clade_frequency_table(
                artifact_root / "clade-frequencies.tsv", clade_frequencies
            ),
            "rf_distance_distribution": write_tree_distance_distribution_table(
                artifact_root / "rf-distance-distribution.tsv", diversity
            ),
            "topology_clusters": write_topology_cluster_table(
                artifact_root / "topology-clusters.tsv", clusters
            ),
            "unstable_taxa": write_tabular_artifact(
                artifact_root / "unstable-taxa.tsv",
                [asdict(row) for row in unstable_taxa.taxa],
            ),
            "unstable_clades": write_unstable_clade_table(
                artifact_root / "unstable-clades.tsv", unstable_clades
            ),
            "clade_credibility_conflicts": write_json_artifact(
                artifact_root / "clade-credibility-conflicts.json",
                scaled_report_note
                if clade_conflicts is None
                else asdict(clade_conflicts),
            ),
            "uncertainty_aware_conclusions": write_json_artifact(
                artifact_root / "uncertainty-aware-conclusions.json",
                scaled_report_note
                if conclusion_summary is None
                else asdict(conclusion_summary),
            ),
            "thinning_sensitivity": write_tabular_artifact(
                artifact_root / "thinning-sensitivity.tsv",
                []
                if thinning_sensitivity is None
                else [asdict(row) for row in thinning_sensitivity.rows],
            ),
            "consensus_threshold_sensitivity": write_tabular_artifact(
                artifact_root / "consensus-threshold-sensitivity.tsv",
                []
                if consensus_sensitivity is None
                else [asdict(row) for row in consensus_sensitivity.rows],
            ),
            "tree_set_benchmark": write_tabular_artifact(
                artifact_root / "tree-set-benchmark.tsv",
                [] if benchmark is None else [asdict(row) for row in benchmark.rows],
            ),
            "topological_diversity": write_json_artifact(
                artifact_root / "topological-diversity.json", asdict(diversity)
            ),
            "topology_multimodality": write_json_artifact(
                artifact_root / "topology-multimodality.json",
                scaled_report_note if multimodality is None else asdict(multimodality),
            ),
            "storage_risk": write_json_artifact(
                artifact_root / "storage-risk.json", asdict(storage_risk)
            ),
            "maturity_gate": write_json_artifact(
                artifact_root / "maturity-gate.json",
                scaled_report_note if maturity is None else asdict(maturity),
            ),
        }
        sections = [
            section(
                "methods-summary-text",
                artifact_paths["methods_summary"].read_text(encoding="utf-8"),
            ),
            section("limitations", limitations),
            section("tree-set-summary", asdict(summary)),
            section(
                "consensus-tree",
                {"newick": dumps_newick(consensus_tree), "report": asdict(consensus)},
            ),
            section(
                "clade-frequencies",
                {
                    "tree_count": clade_frequencies.tree_count,
                    "shared_taxa": clade_frequencies.shared_taxa,
                    "row_count": len(clade_frequencies.clade_frequencies),
                    "truncated_row_count": clade_frequency_truncated,
                    "preview_row_count": min(len(clade_frequency_rows), preview_limit),
                    "preview_rows": preview_report_rows(
                        clade_frequency_rows, limit=preview_limit
                    ),
                    "artifact_path": (
                        artifact_paths["clade_frequencies"]
                        .relative_to(out_path.parent)
                        .as_posix()
                    ),
                },
            ),
            section(
                "rf-distance-distribution",
                {
                    "tree_count": diversity.tree_count,
                    "pair_count": diversity.pair_count,
                    "row_count": len(diversity.rf_distribution),
                    "truncated_row_count": rf_truncated,
                    "preview_row_count": min(len(rf_rows), preview_limit),
                    "preview_rows": preview_report_rows(rf_rows, limit=preview_limit),
                    "artifact_path": (
                        artifact_paths["rf_distance_distribution"]
                        .relative_to(out_path.parent)
                        .as_posix()
                    ),
                },
            ),
            section(
                "topology-clusters",
                {
                    "tree_count": clusters.tree_count,
                    "rooted_topology_count": clusters.rooted_topology_count,
                    "row_count": len(clusters.clusters),
                    "truncated_row_count": cluster_truncated,
                    "preview_row_count": min(len(cluster_rows), preview_limit),
                    "preview_rows": preview_report_rows(
                        cluster_rows, limit=preview_limit
                    ),
                    "artifact_path": (
                        artifact_paths["topology_clusters"]
                        .relative_to(out_path.parent)
                        .as_posix()
                    ),
                },
            ),
            section(
                "topological-diversity",
                {
                    **asdict(diversity),
                    "artifact_path": (
                        artifact_paths["topological_diversity"]
                        .relative_to(out_path.parent)
                        .as_posix()
                    ),
                    "rf_distribution": (
                        f"{len(diversity.rf_distribution)} rows written to linked TSV"
                    ),
                },
            ),
            section(
                "topology-multimodality",
                (
                    {
                        **asdict(multimodality),
                        "artifact_path": (
                            artifact_paths["topology_multimodality"]
                            .relative_to(out_path.parent)
                            .as_posix()
                        ),
                    }
                    if multimodality is not None
                    else {
                        **scaled_report_note,
                        "artifact_path": (
                            artifact_paths["topology_multimodality"]
                            .relative_to(out_path.parent)
                            .as_posix()
                        ),
                    }
                ),
            ),
            section(
                "unstable-taxa",
                {
                    "tree_count": unstable_taxa.tree_count,
                    "row_count": len(unstable_taxa.taxa),
                    "truncated_row_count": unstable_taxa_truncated,
                    "preview_row_count": min(len(unstable_taxa_rows), preview_limit),
                    "preview_rows": preview_report_rows(
                        unstable_taxa_rows, limit=preview_limit
                    ),
                    "artifact_path": (
                        artifact_paths["unstable_taxa"]
                        .relative_to(out_path.parent)
                        .as_posix()
                    ),
                },
            ),
            section(
                "unstable-clades",
                {
                    "tree_count": unstable_clades.tree_count,
                    "row_count": len(unstable_clades.clades),
                    "truncated_row_count": unstable_clade_truncated,
                    "preview_row_count": min(len(unstable_clade_rows), preview_limit),
                    "preview_rows": preview_report_rows(
                        unstable_clade_rows, limit=preview_limit
                    ),
                    "artifact_path": (
                        artifact_paths["unstable_clades"]
                        .relative_to(out_path.parent)
                        .as_posix()
                    ),
                },
            ),
            section(
                "clade-credibility-conflicts",
                {
                    **(
                        {
                            "tree_count": clade_conflicts.tree_count,
                            "credibility_threshold": clade_conflicts.credibility_threshold,
                            "high_credibility_clade_count": clade_conflicts.high_credibility_clade_count,
                            "row_count": len(clade_conflicts.conflicts),
                            "truncated_row_count": conflict_truncated,
                            "preview_row_count": min(len(conflict_rows), preview_limit),
                            "preview_rows": preview_report_rows(
                                conflict_rows, limit=preview_limit
                            ),
                        }
                        if clade_conflicts is not None
                        else scaled_report_note
                    ),
                    "artifact_path": (
                        artifact_paths["clade_credibility_conflicts"]
                        .relative_to(out_path.parent)
                        .as_posix()
                    ),
                },
            ),
            section(
                "uncertainty-aware-conclusions",
                {
                    **(
                        {
                            "tree_count": conclusion_summary.tree_count,
                            "robust_clade_count": conclusion_summary.robust_clade_count,
                            "uncertain_clade_count": conclusion_summary.uncertain_clade_count,
                            "conflicting_clade_count": conclusion_summary.conflicting_clade_count,
                            "robust_rows": preview_report_rows(
                                robust_rows, limit=preview_limit
                            ),
                            "robust_truncated_row_count": robust_truncated,
                            "uncertain_rows": preview_report_rows(
                                uncertain_rows, limit=preview_limit
                            ),
                            "uncertain_truncated_row_count": uncertain_truncated,
                            "conflicting_rows": preview_report_rows(
                                conflicting_rows, limit=preview_limit
                            ),
                            "conflicting_truncated_row_count": conflicting_truncated,
                        }
                        if conclusion_summary is not None
                        else scaled_report_note
                    ),
                    "artifact_path": (
                        artifact_paths["uncertainty_aware_conclusions"]
                        .relative_to(out_path.parent)
                        .as_posix()
                    ),
                },
            ),
            section(
                "storage-risk",
                {
                    **asdict(storage_risk),
                    "artifact_path": (
                        artifact_paths["storage_risk"]
                        .relative_to(out_path.parent)
                        .as_posix()
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
                            "preview_row_count": min(
                                len(consensus_rows), preview_limit
                            ),
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
                            "benchmark_capped": benchmark_tree_count
                            != summary.tree_count
                            or benchmark_taxon_count
                            != max(len(summary.shared_taxa), 2),
                            "row_count": len(benchmark.rows),
                            "truncated_row_count": benchmark_truncated,
                            "preview_row_count": min(
                                len(benchmark_rows), preview_limit
                            ),
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
        core_sections = sections[:11]
        supplemental_sections = sections[11:]
        _current, peak = tracemalloc.get_traced_memory()
        processing = TreeSetProcessingSummary(
            runtime_seconds=round(perf_counter() - started, 6),
            peak_memory_bytes=peak,
            skipped_malformed_tree_count=summary.processing.skipped_malformed_tree_count,
        )
        budget_report = build_tree_set_budget_report(
            budget=budget,
            peak_memory_bytes=processing.peak_memory_bytes,
            truncated_section_names=truncated_sections,
        )
        machine_manifest = {
            "report_kind": "tree-uncertainty",
            "title": title,
            "source_path": str(tree_set_path),
            "input_checksum": sha256(tree_set_path),
            "tree_count": summary.tree_count,
            "rooted_topology_count": summary.rooted_topology_count,
            "processing": asdict(processing),
            "budget": asdict(budget_report),
            "report_mode": "scaled-summary" if scaled_report_mode else "full-review",
            "artifact_root": str(artifact_root),
            "linked_artifact_count": len(artifact_paths),
            "methods_summary_path": artifact_paths["methods_summary"]
            .relative_to(out_path.parent)
            .as_posix(),
            "methods_summary_warning_count": methods_summary_result.warning_count,
            "limitations": limitations,
            "linked_artifacts": {
                name: {
                    "path": path.relative_to(out_path.parent).as_posix(),
                    "byte_count": path.stat().st_size,
                }
                for name, path in artifact_paths.items()
            },
            "sections": [name for name, _ in core_sections],
            "supplemental_sections": [name for name, _ in supplemental_sections],
        }
        summary_metrics = [
            ("tree count", summary.tree_count),
            ("rooted topologies", summary.rooted_topology_count),
            ("report mode", "scaled-summary" if scaled_report_mode else "full-review"),
            ("runtime seconds", processing.runtime_seconds),
            ("peak memory bytes", processing.peak_memory_bytes),
            ("linked artifacts", len(artifact_paths)),
        ]
        artifact_links = [
            (
                name.replace("_", "-"),
                path.relative_to(out_path.parent).as_posix(),
                f"{path.stat().st_size} bytes",
            )
            for name, path in artifact_paths.items()
        ]
        artifact_manifest_path = artifact_root / "tree-uncertainty.manifest.json"
        machine_manifest["artifact_manifest_path"] = artifact_manifest_path.relative_to(
            out_path.parent
        ).as_posix()
        machine_manifest["linked_artifact_count"] = len(artifact_paths) + 1
        machine_manifest["linked_artifacts"]["tree_uncertainty_manifest"] = {
            "path": artifact_manifest_path.relative_to(out_path.parent).as_posix(),
            "byte_count": 0,
        }
        write_json_artifact(artifact_manifest_path, machine_manifest)
        write_html_report(
            title=title,
            sections=sections,
            out_path=out_path,
            embedded_json=machine_manifest,
            summary_metrics=summary_metrics,
            artifact_links=[
                *artifact_links,
                (
                    "tree-uncertainty-manifest",
                    artifact_manifest_path.relative_to(out_path.parent).as_posix(),
                    None,
                ),
            ],
        )
        html_size_bytes = out_path.stat().st_size
        linked_artifact_bytes = sum(
            path.stat().st_size for path in artifact_paths.values()
        )
        manifest_size_bytes = artifact_manifest_path.stat().st_size
        linked_artifact_bytes += manifest_size_bytes
        total_output_bytes = html_size_bytes + linked_artifact_bytes
        machine_manifest["linked_artifacts"]["tree_uncertainty_manifest"] = {
            "path": artifact_manifest_path.relative_to(out_path.parent).as_posix(),
            "byte_count": manifest_size_bytes,
        }
        machine_manifest["html_size_bytes"] = html_size_bytes
        machine_manifest["linked_artifact_bytes"] = linked_artifact_bytes
        machine_manifest["total_output_bytes"] = total_output_bytes
        write_json_artifact(artifact_manifest_path, machine_manifest)
        write_html_report(
            title=title,
            sections=sections,
            out_path=out_path,
            embedded_json=machine_manifest,
            summary_metrics=summary_metrics,
            artifact_links=[
                *artifact_links,
                (
                    "tree-uncertainty-manifest",
                    artifact_manifest_path.relative_to(out_path.parent).as_posix(),
                    f"{artifact_manifest_path.stat().st_size} bytes",
                ),
            ],
        )
        html_size_bytes = out_path.stat().st_size
        total_output_bytes = html_size_bytes + linked_artifact_bytes
        machine_manifest["html_size_bytes"] = html_size_bytes
        machine_manifest["total_output_bytes"] = total_output_bytes
        write_json_artifact(artifact_manifest_path, machine_manifest)
        final_manifest_size_bytes = artifact_manifest_path.stat().st_size
        if final_manifest_size_bytes != manifest_size_bytes:
            linked_artifact_bytes = (
                sum(path.stat().st_size for path in artifact_paths.values())
                + final_manifest_size_bytes
            )
            total_output_bytes = html_size_bytes + linked_artifact_bytes
            machine_manifest["linked_artifacts"]["tree_uncertainty_manifest"] = {
                "path": artifact_manifest_path.relative_to(out_path.parent).as_posix(),
                "byte_count": final_manifest_size_bytes,
            }
            machine_manifest["linked_artifact_bytes"] = linked_artifact_bytes
            machine_manifest["total_output_bytes"] = total_output_bytes
            write_json_artifact(artifact_manifest_path, machine_manifest)
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
