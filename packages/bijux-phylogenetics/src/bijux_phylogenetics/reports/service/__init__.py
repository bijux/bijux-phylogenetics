from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
from time import perf_counter
import tracemalloc

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
    compare_posterior_topological_diversity,
    compare_posterior_tree_sets,
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
from bijux_phylogenetics.trees.uncertainty_methods import (
    build_tree_set_uncertainty_method_report,
)
from bijux_phylogenetics.validation import (
    build_core_workflow_validation_report,
    build_level_one_release_gate_report,
    build_production_scale_readiness_report,
    build_release_truth_report,
)
from .artifacts import (
    preview_report_rows as _preview_report_rows,
    report_sidecar_path as _report_sidecar_path,
    section as _section,
    truncate_report_rows as _truncate_report_rows,
    write_json_artifact as _write_json_artifact,
    write_machine_manifest as _write_machine_manifest,
    write_tabular_artifact as _write_tabular_artifact,
)
from .distance_reports import render_distance_report
from .input_reports import (
    render_alignment_report,
    render_dataset_report,
    render_phylo_inputs_report,
    render_phylogenetics_report,
    render_tree_report,
)
from .ledger import sha256 as _sha256
from .linkage import (
    annotate_tree_against_table,
    summarise_alignment_path,
    write_annotation_report,
)
from .models import (
    AlignmentReportBuildResult,
    DistanceReportBuildResult,
    ProductionScaleReadinessReportBuildResult,
    ReleaseGateReportBuildResult,
    ReleaseTruthReportBuildResult,
    ReportBuildResult,
    ReportInputLedgerEntry,
    TableLinkageReport,
    TaxonReportBuildResult,
    TreeSetComparisonReportBuildResult,
    TreeUncertaintyReportBuildResult,
    WorkflowValidationReportBuildResult,
)
from .summary import distance_method_limitations
from .taxon_reports import render_taxon_report


__all__ = [
    "AlignmentReportBuildResult",
    "DistanceReportBuildResult",
    "ProductionScaleReadinessReportBuildResult",
    "ReleaseGateReportBuildResult",
    "ReleaseTruthReportBuildResult",
    "ReportBuildResult",
    "ReportInputLedgerEntry",
    "TableLinkageReport",
    "TaxonReportBuildResult",
    "TreeSetComparisonReportBuildResult",
    "TreeUncertaintyReportBuildResult",
    "WorkflowValidationReportBuildResult",
    "annotate_tree_against_table",
    "distance_method_limitations",
    "render_alignment_report",
    "render_dataset_report",
    "render_distance_report",
    "render_level_one_release_gate_report",
    "render_phylogenetics_report",
    "render_phylo_inputs_report",
    "render_production_scale_readiness_report",
    "render_release_truth_report",
    "render_taxon_report",
    "render_tree_report",
    "render_tree_set_comparison_report",
    "render_tree_uncertainty_report",
    "render_workflow_validation_report",
    "summarise_alignment_path",
    "write_annotation_report",
]
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
            dict.fromkeys(
                [
                    "consensus support and topology summaries describe the supplied tree set and should not be treated as direct proof of one true history",
                    "alternative rooted modes, unstable taxa, and conflict-prone clades must remain part of interpretation instead of being collapsed into the consensus tree alone",
                    *methods_summary_result.warnings,
                    *(
                        [scaled_report_note["reason"]]
                        if scaled_report_mode
                        else []
                    ),
                ]
            )
        )
        preview_limit = 5
        clade_frequency_rows, clade_frequency_truncated = _truncate_report_rows(
            [asdict(row) for row in clade_frequencies.clade_frequencies],
            limit=budget.max_report_table_rows,
            section_name="clade-frequencies",
            truncated_sections=truncated_sections,
        )
        rf_rows, rf_truncated = _truncate_report_rows(
            [asdict(row) for row in diversity.rf_distribution],
            limit=budget.max_report_table_rows,
            section_name="rf-distance-distribution",
            truncated_sections=truncated_sections,
        )
        cluster_rows, cluster_truncated = _truncate_report_rows(
            [asdict(row) for row in clusters.clusters],
            limit=budget.max_report_table_rows,
            section_name="topology-clusters",
            truncated_sections=truncated_sections,
        )
        unstable_taxa_rows, unstable_taxa_truncated = _truncate_report_rows(
            [asdict(row) for row in unstable_taxa.taxa],
            limit=budget.max_report_table_rows,
            section_name="unstable-taxa",
            truncated_sections=truncated_sections,
        )
        unstable_clade_rows, unstable_clade_truncated = _truncate_report_rows(
            [asdict(row) for row in unstable_clades.clades],
            limit=budget.max_report_table_rows,
            section_name="unstable-clades",
            truncated_sections=truncated_sections,
        )
        conflict_rows, conflict_truncated = _truncate_report_rows(
            []
            if clade_conflicts is None
            else [asdict(row) for row in clade_conflicts.conflicts],
            limit=budget.max_report_table_rows,
            section_name="clade-credibility-conflicts",
            truncated_sections=truncated_sections,
        )
        robust_rows, robust_truncated = _truncate_report_rows(
            []
            if conclusion_summary is None
            else [asdict(row) for row in conclusion_summary.robust_clades],
            limit=budget.max_report_table_rows,
            section_name="uncertainty-aware-conclusions.robust",
            truncated_sections=truncated_sections,
        )
        uncertain_rows, uncertain_truncated = _truncate_report_rows(
            []
            if conclusion_summary is None
            else [asdict(row) for row in conclusion_summary.uncertain_clades],
            limit=budget.max_report_table_rows,
            section_name="uncertainty-aware-conclusions.uncertain",
            truncated_sections=truncated_sections,
        )
        conflicting_rows, conflicting_truncated = _truncate_report_rows(
            []
            if conclusion_summary is None
            else [asdict(row) for row in conclusion_summary.conflicting_clades],
            limit=budget.max_report_table_rows,
            section_name="uncertainty-aware-conclusions.conflicting",
            truncated_sections=truncated_sections,
        )
        thinning_rows, thinning_truncated = _truncate_report_rows(
            []
            if thinning_sensitivity is None
            else [asdict(row) for row in thinning_sensitivity.rows],
            limit=budget.max_report_table_rows,
            section_name="thinning-sensitivity",
            truncated_sections=truncated_sections,
        )
        consensus_rows, consensus_truncated = _truncate_report_rows(
            []
            if consensus_sensitivity is None
            else [asdict(row) for row in consensus_sensitivity.rows],
            limit=budget.max_report_table_rows,
            section_name="consensus-threshold-sensitivity",
            truncated_sections=truncated_sections,
        )
        benchmark_rows, benchmark_truncated = _truncate_report_rows(
            [] if benchmark is None else [asdict(row) for row in benchmark.rows],
            limit=budget.max_report_table_rows,
            section_name="tree-set-benchmark",
            truncated_sections=truncated_sections,
        )
        artifact_paths = {
            "tree_set_summary": _write_json_artifact(
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
            "unstable_taxa": _write_tabular_artifact(
                artifact_root / "unstable-taxa.tsv",
                [asdict(row) for row in unstable_taxa.taxa],
            ),
            "unstable_clades": write_unstable_clade_table(
                artifact_root / "unstable-clades.tsv", unstable_clades
            ),
            "clade_credibility_conflicts": _write_json_artifact(
                artifact_root / "clade-credibility-conflicts.json",
                (
                    scaled_report_note
                    if clade_conflicts is None
                    else asdict(clade_conflicts)
                ),
            ),
            "uncertainty_aware_conclusions": _write_json_artifact(
                artifact_root / "uncertainty-aware-conclusions.json",
                (
                    scaled_report_note
                    if conclusion_summary is None
                    else asdict(conclusion_summary)
                ),
            ),
            "thinning_sensitivity": _write_tabular_artifact(
                artifact_root / "thinning-sensitivity.tsv",
                (
                    []
                    if thinning_sensitivity is None
                    else [asdict(row) for row in thinning_sensitivity.rows]
                ),
            ),
            "consensus_threshold_sensitivity": _write_tabular_artifact(
                artifact_root / "consensus-threshold-sensitivity.tsv",
                (
                    []
                    if consensus_sensitivity is None
                    else [asdict(row) for row in consensus_sensitivity.rows]
                ),
            ),
            "tree_set_benchmark": _write_tabular_artifact(
                artifact_root / "tree-set-benchmark.tsv",
                [] if benchmark is None else [asdict(row) for row in benchmark.rows],
            ),
            "topological_diversity": _write_json_artifact(
                artifact_root / "topological-diversity.json", asdict(diversity)
            ),
            "topology_multimodality": _write_json_artifact(
                artifact_root / "topology-multimodality.json",
                scaled_report_note if multimodality is None else asdict(multimodality),
            ),
            "storage_risk": _write_json_artifact(
                artifact_root / "storage-risk.json", asdict(storage_risk)
            ),
            "maturity_gate": _write_json_artifact(
                artifact_root / "maturity-gate.json",
                scaled_report_note if maturity is None else asdict(maturity),
            ),
        }
        sections = [
            _section(
                "methods-summary-text",
                artifact_paths["methods_summary"].read_text(encoding="utf-8"),
            ),
            _section("limitations", limitations),
            _section("tree-set-summary", asdict(summary)),
            _section(
                "consensus-tree",
                {"newick": dumps_newick(consensus_tree), "report": asdict(consensus)},
            ),
            _section(
                "clade-frequencies",
                {
                    "tree_count": clade_frequencies.tree_count,
                    "shared_taxa": clade_frequencies.shared_taxa,
                    "row_count": len(clade_frequencies.clade_frequencies),
                    "truncated_row_count": clade_frequency_truncated,
                    "preview_row_count": min(len(clade_frequency_rows), preview_limit),
                    "preview_rows": _preview_report_rows(
                        clade_frequency_rows, limit=preview_limit
                    ),
                    "artifact_path": (
                        artifact_paths["clade_frequencies"]
                        .relative_to(out_path.parent)
                        .as_posix()
                    ),
                },
            ),
            _section(
                "rf-distance-distribution",
                {
                    "tree_count": diversity.tree_count,
                    "pair_count": diversity.pair_count,
                    "row_count": len(diversity.rf_distribution),
                    "truncated_row_count": rf_truncated,
                    "preview_row_count": min(len(rf_rows), preview_limit),
                    "preview_rows": _preview_report_rows(rf_rows, limit=preview_limit),
                    "artifact_path": (
                        artifact_paths["rf_distance_distribution"]
                        .relative_to(out_path.parent)
                        .as_posix()
                    ),
                },
            ),
            _section(
                "topology-clusters",
                {
                    "tree_count": clusters.tree_count,
                    "rooted_topology_count": clusters.rooted_topology_count,
                    "row_count": len(clusters.clusters),
                    "truncated_row_count": cluster_truncated,
                    "preview_row_count": min(len(cluster_rows), preview_limit),
                    "preview_rows": _preview_report_rows(
                        cluster_rows, limit=preview_limit
                    ),
                    "artifact_path": (
                        artifact_paths["topology_clusters"]
                        .relative_to(out_path.parent)
                        .as_posix()
                    ),
                },
            ),
            _section(
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
            _section(
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
            _section(
                "unstable-taxa",
                {
                    "tree_count": unstable_taxa.tree_count,
                    "row_count": len(unstable_taxa.taxa),
                    "truncated_row_count": unstable_taxa_truncated,
                    "preview_row_count": min(len(unstable_taxa_rows), preview_limit),
                    "preview_rows": _preview_report_rows(
                        unstable_taxa_rows, limit=preview_limit
                    ),
                    "artifact_path": (
                        artifact_paths["unstable_taxa"]
                        .relative_to(out_path.parent)
                        .as_posix()
                    ),
                },
            ),
            _section(
                "unstable-clades",
                {
                    "tree_count": unstable_clades.tree_count,
                    "row_count": len(unstable_clades.clades),
                    "truncated_row_count": unstable_clade_truncated,
                    "preview_row_count": min(len(unstable_clade_rows), preview_limit),
                    "preview_rows": _preview_report_rows(
                        unstable_clade_rows, limit=preview_limit
                    ),
                    "artifact_path": (
                        artifact_paths["unstable_clades"]
                        .relative_to(out_path.parent)
                        .as_posix()
                    ),
                },
            ),
            _section(
                "clade-credibility-conflicts",
                {
                    **(
                        {
                            "tree_count": clade_conflicts.tree_count,
                            "credibility_threshold": clade_conflicts.credibility_threshold,
                            "high_credibility_clade_count": (
                                clade_conflicts.high_credibility_clade_count
                            ),
                            "row_count": len(clade_conflicts.conflicts),
                            "truncated_row_count": conflict_truncated,
                            "preview_row_count": min(len(conflict_rows), preview_limit),
                            "preview_rows": _preview_report_rows(
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
            _section(
                "uncertainty-aware-conclusions",
                {
                    **(
                        {
                            "tree_count": conclusion_summary.tree_count,
                            "robust_clade_count": conclusion_summary.robust_clade_count,
                            "uncertain_clade_count": conclusion_summary.uncertain_clade_count,
                            "conflicting_clade_count": (
                                conclusion_summary.conflicting_clade_count
                            ),
                            "robust_rows": _preview_report_rows(
                                robust_rows, limit=preview_limit
                            ),
                            "robust_truncated_row_count": robust_truncated,
                            "uncertain_rows": _preview_report_rows(
                                uncertain_rows, limit=preview_limit
                            ),
                            "uncertain_truncated_row_count": uncertain_truncated,
                            "conflicting_rows": _preview_report_rows(
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
            _section(
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
            _section(
                "thinning-sensitivity",
                {
                    **(
                        {
                            "path": str(thinning_sensitivity.path),
                            "original_tree_count": thinning_sensitivity.original_tree_count,
                            "original_rooted_topology_count": (
                                thinning_sensitivity.original_rooted_topology_count
                            ),
                            "original_dominant_topology_frequency": (
                                thinning_sensitivity.original_dominant_topology_frequency
                            ),
                            "warning_count": len(thinning_sensitivity.warnings),
                            "warnings": thinning_sensitivity.warnings,
                            "row_count": len(thinning_sensitivity.rows),
                            "truncated_row_count": thinning_truncated,
                            "preview_row_count": min(len(thinning_rows), preview_limit),
                            "preview_rows": _preview_report_rows(
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
            _section(
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
                            "preview_rows": _preview_report_rows(
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
            _section(
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
                            "preview_rows": _preview_report_rows(
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
            _section(
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
            "input_checksum": _sha256(tree_set_path),
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
        _write_json_artifact(artifact_manifest_path, machine_manifest)
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
        _write_json_artifact(
            artifact_manifest_path,
            machine_manifest,
        )
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
        _write_json_artifact(artifact_manifest_path, machine_manifest)
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
            _write_json_artifact(artifact_manifest_path, machine_manifest)
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
        _section("tree-set-comparison", asdict(comparison)),
        _section("topological-diversity-comparison", asdict(diversity)),
        _section("left-tree-set-summary", asdict(left_summary)),
        _section("right-tree-set-summary", asdict(right_summary)),
        _section("left-topology-clusters", asdict(left_clusters)),
        _section("right-topology-clusters", asdict(right_clusters)),
        _section("left-topology-multimodality", asdict(left_multimodality)),
        _section("right-topology-multimodality", asdict(right_multimodality)),
        _section("left-unstable-taxa", asdict(left_unstable_taxa)),
        _section("right-unstable-taxa", asdict(right_unstable_taxa)),
        _section("left-unstable-clades", asdict(left_unstable_clades)),
        _section("right-unstable-clades", asdict(right_unstable_clades)),
        _section("left-clade-credibility-conflicts", asdict(left_conflicts)),
        _section("right-clade-credibility-conflicts", asdict(right_conflicts)),
        _section("left-uncertainty-aware-conclusions", asdict(left_conclusions)),
        _section("right-uncertainty-aware-conclusions", asdict(right_conclusions)),
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


def render_workflow_validation_report(
    *,
    out_path: Path,
    fixtures_root: Path | None = None,
) -> WorkflowValidationReportBuildResult:
    """Render the Level 1 workflow validation fixture report."""
    validation = build_core_workflow_validation_report(fixtures_root=fixtures_root)
    title = "Bijux Core Workflow Validation Report"
    reviewer_summary = [
        f"fixture checks passed: {validation.passed_fixture_count}/{validation.total_fixture_count}",
        f"validated workflow surfaces: {len(validation.workflows)}",
        f"known failure-gallery cases: {len(validation.failure_gallery)}",
    ]
    sections = [
        _section("reviewer-summary", reviewer_summary),
        _section(
            "validation-overview",
            {
                "total_fixture_count": validation.total_fixture_count,
                "passed_fixture_count": validation.passed_fixture_count,
                "failed_fixture_count": validation.failed_fixture_count,
            },
        ),
        _section("validation-suites", [asdict(suite) for suite in validation.suites]),
        _section("workflow-coverage", [asdict(row) for row in validation.workflows]),
        _section(
            "failure-gallery", [asdict(row) for row in validation.failure_gallery]
        ),
        _section(
            "maturity-classification",
            [asdict(row) for row in validation.maturity_classifications],
        ),
        _section("limitations", validation.limitations),
    ]
    fixture_paths = sorted(
        {
            path
            for suite in validation.suites
            for fixture in suite.fixtures
            for path in fixture.fixture_paths
        }
    )
    machine_manifest = {
        "report_kind": "workflow-validation",
        "title": title,
        "input_paths": [str(path) for path in fixture_paths],
        "input_checksums": {
            str(path): _sha256(path) for path in fixture_paths if path.exists()
        },
        "sections": [name for name, _ in sections],
        "metrics": {
            "total_fixture_count": validation.total_fixture_count,
            "passed_fixture_count": validation.passed_fixture_count,
            "workflow_count": len(validation.workflows),
        },
        "reviewer_summary": reviewer_summary,
        "limitations": validation.limitations,
    }
    machine_manifest_path = _write_machine_manifest(
        _report_sidecar_path(out_path), machine_manifest
    )
    write_html_report(
        title=title,
        sections=sections,
        out_path=out_path,
        embedded_json=machine_manifest,
    )
    return WorkflowValidationReportBuildResult(
        output_path=out_path,
        machine_manifest_path=machine_manifest_path,
        report_kind="workflow-validation",
        title=title,
        validation=validation,
        machine_manifest=machine_manifest,
    )


def render_level_one_release_gate_report(
    *,
    out_path: Path,
    fixtures_root: Path | None = None,
) -> ReleaseGateReportBuildResult:
    """Render the Level 1 release gate for the checked-in workflow fixtures."""
    release_gate = build_level_one_release_gate_report(fixtures_root=fixtures_root)
    title = "Bijux Level 1 Release Gate Report"
    reviewer_summary = [
        f"gate decision: {release_gate.gate.decision}",
        f"dataset readiness: {release_gate.dataset_readiness_decision}",
        f"retained taxa: {len(release_gate.gate.retained_taxa)}, excluded taxa: {len(release_gate.gate.excluded_taxa)}",
    ]
    sections = [
        _section("reviewer-summary", reviewer_summary),
        _section("gate-decision", asdict(release_gate.gate)),
        _section(
            "dataset-readiness",
            {
                "decision": release_gate.dataset_readiness_decision,
                "blockers": release_gate.dataset_blockers,
                "warnings": release_gate.dataset_warnings,
            },
        ),
        _section(
            "taxon-loss-traceability",
            {
                "first_loss_stage": release_gate.taxon_first_loss_stage,
                "exclusion_causes": release_gate.exclusion_causes,
            },
        ),
        _section("workflow-validation", asdict(release_gate.validation)),
        _section("limitations", release_gate.validation.limitations),
    ]
    fixture_paths = sorted(
        {
            path
            for suite in release_gate.validation.suites
            for fixture in suite.fixtures
            for path in fixture.fixture_paths
        }
    )
    machine_manifest = {
        "report_kind": "release-gate",
        "title": title,
        "input_paths": [str(path) for path in fixture_paths],
        "input_checksums": {
            str(path): _sha256(path) for path in fixture_paths if path.exists()
        },
        "sections": [name for name, _ in sections],
        "metrics": {
            "retained_taxa": len(release_gate.gate.retained_taxa),
            "excluded_taxa": len(release_gate.gate.excluded_taxa),
            "blocked_analysis_count": len(release_gate.gate.blocked_analyses),
        },
        "reviewer_summary": reviewer_summary,
        "limitations": release_gate.validation.limitations,
    }
    machine_manifest_path = _write_machine_manifest(
        _report_sidecar_path(out_path), machine_manifest
    )
    write_html_report(
        title=title,
        sections=sections,
        out_path=out_path,
        embedded_json=machine_manifest,
    )
    return ReleaseGateReportBuildResult(
        output_path=out_path,
        machine_manifest_path=machine_manifest_path,
        report_kind="release-gate",
        title=title,
        release_gate=release_gate,
        machine_manifest=machine_manifest,
    )


def render_release_truth_report(
    *,
    out_path: Path,
    test_report_paths: list[Path],
    real_engine_test_report_paths: list[Path],
    fixtures_root: Path | None = None,
    include_extended_parity: bool = False,
    stress_tier: str = "small",
) -> ReleaseTruthReportBuildResult:
    """Render one machine-produced report of the current release truth surface."""
    release_truth = build_release_truth_report(
        test_report_paths=test_report_paths,
        real_engine_test_report_paths=real_engine_test_report_paths,
        fixtures_root=fixtures_root,
        include_extended_parity=include_extended_parity,
        stress_tier=stress_tier,
    )
    title = "Bijux Release Truth Report"
    reviewer_summary = [
        f"total tests: {release_truth.total_tests.passed_tests} passed, {release_truth.total_tests.failed_tests} failed, {release_truth.total_tests.skipped_tests} skipped",
        f"real-engine tests: {release_truth.real_engine_tests.passed_tests} passed, {release_truth.real_engine_tests.failed_tests} failed, {release_truth.real_engine_tests.skipped_tests} skipped",
        f"supported workflows: {len(release_truth.supported_workflows)}, experimental workflows: {len(release_truth.experimental_workflows)}",
        f"flagship datasets: {len(release_truth.flagship_datasets)}, reference parity cases: {release_truth.reference_parity.case_count}, stress workloads: {len(release_truth.stress_suite.observations)}",
        f"release gate decision: {release_truth.release_gate.gate.decision}",
    ]
    sections = [
        _section("reviewer-summary", reviewer_summary),
        _section("total-tests", asdict(release_truth.total_tests)),
        _section("real-engine-tests", asdict(release_truth.real_engine_tests)),
        _section(
            "supported-workflows",
            [asdict(item) for item in release_truth.supported_workflows],
        ),
        _section(
            "experimental-workflows",
            [asdict(item) for item in release_truth.experimental_workflows],
        ),
        _section(
            "advisory-workflows",
            [asdict(item) for item in release_truth.advisory_workflows],
        ),
        _section(
            "parser-only-workflows",
            [asdict(item) for item in release_truth.parser_only_workflows],
        ),
        _section(
            "flagship-datasets",
            [asdict(item) for item in release_truth.flagship_datasets],
        ),
        _section("workflow-validation", asdict(release_truth.workflow_validation)),
        _section("release-gate", asdict(release_truth.release_gate)),
        _section("reference-parity", asdict(release_truth.reference_parity)),
        _section("stress-suite", asdict(release_truth.stress_suite)),
        _section("known-limitations", release_truth.known_limitations),
    ]
    fixture_paths = sorted(
        {
            path
            for suite in release_truth.workflow_validation.suites
            for fixture in suite.fixtures
            for path in fixture.fixture_paths
        }
    )
    input_paths = [
        *test_report_paths,
        *real_engine_test_report_paths,
        *fixture_paths,
    ]
    machine_manifest = {
        "report_kind": "release-truth",
        "title": title,
        "input_paths": [str(path) for path in input_paths],
        "input_checksums": {
            str(path): _sha256(path) for path in input_paths if path.exists()
        },
        "sections": [name for name, _ in sections],
        "metrics": {
            "total_tests": release_truth.total_tests.total_tests,
            "total_tests_passed": release_truth.total_tests.passed_tests,
            "total_tests_failed": release_truth.total_tests.failed_tests,
            "total_tests_skipped": release_truth.total_tests.skipped_tests,
            "real_engine_tests": release_truth.real_engine_tests.total_tests,
            "real_engine_tests_passed": release_truth.real_engine_tests.passed_tests,
            "real_engine_tests_failed": release_truth.real_engine_tests.failed_tests,
            "real_engine_tests_skipped": release_truth.real_engine_tests.skipped_tests,
            "supported_workflow_count": len(release_truth.supported_workflows),
            "experimental_workflow_count": len(release_truth.experimental_workflows),
            "flagship_dataset_count": len(release_truth.flagship_datasets),
            "reference_parity_case_count": release_truth.reference_parity.case_count,
            "stress_workload_count": len(release_truth.stress_suite.observations),
        },
        "reviewer_summary": reviewer_summary,
        "limitations": release_truth.known_limitations,
    }
    machine_manifest_path = _write_machine_manifest(
        _report_sidecar_path(out_path), machine_manifest
    )
    write_html_report(
        title=title,
        sections=sections,
        out_path=out_path,
        embedded_json=machine_manifest,
    )
    return ReleaseTruthReportBuildResult(
        output_path=out_path,
        machine_manifest_path=machine_manifest_path,
        report_kind="release-truth",
        title=title,
        release_truth=release_truth,
        machine_manifest=machine_manifest,
    )


def render_production_scale_readiness_report(
    *,
    out_path: Path,
    replicates: int = 1,
    tree_tip_counts: list[int] | None = None,
    alignment_size_classes: list[tuple[str, int, int]] | None = None,
    tree_set_size_classes: list[tuple[str, int, int]] | None = None,
    stress_tiers: list[str] | None = None,
) -> ProductionScaleReadinessReportBuildResult:
    """Render one reviewer-facing production-scale readiness report from governed benchmark evidence."""
    production_scale_readiness = build_production_scale_readiness_report(
        replicates=replicates,
        tree_tip_counts=tree_tip_counts,
        alignment_size_classes=alignment_size_classes,
        tree_set_size_classes=tree_set_size_classes,
        stress_tiers=stress_tiers,
    )
    title = "Bijux Production-Scale Readiness Report"
    highest_ready_scale_counts = {
        scale: sum(
            1
            for entry in production_scale_readiness.entries
            if entry.highest_ready_scale == scale
        )
        for scale in sorted(
            {
                "below-small",
                *(
                    threshold.scale
                    for threshold in production_scale_readiness.scale_definitions
                ),
            }
        )
    }
    scale_coverage = [
        {
            "scale": threshold.scale,
            "description": threshold.description,
            "minimum_taxa": threshold.minimum_taxa,
            "minimum_sites": threshold.minimum_sites,
            "minimum_tree_count": threshold.minimum_tree_count,
            "minimum_posterior_size": threshold.minimum_posterior_size,
            "ready_workflow_count": sum(
                1
                for entry in production_scale_readiness.entries
                for decision in entry.scale_decisions
                if decision.scale == threshold.scale and decision.ready
            ),
            "ready_workflows": sorted(
                entry.workflow
                for entry in production_scale_readiness.entries
                for decision in entry.scale_decisions
                if decision.scale == threshold.scale and decision.ready
            ),
        }
        for threshold in production_scale_readiness.scale_definitions
    ]
    reviewer_summary = [
        f"workflow count: {len(production_scale_readiness.entries)}",
        "highest ready scale distribution: "
        + ", ".join(
            f"{scale}={count}"
            for scale, count in highest_ready_scale_counts.items()
            if count > 0
        ),
        f"stress tiers: {', '.join(production_scale_readiness.stress_tiers)}",
    ]
    sections = [
        _section("reviewer-summary", reviewer_summary),
        _section(
            "scale-definitions",
            [asdict(item) for item in production_scale_readiness.scale_definitions],
        ),
        _section("scale-coverage", scale_coverage),
        _section(
            "production-scale-readiness",
            [asdict(item) for item in production_scale_readiness.entries],
        ),
        _section("known-limitations", production_scale_readiness.limitations),
    ]
    machine_manifest = {
        "report_kind": "production-scale-readiness",
        "title": title,
        "input_paths": [],
        "input_checksums": {},
        "sections": [name for name, _ in sections],
        "metrics": {
            "goal_id": production_scale_readiness.goal_id,
            "workflow_count": len(production_scale_readiness.entries),
            "replicates": production_scale_readiness.replicates,
            "stress_tier_count": len(production_scale_readiness.stress_tiers),
            "scale_definition_count": len(
                production_scale_readiness.scale_definitions
            ),
            "below_small_workflow_count": highest_ready_scale_counts.get(
                "below-small", 0
            ),
            **{
                f"{threshold.scale}_ready_workflow_count": sum(
                    1
                    for entry in production_scale_readiness.entries
                    for decision in entry.scale_decisions
                    if decision.scale == threshold.scale and decision.ready
                )
                for threshold in production_scale_readiness.scale_definitions
            },
        },
        "reviewer_summary": reviewer_summary,
        "limitations": production_scale_readiness.limitations,
    }
    machine_manifest_path = _write_machine_manifest(
        _report_sidecar_path(out_path),
        machine_manifest,
    )
    write_html_report(
        title=title,
        sections=sections,
        out_path=out_path,
        embedded_json=machine_manifest,
    )
    return ProductionScaleReadinessReportBuildResult(
        output_path=out_path,
        machine_manifest_path=machine_manifest_path,
        report_kind="production-scale-readiness",
        title=title,
        production_scale_readiness=production_scale_readiness,
        machine_manifest=machine_manifest,
    )
