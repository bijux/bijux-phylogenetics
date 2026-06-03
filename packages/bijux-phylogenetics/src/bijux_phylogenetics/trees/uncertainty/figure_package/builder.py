from __future__ import annotations

import json
from pathlib import Path
from time import perf_counter
import tracemalloc

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.render.reproducibility import (
    FigureReproducibilityFilter,
    write_figure_reproducibility_manifest,
)
from bijux_phylogenetics.render.tree_svg import (
    audit_support_label_rendering,
    render_tree_svg,
)

from ...tree_sets import (
    TreeSetProcessingSummary,
    build_tree_set_budget_report,
    build_tree_set_workflow_budget,
    compute_clade_frequency_table,
    compute_consensus_tree,
    enforce_tree_set_tree_budget,
    load_tree_set,
)
from ..instability import (
    detect_unstable_taxa,
    summarize_clade_credibility_conflicts,
    summarize_uncertainty_aware_conclusions,
    write_uncertainty_conclusion_table,
)
from ..methods_text import (
    build_tree_set_uncertainty_method_report,
    write_tree_set_uncertainty_methods_summary_text,
)
from ..topology_diversity import (
    cluster_trees_by_topology,
    detect_posterior_topology_multimodality,
    write_topology_cluster_table,
)
from .contracts import TreeSetUncertaintyFigurePackageResult
from .manifest import build_machine_manifest
from .plots import (
    write_clade_support_plot,
    write_topology_cluster_plot,
    write_unstable_taxa_plot,
)
from .publication_review import (
    build_caption_draft,
    build_legend_entries,
    build_publication_audit,
    build_review_html,
    build_uncertainty_summary_markdown,
    write_caption,
    write_legend_table,
)


def build_tree_set_uncertainty_figure_package(
    tree_set_path: Path,
    *,
    out_dir: Path,
    layout: str = "phylogram",
    plot_row_limit: int = 12,
    max_tree_count: int | None = None,
    max_report_table_rows: int | None = None,
    memory_warning_threshold_bytes: int | None = None,
) -> TreeSetUncertaintyFigurePackageResult:
    """Build a publication-oriented uncertainty figure package from one tree set."""
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
        enforce_tree_set_tree_budget(
            tree_count=summary.tree_count,
            budget=budget,
            workflow_name="tree-set uncertainty figure package",
            source_path=tree_set_path,
        )
        out_dir.mkdir(parents=True, exist_ok=True)
        consensus_tree, consensus = compute_consensus_tree(tree_set_path)
        clade_frequencies = compute_clade_frequency_table(tree_set_path)
        unstable_taxa = detect_unstable_taxa(tree_set_path)
        clusters = cluster_trees_by_topology(tree_set_path)
        multimodality = detect_posterior_topology_multimodality(tree_set_path)
        conflicts = summarize_clade_credibility_conflicts(tree_set_path)
        conclusions = summarize_uncertainty_aware_conclusions(tree_set_path)
        methods_report = build_tree_set_uncertainty_method_report(tree_set_path)

        consensus_tree_path = out_dir / "consensus-tree.nwk"
        consensus_figure_path = out_dir / "consensus-tree.svg"
        clade_support_plot_path = out_dir / "clade-support-plot.svg"
        unstable_taxa_plot_path = out_dir / "unstable-taxa-plot.svg"
        topology_clusters_plot_path = out_dir / "topology-clusters-plot.svg"
        unstable_taxa_table_path = out_dir / "unstable-taxa.tsv"
        topology_clusters_table_path = out_dir / "topology-clusters.tsv"
        uncertainty_conclusions_table_path = out_dir / "uncertainty-conclusions.tsv"
        conclusion_summary_path = out_dir / "uncertainty-summary.md"
        legend_path = out_dir / "figure-legend.tsv"
        caption_path = out_dir / "figure-caption.md"
        methods_summary_path = out_dir / "tree-set-uncertainty-methods-summary.md"
        review_path = out_dir / "uncertainty-review.html"
        manifest_path = out_dir / "uncertainty-package-manifest.json"
        reproducibility_manifest_path = out_dir / "figure-reproducibility.manifest.json"

        write_newick(consensus_tree_path, consensus_tree)
        support_audit = audit_support_label_rendering(consensus_tree_path)
        consensus_render = render_tree_svg(
            consensus_tree_path,
            out_path=consensus_figure_path,
            layout=layout,
            show_support_values=support_audit.validated,
            validated_support_labels=support_audit.labels_by_node,
            support_validation_warnings=support_audit.warnings,
        )
        plotted_clade_support_count = write_clade_support_plot(
            clade_support_plot_path,
            clade_frequencies=clade_frequencies,
            maximum_rows=plot_row_limit,
        )
        write_taxon_rows(
            unstable_taxa_table_path,
            columns=[
                "taxon",
                "unique_placements",
                "dominant_frequency",
                "instability_score",
                "placement_signatures",
            ],
            rows=[
                {
                    "taxon": row.taxon,
                    "unique_placements": str(row.unique_placements),
                    "dominant_frequency": format(row.dominant_frequency, ".15g"),
                    "instability_score": format(row.instability_score, ".15g"),
                    "placement_signatures": "; ".join(
                        f"{placement.signature} ({format(placement.frequency, '.15g')})"
                        for placement in row.placements
                    ),
                }
                for row in unstable_taxa.taxa
            ],
        )
        plotted_unstable_taxon_count = write_unstable_taxa_plot(
            unstable_taxa_plot_path,
            unstable_taxa=unstable_taxa,
            maximum_rows=plot_row_limit,
        )
        write_topology_cluster_table(topology_clusters_table_path, clusters)
        plotted_topology_cluster_count = write_topology_cluster_plot(
            topology_clusters_plot_path,
            topology_clusters=clusters,
            maximum_rows=plot_row_limit,
        )
        write_uncertainty_conclusion_table(
            uncertainty_conclusions_table_path,
            conclusions,
        )
        conclusion_summary_path.write_text(
            build_uncertainty_summary_markdown(
                consensus_newick=consensus.consensus_newick,
                multimodality=multimodality,
                conflicts=conflicts,
                conclusions=conclusions,
            ),
            encoding="utf-8",
        )
        legend_entries = build_legend_entries()
        write_legend_table(legend_path, legend_entries)
        audit = build_publication_audit(
            consensus_render=consensus_render,
            plotted_clade_support_count=plotted_clade_support_count,
            plotted_unstable_taxon_count=plotted_unstable_taxon_count,
            unstable_taxon_count=len(unstable_taxa.taxa),
            plotted_topology_cluster_count=plotted_topology_cluster_count,
            topology_cluster_count=len(clusters.clusters),
            legend_entries=legend_entries,
        )
        caption_draft = build_caption_draft(
            tree_count=summary.tree_count,
            conclusions=conclusions,
            audit=audit,
        )
        write_caption(caption_path, caption_draft)
        methods_summary = write_tree_set_uncertainty_methods_summary_text(
            methods_summary_path,
            methods_report,
        )
        review_path.write_text(
            build_review_html(
                consensus_figure_path=consensus_figure_path,
                clade_support_plot_path=clade_support_plot_path,
                unstable_taxa_plot_path=unstable_taxa_plot_path,
                topology_clusters_plot_path=topology_clusters_plot_path,
                unstable_taxa_table_path=unstable_taxa_table_path,
                topology_clusters_table_path=topology_clusters_table_path,
                uncertainty_conclusions_table_path=uncertainty_conclusions_table_path,
                summary_path=conclusion_summary_path,
                caption_path=caption_path,
                legend_path=legend_path,
                methods_summary_path=methods_summary_path,
                methods_summary_text=methods_summary.text,
                audit=audit,
            ),
            encoding="utf-8",
        )
        _current, peak = tracemalloc.get_traced_memory()
        processing = TreeSetProcessingSummary(
            runtime_seconds=round(perf_counter() - started, 6),
            peak_memory_bytes=peak,
            skipped_malformed_tree_count=summary.processing.skipped_malformed_tree_count,
        )
        budget_report = build_tree_set_budget_report(
            budget=budget,
            peak_memory_bytes=processing.peak_memory_bytes,
        )
        artifact_paths = [
            consensus_tree_path,
            consensus_figure_path,
            clade_support_plot_path,
            unstable_taxa_plot_path,
            topology_clusters_plot_path,
            unstable_taxa_table_path,
            topology_clusters_table_path,
            uncertainty_conclusions_table_path,
            conclusion_summary_path,
            legend_path,
            caption_path,
            methods_summary_path,
            review_path,
        ]
        reproducibility_manifest = write_figure_reproducibility_manifest(
            reproducibility_manifest_path,
            report_kind="tree_set_uncertainty_figure_package",
            input_files=[("tree_set", tree_set_path)],
            generated_figures=[
                ("consensus_tree", consensus_figure_path),
                ("clade_support", clade_support_plot_path),
                ("unstable_taxa", unstable_taxa_plot_path),
                ("topology_clusters", topology_clusters_plot_path),
            ],
            generated_tables=[
                ("unstable_taxa", unstable_taxa_table_path),
                ("topology_clusters", topology_clusters_table_path),
                ("uncertainty_conclusions", uncertainty_conclusions_table_path),
                ("legend", legend_path),
            ],
            filters=[
                *(
                    [
                        FigureReproducibilityFilter(
                            name="max_tree_count",
                            value=str(max_tree_count),
                            detail="enforced upper tree-count budget before uncertainty package processing",
                        )
                    ]
                    if max_tree_count is not None
                    else []
                ),
                *(
                    [
                        FigureReproducibilityFilter(
                            name="max_report_table_rows",
                            value=str(max_report_table_rows),
                            detail="caps report-facing table volume after uncertainty analysis",
                        )
                    ]
                    if max_report_table_rows is not None
                    else []
                ),
            ]
            or None,
            model={
                "kind": "tree_set_uncertainty",
                "name": "consensus-plus-topology-analysis",
            },
            settings={
                "layout": layout,
                "plot_row_limit": plot_row_limit,
                "tree_count": summary.tree_count,
                "memory_warning_threshold_bytes": memory_warning_threshold_bytes,
            },
            linked_artifacts=[
                ("consensus_tree_newick", consensus_tree_path),
                ("uncertainty_summary", conclusion_summary_path),
                ("caption", caption_path),
                ("methods_summary", methods_summary_path),
                ("review", review_path),
            ],
        )
        machine_manifest = build_machine_manifest(
            tree_set_path=tree_set_path,
            artifact_paths=artifact_paths,
            reproducibility_manifest_path=reproducibility_manifest_path,
            reproducibility_manifest=reproducibility_manifest,
            layout=layout,
            plot_row_limit=plot_row_limit,
            summary=summary,
            processing=processing,
            budget_report=budget_report,
            consensus=consensus,
            multimodality=multimodality,
            conflicts=conflicts,
            conclusions=conclusions,
            methods_summary=methods_summary,
            audit=audit,
            methods_summary_path=methods_summary_path,
        )
        manifest_path.write_text(
            json.dumps(machine_manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return TreeSetUncertaintyFigurePackageResult(
            output_dir=out_dir,
            tree_count=summary.tree_count,
            processing=processing,
            budget_report=budget_report,
            consensus_tree_path=consensus_tree_path,
            consensus_figure_path=consensus_figure_path,
            clade_support_plot_path=clade_support_plot_path,
            unstable_taxa_plot_path=unstable_taxa_plot_path,
            topology_clusters_plot_path=topology_clusters_plot_path,
            unstable_taxa_table_path=unstable_taxa_table_path,
            topology_clusters_table_path=topology_clusters_table_path,
            uncertainty_conclusions_table_path=uncertainty_conclusions_table_path,
            conclusion_summary_path=conclusion_summary_path,
            legend_path=legend_path,
            caption_path=caption_path,
            methods_summary_path=methods_summary_path,
            review_path=review_path,
            manifest_path=manifest_path,
            reproducibility_manifest_path=reproducibility_manifest_path,
            consensus_render=consensus_render,
            methods_report=methods_report,
            methods_summary=methods_summary,
            legend_entries=legend_entries,
            caption_draft=caption_draft,
            audit=audit,
            machine_manifest=machine_manifest,
        )
    finally:
        if not started_tracing:
            tracemalloc.stop()
