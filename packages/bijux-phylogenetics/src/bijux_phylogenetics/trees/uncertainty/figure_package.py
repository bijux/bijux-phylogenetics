from __future__ import annotations

from dataclasses import asdict, dataclass
from hashlib import sha256
from html import escape
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
from bijux_phylogenetics.render.svg import (
    TreeRenderResult,
    audit_support_label_rendering,
    render_tree_svg,
)

from ..tree_sets import (
    TreeSetProcessingSummary,
    TreeSetWorkflowBudgetReport,
    build_tree_set_budget_report,
    build_tree_set_workflow_budget,
    compute_clade_frequency_table,
    compute_consensus_tree,
    enforce_tree_set_tree_budget,
    load_tree_set,
)
from .instability import (
    detect_unstable_taxa,
    summarize_clade_credibility_conflicts,
    summarize_uncertainty_aware_conclusions,
    write_uncertainty_conclusion_table,
)
from .methods_text import (
    TreeSetUncertaintyMethodReport,
    TreeSetUncertaintyMethodsSummaryTextResult,
    build_tree_set_uncertainty_method_report,
    write_tree_set_uncertainty_methods_summary_text,
)
from .topology_diversity import (
    cluster_trees_by_topology,
    detect_posterior_topology_multimodality,
    write_topology_cluster_table,
)


@dataclass(frozen=True, slots=True)
class TreeSetUncertaintyLegendEntry:
    """One explicit legend entry for a tree-set uncertainty figure package."""

    surface: str
    label: str
    swatch: str
    detail: str


@dataclass(frozen=True, slots=True)
class TreeSetUncertaintyCaptionDraft:
    """Structured caption draft for a tree-set uncertainty figure package."""

    title: str
    lead_sentence: str
    support_sentence: str
    instability_sentence: str
    cluster_sentence: str
    limitation_sentence: str
    caption_ready: bool


@dataclass(frozen=True, slots=True)
class TreeSetUncertaintyPublicationAudit:
    """Reviewer-facing publication audit for tree-set uncertainty figures."""

    publication_ready: bool
    support_labels_validated: bool
    consensus_visible: bool
    clade_support_visible: bool
    unstable_taxa_visible: bool
    topology_clusters_visible: bool
    legend_complete: bool
    caption_ready: bool
    rendered_support_count: int
    plotted_clade_support_count: int
    plotted_unstable_taxon_count: int
    plotted_topology_cluster_count: int
    unstable_taxon_count: int
    topology_cluster_count: int
    reviewer_summary: list[str]
    limitations: list[str]


@dataclass(slots=True)
class TreeSetUncertaintyFigurePackageResult:
    output_dir: Path
    tree_count: int
    processing: TreeSetProcessingSummary
    budget_report: TreeSetWorkflowBudgetReport
    consensus_tree_path: Path
    consensus_figure_path: Path
    clade_support_plot_path: Path
    unstable_taxa_plot_path: Path
    topology_clusters_plot_path: Path
    unstable_taxa_table_path: Path
    topology_clusters_table_path: Path
    uncertainty_conclusions_table_path: Path
    conclusion_summary_path: Path
    legend_path: Path
    caption_path: Path
    methods_summary_path: Path
    review_path: Path
    manifest_path: Path
    reproducibility_manifest_path: Path
    consensus_render: TreeRenderResult
    methods_report: TreeSetUncertaintyMethodReport
    methods_summary: TreeSetUncertaintyMethodsSummaryTextResult
    legend_entries: list[TreeSetUncertaintyLegendEntry]
    caption_draft: TreeSetUncertaintyCaptionDraft
    audit: TreeSetUncertaintyPublicationAudit
    machine_manifest: dict[str, object]


def _sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_ready(payload: object) -> object:
    return json.loads(json.dumps(payload, default=str))


def _write_horizontal_bar_plot(
    path: Path,
    *,
    title: str,
    rows: list[tuple[str, float, str]],
    empty_message: str,
    fill_color: str,
    maximum_rows: int,
) -> int:
    visible_rows = rows[: max(maximum_rows, 1)]
    width = 920
    height = 110 + 44 * max(len(visible_rows), 1)
    bar_left = 320
    bar_width = 520
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
        f'<text x="24" y="36" font-family="Avenir Next, Segoe UI, sans-serif" font-size="22" fill="#0f766e">{escape(title)}</text>',
    ]
    if not visible_rows:
        lines.extend(
            [
                f'<text x="24" y="86" font-family="Avenir Next, Segoe UI, sans-serif" font-size="16" fill="#1f2937">{escape(empty_message)}</text>',
                "</svg>",
            ]
        )
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return 0
    for index, (label, fraction, note) in enumerate(visible_rows):
        y = 70 + index * 44
        filled_width = round(bar_width * max(0.0, min(fraction, 1.0)), 3)
        lines.extend(
            [
                f'<text x="24" y="{y + 16}" font-family="SFMono-Regular, Consolas, monospace" font-size="14" fill="#1f2937">{escape(label)}</text>',
                f'<text x="24" y="{y + 32}" font-family="Avenir Next, Segoe UI, sans-serif" font-size="12" fill="#476b67">{escape(note)}</text>',
                f'<rect x="{bar_left}" y="{y}" width="{bar_width}" height="22" rx="8" fill="#e2e8f0"/>',
                f'<rect x="{bar_left}" y="{y}" width="{filled_width}" height="22" rx="8" fill="{fill_color}"/>',
                f'<text x="{bar_left + bar_width + 12}" y="{y + 17}" font-family="SFMono-Regular, Consolas, monospace" font-size="13" fill="#1f2937">{format(fraction, ".3f")}</text>',
            ]
        )
    lines.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return len(visible_rows)


def _write_clade_support_plot(
    path: Path,
    *,
    clade_frequencies,
    maximum_rows: int,
) -> int:
    return _write_horizontal_bar_plot(
        path,
        title="Tree-Set Clade Support",
        rows=[
            (
                row.clade,
                row.frequency,
                f"{row.tree_count} / {clade_frequencies.tree_count} trees",
            )
            for row in sorted(
                clade_frequencies.clade_frequencies,
                key=lambda row: (-row.frequency, row.clade),
            )
        ],
        empty_message="No informative clades were available for support plotting.",
        fill_color="#0f766e",
        maximum_rows=maximum_rows,
    )


def _write_unstable_taxa_plot(
    path: Path,
    *,
    unstable_taxa,
    maximum_rows: int,
) -> int:
    return _write_horizontal_bar_plot(
        path,
        title="Unstable Taxa Across Tree Set",
        rows=[
            (
                row.taxon,
                row.instability_score,
                f"{row.unique_placements} placements; dominant frequency {format(row.dominant_frequency, '.3f')}",
            )
            for row in unstable_taxa.taxa
        ],
        empty_message="No unstable taxa were detected; every shared taxon kept one placement signature across the tree set.",
        fill_color="#b45309",
        maximum_rows=maximum_rows,
    )


def _write_topology_cluster_plot(
    path: Path,
    *,
    topology_clusters,
    maximum_rows: int,
) -> int:
    return _write_horizontal_bar_plot(
        path,
        title="Topology Clusters Across Tree Set",
        rows=[
            (
                cluster.rooted_topology_id[:12],
                cluster.frequency,
                f"{cluster.tree_count} trees; representative index {cluster.representative_index}",
            )
            for cluster in topology_clusters.clusters
        ],
        empty_message="No topology clusters were available for plotting.",
        fill_color="#1d4ed8",
        maximum_rows=maximum_rows,
    )


def _build_uncertainty_summary_markdown(
    *,
    consensus_newick: str,
    multimodality,
    conflicts,
    conclusions,
) -> str:
    robust = ", ".join(row.clade for row in conclusions.robust_clades) or "none"
    uncertain = ", ".join(row.clade for row in conclusions.uncertain_clades) or "none"
    conflicting = (
        ", ".join(row.clade for row in conclusions.conflicting_clades) or "none"
    )
    return (
        "# Tree-Set Uncertainty Summary\n\n"
        f"- Consensus tree: `{consensus_newick}`\n"
        f"- Topology modes detected: `{multimodality.mode_count}`\n"
        f"- High-credibility clade conflicts: `{conflicts.conflict_count}`\n"
        f"- Robust clades: `{robust}`\n"
        f"- Uncertain clades: `{uncertain}`\n"
        f"- Conflict-prone clades: `{conflicting}`\n"
    )


def _build_legend_entries() -> list[TreeSetUncertaintyLegendEntry]:
    return [
        TreeSetUncertaintyLegendEntry(
            surface="consensus-tree",
            label="consensus branch support",
            swatch="#0f766e",
            detail="consensus branches render support labels directly on the tree figure",
        ),
        TreeSetUncertaintyLegendEntry(
            surface="clade-support",
            label="clade support frequency",
            swatch="#0f766e",
            detail="bar length shows the fraction of trees retaining each clade",
        ),
        TreeSetUncertaintyLegendEntry(
            surface="unstable-taxa",
            label="taxon instability score",
            swatch="#b45309",
            detail="bar length shows one minus the dominant placement frequency for each taxon",
        ),
        TreeSetUncertaintyLegendEntry(
            surface="topology-clusters",
            label="topology cluster frequency",
            swatch="#1d4ed8",
            detail="bar length shows the fraction of trees assigned to each rooted topology cluster",
        ),
    ]


def _write_legend_table(
    path: Path,
    entries: list[TreeSetUncertaintyLegendEntry],
) -> Path:
    return write_taxon_rows(
        path,
        columns=["surface", "label", "swatch", "detail"],
        rows=[
            {
                "surface": entry.surface,
                "label": entry.label,
                "swatch": entry.swatch,
                "detail": entry.detail,
            }
            for entry in entries
        ],
    )


def _build_publication_audit(
    *,
    consensus_render: TreeRenderResult,
    plotted_clade_support_count: int,
    plotted_unstable_taxon_count: int,
    unstable_taxon_count: int,
    plotted_topology_cluster_count: int,
    topology_cluster_count: int,
    legend_entries: list[TreeSetUncertaintyLegendEntry],
) -> TreeSetUncertaintyPublicationAudit:
    support_labels_validated = (
        consensus_render.support_labels_validated
        and consensus_render.rendered_support_count > 0
    )
    consensus_visible = consensus_render.label_count > 0
    clade_support_visible = plotted_clade_support_count > 0
    unstable_taxa_visible = (
        unstable_taxon_count == 0 or plotted_unstable_taxon_count > 0
    )
    topology_clusters_visible = plotted_topology_cluster_count > 0
    legend_complete = {entry.surface for entry in legend_entries} == {
        "consensus-tree",
        "clade-support",
        "unstable-taxa",
        "topology-clusters",
    }
    caption_ready = (
        consensus_visible
        and clade_support_visible
        and unstable_taxa_visible
        and topology_clusters_visible
    )
    publication_ready = (
        support_labels_validated
        and consensus_visible
        and clade_support_visible
        and unstable_taxa_visible
        and topology_clusters_visible
        and legend_complete
        and caption_ready
    )
    reviewer_summary = [
        f"consensus support labels rendered: {consensus_render.rendered_support_count}",
        f"clade-support rows plotted: {plotted_clade_support_count}",
        f"unstable taxa plotted: {plotted_unstable_taxon_count} of {unstable_taxon_count}",
        f"topology clusters plotted: {plotted_topology_cluster_count} of {topology_cluster_count}",
    ]
    limitations: list[str] = []
    if not support_labels_validated:
        limitations.append(
            "consensus support labels were not validated for direct figure use"
        )
    if not clade_support_visible:
        limitations.append("the package does not currently plot clade support rows")
    if not unstable_taxa_visible:
        limitations.append("the package does not currently show unstable taxa clearly")
    if not topology_clusters_visible:
        limitations.append("the package does not currently show topology clusters")
    if not limitations:
        limitations.append(
            "the current package is safe for publication-oriented tree-set uncertainty review"
        )
    return TreeSetUncertaintyPublicationAudit(
        publication_ready=publication_ready,
        support_labels_validated=support_labels_validated,
        consensus_visible=consensus_visible,
        clade_support_visible=clade_support_visible,
        unstable_taxa_visible=unstable_taxa_visible,
        topology_clusters_visible=topology_clusters_visible,
        legend_complete=legend_complete,
        caption_ready=caption_ready,
        rendered_support_count=consensus_render.rendered_support_count,
        plotted_clade_support_count=plotted_clade_support_count,
        plotted_unstable_taxon_count=plotted_unstable_taxon_count,
        plotted_topology_cluster_count=plotted_topology_cluster_count,
        unstable_taxon_count=unstable_taxon_count,
        topology_cluster_count=topology_cluster_count,
        reviewer_summary=reviewer_summary,
        limitations=limitations,
    )


def _build_caption_draft(
    *,
    tree_count: int,
    conclusions,
    audit: TreeSetUncertaintyPublicationAudit,
) -> TreeSetUncertaintyCaptionDraft:
    return TreeSetUncertaintyCaptionDraft(
        title="Tree-set uncertainty review across consensus, clade support, and topology dispersion",
        lead_sentence=(
            f"This tree-set package summarizes uncertainty across {tree_count} trees by combining one consensus tree with explicit support labels, one clade-support panel, one unstable-taxon panel, and one topology-cluster panel."
        ),
        support_sentence=(
            f"Clade-support plotting keeps {audit.plotted_clade_support_count} high-frequency clades visible, while the consensus figure renders {audit.rendered_support_count} direct support labels."
        ),
        instability_sentence=(
            "The instability panel reports how often shared taxa move between alternative placements instead of leaving that uncertainty buried in tables."
        ),
        cluster_sentence=(
            f"The topology-cluster panel shows {audit.plotted_topology_cluster_count} rooted topology groups, and the conclusion summary distinguishes {conclusions.robust_clade_count} robust, {conclusions.uncertain_clade_count} uncertain, and {conclusions.conflicting_clade_count} conflict-prone clades."
        ),
        limitation_sentence=audit.limitations[0],
        caption_ready=audit.caption_ready,
    )


def _write_caption(path: Path, draft: TreeSetUncertaintyCaptionDraft) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                f"# {draft.title}",
                "",
                draft.lead_sentence,
                draft.support_sentence,
                draft.instability_sentence,
                draft.cluster_sentence,
                draft.limitation_sentence,
                "",
                f"caption_ready: {'true' if draft.caption_ready else 'false'}",
                "",
            ]
        ),
        encoding="utf-8",
    )
    return path


def _build_review_html(
    *,
    consensus_figure_path: Path,
    clade_support_plot_path: Path,
    unstable_taxa_plot_path: Path,
    topology_clusters_plot_path: Path,
    unstable_taxa_table_path: Path,
    topology_clusters_table_path: Path,
    uncertainty_conclusions_table_path: Path,
    summary_path: Path,
    caption_path: Path,
    legend_path: Path,
    methods_summary_path: Path,
    methods_summary_text: str,
    audit: TreeSetUncertaintyPublicationAudit,
) -> str:
    figures = {
        "consensus": consensus_figure_path.read_text(encoding="utf-8"),
        "clade_support": clade_support_plot_path.read_text(encoding="utf-8"),
        "unstable_taxa": unstable_taxa_plot_path.read_text(encoding="utf-8"),
        "topology_clusters": topology_clusters_plot_path.read_text(encoding="utf-8"),
    }
    audit_rows = "".join(
        "<tr><th>" + escape(label) + "</th><td>" + escape(value) + "</td></tr>"
        for label, value in [
            ("publication_ready", str(audit.publication_ready).lower()),
            ("support_labels_validated", str(audit.support_labels_validated).lower()),
            ("clade_support_visible", str(audit.clade_support_visible).lower()),
            ("unstable_taxa_visible", str(audit.unstable_taxa_visible).lower()),
            ("topology_clusters_visible", str(audit.topology_clusters_visible).lower()),
        ]
    )
    limitation_items = "".join(f"<li>{escape(item)}</li>" for item in audit.limitations)
    return "\n".join(
        [
            "<!DOCTYPE html>",
            '<html lang="en">',
            "<head>",
            '  <meta charset="utf-8">',
            "  <title>Bijux Tree-Set Uncertainty Review</title>",
            "  <style>",
            "    body { margin: 0; background: linear-gradient(180deg, #eef6f4 0%, #f8fafc 100%); color: #1b1f24; font: 16px/1.5 'Iowan Old Style', 'Palatino Linotype', serif; }",
            "    main { max-width: 1220px; margin: 0 auto; padding: 24px; }",
            "    h1, h2 { font-family: 'Avenir Next', 'Segoe UI', sans-serif; }",
            "    h1 { color: #0f766e; margin-top: 0; }",
            "    .grid { display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 18px; }",
            "    .panel { background: rgba(255,255,255,0.84); border: 1px solid rgba(15,118,110,0.14); border-radius: 18px; padding: 18px; box-shadow: 0 18px 42px rgba(15,118,110,0.08); }",
            "    .figure-shell svg { width: 100%; height: auto; display: block; }",
            "    table { width: 100%; border-collapse: collapse; }",
            "    th, td { text-align: left; padding: 8px 10px; border-bottom: 1px solid rgba(15,118,110,0.12); vertical-align: top; }",
            "    ul { margin: 8px 0 0 18px; }",
            "    a { color: #0f766e; }",
            "  </style>",
            "</head>",
            "<body>",
            "<main>",
            "  <h1>Bijux Tree-Set Uncertainty Review</h1>",
            "  <p>Publication-oriented uncertainty review for one generic tree set. The package keeps consensus support, clade support, unstable taxa, topology clusters, legend, and caption together instead of splitting the uncertainty signal across unrelated outputs.</p>",
            '  <section class="panel">',
            "    <h2>Publication Audit</h2>",
            f"    <table><tbody>{audit_rows}</tbody></table>",
            "    <ul>" + limitation_items + "</ul>",
            "  </section>",
            '  <section class="grid" style="margin-top: 20px;">',
            '    <section class="panel"><h2>Consensus Tree</h2><div class="figure-shell">'
            + figures["consensus"]
            + "</div></section>",
            '    <section class="panel"><h2>Clade Support</h2><div class="figure-shell">'
            + figures["clade_support"]
            + "</div></section>",
            '    <section class="panel"><h2>Unstable Taxa</h2><div class="figure-shell">'
            + figures["unstable_taxa"]
            + "</div></section>",
            '    <section class="panel"><h2>Topology Clusters</h2><div class="figure-shell">'
            + figures["topology_clusters"]
            + "</div></section>",
            "  </section>",
            '  <section class="panel" style="margin-top: 20px;">',
            "    <h2>Methods Summary</h2>",
            f"    <pre>{escape(methods_summary_text)}</pre>",
            "  </section>",
            '  <section class="panel" style="margin-top: 20px;">',
            "    <h2>Linked Artifacts</h2>",
            "    <ul>",
            f'      <li><a href="{escape(unstable_taxa_table_path.name)}">{escape(unstable_taxa_table_path.name)}</a></li>',
            f'      <li><a href="{escape(topology_clusters_table_path.name)}">{escape(topology_clusters_table_path.name)}</a></li>',
            f'      <li><a href="{escape(uncertainty_conclusions_table_path.name)}">{escape(uncertainty_conclusions_table_path.name)}</a></li>',
            f'      <li><a href="{escape(summary_path.name)}">{escape(summary_path.name)}</a></li>',
            f'      <li><a href="{escape(legend_path.name)}">{escape(legend_path.name)}</a></li>',
            f'      <li><a href="{escape(caption_path.name)}">{escape(caption_path.name)}</a></li>',
            f'      <li><a href="{escape(methods_summary_path.name)}">{escape(methods_summary_path.name)}</a></li>',
            "    </ul>",
            "  </section>",
            "</main>",
            "</body>",
            "</html>",
        ]
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
        plotted_clade_support_count = _write_clade_support_plot(
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
        plotted_unstable_taxon_count = _write_unstable_taxa_plot(
            unstable_taxa_plot_path,
            unstable_taxa=unstable_taxa,
            maximum_rows=plot_row_limit,
        )
        write_topology_cluster_table(topology_clusters_table_path, clusters)
        plotted_topology_cluster_count = _write_topology_cluster_plot(
            topology_clusters_plot_path,
            topology_clusters=clusters,
            maximum_rows=plot_row_limit,
        )
        write_uncertainty_conclusion_table(
            uncertainty_conclusions_table_path,
            conclusions,
        )
        conclusion_summary_path.write_text(
            _build_uncertainty_summary_markdown(
                consensus_newick=consensus.consensus_newick,
                multimodality=multimodality,
                conflicts=conflicts,
                conclusions=conclusions,
            ),
            encoding="utf-8",
        )
        legend_entries = _build_legend_entries()
        _write_legend_table(legend_path, legend_entries)
        audit = _build_publication_audit(
            consensus_render=consensus_render,
            plotted_clade_support_count=plotted_clade_support_count,
            plotted_unstable_taxon_count=plotted_unstable_taxon_count,
            unstable_taxon_count=len(unstable_taxa.taxa),
            plotted_topology_cluster_count=plotted_topology_cluster_count,
            topology_cluster_count=len(clusters.clusters),
            legend_entries=legend_entries,
        )
        caption_draft = _build_caption_draft(
            tree_count=summary.tree_count,
            conclusions=conclusions,
            audit=audit,
        )
        _write_caption(caption_path, caption_draft)
        methods_summary = write_tree_set_uncertainty_methods_summary_text(
            methods_summary_path,
            methods_report,
        )
        review_path.write_text(
            _build_review_html(
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
        machine_manifest = {
            "report_kind": "tree_set_uncertainty_figure_package",
            "source_path": str(tree_set_path),
            "input_checksums": {str(tree_set_path): _sha256(tree_set_path)},
            "output_paths": [str(path) for path in artifact_paths],
            "output_checksums": {str(path): _sha256(path) for path in artifact_paths},
            "reproducibility_manifest_path": str(reproducibility_manifest_path),
            "reproducibility_manifest_checksum": _sha256(reproducibility_manifest_path),
            "reproducibility_manifest": reproducibility_manifest,
            "layout": layout,
            "plot_row_limit": plot_row_limit,
            "processing": asdict(processing),
            "budget": asdict(budget_report),
            "consensus": _json_ready(asdict(consensus)),
            "multimodality": _json_ready(asdict(multimodality)),
            "clade_conflicts": _json_ready(asdict(conflicts)),
            "conclusions": _json_ready(asdict(conclusions)),
            "methods_summary": _json_ready(asdict(methods_summary)),
            "audit": _json_ready(asdict(audit)),
            "outputs": {"methods_summary_path": str(methods_summary_path)},
            "metrics": {"methods_summary_warning_count": methods_summary.warning_count},
            "linked_artifact_count": len(artifact_paths),
        }
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
