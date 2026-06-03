from __future__ import annotations

from html import escape
from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.render.tree_svg import TreeRenderResult

from .contracts import (
    TreeSetUncertaintyCaptionDraft,
    TreeSetUncertaintyLegendEntry,
    TreeSetUncertaintyPublicationAudit,
)


def build_uncertainty_summary_markdown(
    *,
    consensus_newick: str,
    multimodality,
    conflicts,
    conclusions,
) -> str:
    """Build the reviewer-facing uncertainty summary markdown artifact."""
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


def build_legend_entries() -> list[TreeSetUncertaintyLegendEntry]:
    """Build the durable figure legend contract for the package."""
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


def write_legend_table(
    path: Path,
    entries: list[TreeSetUncertaintyLegendEntry],
) -> Path:
    """Write the legend table artifact used by review and publication checks."""
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


def build_publication_audit(
    *,
    consensus_render: TreeRenderResult,
    plotted_clade_support_count: int,
    plotted_unstable_taxon_count: int,
    unstable_taxon_count: int,
    plotted_topology_cluster_count: int,
    topology_cluster_count: int,
    legend_entries: list[TreeSetUncertaintyLegendEntry],
) -> TreeSetUncertaintyPublicationAudit:
    """Assess whether the package is publication-ready for uncertainty review."""
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


def build_caption_draft(
    *,
    tree_count: int,
    conclusions,
    audit: TreeSetUncertaintyPublicationAudit,
) -> TreeSetUncertaintyCaptionDraft:
    """Build the durable caption draft for the figure package."""
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


def write_caption(path: Path, draft: TreeSetUncertaintyCaptionDraft) -> Path:
    """Write the reviewer-facing caption artifact."""
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


def build_review_html(
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
    """Build the reviewer HTML that assembles the package artifacts."""
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
