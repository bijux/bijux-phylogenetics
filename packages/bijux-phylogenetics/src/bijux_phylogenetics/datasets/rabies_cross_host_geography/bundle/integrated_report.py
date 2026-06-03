from __future__ import annotations

from html import escape
from pathlib import Path

from bijux_phylogenetics.biogeography import BiogeographyReportPackageResult
from bijux_phylogenetics.comparative.reporting.analysis_package import (
    ComparativeAnalysisSummaryRow,
    ComparativeInterpretationRow,
)
from bijux_phylogenetics.compare.presentation import ComparisonReportBuildResult
from bijux_phylogenetics.ecology import HostSwitchingReport
from bijux_phylogenetics.engines.inference import FastaToTreeWorkflowReport
from bijux_phylogenetics.trees import BootstrapTreeSetArtifactReport

from ..models import (
    RabiesCrossHostGeographyPanelWorkflowReport,
    RabiesScientificFindingRow,
)
from ..shared import _format_number, _html_list, _support_range_text, _table


def _write_integrated_report(
    path: Path,
    *,
    report: RabiesCrossHostGeographyPanelWorkflowReport,
    workflow_summary_path: Path,
    bootstrap_artifacts: BootstrapTreeSetArtifactReport,
    bootstrap_tree_comparison_report: ComparisonReportBuildResult,
    clade_row_count: int,
    comparative_summary_row: ComparativeAnalysisSummaryRow,
    comparative_interpretation_rows: list[ComparativeInterpretationRow],
    comparative_branch_repair_count: int,
    scientific_finding_rows: list[RabiesScientificFindingRow],
    max_report_table_rows: int | None,
) -> Path:
    support_summary = report.fasta_to_tree.support_summary
    host_summary = report.host_switching.summary
    geography_summary = report.biogeography_report.state_report.summary
    migration_summary = report.biogeography_report.event_report.summary
    bootstrap_summary = bootstrap_artifacts.summary_report
    core_question = (
        "Do the host-associated rabies lineages in this compact panel occupy one "
        "distinct geographic regime while retaining one coherent phylogenetic signal?"
    )
    core_answer = next(
        (
            row.claim
            for row in comparative_interpretation_rows
            if row.topic == "coefficient" and "nominally supported" in row.claim
        ),
        "the comparative layer did not recover a nominally supported host effect",
    )
    html = "\n".join(
        [
            "<!DOCTYPE html>",
            '<html lang="en">',
            "<head>",
            '  <meta charset="utf-8">',
            "  <title>Bijux Rabies Host and Geography Workflow</title>",
            "  <style>",
            "    body { font-family: Georgia, 'Times New Roman', serif; margin: 0; background: linear-gradient(180deg, #f4f1ea 0%, #e7efe7 100%); color: #163222; }",
            "    main { max-width: 1360px; margin: 0 auto; padding: 24px; }",
            "    h1 { margin: 0 0 8px; font-size: 34px; }",
            "    h2 { margin: 0 0 10px; font-size: 22px; }",
            "    p { line-height: 1.55; }",
            "    .cards { display: grid; grid-template-columns: repeat(6, minmax(0, 1fr)); gap: 14px; margin: 18px 0 24px; }",
            "    .card, .panel { background: rgba(255,255,255,0.86); border: 1px solid rgba(22,50,34,0.12); border-radius: 18px; padding: 18px; box-shadow: 0 16px 42px rgba(22,50,34,0.08); }",
            "    .label { color: #5b7466; font-size: 13px; text-transform: uppercase; letter-spacing: 0.04em; }",
            "    .card strong { display: block; font-size: 21px; margin-top: 6px; }",
            "    .grid { display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }",
            "    .full { grid-column: 1 / -1; }",
            "    .figure-shell { overflow: auto; }",
            "    .figure-shell img { width: 100%; height: auto; display: block; }",
            "    table { width: 100%; border-collapse: collapse; margin-top: 10px; font-size: 14px; }",
            "    th, td { border-bottom: 1px solid rgba(22,50,34,0.10); padding: 8px 10px; text-align: left; vertical-align: top; }",
            "    th { color: #365443; }",
            "    ul { margin: 8px 0 0 18px; }",
            "    a { color: #16543a; }",
            "    iframe { width: 100%; min-height: 760px; border: 1px solid rgba(22,50,34,0.12); border-radius: 14px; background: white; }",
            "  </style>",
            "</head>",
            "<body>",
            "<main>",
            "  <h1>Bijux Rabies Host and Geography Workflow</h1>",
            "  <p>Complete end-to-end review for one real rabies nucleoprotein panel. The workflow starts from raw sequences plus combined host and geography metadata, validates the FASTA surface, aligns and trims the panel, infers a bootstrap-supported maximum-likelihood tree, roots that tree on one explicit outgroup, summarizes bootstrap topology uncertainty, extracts clades, reconstructs host and geographic histories, and fits one comparative model over a derived geographic trait.</p>",
            '  <section class="cards">',
            f'    <div class="card"><span class="label">sequences</span><strong>{report.dataset.sequence_count}</strong></div>',
            f'    <div class="card"><span class="label">selected model</span><strong>{escape(report.fasta_to_tree.selected_model)}</strong></div>',
            f'    <div class="card"><span class="label">aligned quality</span><strong>{_format_number(report.aligned_quality.quality_score)}</strong></div>',
            f'    <div class="card"><span class="label">trimmed quality</span><strong>{_format_number(report.trimmed_quality.quality_score)}</strong></div>',
            f'    <div class="card"><span class="label">root host</span><strong>{escape(host_summary.root_host)}</strong></div>',
            f'    <div class="card"><span class="label">root region</span><strong>{escape(geography_summary.root_region)}</strong></div>',
            "  </section>",
            '  <section class="panel">',
            "    <h2>Scientific Question</h2>",
            f"    <p>{escape(core_question)}</p>",
            f"    <p><strong>Working answer:</strong> {escape(core_answer)}. The comparative layer selects {escape(comparative_summary_row.selected_model)} as the better continuous-trait surface, but the residual diagnostics remain cautionary and the panel is intentionally small.</p>",
            _html_list(
                [
                    f"FASTA validation resolved the raw sequence type as {report.fasta_to_tree.sequence_type}.",
                    f"Bootstrap support spans {_support_range_text(support_summary.minimum_support, support_summary.maximum_support)} across the final rooted tree.",
                    f"Host reconstruction inferred {host_summary.host_switch_count} host-switch branches, with {host_summary.certain_host_switch_count} certain and {host_summary.uncertain_host_switch_count} uncertain changes.",
                    f"Geographic reconstruction inferred {migration_summary.event_count} migration events across {geography_summary.changed_branch_count} changed branches.",
                    f"Bootstrap replicate review retained {bootstrap_summary.tree_count} trees across {bootstrap_summary.diversity.rooted_topology_count} rooted topologies.",
                    (
                        "The rooted ML tree versus bootstrap consensus comparison "
                        f"returned rooted RF distance {bootstrap_tree_comparison_report.topology.rooted_robinson_foulds_distance}."
                    ),
                    f"The clade table contains {clade_row_count} node rows and the comparative tree required {comparative_branch_repair_count} explicit branch-length repair(s).",
                    (
                        "Bootstrap review emitted budget warnings: "
                        + "; ".join(bootstrap_artifacts.budget_report.warning_messages)
                    )
                    if bootstrap_artifacts.budget_report.warning_messages
                    else (
                        "Configured workflow budgets covered the bootstrap review "
                        "without tree-count failure or peak-memory warning."
                    ),
                ]
            ),
            "  </section>",
            '  <section class="grid" style="margin-top: 20px;">',
            '    <section class="panel">',
            "      <h2>Sequence-to-Tree Outputs</h2>",
            _html_list(
                [
                    "input validation: input-validation.tsv",
                    "alignment quality: alignment-quality.tsv",
                    "alignment sequence ranking: alignment-sequence-ranking.tsv",
                    "alignment: rabies-cross-host-geography-panel.aln",
                    "trimmed alignment: rabies-cross-host-geography-panel.trimmed.aln",
                    "rooted tree: rabies-cross-host-geography-panel.rooted.tree",
                    "support table: rabies-cross-host-geography-panel.support.tsv",
                    "workflow summary: workflow-summary.tsv",
                    "resource observations: resource-observations.tsv",
                ]
            ),
            _support_table(report.fasta_to_tree, max_rows=max_report_table_rows),
            "    </section>",
            '    <section class="panel">',
            "      <h2>Bootstrap and Clade Review</h2>",
            _html_list(
                [
                    f"bootstrap tree count: {bootstrap_summary.tree_count}",
                    f"rooted topology count: {bootstrap_summary.diversity.rooted_topology_count}",
                    f"unstable branch count: {bootstrap_summary.unstable_branch_count}",
                    f"clade row count: {clade_row_count}",
                    (
                        "see bootstrap-review/ for consensus, clade frequencies, "
                        "instability, distances, topology clusters, and rooted-tree comparison"
                    ),
                ]
            ),
            _table(
                headers=[
                    "tree_count",
                    "rooted_topology_count",
                    "dominant_topology_frequency",
                    "effective_topology_count",
                    "unstable_branch_count",
                ],
                rows=[
                    [
                        str(bootstrap_summary.tree_count),
                        str(bootstrap_summary.diversity.rooted_topology_count),
                        _format_number(
                            bootstrap_summary.diversity.dominant_topology_frequency
                        ),
                        _format_number(
                            bootstrap_summary.diversity.effective_topology_count
                        ),
                        str(bootstrap_summary.unstable_branch_count),
                    ]
                ],
            ),
            "    </section>",
            '    <section class="panel">',
            "      <h2>Host Switching</h2>",
            _html_list(
                [
                    f"workflow trait: {report.dataset.host_trait}",
                    f"root host confidence: {_format_number(host_summary.root_confidence)}",
                    f"host-switch rows: {len(report.host_switching.count_rows)}",
                    "see host-switch-summary.tsv, host-state-nodes.tsv, host-switch-branches.tsv, and host-switch-counts.tsv",
                ]
            ),
            _host_count_table(report.host_switching, max_rows=max_report_table_rows),
            "    </section>",
            '    <section class="panel">',
            "      <h2>Comparative Layer</h2>",
            _html_list(
                [
                    f"formula: {comparative_summary_row.formula}",
                    f"selected trait model: {comparative_summary_row.selected_model}",
                    f"PGLS lambda: {_format_number(comparative_summary_row.pgls_lambda)}",
                    f"PGLS r-squared: {_format_number(comparative_summary_row.pgls_r_squared)}",
                    "see comparative/ for coefficients, model comparison, diagnostics, signal summary, and interpretation tables",
                ]
            ),
            _table(
                headers=["topic", "claim", "evidence"],
                rows=[
                    [row.topic, row.claim, row.evidence]
                    for row in comparative_interpretation_rows[:5]
                ],
                max_rows=max_report_table_rows,
            ),
            "    </section>",
            '    <section class="panel full">',
            "      <h2>Biogeography</h2>",
            '      <p>The bundle includes the detailed biogeography package at <a href="biogeography/biogeography-report.html">biogeography/biogeography-report.html</a> together with the ancestral-region tree SVG and the self-contained geographic map.</p>',
            '      <div class="grid">',
            '        <div class="panel">',
            "          <h2>Ancestral-Region Tree</h2>",
            '          <div class="figure-shell">',
            '            <img src="biogeography/ancestral-region-tree.svg" alt="Ancestral region tree">',
            "          </div>",
            "        </div>",
            '        <div class="panel">',
            "          <h2>Geographic Map</h2>",
            '          <iframe src="biogeography/geographic-region-map.html" title="Geographic region map"></iframe>',
            "        </div>",
            "      </div>",
            _migration_event_table(
                report.biogeography_report,
                max_rows=max_report_table_rows,
            ),
            "    </section>",
            '    <section class="panel full">',
            "      <h2>Scientific Findings Ledger</h2>",
            _table(
                headers=[
                    "finding_id",
                    "question",
                    "claim",
                    "evidence",
                    "caution",
                    "source_artifact",
                ],
                rows=[
                    [
                        row.finding_id,
                        row.question,
                        row.claim,
                        row.evidence,
                        row.caution,
                        row.source_artifact,
                    ]
                    for row in scientific_finding_rows
                ],
                max_rows=max_report_table_rows,
            ),
            "    </section>",
            "  </section>",
            '  <section class="panel" style="margin-top: 20px;">',
            "    <h2>Key Files</h2>",
            _html_list(
                [
                    f'<a href="{workflow_summary_path.name}">{workflow_summary_path.name}</a>',
                    '<a href="resource-observations.tsv">resource-observations.tsv</a>',
                    '<a href="workflow-config-audit.tsv">workflow-config-audit.tsv</a>',
                    '<a href="clade-table.tsv">clade-table.tsv</a>',
                    '<a href="bootstrap-review/bootstrap-review.summary.tsv">bootstrap-review/bootstrap-review.summary.tsv</a>',
                    '<a href="bootstrap-review/rooted-tree-vs-bootstrap-consensus.summary.tsv">bootstrap-review/rooted-tree-vs-bootstrap-consensus.summary.tsv</a>',
                    '<a href="comparative/comparative-report.html">comparative/comparative-report.html</a>',
                    '<a href="comparative/interpretation-table.tsv">comparative/interpretation-table.tsv</a>',
                    '<a href="scientific-findings.tsv">scientific-findings.tsv</a>',
                    '<a href="host-switch-summary.tsv">host-switch-summary.tsv</a>',
                    '<a href="biogeography/event-table.tsv">biogeography/event-table.tsv</a>',
                    '<a href="rabies-cross-host-geography.manifest.json">rabies-cross-host-geography.manifest.json</a>',
                ]
            ),
            "  </section>",
            "</main>",
            "</body>",
            "</html>",
        ]
    )
    path.write_text(html + "\n", encoding="utf-8")
    return path


def _support_table(
    report: FastaToTreeWorkflowReport,
    *,
    max_rows: int | None = None,
) -> str:
    return _table(
        headers=["node", "descendant_taxa", "support", "support_fraction"],
        rows=[
            [
                row.node,
                ", ".join(row.descendant_taxa),
                _format_number(row.support),
                _format_number(row.support_fraction),
            ]
            for row in report.support_rows
        ],
        max_rows=max_rows,
    )


def _host_count_table(
    report: HostSwitchingReport,
    *,
    max_rows: int | None = None,
) -> str:
    return _table(
        headers=[
            "transition",
            "certain_switch_count",
            "uncertain_switch_count",
            "total_switch_count",
        ],
        rows=[
            [
                row.transition,
                str(row.certain_switch_count),
                str(row.uncertain_switch_count),
                str(row.total_switch_count),
            ]
            for row in report.count_rows
        ],
        max_rows=max_rows,
    )


def _migration_event_table(
    report: BiogeographyReportPackageResult,
    *,
    max_rows: int | None = None,
) -> str:
    return _table(
        headers=[
            "branch_id",
            "source_region",
            "target_region",
            "support",
            "midpoint_depth",
        ],
        rows=[
            [
                row.branch_id,
                row.source_region,
                row.target_region,
                _format_number(row.support),
                _format_number(row.midpoint_depth),
            ]
            for row in report.event_report.event_rows
        ],
        max_rows=max_rows,
    )
