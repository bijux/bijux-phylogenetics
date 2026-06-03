from __future__ import annotations

from html import escape
from pathlib import Path

from ..models import (
    _FLAGSHIP_QUESTION,
    RabiesCrossHostGeographyPanelDataset,
    RabiesCrossHostGeographyPanelExportResult,
    RabiesCrossHostGeographyPanelWorkflowBundle,
    RabiesCrossHostGeographyPanelWorkflowConfig,
)
from ..shared import _html_list


def _write_demo_overview_html(
    path: Path,
    *,
    dataset: RabiesCrossHostGeographyPanelDataset,
    dataset_export: RabiesCrossHostGeographyPanelExportResult,
    workflow_bundle: RabiesCrossHostGeographyPanelWorkflowBundle,
    config: RabiesCrossHostGeographyPanelWorkflowConfig,
    short_answer: str,
    artifact_inventory_path: Path,
    reproducibility_checklist_path: Path,
) -> Path:
    html = "\n".join(
        [
            "<!DOCTYPE html>",
            '<html lang="en">',
            "<head>",
            '  <meta charset="utf-8">',
            "  <title>Rabies Cross-Host Geography Package Overview</title>",
            "  <style>",
            "    body { font-family: Georgia, 'Times New Roman', serif; margin: 0; background: linear-gradient(180deg, #f4f1ea 0%, #eef5ef 100%); color: #173024; }",
            "    main { max-width: 1040px; margin: 0 auto; padding: 28px; }",
            "    h1, h2 { margin: 0 0 10px; }",
            "    p { line-height: 1.6; }",
            "    .panel { background: rgba(255,255,255,0.9); border: 1px solid rgba(23,48,36,0.12); border-radius: 18px; padding: 18px; margin-top: 18px; box-shadow: 0 14px 36px rgba(23,48,36,0.08); }",
            "    .cards { display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 14px; margin-top: 18px; }",
            "    .card { background: rgba(255,255,255,0.9); border: 1px solid rgba(23,48,36,0.12); border-radius: 16px; padding: 16px; }",
            "    .label { color: #5f7469; font-size: 13px; text-transform: uppercase; letter-spacing: 0.04em; }",
            "    .value { display: block; font-size: 21px; margin-top: 6px; }",
            "    ul { margin: 8px 0 0 18px; }",
            "    code { background: rgba(23,48,36,0.06); padding: 0 4px; border-radius: 4px; }",
            "    a { color: #23523b; }",
            "  </style>",
            "</head>",
            "<body>",
            "<main>",
            "  <h1>Rabies Cross-Host Geography Package</h1>",
            "  <p>This public flagship workflow starts from one real rabies nucleoprotein FASTA panel plus grouped host and geographic metadata, reruns the full owned sequence-to-tree workflow, and then connects the rooted tree to host-switching, geographic transition, migration, and comparative review surfaces.</p>",
            '  <section class="panel">',
            "    <h2>Biological Question</h2>",
            f"    <p>{escape(_FLAGSHIP_QUESTION)}</p>",
            '    <h2 style="margin-top: 16px;">Short Answer</h2>',
            f"    <p>{escape(short_answer)}</p>",
            "  </section>",
            '  <section class="cards">',
            f'    <div class="card"><span class="label">dataset id</span><span class="value">{escape(dataset.dataset_id)}</span></div>',
            f'    <div class="card"><span class="label">sequence count</span><span class="value">{dataset.sequence_count}</span></div>',
            f'    <div class="card"><span class="label">selected model</span><span class="value">{escape(workflow_bundle.selected_model)}</span></div>',
            f'    <div class="card"><span class="label">comparative model</span><span class="value">{escape(workflow_bundle.comparative_selected_model)}</span></div>',
            "  </section>",
            '  <section class="panel">',
            "    <h2>Reproducibility Surface</h2>",
            _html_list(
                [
                    f'workflow config: <a href="dataset/{dataset_export.workflow_config_path.name}">dataset/{dataset_export.workflow_config_path.name}</a>',
                    f'source accession ledger: <a href="dataset/{dataset_export.accession_table_path.name}">dataset/{dataset_export.accession_table_path.name}</a>',
                    f'final workflow manifest: <a href="workflow/{workflow_bundle.final_manifest_path.name}">workflow/{workflow_bundle.final_manifest_path.name}</a>',
                    f'package artifact inventory: <a href="{artifact_inventory_path.name}">{artifact_inventory_path.name}</a>',
                    f'package reproducibility checklist: <a href="{reproducibility_checklist_path.name}">{reproducibility_checklist_path.name}</a>',
                    f'package manifest: <a href="{path.name.replace("-overview.html", "-package.manifest.json")}">{path.name.replace("-overview.html", "-package.manifest.json")}</a>',
                ]
            ),
            "  </section>",
            '  <section class="panel">',
            "    <h2>Key Outputs</h2>",
            _html_list(
                [
                    f'<a href="workflow/{workflow_bundle.final_report_path.name}">workflow/{workflow_bundle.final_report_path.name}</a>',
                    f'<a href="workflow/{workflow_bundle.workflow_summary_path.name}">workflow/{workflow_bundle.workflow_summary_path.name}</a>',
                    f'<a href="workflow/bootstrap-review/{workflow_bundle.bootstrap_tree_comparison_summary_path.name}">workflow/bootstrap-review/{workflow_bundle.bootstrap_tree_comparison_summary_path.name}</a>',
                    f'<a href="workflow/comparative/{workflow_bundle.comparative_report_path.name}">workflow/comparative/{workflow_bundle.comparative_report_path.name}</a>',
                    f'<a href="workflow/conclusion-stability/{workflow_bundle.conclusion_stability_report_path.name}">workflow/conclusion-stability/{workflow_bundle.conclusion_stability_report_path.name}</a>',
                    f'<a href="workflow/{workflow_bundle.scientific_findings_path.name}">workflow/{workflow_bundle.scientific_findings_path.name}</a>',
                    f'<a href="dataset/{dataset_export.workflow_config_path.name}">dataset/{dataset_export.workflow_config_path.name}</a>',
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
