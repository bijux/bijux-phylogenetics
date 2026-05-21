from __future__ import annotations

from pathlib import Path

from ..models import (
    _FLAGSHIP_QUESTION,
    RabiesCrossHostGeographyPanelDataset,
    RabiesCrossHostGeographyPanelWorkflowBundle,
    RabiesCrossHostGeographyPanelWorkflowConfig,
)


def _write_overview(
    path: Path,
    *,
    dataset: RabiesCrossHostGeographyPanelDataset,
    workflow_bundle: RabiesCrossHostGeographyPanelWorkflowBundle,
    config: RabiesCrossHostGeographyPanelWorkflowConfig,
    short_answer: str,
    artifact_inventory_path: Path,
    reproducibility_checklist_path: Path,
) -> Path:
    lines = [
        "# Rabies Cross-Host Geography Demo",
        "",
        f"- dataset id: `{dataset.dataset_id}`",
        f"- sequence count: `{dataset.sequence_count}`",
        f"- workflow config: `{config.config_path.name}`",
        f"- biological question: {_FLAGSHIP_QUESTION}",
        f"- short answer: {short_answer}",
        f"- host workflow trait: `{dataset.host_trait}`",
        f"- geography workflow trait: `{dataset.geography_trait}`",
        f"- comparative formula: `{workflow_bundle.comparative_formula}`",
        "",
        "Generated outputs:",
        "",
        "- source accession ledger: `dataset/source-accessions.tsv`",
        f"- workflow summary: `{workflow_bundle.workflow_summary_path.name}`",
        f"- resource observations: `{workflow_bundle.resource_observations_path.name}`",
        f"- clade table: `{workflow_bundle.clade_table_path.name}`",
        f"- bootstrap review: `bootstrap-review/{workflow_bundle.bootstrap_summary_path.name}`",
        (
            "- rooted-versus-consensus comparison: "
            f"`bootstrap-review/{workflow_bundle.bootstrap_tree_comparison_summary_path.name}`"
        ),
        f"- comparative report: `comparative/{workflow_bundle.comparative_report_path.name}`",
        (
            "- conclusion stability report: "
            f"`conclusion-stability/{workflow_bundle.conclusion_stability_report_path.name}`"
        ),
        f"- final report: `{workflow_bundle.final_report_path.name}`",
        f"- package artifact inventory: `{artifact_inventory_path.name}`",
        f"- package reproducibility checklist: `{reproducibility_checklist_path.name}`",
        "- package overview html: `rabies-cross-host-geography-overview.html`",
        "- package manifest: `rabies-cross-host-geography-package.manifest.json`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _build_flagship_answer_summary(
    workflow_bundle: RabiesCrossHostGeographyPanelWorkflowBundle,
) -> str:
    return (
        "The rooted panel remains anchored in "
        f"{workflow_bundle.root_host} and {workflow_bundle.root_region}, and "
        "`host_group[canid]` shows a nominally supported positive longitude "
        "association under the selected comparative model, but the inference "
        "remains cautionary because the panel is intentionally compact."
    )
