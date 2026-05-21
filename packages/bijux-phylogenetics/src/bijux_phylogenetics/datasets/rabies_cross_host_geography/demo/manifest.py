from __future__ import annotations

import json
from pathlib import Path

from ..models import (
    _FLAGSHIP_QUESTION,
    RabiesCrossHostGeographyPanelDataset,
    RabiesCrossHostGeographyPanelExportResult,
    RabiesCrossHostGeographyPanelWorkflowBundle,
    RabiesCrossHostGeographyPanelWorkflowConfig,
)
from ..shared import _checksum


def _package_inventory_counts(
    inventory_rows: list[dict[str, str]],
) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in inventory_rows:
        section = row["section"]
        counts[section] = counts.get(section, 0) + 1
    return counts


def _write_demo_package_manifest(
    path: Path,
    *,
    dataset: RabiesCrossHostGeographyPanelDataset,
    dataset_export: RabiesCrossHostGeographyPanelExportResult,
    workflow_bundle: RabiesCrossHostGeographyPanelWorkflowBundle,
    config: RabiesCrossHostGeographyPanelWorkflowConfig,
    short_answer: str,
    artifact_inventory_path: Path,
    artifact_inventory_rows: list[dict[str, str]],
    reproducibility_checklist_path: Path,
    checklist_rows: list[dict[str, str]],
) -> Path:
    inventory_counts = _package_inventory_counts(artifact_inventory_rows)
    blocked_check_count = len(
        [row for row in checklist_rows if row["status"] == "blocked"]
    )
    risk_check_count = len([row for row in checklist_rows if row["status"] == "risk"])
    payload = {
        "report_kind": "rabies_cross_host_geography_package",
        "dataset_id": dataset.dataset_id,
        "label": dataset.label,
        "biological_question": _FLAGSHIP_QUESTION,
        "short_answer": short_answer,
        "package_files": {
            "overview_markdown": {
                "path": "overview.md",
                "checksum": _checksum(path.parent / "overview.md"),
            },
            "overview_html": {
                "path": "rabies-cross-host-geography-overview.html",
                "checksum": _checksum(
                    path.parent / "rabies-cross-host-geography-overview.html"
                ),
            },
            "artifact_inventory": {
                "path": artifact_inventory_path.name,
                "checksum": _checksum(artifact_inventory_path),
                "artifact_count": len(artifact_inventory_rows),
                "section_counts": inventory_counts,
            },
            "reproducibility_checklist": {
                "path": reproducibility_checklist_path.name,
                "checksum": _checksum(reproducibility_checklist_path),
                "item_count": len(checklist_rows),
                "blocked_count": blocked_check_count,
                "risk_count": risk_check_count,
            },
        },
        "config": {
            "path": f"dataset/{dataset_export.workflow_config_path.name}",
            "checksum": _checksum(dataset_export.workflow_config_path),
            "workflow_prefix": config.workflow_prefix,
            "alignment_mode": config.alignment_mode,
            "trimming_mode": config.trimming_mode,
            "trim_gap_threshold": config.trim_gap_threshold,
            "bootstrap_consensus_threshold": config.bootstrap_consensus_threshold,
            "bootstrap_robust_support_threshold": (
                config.bootstrap_robust_support_threshold
            ),
            "comparative_formula": config.comparative_formula,
            "comparative_response": config.comparative_response,
            "comparative_branch_length_floor": (config.comparative_branch_length_floor),
            "timeout_seconds": config.timeout_seconds,
            "max_bootstrap_tree_count": config.max_bootstrap_tree_count,
            "max_report_table_rows": config.max_report_table_rows,
            "memory_warning_threshold_bytes": config.memory_warning_threshold_bytes,
        },
        "dataset_files": {
            "readme": {
                "path": f"dataset/{dataset_export.readme_path.name}",
                "checksum": _checksum(dataset_export.readme_path),
            },
            "sequences": {
                "path": f"dataset/{dataset_export.sequences_path.name}",
                "checksum": _checksum(dataset_export.sequences_path),
            },
            "metadata": {
                "path": f"dataset/{dataset_export.metadata_path.name}",
                "checksum": _checksum(dataset_export.metadata_path),
            },
            "centroids": {
                "path": f"dataset/{dataset_export.centroids_path.name}",
                "checksum": _checksum(dataset_export.centroids_path),
            },
            "source_accessions": {
                "path": f"dataset/{dataset_export.accession_table_path.name}",
                "checksum": _checksum(dataset_export.accession_table_path),
            },
        },
        "workflow_files": {
            "final_report": {
                "path": f"workflow/{workflow_bundle.final_report_path.name}",
                "checksum": _checksum(workflow_bundle.final_report_path),
            },
            "workflow_log": {
                "path": f"workflow/{workflow_bundle.log_path.name}",
                "checksum": _checksum(workflow_bundle.log_path),
            },
            "workflow_summary": {
                "path": f"workflow/{workflow_bundle.workflow_summary_path.name}",
                "checksum": _checksum(workflow_bundle.workflow_summary_path),
            },
            "resource_observations": {
                "path": f"workflow/{workflow_bundle.resource_observations_path.name}",
                "checksum": _checksum(workflow_bundle.resource_observations_path),
            },
            "final_manifest": {
                "path": f"workflow/{workflow_bundle.final_manifest_path.name}",
                "checksum": _checksum(workflow_bundle.final_manifest_path),
            },
            "rooted_tree": {
                "path": f"workflow/{workflow_bundle.tree_path.name}",
                "checksum": _checksum(workflow_bundle.tree_path),
            },
            "rooting_report": {
                "path": f"workflow/{workflow_bundle.rooting_report_path.name}",
                "checksum": _checksum(workflow_bundle.rooting_report_path),
            },
            "model_table": {
                "path": f"workflow/{workflow_bundle.model_table_path.name}",
                "checksum": _checksum(workflow_bundle.model_table_path),
            },
            "support_table": {
                "path": f"workflow/{workflow_bundle.support_table_path.name}",
                "checksum": _checksum(workflow_bundle.support_table_path),
            },
            "bootstrap_summary": {
                "path": (
                    "workflow/bootstrap-review/"
                    f"{workflow_bundle.bootstrap_summary_path.name}"
                ),
                "checksum": _checksum(workflow_bundle.bootstrap_summary_path),
            },
            "bootstrap_tree_comparison_summary": {
                "path": (
                    "workflow/bootstrap-review/"
                    f"{workflow_bundle.bootstrap_tree_comparison_summary_path.name}"
                ),
                "checksum": _checksum(
                    workflow_bundle.bootstrap_tree_comparison_summary_path
                ),
            },
            "host_switch_summary": {
                "path": f"workflow/{workflow_bundle.host_switch_summary_path.name}",
                "checksum": _checksum(workflow_bundle.host_switch_summary_path),
            },
            "biogeography_report": {
                "path": (
                    "workflow/biogeography/"
                    f"{workflow_bundle.biogeography_report_path.name}"
                ),
                "checksum": _checksum(workflow_bundle.biogeography_report_path),
            },
            "comparative_report": {
                "path": (
                    "workflow/comparative/"
                    f"{workflow_bundle.comparative_report_path.name}"
                ),
                "checksum": _checksum(workflow_bundle.comparative_report_path),
            },
            "conclusion_stability_report": {
                "path": (
                    "workflow/conclusion-stability/"
                    f"{workflow_bundle.conclusion_stability_report_path.name}"
                ),
                "checksum": _checksum(workflow_bundle.conclusion_stability_report_path),
            },
            "scientific_findings": {
                "path": f"workflow/{workflow_bundle.scientific_findings_path.name}",
                "checksum": _checksum(workflow_bundle.scientific_findings_path),
            },
        },
        "metrics": {
            "sequence_count": dataset.sequence_count,
            "selected_model": workflow_bundle.selected_model,
            "root_host": workflow_bundle.root_host,
            "root_region": workflow_bundle.root_region,
            "bootstrap_tree_count": workflow_bundle.bootstrap_tree_count,
            "workflow_runtime_seconds": workflow_bundle.workflow_runtime_seconds,
            "bootstrap_review_runtime_seconds": (
                workflow_bundle.bootstrap_review_runtime_seconds
            ),
            "bootstrap_review_peak_memory_bytes": (
                workflow_bundle.bootstrap_review_peak_memory_bytes
            ),
            "budget_warning_count": workflow_bundle.budget_warning_count,
            "host_switch_count": workflow_bundle.host_switch_count,
            "migration_event_count": workflow_bundle.migration_event_count,
            "comparative_selected_model": workflow_bundle.comparative_selected_model,
            "comparative_pgls_lambda": workflow_bundle.comparative_pgls_lambda,
            "comparative_pgls_r_squared": workflow_bundle.comparative_pgls_r_squared,
            "comparative_branch_repair_count": (
                workflow_bundle.comparative_branch_repair_count
            ),
            "conclusion_stable_count": workflow_bundle.conclusion_stable_count,
            "conclusion_weak_count": workflow_bundle.conclusion_weak_count,
            "conclusion_unstable_count": workflow_bundle.conclusion_unstable_count,
            "scientific_finding_count": workflow_bundle.scientific_finding_count,
        },
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path
