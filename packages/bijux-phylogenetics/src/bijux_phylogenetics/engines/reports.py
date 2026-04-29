from __future__ import annotations

import json
from dataclasses import dataclass
from dataclasses import asdict
from pathlib import Path

from bijux_phylogenetics.render.html import write_html_report

from .common import load_engine_manifest
from .validation import validate_model_selection_against_engine_outputs


@dataclass(slots=True)
class InferenceWorkflowReportBuildResult:
    output_path: Path
    report_kind: str
    title: str
    manifest_path: Path
    engine_name: str
    workflow: str
    warning_count: int
    machine_manifest: dict[str, object]


@dataclass(slots=True)
class ModelSelectionLimitationsReport:
    manifest_path: Path
    selected_model: str | None
    validation_issues: list[str]
    limitations: list[str]
    interpretation_limits: list[str]


@dataclass(slots=True)
class ModelSelectionLimitationsReportBuildResult:
    output_path: Path
    report_kind: str
    title: str
    manifest_path: Path
    selected_model: str | None
    machine_manifest: dict[str, object]


def render_inference_workflow_report(*, manifest_path: Path, out_path: Path) -> InferenceWorkflowReportBuildResult:
    """Render a deterministic HTML report for one engine workflow manifest."""
    manifest = load_engine_manifest(manifest_path)
    title = f"Bijux Inference Workflow Report: {manifest['workflow']}"
    sections = [
        ("workflow-summary", json.dumps(
            {
                "engine_name": manifest["engine_name"],
                "workflow": manifest["workflow"],
                "resumed": manifest.get("resumed", False),
                "selected_model": manifest.get("selected_model"),
                "notes": manifest.get("notes", []),
            },
            indent=2,
            sort_keys=True,
        )),
        ("inputs", json.dumps(
            {
                "input_paths": manifest["input_paths"],
                "input_checksums": manifest.get("input_checksums", {}),
            },
            indent=2,
            sort_keys=True,
        )),
        ("outputs", json.dumps(
            {
                "output_paths": manifest["output_paths"],
                "output_checksums": manifest.get("output_checksums", {}),
            },
            indent=2,
            sort_keys=True,
        )),
        ("engine-run", json.dumps(manifest["run"], indent=2, sort_keys=True)),
    ]
    machine_manifest = {
        "report_kind": "inference-workflow",
        "title": title,
        "manifest_path": str(manifest_path),
        "engine_name": manifest["engine_name"],
        "workflow": manifest["workflow"],
        "warning_count": len(manifest["run"].get("warning_lines", [])),
        "sections": [name for name, _ in sections],
    }
    write_html_report(title=title, sections=sections, out_path=out_path, embedded_json=machine_manifest)
    return InferenceWorkflowReportBuildResult(
        output_path=out_path,
        report_kind="inference-workflow",
        title=title,
        manifest_path=manifest_path,
        engine_name=str(manifest["engine_name"]),
        workflow=str(manifest["workflow"]),
        warning_count=len(manifest["run"].get("warning_lines", [])),
        machine_manifest=machine_manifest,
    )


def build_model_selection_limitations_report(manifest_path: Path) -> ModelSelectionLimitationsReport:
    """Describe why one selected substitution model should not be treated as biological truth."""
    manifest = load_engine_manifest(manifest_path)
    validation = validate_model_selection_against_engine_outputs(manifest_path)
    selected_model = manifest.get("selected_model")
    limitations = [
        "the selected model is the best fit among the engine's tested candidates, not proof of the true evolutionary process",
        "closely scoring models can yield similar likelihood support while implying different branch lengths or uncertainty",
        "model fit reflects the provided alignment, trimming policy, and missing-data pattern as much as the underlying biology",
        "an apparently optimal model can still be inappropriate if the alignment is compositionally biased, saturated, or poorly aligned",
    ]
    interpretation_limits = [
        "do not describe the selected model as the biologically correct mutation process without additional validation",
        "treat downstream topology or branch-length differences across plausible models as sensitivity, not noise",
        "report whether alternative near-best models were considered when scientific conclusions depend on branch support or branch lengths",
    ]
    return ModelSelectionLimitationsReport(
        manifest_path=manifest_path,
        selected_model=None if selected_model is None else str(selected_model),
        validation_issues=validation.issues,
        limitations=limitations,
        interpretation_limits=interpretation_limits,
    )


def render_model_selection_limitations_report(
    *,
    manifest_path: Path,
    out_path: Path,
) -> ModelSelectionLimitationsReportBuildResult:
    """Render a reviewer-facing model-selection limitations report."""
    report = build_model_selection_limitations_report(manifest_path)
    title = "Bijux Model Selection Limitations Report"
    sections = [
        ("model-selection-summary", json.dumps(asdict(report), default=str, indent=2, sort_keys=True)),
    ]
    machine_manifest = {
        "report_kind": "model-selection-limitations",
        "title": title,
        "manifest_path": str(manifest_path),
        "selected_model": report.selected_model,
        "sections": [name for name, _ in sections],
    }
    write_html_report(title=title, sections=sections, out_path=out_path, embedded_json=machine_manifest)
    return ModelSelectionLimitationsReportBuildResult(
        output_path=out_path,
        report_kind="model-selection-limitations",
        title=title,
        manifest_path=manifest_path,
        selected_model=report.selected_model,
        machine_manifest=machine_manifest,
    )
