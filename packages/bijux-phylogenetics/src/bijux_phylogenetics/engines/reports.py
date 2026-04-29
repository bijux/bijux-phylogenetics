from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.render.html import write_html_report

from .common import load_engine_manifest


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
