from __future__ import annotations

from dataclasses import dataclass
from datetime import UTC, datetime
import hashlib
from html import escape
import json
from pathlib import Path
import shutil
from typing import Any

from bijux_phylogenetics.runtime.errors import EngineWorkflowError

from .common import load_engine_manifest
from .reports import render_inference_workflow_report

__all__ = [
    "WorkflowResultBundleExtraInput",
    "WorkflowResultBundleFile",
    "WorkflowResultBundleIssue",
    "WorkflowResultBundleReport",
    "WorkflowResultBundleValidationReport",
    "export_workflow_result_bundle",
    "validate_workflow_result_bundle",
]

_BUNDLE_MANIFEST_NAME = "bundle.manifest.json"
_WORKFLOW_REPORT_NAME = "workflow-report.html"
_WORKFLOW_CONFIG_NAME = "workflow-config.json"
_WORKFLOW_RERUN_NAME = "workflow-rerun.json"


@dataclass(frozen=True, slots=True)
class WorkflowResultBundleExtraInput:
    label: str
    source_path: Path


@dataclass(slots=True)
class WorkflowResultBundleFile:
    role: str
    label: str
    relative_path: Path
    sha256: str
    size_bytes: int
    source_path: str | None


@dataclass(slots=True)
class WorkflowResultBundleReport:
    bundle_root: Path
    workflow: str
    source_manifest_path: Path
    bundle_manifest_path: Path
    workflow_manifest_path: Path
    config_path: Path
    rerun_path: Path
    report_path: Path
    copied_input_count: int
    copied_output_count: int
    copied_step_manifest_count: int
    copied_step_output_count: int
    copied_report_count: int
    input_hash_count: int
    includes_input_files: bool
    missing_input_paths: list[Path]
    file_count: int
    files: list[WorkflowResultBundleFile]
    notes: list[str]


@dataclass(slots=True)
class WorkflowResultBundleIssue:
    kind: str
    label: str
    detail: str
    relative_path: Path | None = None


@dataclass(slots=True)
class WorkflowResultBundleValidationReport:
    bundle_root: Path
    workflow: str | None
    valid: bool
    file_count: int
    issues: list[WorkflowResultBundleIssue]


def export_workflow_result_bundle(
    manifest_path: Path,
    *,
    bundle_root: Path,
    config_payload: dict[str, object] | None = None,
    extra_inputs: list[WorkflowResultBundleExtraInput] | None = None,
    extra_notes: list[str] | None = None,
) -> WorkflowResultBundleReport:
    """Export one portable workflow-result bundle rooted on one workflow manifest."""
    payload = load_engine_manifest(manifest_path)
    workflow = _payload_workflow(payload)
    if bundle_root.exists():
        shutil.rmtree(bundle_root)
    (bundle_root / "inputs").mkdir(parents=True, exist_ok=True)
    (bundle_root / "outputs" / "final").mkdir(parents=True, exist_ok=True)
    (bundle_root / "outputs" / "engine-artifacts").mkdir(parents=True, exist_ok=True)
    (bundle_root / "manifests" / "steps").mkdir(parents=True, exist_ok=True)
    (bundle_root / "reports" / "steps").mkdir(parents=True, exist_ok=True)

    files: list[WorkflowResultBundleFile] = []
    missing_input_paths: list[Path] = []
    notes: list[str] = []
    copied_input_count = 0
    copied_output_count = 0
    copied_step_manifest_count = 0
    copied_step_output_count = 0
    copied_report_count = 0

    workflow_manifest_path = bundle_root / "manifests" / "workflow.manifest.json"
    _copy_file(Path(manifest_path), workflow_manifest_path)
    files.append(
        _record_bundle_file(
            role="workflow_manifest",
            label="workflow_manifest",
            bundle_root=bundle_root,
            path=workflow_manifest_path,
            source_path=Path(manifest_path),
        )
    )

    run_manifest_path = _maybe_path(payload.get("run_manifest_path"))
    if run_manifest_path is not None:
        if not run_manifest_path.exists():
            raise EngineWorkflowError(
                f"workflow bundle source is missing the run manifest: {run_manifest_path}",
                code="workflow_bundle_missing_run_manifest",
                details={"path": str(run_manifest_path)},
            )
        bundled_run_manifest = bundle_root / "manifests" / "workflow.run.json"
        _copy_file(run_manifest_path, bundled_run_manifest)
        files.append(
            _record_bundle_file(
                role="run_manifest",
                label="workflow_run_manifest",
                bundle_root=bundle_root,
                path=bundled_run_manifest,
                source_path=run_manifest_path,
            )
        )

    config_path = bundle_root / _WORKFLOW_CONFIG_NAME
    bundle_config_payload = (
        dict(payload.get("config", {}))
        if config_payload is None
        else dict(config_payload)
    )
    _write_json(config_path, bundle_config_payload)
    files.append(
        _record_bundle_file(
            role="workflow_config",
            label="workflow_config",
            bundle_root=bundle_root,
            path=config_path,
            source_path=None,
        )
    )

    rerun_path = bundle_root / _WORKFLOW_RERUN_NAME
    rerun_payload = _build_bundle_rerun_payload(payload, bundle_root=bundle_root)
    _write_json(rerun_path, rerun_payload)
    files.append(
        _record_bundle_file(
            role="workflow_rerun",
            label="workflow_rerun",
            bundle_root=bundle_root,
            path=rerun_path,
            source_path=None,
        )
    )

    input_entries: list[dict[str, str]] = []
    for index, input_path in enumerate(_recorded_input_paths(payload), start=1):
        label = _input_label(index=index, path=input_path)
        if input_path.exists():
            destination = bundle_root / "inputs" / label
            _copy_file(input_path, destination)
            files.append(
                _record_bundle_file(
                    role="input_file",
                    label=label,
                    bundle_root=bundle_root,
                    path=destination,
                    source_path=input_path,
                )
            )
            copied_input_count += 1
            input_entries.append(
                {
                    "label": label,
                    "relative_path": destination.relative_to(bundle_root).as_posix(),
                    "source_path": str(input_path),
                }
            )
            continue
        missing_input_paths.append(input_path)
        input_entries.append(
            {
                "label": label,
                "source_path": str(input_path),
            }
        )

    prepared_input_path = _maybe_path(payload.get("prepared_input_path"))
    if prepared_input_path is not None and prepared_input_path.exists():
        prepared_label = _prepared_input_label(prepared_input_path)
        destination = bundle_root / "inputs" / prepared_label
        if not destination.exists():
            _copy_file(prepared_input_path, destination)
            files.append(
                _record_bundle_file(
                    role="input_file",
                    label=prepared_label,
                    bundle_root=bundle_root,
                    path=destination,
                    source_path=prepared_input_path,
                )
            )
            copied_input_count += 1
        input_entries.append(
            {
                "label": prepared_label,
                "relative_path": destination.relative_to(bundle_root).as_posix(),
                "source_path": str(prepared_input_path),
            }
        )

    for extra_input in [] if extra_inputs is None else extra_inputs:
        if not extra_input.source_path.exists():
            missing_input_paths.append(extra_input.source_path)
            input_entries.append(
                {
                    "label": extra_input.label,
                    "source_path": str(extra_input.source_path),
                }
            )
            continue
        destination = bundle_root / "inputs" / extra_input.label
        if destination.exists():
            raise EngineWorkflowError(
                "workflow bundle extra input label would overwrite an existing bundled file",
                code="workflow_bundle_duplicate_extra_input",
                details={
                    "label": extra_input.label,
                    "destination": str(destination),
                },
            )
        _copy_file(extra_input.source_path, destination)
        files.append(
            _record_bundle_file(
                role="input_file",
                label=extra_input.label,
                bundle_root=bundle_root,
                path=destination,
                source_path=extra_input.source_path,
            )
        )
        copied_input_count += 1
        input_entries.append(
            {
                "label": extra_input.label,
                "relative_path": destination.relative_to(bundle_root).as_posix(),
                "source_path": str(extra_input.source_path),
            }
        )

    output_entries: dict[str, str] = {}
    for label, source_path in _required_output_paths(payload).items():
        if (
            label == "manifest"
            and source_path.resolve() == Path(manifest_path).resolve()
        ):
            continue
        destination = (
            bundle_root
            / "outputs"
            / "final"
            / _output_filename(label=label, source_path=source_path)
        )
        _copy_file(source_path, destination)
        files.append(
            _record_bundle_file(
                role="workflow_output",
                label=label,
                bundle_root=bundle_root,
                path=destination,
                source_path=source_path,
            )
        )
        copied_output_count += 1
        output_entries[label] = destination.relative_to(bundle_root).as_posix()

    step_manifest_entries: dict[str, str] = {}
    step_output_entries: dict[str, dict[str, str]] = {}
    for step_label, source_manifest_path in _step_manifest_paths(payload).items():
        destination = (
            bundle_root / "manifests" / "steps" / f"{step_label}.manifest.json"
        )
        _copy_file(source_manifest_path, destination)
        files.append(
            _record_bundle_file(
                role="step_manifest",
                label=step_label,
                bundle_root=bundle_root,
                path=destination,
                source_path=source_manifest_path,
            )
        )
        copied_step_manifest_count += 1
        step_manifest_entries[step_label] = destination.relative_to(
            bundle_root
        ).as_posix()

        step_payload = load_engine_manifest(source_manifest_path)
        step_output_entries[step_label] = {}
        for output_label, output_path in _required_output_paths(step_payload).items():
            step_destination = (
                bundle_root
                / "outputs"
                / "engine-artifacts"
                / step_label
                / _output_filename(label=output_label, source_path=output_path)
            )
            _copy_file(output_path, step_destination)
            files.append(
                _record_bundle_file(
                    role="step_output",
                    label=f"{step_label}:{output_label}",
                    bundle_root=bundle_root,
                    path=step_destination,
                    source_path=output_path,
                )
            )
            copied_step_output_count += 1
            step_output_entries[step_label][output_label] = (
                step_destination.relative_to(bundle_root).as_posix()
            )

        step_report_path = (
            bundle_root / "reports" / "steps" / f"{step_label}-report.html"
        )
        render_inference_workflow_report(
            manifest_path=source_manifest_path,
            out_path=step_report_path,
        )
        files.append(
            _record_bundle_file(
                role="step_report",
                label=step_label,
                bundle_root=bundle_root,
                path=step_report_path,
                source_path=None,
            )
        )
        copied_report_count += 1

    if missing_input_paths:
        notes.append(
            "bundle exported with missing source input files; original manifest input checksums remain recorded for reproducibility"
        )
    if extra_notes:
        notes.extend(extra_notes)

    report_path = bundle_root / "reports" / _WORKFLOW_REPORT_NAME
    _write_bundle_report(
        report_path,
        workflow=workflow,
        source_manifest_path=Path(manifest_path),
        input_entries=input_entries,
        output_entries=output_entries,
        step_manifest_entries=step_manifest_entries,
        notes=notes,
        stage_fingerprints=payload.get("stage_fingerprints"),
        embedded_manifest={
            "workflow": workflow,
            "config": bundle_config_payload,
            "input_checksums": payload.get("input_checksums", {}),
            "output_checksums": payload.get("output_checksums", {}),
        },
    )
    files.append(
        _record_bundle_file(
            role="workflow_report",
            label="workflow_report",
            bundle_root=bundle_root,
            path=report_path,
            source_path=None,
        )
    )
    copied_report_count += 1

    readme_path = bundle_root / "README.md"
    readme_path.write_text(
        _render_bundle_readme(
            workflow=workflow,
            source_manifest_path=Path(manifest_path),
            includes_input_files=bool(input_entries),
            missing_input_paths=missing_input_paths,
            output_labels=sorted(output_entries),
            step_labels=sorted(step_manifest_entries),
        ),
        encoding="utf-8",
    )
    files.append(
        _record_bundle_file(
            role="bundle_readme",
            label="bundle_readme",
            bundle_root=bundle_root,
            path=readme_path,
            source_path=None,
        )
    )

    bundle_manifest_path = bundle_root / _BUNDLE_MANIFEST_NAME
    bundle_manifest = {
        "workflow": workflow,
        "source_manifest_path": str(Path(manifest_path)),
        "created_at_utc": datetime.now(UTC)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z"),
        "bundle_root": str(bundle_root),
        "workflow_manifest": workflow_manifest_path.relative_to(bundle_root).as_posix(),
        "workflow_config": config_path.relative_to(bundle_root).as_posix(),
        "workflow_rerun": rerun_path.relative_to(bundle_root).as_posix(),
        "workflow_report": report_path.relative_to(bundle_root).as_posix(),
        "workflow_input_checksums": dict(payload.get("input_checksums", {})),
        "workflow_output_checksums": dict(payload.get("output_checksums", {})),
        "input_files": input_entries,
        "missing_input_paths": [str(path) for path in missing_input_paths],
        "workflow_outputs": output_entries,
        "step_manifests": step_manifest_entries,
        "step_outputs": step_output_entries,
        "required_output_labels": sorted(output_entries),
        "required_step_manifest_labels": sorted(step_manifest_entries),
        "files": [
            {
                "role": file.role,
                "label": file.label,
                "relative_path": file.relative_path.as_posix(),
                "sha256": file.sha256,
                "size_bytes": file.size_bytes,
                "source_path": file.source_path,
            }
            for file in files
        ],
        "notes": notes,
    }
    _write_json(bundle_manifest_path, bundle_manifest)

    return WorkflowResultBundleReport(
        bundle_root=bundle_root,
        workflow=workflow,
        source_manifest_path=Path(manifest_path),
        bundle_manifest_path=bundle_manifest_path,
        workflow_manifest_path=workflow_manifest_path,
        config_path=config_path,
        rerun_path=rerun_path,
        report_path=report_path,
        copied_input_count=copied_input_count,
        copied_output_count=copied_output_count,
        copied_step_manifest_count=copied_step_manifest_count,
        copied_step_output_count=copied_step_output_count,
        copied_report_count=copied_report_count,
        input_hash_count=len(dict(payload.get("input_checksums", {}))),
        includes_input_files=copied_input_count > 0,
        missing_input_paths=missing_input_paths,
        file_count=len(files),
        files=files,
        notes=notes,
    )


def validate_workflow_result_bundle(
    bundle_root: Path,
) -> WorkflowResultBundleValidationReport:
    """Validate one workflow-result bundle for checksum integrity and completeness."""
    issues: list[WorkflowResultBundleIssue] = []
    bundle_manifest_path = bundle_root / _BUNDLE_MANIFEST_NAME
    if not bundle_manifest_path.exists():
        return WorkflowResultBundleValidationReport(
            bundle_root=bundle_root,
            workflow=None,
            valid=False,
            file_count=0,
            issues=[
                WorkflowResultBundleIssue(
                    kind="missing-bundle-manifest",
                    label="bundle_manifest",
                    detail=f"bundle manifest not found: {bundle_manifest_path}",
                    relative_path=Path(_BUNDLE_MANIFEST_NAME),
                )
            ],
        )
    payload = json.loads(bundle_manifest_path.read_text(encoding="utf-8"))
    workflow = (
        str(payload.get("workflow")) if payload.get("workflow") is not None else None
    )
    file_entries = list(payload.get("files", []))
    for entry in file_entries:
        relative_path = Path(str(entry["relative_path"]))
        target = bundle_root / relative_path
        if not target.exists():
            issues.append(
                WorkflowResultBundleIssue(
                    kind="missing-file",
                    label=str(entry["label"]),
                    detail=f"bundle file is missing: {relative_path.as_posix()}",
                    relative_path=relative_path,
                )
            )
            continue
        observed = _sha256(target)
        if observed != str(entry["sha256"]):
            issues.append(
                WorkflowResultBundleIssue(
                    kind="checksum-mismatch",
                    label=str(entry["label"]),
                    detail=(
                        f"expected {entry['sha256']} but observed {observed} for "
                        f"{relative_path.as_posix()}"
                    ),
                    relative_path=relative_path,
                )
            )

    for label_name, relative_key in (
        ("workflow_manifest", "workflow_manifest"),
        ("workflow_config", "workflow_config"),
        ("workflow_rerun", "workflow_rerun"),
        ("workflow_report", "workflow_report"),
    ):
        _require_payload_path(
            issues,
            bundle_root=bundle_root,
            payload=payload,
            payload_key=relative_key,
            label=label_name,
        )

    if not list(payload.get("input_files", [])) and not dict(
        payload.get("workflow_input_checksums", {})
    ):
        issues.append(
            WorkflowResultBundleIssue(
                kind="missing-input-evidence",
                label="inputs",
                detail="bundle is missing both copied input files and recorded input checksums",
            )
        )

    workflow_output_labels = dict(payload.get("workflow_outputs", {}))
    if not workflow_output_labels:
        issues.append(
            WorkflowResultBundleIssue(
                kind="missing-workflow-outputs",
                label="workflow_outputs",
                detail="bundle does not record any copied workflow outputs",
            )
        )
    for label, relative_path_text in workflow_output_labels.items():
        target = bundle_root / Path(relative_path_text)
        if not target.exists():
            issues.append(
                WorkflowResultBundleIssue(
                    kind="missing-required-output",
                    label=label,
                    detail=f"expected workflow output is missing: {relative_path_text}",
                    relative_path=Path(relative_path_text),
                )
            )

    step_manifest_labels = dict(payload.get("step_manifests", {}))
    for label, relative_path_text in step_manifest_labels.items():
        target = bundle_root / Path(relative_path_text)
        if not target.exists():
            issues.append(
                WorkflowResultBundleIssue(
                    kind="missing-step-manifest",
                    label=label,
                    detail=f"expected step manifest is missing: {relative_path_text}",
                    relative_path=Path(relative_path_text),
                )
            )

    step_output_labels = {
        str(step_label): dict(step_outputs)
        for step_label, step_outputs in dict(payload.get("step_outputs", {})).items()
    }
    for step_label, outputs in step_output_labels.items():
        for output_label, relative_path_text in outputs.items():
            target = bundle_root / Path(relative_path_text)
            if not target.exists():
                issues.append(
                    WorkflowResultBundleIssue(
                        kind="missing-step-output",
                        label=f"{step_label}:{output_label}",
                        detail=f"expected step output is missing: {relative_path_text}",
                        relative_path=Path(relative_path_text),
                    )
                )

    if workflow == "fasta-to-tree":
        _validate_fasta_to_tree_bundle(payload, bundle_root=bundle_root, issues=issues)

    return WorkflowResultBundleValidationReport(
        bundle_root=bundle_root,
        workflow=workflow,
        valid=not issues,
        file_count=len(file_entries),
        issues=issues,
    )


def _payload_workflow(payload: dict[str, Any]) -> str:
    workflow = payload.get("workflow")
    if workflow is None:
        raise EngineWorkflowError(
            "workflow bundle requires a workflow identifier",
            code="workflow_bundle_missing_workflow",
        )
    return str(workflow)


def _recorded_input_paths(payload: dict[str, Any]) -> list[Path]:
    if "input_paths" in payload:
        return [Path(path) for path in payload["input_paths"]]
    if "input_path" in payload:
        return [Path(payload["input_path"])]
    raise EngineWorkflowError(
        "workflow bundle requires recorded input paths",
        code="workflow_bundle_missing_inputs",
    )


def _required_output_paths(payload: dict[str, Any]) -> dict[str, Path]:
    output_paths = {
        str(label): Path(path)
        for label, path in dict(payload.get("output_paths", {})).items()
    }
    if not output_paths:
        raise EngineWorkflowError(
            "workflow bundle requires recorded output paths",
            code="workflow_bundle_missing_outputs",
        )
    missing = {label: path for label, path in output_paths.items() if not path.exists()}
    if missing:
        missing_payload = {label: str(path) for label, path in missing.items()}
        raise EngineWorkflowError(
            "workflow bundle source is missing one or more declared outputs",
            code="workflow_bundle_missing_output",
            details={"missing_outputs": missing_payload},
        )
    return output_paths


def _step_manifest_paths(payload: dict[str, Any]) -> dict[str, Path]:
    step_manifests = {
        str(label): Path(path)
        for label, path in dict(payload.get("step_manifests", {})).items()
    }
    missing = {
        label: path for label, path in step_manifests.items() if not path.exists()
    }
    if missing:
        missing_payload = {label: str(path) for label, path in missing.items()}
        raise EngineWorkflowError(
            "workflow bundle source is missing one or more step manifests",
            code="workflow_bundle_missing_step_manifest",
            details={"missing_step_manifests": missing_payload},
        )
    return step_manifests


def _build_bundle_rerun_payload(
    payload: dict[str, Any], *, bundle_root: Path
) -> dict[str, object]:
    return {
        "workflow": _payload_workflow(payload),
        "config": dict(payload.get("config", {})),
        "engine_versions": dict(payload.get("engine_versions", {})),
        "bundle_local_inputs": [
            {
                "source_path": str(path),
                "relative_path": (
                    bundle_root / "inputs" / _input_label(index=index, path=path)
                )
                .relative_to(bundle_root)
                .as_posix(),
            }
            for index, path in enumerate(_recorded_input_paths(payload), start=1)
            if path.exists()
        ],
        "input_checksums": dict(payload.get("input_checksums", {})),
        "notes": [
            "Use the bundled input files together with this config to rerun the workflow in a new output directory.",
            "The copied workflow manifest preserves original source paths for provenance; use this rerun ledger for bundle-local execution.",
        ],
    }


def _write_bundle_report(
    path: Path,
    *,
    workflow: str,
    source_manifest_path: Path,
    input_entries: list[dict[str, str]],
    output_entries: dict[str, str],
    step_manifest_entries: dict[str, str],
    notes: list[str],
    stage_fingerprints: object,
    embedded_manifest: dict[str, object],
) -> Path:
    input_lines = [
        f"{entry['label']}: {entry.get('relative_path', '[checksum only]')}"
        for entry in input_entries
    ]
    sections = [
        ("workflow", f"workflow: {workflow}\nsource_manifest: {source_manifest_path}"),
        ("inputs", "\n".join(input_lines) if input_lines else "none copied"),
        (
            "outputs",
            "\n".join(
                f"{label}: {path_text}" for label, path_text in output_entries.items()
            ),
        ),
        (
            "step-manifests",
            "\n".join(
                f"{label}: {path_text}"
                for label, path_text in step_manifest_entries.items()
            )
            if step_manifest_entries
            else "none",
        ),
    ]
    if stage_fingerprints is not None:
        sections.append(
            (
                "stage-fingerprints",
                json.dumps(stage_fingerprints, indent=2, sort_keys=True),
            )
        )
    if notes:
        sections.append(("notes", "\n".join(f"- {note}" for note in notes)))
    body = "\n".join(
        f"<section><h2>{escape(title)}</h2><pre>{escape(content)}</pre></section>"
        for title, content in sections
    )
    embedded_json = escape(
        json.dumps(embedded_manifest, indent=2, sort_keys=True).replace("</", "<\\/")
    )
    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{escape(f"Bijux Workflow Result Bundle: {workflow}")}</title>
  <style>
    :root {{
      color-scheme: light;
      --ink: #1b1f24;
      --bg: #f8fafc;
      --panel: #ffffff;
      --accent: #0f766e;
      --rule: #d6dee8;
      --mono: "SFMono-Regular", "SF Mono", Consolas, monospace;
    }}
    body {{
      margin: 0;
      padding: 2rem;
      background: linear-gradient(180deg, #eef6f4 0%, var(--bg) 100%);
      color: var(--ink);
      font: 16px/1.5 "Iowan Old Style", "Palatino Linotype", serif;
    }}
    main {{
      max-width: 960px;
      margin: 0 auto;
      background: var(--panel);
      border: 1px solid var(--rule);
      border-radius: 18px;
      padding: 2rem;
      box-shadow: 0 24px 80px rgba(15, 118, 110, 0.08);
    }}
    h1, h2 {{
      font-family: "Avenir Next", "Segoe UI", sans-serif;
      letter-spacing: 0.02em;
    }}
    h1 {{
      margin-top: 0;
      color: var(--accent);
    }}
    section + section {{
      margin-top: 1.5rem;
      padding-top: 1.5rem;
      border-top: 1px solid var(--rule);
    }}
    pre {{
      margin: 0;
      white-space: pre-wrap;
      word-break: break-word;
      font-family: var(--mono);
      background: #f6f8fa;
      border-radius: 12px;
      padding: 1rem;
    }}
  </style>
</head>
<body>
  <main>
    <h1>{escape(f"Bijux Workflow Result Bundle: {workflow}")}</h1>
    <script id="bijux-report-manifest" type="application/json">{embedded_json}</script>
    {body}
  </main>
</body>
</html>
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(html, encoding="utf-8")
    return path


def _render_bundle_readme(
    *,
    workflow: str,
    source_manifest_path: Path,
    includes_input_files: bool,
    missing_input_paths: list[Path],
    output_labels: list[str],
    step_labels: list[str],
) -> str:
    lines = [
        "# Bijux Workflow Result Bundle",
        "",
        f"Workflow: `{workflow}`",
        f"Source manifest: `{source_manifest_path}`",
        "",
        "## Contents",
        "",
        "- `bundle.manifest.json`: bundle inventory, checksums, and required workflow entries.",
        "- `workflow-config.json`: extracted workflow config for reruns and review.",
        "- `workflow-rerun.json`: bundle-local rerun ledger using the copied input files.",
        "- `reports/workflow-report.html`: reviewer-facing summary for the bundled workflow.",
        "- `manifests/`: copied workflow manifest plus any step manifests.",
        "- `inputs/`: copied workflow input files when they were still available at export time.",
        "- `outputs/final/`: copied reviewer-facing workflow outputs.",
        "- `outputs/engine-artifacts/`: copied native engine artifacts declared by the step manifests.",
        "",
        "## Review Notes",
        "",
        f"- copied input files present: `{includes_input_files}`",
        f"- missing source inputs at export time: `{len(missing_input_paths)}`",
        f"- final output labels: `{', '.join(output_labels)}`",
        f"- step manifests: `{', '.join(step_labels) if step_labels else 'none'}`",
        "",
    ]
    if missing_input_paths:
        lines.append("Missing input paths:")
        lines.extend(f"- `{path}`" for path in missing_input_paths)
        lines.append("")
    return "\n".join(lines)


def _validate_fasta_to_tree_bundle(
    payload: dict[str, Any],
    *,
    bundle_root: Path,
    issues: list[WorkflowResultBundleIssue],
) -> None:
    required_output_labels = {
        "alignment",
        "log",
        "model_table",
        "run_manifest",
        "support_table",
        "tree",
        "trimmed_alignment",
    }
    actual_output_labels = set(dict(payload.get("workflow_outputs", {})))
    missing_output_labels = sorted(required_output_labels - actual_output_labels)
    for label in missing_output_labels:
        issues.append(
            WorkflowResultBundleIssue(
                kind="missing-required-output-label",
                label=label,
                detail="fasta-to-tree bundle is missing one required final output entry",
            )
        )

    required_step_labels = {
        "alignment",
        "bootstrap_support",
        "maximum_likelihood",
        "model_selection",
        "trimming",
    }
    actual_step_labels = set(dict(payload.get("step_manifests", {})))
    missing_step_labels = sorted(required_step_labels - actual_step_labels)
    for label in missing_step_labels:
        issues.append(
            WorkflowResultBundleIssue(
                kind="missing-required-step-manifest-label",
                label=label,
                detail="fasta-to-tree bundle is missing one required step manifest entry",
            )
        )

    report_path = bundle_root / Path(str(payload["workflow_report"]))
    if not report_path.exists():
        issues.append(
            WorkflowResultBundleIssue(
                kind="missing-workflow-report",
                label="workflow_report",
                detail="fasta-to-tree bundle is missing its reviewer-facing report",
                relative_path=Path(str(payload["workflow_report"])),
            )
        )


def _require_payload_path(
    issues: list[WorkflowResultBundleIssue],
    *,
    bundle_root: Path,
    payload: dict[str, Any],
    payload_key: str,
    label: str,
) -> None:
    value = payload.get(payload_key)
    if value is None:
        issues.append(
            WorkflowResultBundleIssue(
                kind="missing-required-entry",
                label=label,
                detail=f"bundle manifest is missing required key: {payload_key}",
            )
        )
        return
    path = bundle_root / Path(str(value))
    if not path.exists():
        issues.append(
            WorkflowResultBundleIssue(
                kind="missing-required-file",
                label=label,
                detail=f"required bundle file is missing: {value}",
                relative_path=Path(str(value)),
            )
        )


def _copy_file(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def _record_bundle_file(
    *,
    role: str,
    label: str,
    bundle_root: Path,
    path: Path,
    source_path: Path | None,
) -> WorkflowResultBundleFile:
    relative_path = path.relative_to(bundle_root)
    return WorkflowResultBundleFile(
        role=role,
        label=label,
        relative_path=relative_path,
        sha256=_sha256(path),
        size_bytes=path.stat().st_size,
        source_path=None if source_path is None else str(source_path),
    )


def _maybe_path(value: object) -> Path | None:
    if value is None:
        return None
    return Path(str(value))


def _input_label(*, index: int, path: Path) -> str:
    return f"{index:02d}-{path.name}"


def _prepared_input_label(path: Path) -> str:
    return f"prepared-{path.name}"


def _output_filename(*, label: str, source_path: Path) -> str:
    return (
        f"{label}{source_path.suffix}"
        if source_path.suffix
        else f"{label}-{source_path.name}"
    )


def _write_json(path: Path, payload: dict[str, object]) -> Path:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()
