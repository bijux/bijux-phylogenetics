from __future__ import annotations

from datetime import UTC, datetime
import json
from pathlib import Path
import shutil
from typing import Any

from bijux_phylogenetics.runtime.errors import EngineWorkflowError

from ..common import load_engine_manifest
from ..presentation import render_inference_workflow_report
from .workflow_result_bundle import (
    BUNDLE_MANIFEST_NAME,
    WorkflowResultBundleExtraInput,
    WorkflowResultBundleFile,
    WorkflowResultBundleIssue,
    WorkflowResultBundleReport,
    WorkflowResultBundleValidationReport,
    WORKFLOW_CONFIG_NAME,
    WORKFLOW_REPORT_NAME,
    WORKFLOW_RERUN_NAME,
    build_bundle_rerun_payload,
    copy_bundle_file,
    input_label,
    maybe_path,
    output_filename,
    payload_workflow,
    prepared_input_label,
    recorded_input_paths,
    record_bundle_file,
    render_bundle_readme,
    required_output_paths,
    sha256_file,
    step_manifest_paths,
    write_bundle_report,
    write_bundle_json,
)

__all__ = [
    "WorkflowResultBundleExtraInput",
    "WorkflowResultBundleFile",
    "WorkflowResultBundleIssue",
    "WorkflowResultBundleReport",
    "WorkflowResultBundleValidationReport",
    "export_workflow_result_bundle",
    "validate_workflow_result_bundle",
]

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
    workflow = payload_workflow(payload)
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
    copy_bundle_file(Path(manifest_path), workflow_manifest_path)
    files.append(
        record_bundle_file(
            role="workflow_manifest",
            label="workflow_manifest",
            bundle_root=bundle_root,
            path=workflow_manifest_path,
            source_path=Path(manifest_path),
        )
    )

    run_manifest_path = maybe_path(payload.get("run_manifest_path"))
    if run_manifest_path is not None:
        if not run_manifest_path.exists():
            raise EngineWorkflowError(
                f"workflow bundle source is missing the run manifest: {run_manifest_path}",
                code="workflow_bundle_missing_run_manifest",
                details={"path": str(run_manifest_path)},
            )
        bundled_run_manifest = bundle_root / "manifests" / "workflow.run.json"
        copy_bundle_file(run_manifest_path, bundled_run_manifest)
        files.append(
            record_bundle_file(
                role="run_manifest",
                label="workflow_run_manifest",
                bundle_root=bundle_root,
                path=bundled_run_manifest,
                source_path=run_manifest_path,
            )
        )

    config_path = bundle_root / WORKFLOW_CONFIG_NAME
    bundle_config_payload = (
        dict(payload.get("config", {}))
        if config_payload is None
        else dict(config_payload)
    )
    write_bundle_json(config_path, bundle_config_payload)
    files.append(
        record_bundle_file(
            role="workflow_config",
            label="workflow_config",
            bundle_root=bundle_root,
            path=config_path,
            source_path=None,
        )
    )

    rerun_path = bundle_root / WORKFLOW_RERUN_NAME
    rerun_payload = build_bundle_rerun_payload(payload, bundle_root=bundle_root)
    write_bundle_json(rerun_path, rerun_payload)
    files.append(
        record_bundle_file(
            role="workflow_rerun",
            label="workflow_rerun",
            bundle_root=bundle_root,
            path=rerun_path,
            source_path=None,
        )
    )

    input_entries: list[dict[str, str]] = []
    for index, input_path in enumerate(recorded_input_paths(payload), start=1):
        label = input_label(index=index, path=input_path)
        if input_path.exists():
            destination = bundle_root / "inputs" / label
            copy_bundle_file(input_path, destination)
            files.append(
                record_bundle_file(
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

    prepared_input_path = maybe_path(payload.get("prepared_input_path"))
    if prepared_input_path is not None and prepared_input_path.exists():
        prepared_label = prepared_input_label(prepared_input_path)
        destination = bundle_root / "inputs" / prepared_label
        if not destination.exists():
            copy_bundle_file(prepared_input_path, destination)
            files.append(
                record_bundle_file(
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
        copy_bundle_file(extra_input.source_path, destination)
        files.append(
            record_bundle_file(
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
    for label, source_path in required_output_paths(payload).items():
        if (
            label == "manifest"
            and source_path.resolve() == Path(manifest_path).resolve()
        ):
            continue
        destination = (
            bundle_root
            / "outputs"
            / "final"
            / output_filename(label=label, source_path=source_path)
        )
        copy_bundle_file(source_path, destination)
        files.append(
            record_bundle_file(
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
    for step_label, source_manifest_path in step_manifest_paths(payload).items():
        destination = (
            bundle_root / "manifests" / "steps" / f"{step_label}.manifest.json"
        )
        copy_bundle_file(source_manifest_path, destination)
        files.append(
            record_bundle_file(
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
        for output_label, output_path in required_output_paths(step_payload).items():
            step_destination = (
                bundle_root
                / "outputs"
                / "engine-artifacts"
                / step_label
                / output_filename(label=output_label, source_path=output_path)
            )
            copy_bundle_file(output_path, step_destination)
            files.append(
                record_bundle_file(
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
            record_bundle_file(
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

    report_path = bundle_root / "reports" / WORKFLOW_REPORT_NAME
    write_bundle_report(
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
        record_bundle_file(
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
        render_bundle_readme(
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
        record_bundle_file(
            role="bundle_readme",
            label="bundle_readme",
            bundle_root=bundle_root,
            path=readme_path,
            source_path=None,
        )
    )

    bundle_manifest_path = bundle_root / BUNDLE_MANIFEST_NAME
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
    write_bundle_json(bundle_manifest_path, bundle_manifest)

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
    bundle_manifest_path = bundle_root / BUNDLE_MANIFEST_NAME
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
                    relative_path=Path(BUNDLE_MANIFEST_NAME),
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
        observed = sha256_file(target)
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
