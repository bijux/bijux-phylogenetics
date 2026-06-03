from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .bundle_files import sha256_file
from .contracts import (
    WorkflowResultBundleIssue,
    WorkflowResultBundleValidationReport,
)
from .layout import BUNDLE_MANIFEST_NAME


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
        require_payload_path(
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
        validate_fasta_to_tree_bundle(payload, bundle_root=bundle_root, issues=issues)

    return WorkflowResultBundleValidationReport(
        bundle_root=bundle_root,
        workflow=workflow,
        valid=not issues,
        file_count=len(file_entries),
        issues=issues,
    )


def validate_fasta_to_tree_bundle(
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


def require_payload_path(
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
