from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


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
