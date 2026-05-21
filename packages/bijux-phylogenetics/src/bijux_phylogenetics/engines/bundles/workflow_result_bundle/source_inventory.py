from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.runtime.errors import EngineWorkflowError

from .bundle_files import input_label


def payload_workflow(payload: dict[str, Any]) -> str:
    workflow = payload.get("workflow")
    if workflow is None:
        raise EngineWorkflowError(
            "workflow bundle requires a workflow identifier",
            code="workflow_bundle_missing_workflow",
        )
    return str(workflow)


def recorded_input_paths(payload: dict[str, Any]) -> list[Path]:
    if "input_paths" in payload:
        return [Path(path) for path in payload["input_paths"]]
    if "input_path" in payload:
        return [Path(payload["input_path"])]
    raise EngineWorkflowError(
        "workflow bundle requires recorded input paths",
        code="workflow_bundle_missing_inputs",
    )


def required_output_paths(payload: dict[str, Any]) -> dict[str, Path]:
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


def step_manifest_paths(payload: dict[str, Any]) -> dict[str, Path]:
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


def build_bundle_rerun_payload(
    payload: dict[str, Any], *, bundle_root: Path
) -> dict[str, object]:
    return {
        "workflow": payload_workflow(payload),
        "config": dict(payload.get("config", {})),
        "engine_versions": dict(payload.get("engine_versions", {})),
        "bundle_local_inputs": [
            {
                "source_path": str(path),
                "relative_path": (
                    bundle_root / "inputs" / input_label(index=index, path=path)
                )
                .relative_to(bundle_root)
                .as_posix(),
            }
            for index, path in enumerate(recorded_input_paths(payload), start=1)
            if path.exists()
        ],
        "input_checksums": dict(payload.get("input_checksums", {})),
        "notes": [
            "Use the bundled input files together with this config to rerun the workflow in a new output directory.",
            "The copied workflow manifest preserves original source paths for provenance; use this rerun ledger for bundle-local execution.",
        ],
    }
