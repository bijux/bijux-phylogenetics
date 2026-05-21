from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.runtime.errors import EngineWorkflowError

from ...common import load_engine_manifest
from .contracts import ManifestReplayReport
from .manifest_policy import (
    collect_engine_version_drift,
    collect_input_drift,
    default_replay_out_dir,
    payload_workflow,
)
from .output_comparison import compare_outputs
from .workflow_execution import (
    replay_composite_workflow,
    replay_engine_workflow,
)


def replay_workflow_manifest(
    manifest_path: Path,
    *,
    out_dir: Path | None = None,
    executables: dict[str, str | Path | None] | None = None,
) -> ManifestReplayReport:
    payload = load_engine_manifest(manifest_path)
    workflow = payload_workflow(payload)
    replay_out_dir = (
        default_replay_out_dir(manifest_path) if out_dir is None else out_dir
    )
    replay_out_dir.mkdir(parents=True, exist_ok=True)
    executable_overrides = {} if executables is None else dict(executables)

    input_drift = collect_input_drift(payload)
    input_drift_detected = any(not drift.matched for drift in input_drift)
    if input_drift_detected:
        raise EngineWorkflowError(
            "manifest replay detected changed inputs and refused to rerun",
            code="manifest_replay_input_changed",
            details={
                "manifest_path": str(manifest_path),
                "workflow": workflow,
                "changed_inputs": [
                    drift.label for drift in input_drift if not drift.matched
                ],
            },
        )

    engine_version_drift = collect_engine_version_drift(
        payload,
        executables=executable_overrides,
    )
    if "run" in payload:
        replay_report = replay_engine_workflow(
            payload,
            replay_out_dir=replay_out_dir,
            executables=executable_overrides,
        )
    else:
        replay_report = replay_composite_workflow(
            payload,
            replay_out_dir=replay_out_dir,
            executables=executable_overrides,
        )
    comparisons = compare_outputs(payload, replay_report)
    return ManifestReplayReport(
        manifest_path=manifest_path,
        workflow=workflow,
        replay_out_dir=replay_out_dir,
        replay_manifest_path=Path(replay_report.manifest_path),
        input_drift=input_drift,
        engine_version_drift=engine_version_drift,
        comparisons=comparisons,
        input_drift_detected=False,
        engine_version_drift_detected=any(
            not drift.matched for drift in engine_version_drift
        ),
        outputs_equivalent=all(
            comparison.status in {"exact", "equivalent"} for comparison in comparisons
        ),
        notes=[
            "replay reran the governed workflow from its recorded manifest inputs and configuration"
        ],
    )
