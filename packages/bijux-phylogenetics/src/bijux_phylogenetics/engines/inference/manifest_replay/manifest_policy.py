from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.runtime.errors import EngineWorkflowError

from ...common import build_file_checksums, load_engine_manifest, read_engine_version
from .contracts import ManifestReplayDrift

ENGINE_VERSION_ARGS: dict[str, tuple[str, ...]] = {
    "mafft": ("--version",),
    "trimal": ("--version",),
    "iqtree": ("--version",),
    "fasttree": ("-help",),
    "mrbayes": ("-v",),
    "beast": ("-version",),
}


def payload_workflow(payload: dict[str, Any]) -> str:
    workflow = payload.get("workflow")
    if workflow is None:
        raise EngineWorkflowError(
            "manifest replay requires a workflow identifier",
            code="manifest_replay_missing_workflow",
        )
    return str(workflow)


def path_map(values: dict[str, Any]) -> dict[str, Path]:
    return {str(key): Path(value) for key, value in values.items()}


def recorded_input_paths(payload: dict[str, Any]) -> list[Path]:
    if "input_paths" in payload:
        return [Path(path) for path in payload["input_paths"]]
    return [Path(payload["input_path"])]


def default_replay_out_dir(manifest_path: Path) -> Path:
    stem = manifest_path.name.removesuffix(".manifest.json")
    return manifest_path.parent / "replay" / stem


def engine_key_from_name(engine_name: str) -> str:
    normalized = engine_name.strip().lower()
    if normalized in {"mafft", "trimal", "iqtree", "fasttree", "mrbayes", "beast"}:
        return normalized
    raise EngineWorkflowError(
        f"unsupported replay engine name: {engine_name}",
        code="manifest_replay_unknown_engine",
        details={"engine_name": engine_name},
    )


def version_drift_for_engine_manifest(
    payload: dict[str, Any],
    *,
    executables: dict[str, str | Path | None],
) -> list[ManifestReplayDrift]:
    run_payload = dict(payload["run"])
    version_payload = dict(run_payload["version"])
    engine_key = engine_key_from_name(str(payload["engine_name"]))
    executable = executables.get(engine_key) or str(run_payload["executable"])
    version = read_engine_version(
        str(payload["engine_name"]),
        executable,
        version_args=ENGINE_VERSION_ARGS[engine_key],
        timeout_seconds=(
            None
            if run_payload.get("timeout_seconds") is None
            else float(run_payload["timeout_seconds"])
        ),
    )
    expected = str(version_payload["text"])
    observed = version.text
    return [
        ManifestReplayDrift(
            kind="engine-version",
            label=str(payload["engine_name"]),
            expected=expected,
            observed=observed,
            matched=expected == observed,
        )
    ]


def version_drift_for_large_alignment(
    payload: dict[str, Any],
    *,
    executables: dict[str, str | Path | None],
) -> list[ManifestReplayDrift]:
    engine_key = engine_key_from_name(str(payload["engine_name"]))
    command = [str(item) for item in payload["command"]]
    executable = executables.get(engine_key) or command[0]
    version = read_engine_version(
        str(payload["engine_name"]),
        executable,
        version_args=ENGINE_VERSION_ARGS[engine_key],
    )
    expected = str(payload["engine_version_text"])
    observed = version.text
    return [
        ManifestReplayDrift(
            kind="engine-version",
            label=str(payload["engine_name"]),
            expected=expected,
            observed=observed,
            matched=expected == observed,
        )
    ]


def collect_engine_version_drift(
    payload: dict[str, Any],
    *,
    executables: dict[str, str | Path | None],
) -> list[ManifestReplayDrift]:
    if "run" in payload:
        return version_drift_for_engine_manifest(payload, executables=executables)
    if "step_manifests" in payload:
        drifts: list[ManifestReplayDrift] = []
        for step_manifest_path in path_map(dict(payload["step_manifests"])).values():
            drifts.extend(
                version_drift_for_engine_manifest(
                    load_engine_manifest(step_manifest_path),
                    executables=executables,
                )
            )
        return drifts
    if "engine_version_text" in payload and "command" in payload:
        return version_drift_for_large_alignment(payload, executables=executables)
    raise EngineWorkflowError(
        "manifest replay could not determine engine versions from the manifest",
        code="manifest_replay_missing_engine_versions",
    )


def collect_input_drift(payload: dict[str, Any]) -> list[ManifestReplayDrift]:
    recorded = {
        str(key): str(value)
        for key, value in dict(payload.get("input_checksums", {})).items()
    }
    observed = build_file_checksums([Path(path) for path in recorded])
    drifts: list[ManifestReplayDrift] = []
    for path_text, expected in recorded.items():
        current = observed.get(path_text)
        drifts.append(
            ManifestReplayDrift(
                kind="input-checksum",
                label=path_text,
                expected=expected,
                observed=current,
                matched=current == expected,
            )
        )
    return drifts
