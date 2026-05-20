from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.evidence.provenance.bundles import (
    ArtifactBundleReport,
    ArtifactBundleValidationReport,
    bundle_artifact_files,
    validate_artifact_bundle,
)

from ..common import load_engine_manifest


@dataclass(slots=True)
class InferenceEvidenceBundleResult:
    bundle_root: Path
    workflow_count: int
    input_file_count: int
    output_file_count: int
    manifest_file_count: int
    bundle: ArtifactBundleReport
    validation: ArtifactBundleValidationReport


def bundle_inference_workflow_evidence(
    manifest_paths: list[Path],
    *,
    bundle_root: Path,
) -> InferenceEvidenceBundleResult:
    """Bundle one or more inference workflow manifests together with their declared inputs and outputs."""
    if not manifest_paths:
        raise ValueError("at least one inference manifest is required")
    normalized_manifests = [path.resolve() for path in manifest_paths]
    input_paths: list[Path] = []
    output_paths: list[Path] = []
    for manifest_path in normalized_manifests:
        manifest = load_engine_manifest(manifest_path)
        input_paths.extend(Path(path) for path in manifest["input_paths"])
        output_paths.extend(
            Path(path) for path in dict(manifest["output_paths"]).values()
        )
    unique_inputs = _deduplicate_existing_paths(input_paths)
    unique_outputs = _deduplicate_existing_paths([*output_paths, *normalized_manifests])
    bundle = bundle_artifact_files(unique_inputs, unique_outputs, bundle_root)
    validation = validate_artifact_bundle(bundle_root)
    return InferenceEvidenceBundleResult(
        bundle_root=bundle_root,
        workflow_count=len(normalized_manifests),
        input_file_count=len(unique_inputs),
        output_file_count=len(_deduplicate_existing_paths(output_paths)),
        manifest_file_count=len(normalized_manifests),
        bundle=bundle,
        validation=validation,
    )


def _deduplicate_existing_paths(paths: list[Path]) -> list[Path]:
    seen: set[Path] = set()
    deduplicated: list[Path] = []
    for path in paths:
        resolved = path.resolve()
        if not resolved.exists() or resolved in seen:
            continue
        seen.add(resolved)
        deduplicated.append(resolved)
    return deduplicated
