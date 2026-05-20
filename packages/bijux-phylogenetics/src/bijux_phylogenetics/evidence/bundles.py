from __future__ import annotations

from pathlib import Path

from .provenance.bundles import (
    ArtifactBundleReport as EvidenceBundleReport,
)
from .provenance.bundles import (
    ArtifactBundleValidationReport as EvidenceValidationReport,
)
from .provenance.bundles import (
    bundle_artifact_directories,
    bundle_artifact_files,
    validate_artifact_bundle,
)


def bundle_directory(
    input_roots: list[Path], output_roots: list[Path], bundle_root: Path
) -> EvidenceBundleReport:
    """Copy explicit input and output directories into a checksummed evidence bundle."""
    return bundle_artifact_directories(input_roots, output_roots, bundle_root)


def bundle_file_paths(
    input_paths: list[Path], output_paths: list[Path], bundle_root: Path
) -> EvidenceBundleReport:
    """Copy explicit input and output files into a checksummed evidence bundle."""
    return bundle_artifact_files(input_paths, output_paths, bundle_root)


def validate_bundle(bundle_root: Path) -> EvidenceValidationReport:
    """Validate an evidence bundle against its manifest checksums."""
    return validate_artifact_bundle(bundle_root)
