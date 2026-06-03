from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.evidence.provenance.bundles import (
    bundle_artifact_files,
    validate_artifact_bundle,
)


@dataclass(slots=True)
class BayesianEvidencePackageReport:
    bundle_root: Path
    file_count: int
    input_count: int
    config_count: int
    tree_count: int
    log_count: int
    diagnostic_count: int
    report_count: int
    valid: bool


def build_bayesian_evidence_package(
    *,
    bundle_root: Path,
    input_paths: list[Path],
    config_paths: list[Path],
    tree_paths: list[Path],
    log_paths: list[Path],
    diagnostic_paths: list[Path],
    report_paths: list[Path],
) -> BayesianEvidencePackageReport:
    """Bundle Bayesian inputs and outputs into one checksummed evidence package."""
    if not input_paths:
        raise ValueError("bayesian evidence package requires at least one input path")
    if not config_paths:
        raise ValueError(
            "bayesian evidence package requires at least one config artifact"
        )
    if not tree_paths:
        raise ValueError(
            "bayesian evidence package requires at least one tree artifact"
        )
    if not log_paths:
        raise ValueError(
            "bayesian evidence package requires at least one log or trace artifact"
        )
    if not diagnostic_paths:
        raise ValueError(
            "bayesian evidence package requires at least one diagnostic artifact"
        )
    if not report_paths:
        raise ValueError(
            "bayesian evidence package requires at least one rendered report artifact"
        )

    bundle = bundle_artifact_files(
        input_paths=input_paths,
        output_paths=[
            *config_paths,
            *tree_paths,
            *log_paths,
            *diagnostic_paths,
            *report_paths,
        ],
        bundle_root=bundle_root,
    )
    validation = validate_artifact_bundle(bundle_root)
    return BayesianEvidencePackageReport(
        bundle_root=bundle_root,
        file_count=bundle.file_count,
        input_count=len(input_paths),
        config_count=len(config_paths),
        tree_count=len(tree_paths),
        log_count=len(log_paths),
        diagnostic_count=len(diagnostic_paths),
        report_count=len(report_paths),
        valid=validation.valid,
    )
