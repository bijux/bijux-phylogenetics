from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ..audit import (
    audit_rabies_method_sensitivity_workflow_bundle,
    write_rabies_method_sensitivity_reproducibility_audit_json,
    write_rabies_method_sensitivity_reproducibility_checks_table,
    write_rabies_method_sensitivity_variant_audit_table,
)
from ..models import RabiesMethodSensitivityPanelWorkflowReport


@dataclass(slots=True)
class RabiesMethodSensitivityReproducibilityArtifacts:
    """Owned reproducibility review artifacts for the workflow bundle."""

    report: object
    checks_path: Path
    variant_audit_path: Path
    audit_path: Path


def _write_reproducibility_artifacts(
    output_root: Path,
    report: RabiesMethodSensitivityPanelWorkflowReport,
) -> RabiesMethodSensitivityReproducibilityArtifacts:
    reproducibility_report = audit_rabies_method_sensitivity_workflow_bundle(
        output_root,
        sequences_path=report.dataset.sequences_path,
        metadata_path=report.dataset.metadata_path,
    )
    checks_path = write_rabies_method_sensitivity_reproducibility_checks_table(
        output_root / "reproducibility-checks.tsv",
        reproducibility_report,
    )
    variant_audit_path = write_rabies_method_sensitivity_variant_audit_table(
        output_root / "reproducibility-variants.tsv",
        reproducibility_report,
    )
    audit_path = write_rabies_method_sensitivity_reproducibility_audit_json(
        output_root / "reproducibility-audit.json",
        reproducibility_report,
    )
    return RabiesMethodSensitivityReproducibilityArtifacts(
        report=reproducibility_report,
        checks_path=checks_path,
        variant_audit_path=variant_audit_path,
        audit_path=audit_path,
    )
