from __future__ import annotations

from pathlib import Path
import shutil

from ..models import (
    RabiesMethodSensitivityPanelWorkflowBundle,
    RabiesMethodSensitivityPanelWorkflowReport,
)
from .finalization import _build_workflow_bundle
from .manifest_artifacts import _write_bundle_manifest_artifacts
from .report_artifacts import _write_bundle_report_artifacts
from .reproducibility_artifacts import _write_reproducibility_artifacts
from .shared import (
    _format_float,
    _format_optional_bool,
    _format_optional_float,
    _write_tsv,
)
from .slurm_artifacts import _write_slurm_bundle_artifacts
from .variant_artifacts import (
    _copy_output,
    _copy_task_logs,
    _write_rooting_summary_table,
    _write_variant_outputs,
)
from .workflow_artifacts import _write_workflow_bundle_artifacts


def write_rabies_method_sensitivity_panel_workflow_bundle(
    output_root: Path,
    report: RabiesMethodSensitivityPanelWorkflowReport,
) -> RabiesMethodSensitivityPanelWorkflowBundle:
    """Write the governed reviewer-facing bundle for the method-sensitivity workflow."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    workflow_artifacts = _write_workflow_bundle_artifacts(
        output_root,
        report,
    )
    slurm_artifacts = _write_slurm_bundle_artifacts(
        output_root,
        report,
        execution_record_path=workflow_artifacts.execution_record_path,
    )
    manifest_artifacts = _write_bundle_manifest_artifacts(
        output_root,
        report,
        workflow_artifacts,
        slurm_artifacts,
    )
    reproducibility_artifacts = _write_reproducibility_artifacts(output_root, report)
    report_artifacts = _write_bundle_report_artifacts(
        output_root,
        report,
        workflow_artifacts,
        manifest_artifacts,
        slurm_artifacts,
        reproducibility_artifacts,
    )
    return _build_workflow_bundle(
        output_root,
        report,
        workflow_artifacts,
        manifest_artifacts,
        slurm_artifacts,
        reproducibility_artifacts,
        report_artifacts,
    )
