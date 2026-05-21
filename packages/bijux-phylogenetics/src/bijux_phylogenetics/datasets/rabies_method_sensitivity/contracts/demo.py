from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .bundle import RabiesMethodSensitivityPanelWorkflowBundle
from .dataset import (
    RabiesMethodSensitivityPanelDataset,
    RabiesMethodSensitivityPanelExportResult,
)


@dataclass(slots=True)
class RabiesMethodSensitivityPanelDemoResult:
    """Dataset export plus workflow outputs for the public method-sensitivity demo."""

    output_root: Path
    dataset: RabiesMethodSensitivityPanelDataset
    dataset_export: RabiesMethodSensitivityPanelExportResult
    workflow_bundle: RabiesMethodSensitivityPanelWorkflowBundle
    overview_path: Path
