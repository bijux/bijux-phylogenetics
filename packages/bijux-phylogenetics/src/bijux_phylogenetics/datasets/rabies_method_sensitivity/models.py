from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .contracts.bundle import RabiesMethodSensitivityPanelWorkflowBundle
from .contracts.dataset import (
    RabiesMethodSensitivityPanelDataset,
    RabiesMethodSensitivityPanelExportResult,
    RabiesMethodSensitivityVariant,
)
from .contracts.workflow import (
    RabiesMethodSensitivityCladeRow,
    RabiesMethodSensitivityConclusionRow,
    RabiesMethodSensitivityPanelWorkflowReport,
    RabiesMethodSensitivityPreprocessingComparisonRow,
    RabiesMethodSensitivityTaskRecord,
    RabiesMethodSensitivityVariantRun,
)

__all__ = [
    "RabiesMethodSensitivityCladeRow",
    "RabiesMethodSensitivityConclusionRow",
    "RabiesMethodSensitivityPanelDataset",
    "RabiesMethodSensitivityPanelDemoResult",
    "RabiesMethodSensitivityPanelExportResult",
    "RabiesMethodSensitivityPanelWorkflowBundle",
    "RabiesMethodSensitivityPanelWorkflowReport",
    "RabiesMethodSensitivityPreprocessingComparisonRow",
    "RabiesMethodSensitivityTaskRecord",
    "RabiesMethodSensitivityVariant",
    "RabiesMethodSensitivityVariantRun",
]
@dataclass(slots=True)
class RabiesMethodSensitivityPanelDemoResult:
    """Dataset export plus workflow outputs for the public method-sensitivity demo."""

    output_root: Path
    dataset: RabiesMethodSensitivityPanelDataset
    dataset_export: RabiesMethodSensitivityPanelExportResult
    workflow_bundle: RabiesMethodSensitivityPanelWorkflowBundle
    overview_path: Path
