from __future__ import annotations

from pathlib import Path

from .cleanup import build_catarrhine_data_quality_stress_panel_workflow_report
from .models import CatarrhineDataQualityStressPanelWorkflowReport
from .panel import load_catarrhine_data_quality_stress_panel_dataset


def run_catarrhine_data_quality_stress_panel_workflow(
    out_dir: Path,
) -> CatarrhineDataQualityStressPanelWorkflowReport:
    """Audit the raw stress panel and write one cleaned comparative subset."""
    dataset = load_catarrhine_data_quality_stress_panel_dataset()
    return build_catarrhine_data_quality_stress_panel_workflow_report(dataset, out_dir)
