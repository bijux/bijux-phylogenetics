from __future__ import annotations

from bijux_phylogenetics.comparative.continuous_mode_recovery import (
    run_continuous_mode_recovery,
)

from .models import ContinuousModeRecoveryPanelWorkflowReport
from .panel import load_continuous_mode_recovery_panel_dataset
from .scenarios import load_continuous_mode_recovery_panel_scenarios


def run_continuous_mode_recovery_panel_workflow() -> (
    ContinuousModeRecoveryPanelWorkflowReport
):
    """Run the governed recovery workflow over the packaged continuous-mode panel."""
    dataset = load_continuous_mode_recovery_panel_dataset()
    recovery_report = run_continuous_mode_recovery(
        dataset.default_tree_path,
        load_continuous_mode_recovery_panel_scenarios(dataset),
    )
    return ContinuousModeRecoveryPanelWorkflowReport(
        dataset=dataset,
        recovery_report=recovery_report,
    )
