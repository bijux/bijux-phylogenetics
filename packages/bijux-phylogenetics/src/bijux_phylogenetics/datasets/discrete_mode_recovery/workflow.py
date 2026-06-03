from __future__ import annotations

from bijux_phylogenetics.comparative.discrete_mode_recovery import (
    run_discrete_mode_recovery,
)

from .models import DiscreteModeRecoveryPanelWorkflowReport
from .panel import load_discrete_mode_recovery_panel_dataset
from .scenarios import load_discrete_mode_recovery_panel_scenarios


def run_discrete_mode_recovery_panel_workflow() -> (
    DiscreteModeRecoveryPanelWorkflowReport
):
    """Run the governed recovery workflow over the packaged discrete-mode panel."""
    dataset = load_discrete_mode_recovery_panel_dataset()
    recovery_report = run_discrete_mode_recovery(
        dataset.default_tree_path,
        load_discrete_mode_recovery_panel_scenarios(dataset),
    )
    return DiscreteModeRecoveryPanelWorkflowReport(
        dataset=dataset,
        recovery_report=recovery_report,
    )
