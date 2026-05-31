from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.datasets.continuous_mode_recovery import (
    load_continuous_mode_recovery_panel_dataset,
)
from bijux_phylogenetics.datasets.discrete_mode_recovery import (
    load_discrete_mode_recovery_panel_dataset,
)
from bijux_phylogenetics.datasets.known_answer_reference import (
    load_known_answer_reference_dataset,
)

from .models import MacroevolutionRecoverySuiteDataset

DATASET_ID = "macroevolution_recovery_suite"
DATASET_LABEL = "Macroevolution recovery suite"


def load_macroevolution_recovery_suite_dataset() -> MacroevolutionRecoverySuiteDataset:
    """Expose the governed macroevolution recovery suite as one public surface."""
    continuous_panel = load_continuous_mode_recovery_panel_dataset()
    discrete_panel = load_discrete_mode_recovery_panel_dataset()
    known_answer_panel = load_known_answer_reference_dataset()
    dataset_root = (
        Path(__file__).resolve().parents[2]
        / "resources"
        / "datasets"
        / "simulation"
        / DATASET_ID
    )
    return MacroevolutionRecoverySuiteDataset(
        dataset_id=DATASET_ID,
        label=DATASET_LABEL,
        dataset_root=dataset_root,
        reference_output_root=dataset_root / "expected",
        continuous_panel=continuous_panel,
        discrete_panel=discrete_panel,
        known_answer_panel=known_answer_panel,
        component_count=3,
        geiger_component_count=2,
        max_taxon_count=max(
            continuous_panel.taxon_count,
            discrete_panel.taxon_count,
            known_answer_panel.taxon_count,
        ),
        total_recovery_case_count=(
            continuous_panel.case_count + discrete_panel.case_count + 11
        ),
        geiger_recovery_case_count=(
            continuous_panel.case_count + discrete_panel.case_count
        ),
        truth_threshold_row_count=11,
        source_summary=(
            "Unified governed recovery suite over the continuous, discrete, and "
            "known-answer macroevolution panels."
        ),
    )
