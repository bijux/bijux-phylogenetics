from __future__ import annotations

from pathlib import Path
import shutil

from bijux_phylogenetics.datasets.continuous_mode_recovery import (
    export_continuous_mode_recovery_panel_dataset,
)
from bijux_phylogenetics.datasets.discrete_mode_recovery import (
    export_discrete_mode_recovery_panel_dataset,
)
from bijux_phylogenetics.datasets.known_answer_reference import (
    export_known_answer_reference_dataset,
)

from .models import MacroevolutionRecoverySuiteExportResult
from .panel import load_macroevolution_recovery_suite_dataset


def export_macroevolution_recovery_suite_dataset(
    destination: Path,
) -> MacroevolutionRecoverySuiteExportResult:
    """Copy the packaged macroevolution recovery suite and its component panels."""
    dataset = load_macroevolution_recovery_suite_dataset()
    if destination.exists():
        shutil.rmtree(destination)
    destination.mkdir(parents=True, exist_ok=True)
    readme_path = shutil.copy2(
        dataset.dataset_root / "README.md",
        destination / "README.md",
    )
    component_root = destination / "components"
    component_root.mkdir(parents=True, exist_ok=True)
    continuous_panel_export = export_continuous_mode_recovery_panel_dataset(
        component_root / "continuous-mode-recovery-panel"
    )
    discrete_panel_export = export_discrete_mode_recovery_panel_dataset(
        component_root / "discrete-mode-recovery-panel"
    )
    known_answer_panel_export = export_known_answer_reference_dataset(
        component_root / "known-answer-reference-panel"
    )
    expected_output_root = destination / "expected"
    shutil.copytree(dataset.reference_output_root, expected_output_root)
    return MacroevolutionRecoverySuiteExportResult(
        output_root=destination,
        readme_path=Path(readme_path),
        component_root=component_root,
        continuous_panel_export=continuous_panel_export,
        discrete_panel_export=discrete_panel_export,
        known_answer_panel_export=known_answer_panel_export,
        expected_output_root=expected_output_root,
    )
