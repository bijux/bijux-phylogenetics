from __future__ import annotations

from pathlib import Path
import shutil

from .models import ContinuousModeRecoveryPanelExportResult
from .panel import DEFAULT_TREE_FILE, load_continuous_mode_recovery_panel_dataset


def export_continuous_mode_recovery_panel_dataset(
    destination: Path,
) -> ContinuousModeRecoveryPanelExportResult:
    """Copy the packaged continuous-mode recovery panel and reference outputs."""
    dataset = load_continuous_mode_recovery_panel_dataset()
    if destination.exists():
        shutil.rmtree(destination)
    destination.mkdir(parents=True, exist_ok=True)
    readme_path = shutil.copy2(
        dataset.dataset_root / "README.md",
        destination / "README.md",
    )
    reference_tree_root = destination / "trees"
    shutil.copytree(dataset.dataset_root / "trees", reference_tree_root)
    simulation_cases_path = shutil.copy2(
        dataset.simulation_cases_path,
        destination / "simulation-cases.tsv",
    )
    expected_output_root = destination / "expected"
    shutil.copytree(dataset.reference_output_root, expected_output_root)
    return ContinuousModeRecoveryPanelExportResult(
        output_root=destination,
        readme_path=Path(readme_path),
        default_tree_path=reference_tree_root / Path(DEFAULT_TREE_FILE).name,
        reference_tree_root=reference_tree_root,
        simulation_cases_path=Path(simulation_cases_path),
        expected_output_root=expected_output_root,
    )
