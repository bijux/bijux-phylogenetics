from __future__ import annotations

from pathlib import Path
import shutil

from ..models import RabiesMethodSensitivityPanelExportResult
from .panel_dataset import load_rabies_method_sensitivity_panel_dataset


def export_rabies_method_sensitivity_panel_dataset(
    destination: Path,
) -> RabiesMethodSensitivityPanelExportResult:
    """Copy the packaged rabies method-sensitivity dataset and expected outputs."""
    dataset = load_rabies_method_sensitivity_panel_dataset()
    if destination.exists():
        shutil.rmtree(destination)
    destination.mkdir(parents=True, exist_ok=True)
    readme_path = Path(shutil.copy2(dataset.readme_path, destination / "README.md"))
    config_path = Path(
        shutil.copy2(dataset.config_path, destination / "workflow-config.json")
    )
    sequences_path = Path(
        shutil.copy2(dataset.sequences_path, destination / "sequences.fasta")
    )
    metadata_path = Path(
        shutil.copy2(dataset.metadata_path, destination / "metadata.csv")
    )
    expected_output_root = destination / "expected"
    shutil.copytree(dataset.reference_output_root, expected_output_root)
    return RabiesMethodSensitivityPanelExportResult(
        output_root=destination,
        readme_path=readme_path,
        config_path=config_path,
        sequences_path=sequences_path,
        metadata_path=metadata_path,
        expected_output_root=expected_output_root,
    )
