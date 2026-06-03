from __future__ import annotations

from pathlib import Path
import shutil

from .models import InfluenzaAHAReferenceDatasetExportResult
from .panel import load_influenza_a_ha_reference_dataset


def export_influenza_a_ha_reference_dataset(
    destination: Path,
) -> InfluenzaAHAReferenceDatasetExportResult:
    """Copy the packaged influenza A HA dataset and reference outputs to one directory."""
    dataset = load_influenza_a_ha_reference_dataset()
    if destination.exists():
        shutil.rmtree(destination)
    destination.mkdir(parents=True, exist_ok=True)
    readme_path = shutil.copy2(
        dataset.dataset_root / "README.md",
        destination / "README.md",
    )
    sequences_path = shutil.copy2(
        dataset.sequences_path,
        destination / "sequences.fasta",
    )
    expected_output_root = destination / "expected"
    shutil.copytree(dataset.reference_output_root, expected_output_root)
    return InfluenzaAHAReferenceDatasetExportResult(
        output_root=destination,
        readme_path=Path(readme_path),
        sequences_path=Path(sequences_path),
        expected_output_root=expected_output_root,
    )
