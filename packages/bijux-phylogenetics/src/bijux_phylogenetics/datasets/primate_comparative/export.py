from __future__ import annotations

from pathlib import Path
import shutil

from .models import PrimateComparativeDatasetExportResult
from .panel import load_primate_comparative_dataset


def export_primate_comparative_dataset(
    destination: Path,
) -> PrimateComparativeDatasetExportResult:
    """Copy the packaged primate dataset and reference outputs to one directory."""
    dataset = load_primate_comparative_dataset()
    if destination.exists():
        shutil.rmtree(destination)
    destination.mkdir(parents=True, exist_ok=True)
    readme_path = shutil.copy2(
        dataset.dataset_root / "README.md",
        destination / "README.md",
    )
    tree_path = shutil.copy2(dataset.tree_path, destination / "tree.nwk")
    traits_path = shutil.copy2(dataset.traits_path, destination / "traits.csv")
    expected_output_root = destination / "expected"
    shutil.copytree(dataset.reference_output_root, expected_output_root)
    return PrimateComparativeDatasetExportResult(
        output_root=destination,
        readme_path=Path(readme_path),
        tree_path=Path(tree_path),
        traits_path=Path(traits_path),
        expected_output_root=expected_output_root,
    )
