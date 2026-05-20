from __future__ import annotations

from pathlib import Path
import shutil

from .models import GnathostomeOrthologProteinBenchmarkExportResult
from .panel import load_gnathostome_ortholog_protein_benchmark_dataset


def export_gnathostome_ortholog_protein_benchmark_dataset(
    destination: Path,
) -> GnathostomeOrthologProteinBenchmarkExportResult:
    """Copy the packaged protein benchmark dataset and governed outputs."""
    dataset = load_gnathostome_ortholog_protein_benchmark_dataset()
    if destination.exists():
        shutil.rmtree(destination)
    destination.mkdir(parents=True, exist_ok=True)
    readme_path = Path(shutil.copy2(dataset.readme_path, destination / "README.md"))
    sequences_path = Path(
        shutil.copy2(dataset.sequences_path, destination / "sequences.fasta")
    )
    expected_output_root = destination / "expected"
    shutil.copytree(dataset.reference_output_root, expected_output_root)
    return GnathostomeOrthologProteinBenchmarkExportResult(
        output_root=destination,
        readme_path=readme_path,
        sequences_path=sequences_path,
        expected_output_root=expected_output_root,
    )
