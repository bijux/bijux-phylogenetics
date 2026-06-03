from __future__ import annotations

from pathlib import Path
import shutil

from .models import CatarrhineDataQualityStressPanelExportResult
from .panel import (
    RAW_ALIGNMENT_NAME,
    RAW_CODING_SEQUENCE_NAME,
    RAW_SEQUENCE_INPUT_NAME,
    RAW_TRAIT_MISMATCH_NAME,
    RAW_TRAITS_NAME,
    RAW_TREE_NAME,
    load_catarrhine_data_quality_stress_panel_dataset,
)


def export_catarrhine_data_quality_stress_panel_dataset(
    destination: Path,
) -> CatarrhineDataQualityStressPanelExportResult:
    """Copy the packaged raw stress inputs and governed expected outputs."""
    dataset = load_catarrhine_data_quality_stress_panel_dataset()
    if destination.exists():
        shutil.rmtree(destination)
    destination.mkdir(parents=True, exist_ok=True)
    readme_path = Path(shutil.copy2(dataset.readme_path, destination / "README.md"))
    raw_root = destination / "raw"
    raw_root.mkdir(parents=True, exist_ok=True)
    raw_alignment_path = Path(
        shutil.copy2(dataset.raw_alignment_path, raw_root / RAW_ALIGNMENT_NAME)
    )
    raw_sequence_input_path = Path(
        shutil.copy2(
            dataset.raw_sequence_input_path,
            raw_root / RAW_SEQUENCE_INPUT_NAME,
        )
    )
    raw_coding_sequences_path = Path(
        shutil.copy2(
            dataset.raw_coding_sequences_path,
            raw_root / RAW_CODING_SEQUENCE_NAME,
        )
    )
    raw_tree_path = Path(shutil.copy2(dataset.raw_tree_path, raw_root / RAW_TREE_NAME))
    raw_traits_path = Path(
        shutil.copy2(dataset.raw_traits_path, raw_root / RAW_TRAITS_NAME)
    )
    raw_trait_mismatch_path = Path(
        shutil.copy2(
            dataset.raw_trait_mismatch_path,
            raw_root / RAW_TRAIT_MISMATCH_NAME,
        )
    )
    expected_output_root = destination / "expected"
    shutil.copytree(dataset.reference_output_root, expected_output_root)
    return CatarrhineDataQualityStressPanelExportResult(
        output_root=destination,
        readme_path=readme_path,
        raw_alignment_path=raw_alignment_path,
        raw_tree_path=raw_tree_path,
        raw_traits_path=raw_traits_path,
        raw_sequence_input_path=raw_sequence_input_path,
        raw_coding_sequences_path=raw_coding_sequences_path,
        raw_trait_mismatch_path=raw_trait_mismatch_path,
        expected_output_root=expected_output_root,
    )
