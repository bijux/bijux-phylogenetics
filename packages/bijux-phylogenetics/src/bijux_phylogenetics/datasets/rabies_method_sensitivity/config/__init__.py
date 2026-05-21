from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import shutil

from ..models import (
    RabiesMethodSensitivityPanelDataset,
    RabiesMethodSensitivityPanelExportResult,
    RabiesMethodSensitivityVariant,
)
from .panel_dataset import load_rabies_method_sensitivity_panel_dataset

__all__ = [
    "export_rabies_method_sensitivity_panel_dataset",
    "load_rabies_method_sensitivity_panel_dataset",
]


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


def _resolve_selected_variant_dataset(
    dataset: RabiesMethodSensitivityPanelDataset,
    *,
    variant_ids: tuple[str, ...] | None,
) -> RabiesMethodSensitivityPanelDataset:
    """Return either the full dataset or an explicit variant-scoped subset."""
    if variant_ids is None:
        return dataset
    if not variant_ids:
        raise ValueError("variant_ids must not be empty when provided")
    variants_by_id = {variant.variant_id: variant for variant in dataset.variants}
    selected_variants: list[RabiesMethodSensitivityVariant] = []
    seen_variant_ids: set[str] = set()
    for variant_id in variant_ids:
        if variant_id in seen_variant_ids:
            raise ValueError(f"duplicate variant_id requested: {variant_id}")
        seen_variant_ids.add(variant_id)
        variant = variants_by_id.get(variant_id)
        if variant is None:
            known = ", ".join(sorted(variants_by_id))
            raise ValueError(
                f"unknown variant_id '{variant_id}'; known variants: {known}"
            )
        selected_variants.append(variant)
    return replace(dataset, variants=tuple(selected_variants))
