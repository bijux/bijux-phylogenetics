from __future__ import annotations

from dataclasses import replace
import json
from pathlib import Path
import shutil

from bijux_phylogenetics.io.fasta.records import validate_fasta_input

from .models import (
    RabiesMethodSensitivityPanelDataset,
    RabiesMethodSensitivityPanelExportResult,
    RabiesMethodSensitivityVariant,
)
from .shared import (
    _DATASET_ID,
    _DATASET_LABEL,
    _SEQUENCE_TYPE,
    _SOURCE_ACCESSIONS,
    _WORKFLOW_PREFIX,
    _resource_root,
)

__all__ = [
    "export_rabies_method_sensitivity_panel_dataset",
    "load_rabies_method_sensitivity_panel_dataset",
]


def load_rabies_method_sensitivity_panel_dataset() -> (
    RabiesMethodSensitivityPanelDataset
):
    """Expose the packaged rabies method-sensitivity panel as a first-class surface."""
    dataset_root = _resource_root()
    config = json.loads(
        (dataset_root / "workflow-config.json").read_text(encoding="utf-8")
    )
    sequences_path = dataset_root / "sequences.fasta"
    metadata_path = dataset_root / "metadata.csv"
    validation = validate_fasta_input(sequences_path, sequence_type=_SEQUENCE_TYPE)
    variants = tuple(
        RabiesMethodSensitivityVariant(
            variant_id=str(item["variant_id"]),
            label=str(item["label"]),
            alignment_mode=str(item["alignment_mode"]),
            trimming_mode=str(item["trimming_mode"]),
            trim_gap_threshold=float(item["trim_gap_threshold"]),
        )
        for item in config["variants"]
    )
    return RabiesMethodSensitivityPanelDataset(
        dataset_id=_DATASET_ID,
        label=_DATASET_LABEL,
        dataset_root=dataset_root,
        readme_path=dataset_root / "README.md",
        config_path=dataset_root / "workflow-config.json",
        sequences_path=sequences_path,
        metadata_path=metadata_path,
        reference_output_root=dataset_root / "expected",
        taxon_count=validation.summary.sequence_count,
        sequence_type=_SEQUENCE_TYPE,
        workflow_prefix=_WORKFLOW_PREFIX,
        outgroup_taxa=tuple(str(value) for value in config["outgroup_taxa"]),
        iqtree_seed=int(config["iqtree_seed"]),
        iqtree_threads=int(config["iqtree_threads"]),
        bootstrap_replicates=int(config["bootstrap_replicates"]),
        parallel_workers=int(config.get("parallel_workers", 1)),
        source_accessions=_SOURCE_ACCESSIONS,
        variants=variants,
        source_summary=(
            "Real rabies virus nucleoprotein sequences packaged with grouped host "
            "and geography metadata for one owned sensitivity workflow that checks "
            "how alignment, trimming, and inference-engine choices change or "
            "preserve rooted biological conclusions."
        ),
    )


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
