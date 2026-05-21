from __future__ import annotations

import json

from bijux_phylogenetics.io.fasta.records import validate_fasta_input

from ..models import RabiesMethodSensitivityPanelDataset, RabiesMethodSensitivityVariant
from ..shared import (
    _DATASET_ID,
    _DATASET_LABEL,
    _SEQUENCE_TYPE,
    _SOURCE_ACCESSIONS,
    _WORKFLOW_PREFIX,
    _resource_root,
)


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
