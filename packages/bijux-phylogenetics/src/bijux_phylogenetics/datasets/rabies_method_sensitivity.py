from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
import shutil

from bijux_phylogenetics.io.fasta import validate_fasta_input

_DATASET_ID = "rabies_method_sensitivity_panel"
_DATASET_LABEL = "Rabies method-sensitivity panel"
_SEQUENCE_TYPE = "dna"
_WORKFLOW_PREFIX = "rabies-method-sensitivity-panel"
_SOURCE_ACCESSIONS = (
    "MG458305",
    "MG458304",
    "PV641713",
    "PX845689",
    "OQ693985",
    "PX845683",
    "PX845681",
    "PX845678",
    "PX845676",
)


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivityVariant:
    """One declared method combination in the packaged sensitivity matrix."""

    variant_id: str
    label: str
    alignment_mode: str
    trimming_mode: str
    trim_gap_threshold: float


@dataclass(slots=True)
class RabiesMethodSensitivityPanelDataset:
    """Packaged rabies dataset for method-sensitivity workflow review."""

    dataset_id: str
    label: str
    dataset_root: Path
    readme_path: Path
    config_path: Path
    sequences_path: Path
    metadata_path: Path
    reference_output_root: Path
    taxon_count: int
    sequence_type: str
    workflow_prefix: str
    outgroup_taxa: tuple[str, ...]
    iqtree_seed: int
    iqtree_threads: int
    bootstrap_replicates: int
    source_accessions: tuple[str, ...]
    variants: tuple[RabiesMethodSensitivityVariant, ...]
    source_summary: str


@dataclass(slots=True)
class RabiesMethodSensitivityPanelExportResult:
    """Materialized copy of the packaged rabies method-sensitivity dataset."""

    output_root: Path
    readme_path: Path
    config_path: Path
    sequences_path: Path
    metadata_path: Path
    expected_output_root: Path


def load_rabies_method_sensitivity_panel_dataset() -> (
    RabiesMethodSensitivityPanelDataset
):
    """Expose the packaged rabies method-sensitivity panel as a first-class surface."""
    dataset_root = _resource_root()
    config = json.loads((dataset_root / "workflow-config.json").read_text(encoding="utf-8"))
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
    metadata_path = Path(shutil.copy2(dataset.metadata_path, destination / "metadata.csv"))
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


def _resource_root() -> Path:
    return (
        Path(__file__).resolve().parent.parent
        / "resources"
        / "datasets"
        / "pathogens"
        / _DATASET_ID
    )
