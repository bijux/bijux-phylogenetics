from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.io.fasta.records import validate_fasta_input

from .models import InfluenzaAHAReferenceDataset

DATASET_ID = "influenza_a_ha_reference_panel"
DATASET_LABEL = "Influenza A hemagglutinin reference panel"
SEQUENCE_TYPE = "dna"
WORKFLOW_PREFIX = "influenza-a-ha-reference-panel"
IQTREE_SEED = 1
IQTREE_THREADS = 1
BOOTSTRAP_REPLICATES = 1000
SOURCE_ACCESSIONS = (
    "NC_002017.1",
    "CY033655.1",
    "CY046787.1",
    "NC_007366.1",
    "NC_007374.1",
    "AY653200.1",
)


def load_influenza_a_ha_reference_dataset() -> InfluenzaAHAReferenceDataset:
    """Expose the packaged influenza A HA panel as a first-class runtime surface."""
    dataset_root = _resource_root()
    sequences_path = dataset_root / "sequences.fasta"
    validation = validate_fasta_input(sequences_path, sequence_type=SEQUENCE_TYPE)
    return InfluenzaAHAReferenceDataset(
        dataset_id=DATASET_ID,
        label=DATASET_LABEL,
        dataset_root=dataset_root,
        sequences_path=sequences_path,
        reference_output_root=dataset_root / "expected",
        sequence_count=validation.summary.sequence_count,
        sequence_type=SEQUENCE_TYPE,
        workflow_prefix=WORKFLOW_PREFIX,
        iqtree_seed=IQTREE_SEED,
        iqtree_threads=IQTREE_THREADS,
        bootstrap_replicates=BOOTSTRAP_REPLICATES,
        source_accessions=SOURCE_ACCESSIONS,
        source_summary=(
            "Published influenza A hemagglutinin segment-4 panel assembled from "
            "stable NCBI GenBank and RefSeq accessions spanning H1N1, H2N2, "
            "H3N2, and H5N1 lineages."
        ),
    )


def _resource_root() -> Path:
    return (
        Path(__file__).resolve().parent.parent.parent
        / "resources"
        / "datasets"
        / "viruses"
        / DATASET_ID
    )
