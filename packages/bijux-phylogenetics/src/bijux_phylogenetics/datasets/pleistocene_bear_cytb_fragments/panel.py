from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.io.fasta.records import validate_fasta_input

from .models import PleistoceneBearCytbFragmentDataset

DATASET_ID = "pleistocene_bear_cytb_fragments"
DATASET_LABEL = "Pleistocene bear CYTB fragment panel"
SEQUENCE_TYPE = "dna"
WORKFLOW_PREFIX = "pleistocene-bear-cytb-fragments"
IQTREE_SEED = 1
IQTREE_THREADS = 1
BOOTSTRAP_REPLICATES = 1000
SITE_MISSINGNESS_THRESHOLD = 0.15
SEQUENCE_MISSINGNESS_THRESHOLD = 0.15
DEGRADED_SEQUENCE_IDS = (
    "cave_bear_ud1838_fragment",
    "cave_bear_wk01_fragment",
)
SOURCE_ACCESSIONS = (
    "OQ318974.1",
    "NC_003428.1",
    "OQ318956.1",
    "KX641337.1",
    "KX641335.1",
)


def load_pleistocene_bear_cytb_fragment_dataset() -> PleistoceneBearCytbFragmentDataset:
    """Expose the packaged bear fragment panel as a first-class runtime surface."""
    dataset_root = _resource_root()
    sequences_path = dataset_root / "sequences.fasta"
    validation = validate_fasta_input(sequences_path, sequence_type=SEQUENCE_TYPE)
    return PleistoceneBearCytbFragmentDataset(
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
        site_missingness_threshold=SITE_MISSINGNESS_THRESHOLD,
        sequence_missingness_threshold=SEQUENCE_MISSINGNESS_THRESHOLD,
        degraded_sequence_ids=DEGRADED_SEQUENCE_IDS,
        source_accessions=SOURCE_ACCESSIONS,
        source_summary=(
            "Modern bear CYTB references paired with real ancient cave-bear CYTB "
            "sequences reduced to short fragment-style inputs with explicit "
            "internal missing-data blocks."
        ),
    )


def _resource_root() -> Path:
    return (
        Path(__file__).resolve().parent.parent.parent
        / "resources"
        / "datasets"
        / "ancient_dna"
        / DATASET_ID
    )
