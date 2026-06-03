from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.io.fasta import load_fasta_records
from bijux_phylogenetics.io.fasta.records import validate_fasta_input

from .models import GnathostomeOrthologProteinBenchmarkDataset

DATASET_ID = "gnathostome_ortholog_protein_benchmark"
DATASET_LABEL = "Gnathostome ortholog protein benchmark"
SEQUENCE_TYPE = "protein"
WORKFLOW_PREFIX = "gnathostome-ortholog-protein-benchmark"
IQTREE_SEED = 1
IQTREE_THREADS = 1
BOOTSTRAP_REPLICATES = 1000
SOURCE_REFERENCE = "trimAl governed reference corpus example.009.AA"
SOURCE_TRANSFORMATION = (
    "Removed alignment gap characters and placeholder missing-data marks from "
    "the aligned reference FASTA to recover the raw protein inputs used for "
    "the packaged end-to-end benchmark."
)


def load_gnathostome_ortholog_protein_benchmark_dataset() -> (
    GnathostomeOrthologProteinBenchmarkDataset
):
    """Expose the packaged gnathostome protein panel as a first-class surface."""
    dataset_root = _resource_root()
    sequences_path = dataset_root / "sequences.fasta"
    validation = validate_fasta_input(sequences_path, sequence_type=SEQUENCE_TYPE)
    sequence_lengths = [
        len(record.sequence) for record in load_fasta_records(sequences_path)
    ]
    return GnathostomeOrthologProteinBenchmarkDataset(
        dataset_id=DATASET_ID,
        label=DATASET_LABEL,
        dataset_root=dataset_root,
        readme_path=dataset_root / "README.md",
        sequences_path=sequences_path,
        reference_output_root=dataset_root / "expected",
        sequence_count=validation.summary.sequence_count,
        sequence_type=SEQUENCE_TYPE,
        workflow_prefix=WORKFLOW_PREFIX,
        iqtree_seed=IQTREE_SEED,
        iqtree_threads=IQTREE_THREADS,
        bootstrap_replicates=BOOTSTRAP_REPLICATES,
        source_reference=SOURCE_REFERENCE,
        source_transformation=SOURCE_TRANSFORMATION,
        minimum_sequence_length=min(sequence_lengths),
        maximum_sequence_length=max(sequence_lengths),
        source_summary=(
            "Nine ungapped gnathostome ortholog proteins packaged as one public "
            "amino-acid benchmark for MAFFT alignment, trimAl trimming, IQ-TREE "
            "protein model selection, maximum-likelihood inference, and "
            "bootstrap support review."
        ),
    )


def _resource_root() -> Path:
    return (
        Path(__file__).resolve().parent.parent.parent
        / "resources"
        / "datasets"
        / "vertebrates"
        / DATASET_ID
    )
