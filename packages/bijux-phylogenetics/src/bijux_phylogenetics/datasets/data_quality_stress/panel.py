from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.io.fasta.records import validate_fasta_input

from .models import CatarrhineDataQualityStressPanelDataset
from .traits import load_permissive_trait_rows

DATASET_ID = "catarrhine_data_quality_stress_panel"
DATASET_LABEL = "Catarrhine data quality stress panel"
RAW_ALIGNMENT_NAME = "alignment.fasta"
RAW_SEQUENCE_INPUT_NAME = "sequence-input.fasta"
RAW_CODING_SEQUENCE_NAME = "coding-sequences.fasta"
RAW_TREE_NAME = "tree.nwk"
RAW_TRAITS_NAME = "traits.csv"
RAW_TRAIT_MISMATCH_NAME = "traits-mismatch.csv"
SEQUENCE_TYPE = "dna"
WORKFLOW_REQUIRED_TRAITS = ("body_mass_g", "gestation_days")
TREE_BRANCH_FLOOR = 1e-6


def load_catarrhine_data_quality_stress_panel_dataset() -> (
    CatarrhineDataQualityStressPanelDataset
):
    """Expose the packaged catarrhine stress panel as a first-class dataset surface."""
    dataset_root = _resource_root()
    raw_root = dataset_root / "raw"
    raw_alignment_path = raw_root / RAW_ALIGNMENT_NAME
    raw_sequence_input_path = raw_root / RAW_SEQUENCE_INPUT_NAME
    raw_coding_sequences_path = raw_root / RAW_CODING_SEQUENCE_NAME
    raw_tree_path = raw_root / RAW_TREE_NAME
    raw_traits_path = raw_root / RAW_TRAITS_NAME
    raw_trait_mismatch_path = raw_root / RAW_TRAIT_MISMATCH_NAME
    taxon_count = validate_fasta_input(
        raw_alignment_path, sequence_type=SEQUENCE_TYPE
    ).summary.sequence_count
    raw_trait_row_count = len(load_permissive_trait_rows(raw_traits_path))
    return CatarrhineDataQualityStressPanelDataset(
        dataset_id=DATASET_ID,
        label=DATASET_LABEL,
        dataset_root=dataset_root,
        readme_path=dataset_root / "README.md",
        raw_alignment_path=raw_alignment_path,
        raw_tree_path=raw_tree_path,
        raw_traits_path=raw_traits_path,
        raw_sequence_input_path=raw_sequence_input_path,
        raw_coding_sequences_path=raw_coding_sequences_path,
        raw_trait_mismatch_path=raw_trait_mismatch_path,
        reference_output_root=dataset_root / "expected",
        taxon_count=taxon_count,
        raw_trait_row_count=raw_trait_row_count,
        required_traits=WORKFLOW_REQUIRED_TRAITS,
        sequence_type=SEQUENCE_TYPE,
        source_summary=(
            "Real catarrhine mitochondrial sequence and topology material packaged "
            "with deliberate data-quality defects so duplicate FASTA identifiers, "
            "illegal or empty sequences, coding frame and stop-codon failures, "
            "trait duplication, tree-trait mismatch, sequence outliers, and "
            "branch-length pathologies remain visible and reviewable."
        ),
    )


def _resource_root() -> Path:
    return (
        Path(__file__).resolve().parent.parent.parent
        / "resources"
        / "datasets"
        / "stress"
        / DATASET_ID
    )
