from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.diagnostics.validation import (
    TreeInspectionReport,
    TreeValidationReport,
    inspect_tree_path,
    validate_tree_path,
)
from bijux_phylogenetics.io.fasta import AlignmentRecord, write_fasta_alignment
from bijux_phylogenetics.io.fasta.cleaning import detect_composition_outlier_sequences
from bijux_phylogenetics.io.fasta.coding import prepare_coding_sequences_for_alignment
from bijux_phylogenetics.io.fasta.records import (
    FastaInputValidationReport,
    FastaRepairReport,
    SequenceCompositionOutlier,
    SequenceLengthOutlier,
    detect_sequence_length_outliers,
    repair_fasta_input,
    validate_fasta_input,
)
from bijux_phylogenetics.phylo.alignment import CodingSequencePreparationReport

from ..models import CatarrhineDataQualityStressPanelDataset


@dataclass(slots=True)
class SequenceCleanupSurfaces:
    sequence_type: str
    raw_sequence_input_validation: FastaInputValidationReport
    raw_sequence_input_repair: FastaRepairReport
    raw_sequence_length_outliers: list[SequenceLengthOutlier]
    repaired_sequence_input_validation: FastaInputValidationReport
    coding_sequence_preparation: CodingSequencePreparationReport
    raw_alignment_validation: FastaInputValidationReport
    sequence_outliers: list[SequenceCompositionOutlier]
    raw_tree_inspection: TreeInspectionReport
    raw_tree_validation: TreeValidationReport
    repaired_sequence_input_path: Path
    prepared_coding_sequences_path: Path
    sequence_outlier_taxa: list[str]
    tree_outlier_taxa: list[str]

    def validate_cleaned_alignment(
        self, cleaned_alignment_path: Path
    ) -> FastaInputValidationReport:
        return validate_fasta_input(
            cleaned_alignment_path,
            sequence_type=self.sequence_type,
        )

    def validate_cleaned_tree(self, cleaned_tree_path: Path) -> TreeValidationReport:
        return validate_tree_path(cleaned_tree_path)


def prepare_sequence_surfaces(
    *,
    dataset: CatarrhineDataQualityStressPanelDataset,
    assembled_root: Path,
) -> SequenceCleanupSurfaces:
    raw_sequence_input_validation = validate_fasta_input(
        dataset.raw_sequence_input_path,
        sequence_type=dataset.sequence_type,
    )
    repaired_sequence_records, raw_sequence_input_repair = repair_fasta_input(
        dataset.raw_sequence_input_path,
        sequence_type=dataset.sequence_type,
        normalize_identifiers=True,
        remove_invalid_records=True,
    )
    provisional_repaired_sequence_path = write_fasta_alignment(
        assembled_root / "repaired-sequence-input.provisional.fasta",
        repaired_sequence_records,
    )
    raw_sequence_length_outliers = detect_sequence_length_outliers(
        provisional_repaired_sequence_path
    )
    repaired_sequence_input_path = write_fasta_alignment(
        assembled_root / "repaired-sequence-input.fasta",
        _retained_repaired_sequence_records(
            repaired_sequence_records,
            raw_sequence_length_outliers,
        ),
    )
    repaired_sequence_input_validation = validate_fasta_input(
        repaired_sequence_input_path,
        sequence_type=dataset.sequence_type,
    )
    prepared_coding_sequences, coding_sequence_preparation = (
        prepare_coding_sequences_for_alignment(
            dataset.raw_coding_sequences_path,
            sequence_type=dataset.sequence_type,
        )
    )
    prepared_coding_sequences_path = write_fasta_alignment(
        assembled_root / "prepared-coding-sequences.fasta",
        prepared_coding_sequences,
    )
    raw_alignment_validation = validate_fasta_input(
        dataset.raw_alignment_path,
        sequence_type=dataset.sequence_type,
    )
    sequence_outliers = detect_composition_outlier_sequences(dataset.raw_alignment_path)
    raw_tree_inspection = inspect_tree_path(dataset.raw_tree_path)
    raw_tree_validation = validate_tree_path(
        dataset.raw_tree_path,
        allow_negative_branch_lengths=True,
    )
    return SequenceCleanupSurfaces(
        sequence_type=dataset.sequence_type,
        raw_sequence_input_validation=raw_sequence_input_validation,
        raw_sequence_input_repair=raw_sequence_input_repair,
        raw_sequence_length_outliers=raw_sequence_length_outliers,
        repaired_sequence_input_validation=repaired_sequence_input_validation,
        coding_sequence_preparation=coding_sequence_preparation,
        raw_alignment_validation=raw_alignment_validation,
        sequence_outliers=sequence_outliers,
        raw_tree_inspection=raw_tree_inspection,
        raw_tree_validation=raw_tree_validation,
        repaired_sequence_input_path=repaired_sequence_input_path,
        prepared_coding_sequences_path=prepared_coding_sequences_path,
        sequence_outlier_taxa=sorted(
            {row.identifier for row in sequence_outliers if row.identifier}
        ),
        tree_outlier_taxa=sorted(raw_tree_inspection.long_branch_taxa),
    )


def _retained_repaired_sequence_records(
    repaired_sequence_records: list[AlignmentRecord],
    raw_sequence_length_outliers: list[SequenceLengthOutlier],
) -> list[AlignmentRecord]:
    excluded_identifiers = {row.identifier for row in raw_sequence_length_outliers}
    return [
        record
        for record in repaired_sequence_records
        if record.identifier not in excluded_identifiers
    ]
