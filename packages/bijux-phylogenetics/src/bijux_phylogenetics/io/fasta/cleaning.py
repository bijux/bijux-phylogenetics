from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import load_taxon_table
from bijux_phylogenetics.phylo.alignment import (
    AlignmentCleaningReport,
    AlignmentComparisonReport,
    AlignmentCompositionShift,
    AlignmentFilterProfile,
    AlignmentGroupRetention,
    AlignmentRecord,
    AlignmentSignalWarning,
    AlignmentSummary,
    AlignmentTrimReport,
    DuplicateSequenceGroup,
    NearDuplicateSequencePair,
    PairwiseSequenceIdentity,
    RemovedAlignmentSequence,
    SequenceCompositionOutlier,
    SequenceIdentityMatrix,
    TrimmedAlignmentColumn,
)

from .core.alignment_io import load_fasta_alignment
from .core.character_policy import (
    _GAP_CHARACTERS,
    _is_explicit_missing,
    _validate_fraction_threshold,
)
from .records import (
    _detect_composition_outlier_sequences_records,
    _detect_identical_duplicate_sequences_records,
    _detect_near_duplicate_sequences_records,
    _pairwise_identity,
    detect_sequences_with_excessive_missing_data,
    detect_sites_with_excessive_missing_data,
    summarise_fasta,
    summarise_records_as_alignment_summary,
)

_ALIGNMENT_FILTER_PROFILES: dict[str, AlignmentFilterProfile] = {
    "conservative": AlignmentFilterProfile(
        name="conservative",
        remove_all_gap_sites=True,
        remove_all_missing_sites=True,
        site_missingness_threshold=0.8,
        sequence_missingness_threshold=0.6,
        preserve_codon_structure=False,
        note="retain most data while removing only the most uncertainty-heavy sites and sequences",
    ),
    "moderate": AlignmentFilterProfile(
        name="moderate",
        remove_all_gap_sites=True,
        remove_all_missing_sites=True,
        site_missingness_threshold=0.5,
        sequence_missingness_threshold=0.4,
        preserve_codon_structure=False,
        note="balanced trimming for routine phylogenetic alignment cleaning",
    ),
    "aggressive": AlignmentFilterProfile(
        name="aggressive",
        remove_all_gap_sites=True,
        remove_all_missing_sites=True,
        site_missingness_threshold=0.3,
        sequence_missingness_threshold=0.25,
        preserve_codon_structure=False,
        note="remove strongly uncertainty-heavy sites and weak sequences even at the cost of data loss",
    ),
    "coding-safe": AlignmentFilterProfile(
        name="coding-safe",
        remove_all_gap_sites=True,
        remove_all_missing_sites=True,
        site_missingness_threshold=0.5,
        sequence_missingness_threshold=0.35,
        preserve_codon_structure=True,
        note="trim codon-aligned nucleotide data while preserving codon phase by removing full codons together",
    ),
    "phylogenomics-scale": AlignmentFilterProfile(
        name="phylogenomics-scale",
        remove_all_gap_sites=True,
        remove_all_missing_sites=True,
        site_missingness_threshold=0.7,
        sequence_missingness_threshold=0.7,
        preserve_codon_structure=False,
        note="favor taxon retention for large concatenated matrices while removing empty or highly degraded regions",
    ),
}


def _trim_columns(
    records: list[AlignmentRecord],
    *,
    keep_positions: list[int],
) -> list[AlignmentRecord]:
    return [
        AlignmentRecord(
            identifier=record.identifier,
            sequence="".join(record.sequence[index] for index in keep_positions),
        )
        for record in records
    ]


def _expand_removed_positions_to_groups(
    removed_positions: set[int],
    *,
    alignment_length: int,
    group_size: int,
) -> set[int]:
    if group_size <= 1:
        return set(removed_positions)
    expanded: set[int] = set()
    for position in removed_positions:
        group_start = ((position - 1) // group_size) * group_size + 1
        for grouped_position in range(
            group_start,
            min(group_start + group_size, alignment_length + 1),
        ):
            expanded.add(grouped_position)
    return expanded


def _composition_for_comparison(summary: AlignmentSummary) -> dict[str, float]:
    if summary.nucleotide_composition:
        return summary.nucleotide_composition
    return summary.amino_acid_composition


def _signal_warnings_for_cleaning(
    original: AlignmentSummary,
    cleaned: AlignmentSummary,
) -> list[AlignmentSignalWarning]:
    warnings: list[AlignmentSignalWarning] = []
    if original.parsimony_informative_site_count > 0:
        retained_fraction = (
            cleaned.parsimony_informative_site_count
            / original.parsimony_informative_site_count
        )
        if retained_fraction < 0.5:
            warnings.append(
                AlignmentSignalWarning(
                    code="informative-site-loss",
                    message="cleaning retained less than half of the original parsimony-informative sites",
                )
            )
    if cleaned.sequence_count < max(2, original.sequence_count // 2):
        warnings.append(
            AlignmentSignalWarning(
                code="taxon-coverage-loss",
                message="cleaning removed more than half of the original taxa",
            )
        )
    if cleaned.variable_site_count == 0 and original.variable_site_count > 0:
        warnings.append(
            AlignmentSignalWarning(
                code="variable-site-collapse",
                message="cleaning removed all variable sites from the alignment",
            )
        )
    return warnings


def _group_retention_after_cleaning(
    original_taxa: list[str],
    retained_taxa: list[str],
    *,
    table_path: Path | None,
    group_columns: list[str] | None,
) -> list[AlignmentGroupRetention]:
    if table_path is None or not group_columns:
        return []
    table = load_taxon_table(table_path)
    rows_by_taxon = {row[table.taxon_column]: row for row in table.rows}
    original_taxa_set = set(original_taxa)
    retained_taxa_set = set(retained_taxa)
    reports: list[AlignmentGroupRetention] = []
    for column in group_columns:
        if column not in table.columns:
            raise ValueError(f"group column '{column}' is not present in {table_path}")
        values = sorted(
            {
                rows_by_taxon[taxon][column]
                for taxon in original_taxa
                if taxon in rows_by_taxon and rows_by_taxon[taxon][column]
            }
        )
        for value in values:
            original_count = sum(
                1
                for taxon in original_taxa_set
                if taxon in rows_by_taxon and rows_by_taxon[taxon][column] == value
            )
            retained_count = sum(
                1
                for taxon in retained_taxa_set
                if taxon in rows_by_taxon and rows_by_taxon[taxon][column] == value
            )
            removed_count = original_count - retained_count
            reports.append(
                AlignmentGroupRetention(
                    column=column,
                    value=value,
                    original_count=original_count,
                    retained_count=retained_count,
                    removed_count=removed_count,
                    removed_fraction=0.0
                    if original_count == 0
                    else round(removed_count / original_count, 15),
                )
            )
    return reports


def list_alignment_filter_profiles() -> list[AlignmentFilterProfile]:
    """Return the supported named alignment-cleaning profiles."""
    return list(_ALIGNMENT_FILTER_PROFILES.values())


def get_alignment_filter_profile(name: str) -> AlignmentFilterProfile:
    """Resolve one named alignment-cleaning profile."""
    try:
        return _ALIGNMENT_FILTER_PROFILES[name]
    except KeyError as error:
        available = ", ".join(sorted(_ALIGNMENT_FILTER_PROFILES))
        raise ValueError(
            f"unknown alignment filter profile '{name}', expected one of: {available}"
        ) from error


def detect_composition_outlier_sequences(
    path: Path,
    *,
    deviation_threshold: float = 0.25,
) -> list[SequenceCompositionOutlier]:
    """Detect sequences with unusually deviant GC or amino-acid composition."""
    records = load_fasta_alignment(path)
    return _detect_composition_outlier_sequences_records(
        records,
        deviation_threshold=deviation_threshold,
    )


def detect_identical_duplicate_sequences(path: Path) -> list[DuplicateSequenceGroup]:
    """Group sequences that are exactly identical over the full aligned string."""
    records = load_fasta_alignment(path)
    return _detect_identical_duplicate_sequences_records(records)


def detect_near_duplicate_sequences(
    path: Path,
    *,
    identity_threshold: float,
) -> list[NearDuplicateSequencePair]:
    """Return sequence pairs above the given identity threshold."""
    _validate_fraction_threshold(identity_threshold)
    records = load_fasta_alignment(path)
    return _detect_near_duplicate_sequences_records(
        records,
        identity_threshold=identity_threshold,
    )


def remove_all_gap_columns(
    path: Path,
) -> tuple[list[AlignmentRecord], AlignmentTrimReport]:
    """Remove columns composed entirely of gap characters."""
    records = load_fasta_alignment(path)
    summary = summarise_fasta(path)
    removed_positions = set(summary.all_gap_columns)
    keep_positions = [
        index
        for index in range(summary.alignment_length)
        if (index + 1) not in removed_positions
    ]
    trimmed_records = _trim_columns(records, keep_positions=keep_positions)
    return trimmed_records, AlignmentTrimReport(
        path=path,
        original_sequence_count=summary.sequence_count,
        trimmed_sequence_count=len(trimmed_records),
        original_alignment_length=summary.alignment_length,
        trimmed_alignment_length=len(keep_positions),
        removed_columns=[
            TrimmedAlignmentColumn(position=position, reason="all-gap")
            for position in summary.all_gap_columns
        ],
        removed_sequences=[],
    )


def remove_all_missing_columns(
    path: Path,
) -> tuple[list[AlignmentRecord], AlignmentTrimReport]:
    """Remove columns composed entirely of missing-data symbols."""
    records = load_fasta_alignment(path)
    summary = summarise_fasta(path)
    removed_positions = set(summary.all_missing_columns)
    keep_positions = [
        index
        for index in range(summary.alignment_length)
        if (index + 1) not in removed_positions
    ]
    trimmed_records = _trim_columns(records, keep_positions=keep_positions)
    return trimmed_records, AlignmentTrimReport(
        path=path,
        original_sequence_count=summary.sequence_count,
        trimmed_sequence_count=len(trimmed_records),
        original_alignment_length=summary.alignment_length,
        trimmed_alignment_length=len(keep_positions),
        removed_columns=[
            TrimmedAlignmentColumn(position=position, reason="all-missing")
            for position in summary.all_missing_columns
        ],
        removed_sequences=[],
    )


def trim_columns_above_missingness_threshold(
    path: Path,
    *,
    threshold: float,
) -> tuple[list[AlignmentRecord], AlignmentTrimReport]:
    """Remove columns whose missing-data fraction exceeds the given threshold."""
    _validate_fraction_threshold(threshold)
    records = load_fasta_alignment(path)
    summary = summarise_fasta(path)
    removed_positions = {
        row.position
        for row in detect_sites_with_excessive_missing_data(path, threshold=threshold)
    }
    keep_positions = [
        index
        for index in range(summary.alignment_length)
        if (index + 1) not in removed_positions
    ]
    trimmed_records = _trim_columns(records, keep_positions=keep_positions)
    return trimmed_records, AlignmentTrimReport(
        path=path,
        original_sequence_count=summary.sequence_count,
        trimmed_sequence_count=len(trimmed_records),
        original_alignment_length=summary.alignment_length,
        trimmed_alignment_length=len(keep_positions),
        removed_columns=[
            TrimmedAlignmentColumn(position=position, reason="missingness-threshold")
            for position in sorted(removed_positions)
        ],
        removed_sequences=[],
    )


def remove_sequences_above_missingness_threshold(
    path: Path,
    *,
    threshold: float,
) -> tuple[list[AlignmentRecord], AlignmentTrimReport]:
    """Remove sequences whose missing-data fraction exceeds the given threshold."""
    _validate_fraction_threshold(threshold)
    records = load_fasta_alignment(path)
    summary = summarise_fasta(path)
    excessive = {
        row.identifier: row.missing_fraction
        for row in detect_sequences_with_excessive_missing_data(
            path, threshold=threshold
        )
    }
    trimmed_records = [
        record for record in records if record.identifier not in excessive
    ]
    removed_sequences = [
        RemovedAlignmentSequence(
            identifier=record.identifier,
            missing_fraction=excessive[record.identifier],
            reason="missingness-threshold",
        )
        for record in records
        if record.identifier in excessive
    ]
    return trimmed_records, AlignmentTrimReport(
        path=path,
        original_sequence_count=summary.sequence_count,
        trimmed_sequence_count=len(trimmed_records),
        original_alignment_length=summary.alignment_length,
        trimmed_alignment_length=summary.alignment_length,
        removed_columns=[],
        removed_sequences=removed_sequences,
    )


def trim_alignment(
    path: Path,
    *,
    remove_all_gap_sites: bool = True,
    remove_all_missing_sites: bool = True,
    site_missingness_threshold: float | None = None,
    sequence_missingness_threshold: float | None = None,
    preserve_codon_structure: bool = False,
) -> tuple[list[AlignmentRecord], AlignmentTrimReport]:
    """Trim an alignment with explicit deterministic transform reporting."""
    if site_missingness_threshold is not None:
        _validate_fraction_threshold(site_missingness_threshold)
    if sequence_missingness_threshold is not None:
        _validate_fraction_threshold(sequence_missingness_threshold)

    original_records = load_fasta_alignment(path)
    summary = summarise_fasta(path)
    records = list(original_records)
    removed_columns: list[TrimmedAlignmentColumn] = []
    removed_sequences: list[RemovedAlignmentSequence] = []

    if sequence_missingness_threshold is not None:
        excessive = {
            row.identifier: row.missing_fraction
            for row in detect_sequences_with_excessive_missing_data(
                path,
                threshold=sequence_missingness_threshold,
            )
        }
        removed_sequences.extend(
            RemovedAlignmentSequence(
                identifier=record.identifier,
                missing_fraction=excessive[record.identifier],
                reason="missingness-threshold",
            )
            for record in records
            if record.identifier in excessive
        )
        records = [record for record in records if record.identifier not in excessive]

    removed_positions: set[int] = set()
    if remove_all_gap_sites:
        removed_positions.update(summary.all_gap_columns)
        removed_columns.extend(
            TrimmedAlignmentColumn(position=position, reason="all-gap")
            for position in summary.all_gap_columns
        )
    if remove_all_missing_sites:
        removed_positions.update(
            position
            for position in summary.all_missing_columns
            if position not in removed_positions
        )
        removed_columns.extend(
            TrimmedAlignmentColumn(position=position, reason="all-missing")
            for position in summary.all_missing_columns
            if position not in summary.all_gap_columns or not remove_all_gap_sites
        )
    if site_missingness_threshold is not None:
        excessive_sites = [
            row.position
            for row in detect_sites_with_excessive_missing_data(
                path,
                threshold=site_missingness_threshold,
            )
            if row.position not in removed_positions
        ]
        removed_positions.update(excessive_sites)
        removed_columns.extend(
            TrimmedAlignmentColumn(position=position, reason="missingness-threshold")
            for position in excessive_sites
        )

    if preserve_codon_structure and summary.inferred_alphabet in {"dna", "rna"}:
        expanded_positions = _expand_removed_positions_to_groups(
            removed_positions,
            alignment_length=summary.alignment_length,
            group_size=3,
        )
        codon_phase_positions = sorted(expanded_positions - removed_positions)
        removed_positions = expanded_positions
        removed_columns.extend(
            TrimmedAlignmentColumn(position=position, reason="codon-phase-preservation")
            for position in codon_phase_positions
        )

    keep_positions = [
        index
        for index in range(summary.alignment_length)
        if (index + 1) not in removed_positions
    ]
    trimmed_records = _trim_columns(records, keep_positions=keep_positions)
    removed_columns.sort(key=lambda item: item.position)
    removed_sequences.sort(key=lambda item: item.identifier)
    return trimmed_records, AlignmentTrimReport(
        path=path,
        original_sequence_count=summary.sequence_count,
        trimmed_sequence_count=len(trimmed_records),
        original_alignment_length=summary.alignment_length,
        trimmed_alignment_length=len(keep_positions),
        removed_columns=removed_columns,
        removed_sequences=removed_sequences,
    )


def compare_alignment_versions(
    left_path: Path, right_path: Path
) -> AlignmentComparisonReport:
    """Compare two alignments over taxa, sites, uncertainty burden, and composition."""
    left = summarise_fasta(left_path)
    right = summarise_fasta(right_path)
    composition_shifts: list[AlignmentCompositionShift] = []
    left_composition = _composition_for_comparison(left)
    right_composition = _composition_for_comparison(right)
    for component in sorted(set(left_composition) | set(right_composition)):
        before = left_composition.get(component, 0.0)
        after = right_composition.get(component, 0.0)
        composition_shifts.append(
            AlignmentCompositionShift(
                component=component,
                before=before,
                after=after,
                delta=round(after - before, 15),
            )
        )

    warnings: list[str] = []
    if left.inferred_alphabet != right.inferred_alphabet:
        warnings.append("alignment versions infer different alphabets")
    if set(left.ids) != set(right.ids):
        warnings.append("alignment versions contain different taxon sets")

    return AlignmentComparisonReport(
        left_path=left_path,
        right_path=right_path,
        shared_taxa=sorted(set(left.ids) & set(right.ids)),
        left_only_taxa=sorted(set(left.ids) - set(right.ids)),
        right_only_taxa=sorted(set(right.ids) - set(left.ids)),
        left_alignment_length=left.alignment_length,
        right_alignment_length=right.alignment_length,
        left_missing_data_fraction=left.missing_data_fraction,
        right_missing_data_fraction=right.missing_data_fraction,
        left_gap_fraction=left.gap_fraction,
        right_gap_fraction=right.gap_fraction,
        left_variable_site_count=left.variable_site_count,
        right_variable_site_count=right.variable_site_count,
        left_parsimony_informative_site_count=left.parsimony_informative_site_count,
        right_parsimony_informative_site_count=right.parsimony_informative_site_count,
        composition_shifts=composition_shifts,
        warnings=warnings,
    )


def compare_alignment_summaries(
    left_path: Path,
    left: AlignmentSummary,
    right: AlignmentSummary,
) -> AlignmentComparisonReport:
    """Compare two already-computed alignment summaries."""
    composition_shifts: list[AlignmentCompositionShift] = []
    left_composition = _composition_for_comparison(left)
    right_composition = _composition_for_comparison(right)
    for component in sorted(set(left_composition) | set(right_composition)):
        before = left_composition.get(component, 0.0)
        after = right_composition.get(component, 0.0)
        composition_shifts.append(
            AlignmentCompositionShift(
                component=component,
                before=before,
                after=after,
                delta=round(after - before, 15),
            )
        )
    warnings: list[str] = []
    if left.inferred_alphabet != right.inferred_alphabet:
        warnings.append("alignment versions infer different alphabets")
    if set(left.ids) != set(right.ids):
        warnings.append("alignment versions contain different taxon sets")
    return AlignmentComparisonReport(
        left_path=left_path,
        right_path=right.path,
        shared_taxa=sorted(set(left.ids) & set(right.ids)),
        left_only_taxa=sorted(set(left.ids) - set(right.ids)),
        right_only_taxa=sorted(set(right.ids) - set(left.ids)),
        left_alignment_length=left.alignment_length,
        right_alignment_length=right.alignment_length,
        left_missing_data_fraction=left.missing_data_fraction,
        right_missing_data_fraction=right.missing_data_fraction,
        left_gap_fraction=left.gap_fraction,
        right_gap_fraction=right.gap_fraction,
        left_variable_site_count=left.variable_site_count,
        right_variable_site_count=right.variable_site_count,
        left_parsimony_informative_site_count=left.parsimony_informative_site_count,
        right_parsimony_informative_site_count=right.parsimony_informative_site_count,
        composition_shifts=composition_shifts,
        warnings=warnings,
    )


def clean_alignment_with_profile(
    path: Path,
    *,
    profile_name: str,
    group_table_path: Path | None = None,
    group_columns: list[str] | None = None,
) -> tuple[list[AlignmentRecord], AlignmentCleaningReport]:
    """Clean an alignment using one named profile and compare the result against the original."""
    profile = get_alignment_filter_profile(profile_name)
    cleaned_records, trim_report = trim_alignment(
        path,
        remove_all_gap_sites=profile.remove_all_gap_sites,
        remove_all_missing_sites=profile.remove_all_missing_sites,
        site_missingness_threshold=profile.site_missingness_threshold,
        sequence_missingness_threshold=profile.sequence_missingness_threshold,
        preserve_codon_structure=profile.preserve_codon_structure,
    )

    original = summarise_fasta(path)
    cleaned_summary = summarise_records_as_alignment_summary(
        path=path, records=cleaned_records
    )
    comparison = compare_alignment_summaries(path, original, cleaned_summary)
    signal_warnings = _signal_warnings_for_cleaning(original, cleaned_summary)
    group_retention = _group_retention_after_cleaning(
        original.ids,
        cleaned_summary.ids,
        table_path=group_table_path,
        group_columns=group_columns,
    )
    warnings = [warning.message for warning in signal_warnings]
    if any(
        group.removed_fraction >= 0.5 and group.original_count >= 1
        for group in group_retention
    ):
        warnings.append(
            "cleaning removed most taxa from one or more metadata or trait groups"
        )
    return cleaned_records, AlignmentCleaningReport(
        profile=profile,
        trim=trim_report,
        comparison=comparison,
        signal_warnings=signal_warnings,
        group_retention=group_retention,
        warnings=warnings,
    )


def compute_pairwise_sequence_identity_matrix(path: Path) -> SequenceIdentityMatrix:
    """Compute a deterministic pairwise sequence identity matrix."""
    records = load_fasta_alignment(path)
    pairs: list[PairwiseSequenceIdentity] = []
    for left_index, left in enumerate(records):
        for right_index, right in enumerate(records):
            if right_index < left_index:
                continue
            if left_index == right_index:
                pairs.append(
                    PairwiseSequenceIdentity(
                        left_identifier=left.identifier,
                        right_identifier=right.identifier,
                        identity=1.0,
                        comparable_sites=len(
                            [
                                residue
                                for residue in left.sequence
                                if residue not in _GAP_CHARACTERS
                                and not _is_explicit_missing(residue)
                            ]
                        ),
                    )
                )
                continue
            identity, comparable_sites = _pairwise_identity(
                left.sequence, right.sequence
            )
            pairs.append(
                PairwiseSequenceIdentity(
                    left_identifier=left.identifier,
                    right_identifier=right.identifier,
                    identity=identity if comparable_sites > 0 else None,
                    comparable_sites=comparable_sites,
                )
            )
    return SequenceIdentityMatrix(
        path=path,
        identifiers=[record.identifier for record in records],
        pairs=pairs,
    )


def write_sequence_identity_matrix(path: Path, report: SequenceIdentityMatrix) -> Path:
    """Write a pairwise sequence identity matrix as a deterministic TSV."""
    rows = {
        (pair.left_identifier, pair.right_identifier): pair for pair in report.pairs
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["left_identifier\tright_identifier\tidentity\tcomparable_sites"]
    for left in report.identifiers:
        for right in report.identifiers:
            pair = rows.get((left, right)) or rows.get((right, left))
            if pair is None:
                continue
            identity = "" if pair.identity is None else format(pair.identity, ".15g")
            lines.append(f"{left}\t{right}\t{identity}\t{pair.comparable_sites}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
