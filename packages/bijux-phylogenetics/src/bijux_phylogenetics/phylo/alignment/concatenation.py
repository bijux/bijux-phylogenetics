from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import NamedTuple

from bijux_phylogenetics.phylo.alignment.models import (
    AlignmentAlphabet,
    AlignmentRecord,
)
from bijux_phylogenetics.phylo.alignment.occupancy import (
    LocusCoverageRow,
    LocusOccupancyCell,
    LocusOccupancyReport,
    TaxonCoverageRow,
)
from bijux_phylogenetics.phylo.alignment.partitions import (
    LocusPartition,
    LocusSegment,
    normalize_partition_data_type,
)
from bijux_phylogenetics.runtime.errors import (
    InvalidAlignmentError,
    InvalidPartitionError,
)

__all__ = [
    "ConcatenatedAlignmentLocusRow",
    "ConcatenatedAlignmentReport",
    "concatenate_locus_alignments",
]

_DEFAULT_MISSING_SYMBOL = "?"


@dataclass(frozen=True, slots=True)
class ConcatenatedAlignmentLocusRow:
    """One input locus summarized inside a concatenated supermatrix assembly."""

    locus_name: str
    alignment_path: Path
    inferred_alphabet: AlignmentAlphabet
    data_type: str | None
    taxon_count: int
    alignment_length: int
    missing_taxa: list[str]


@dataclass(slots=True)
class ConcatenatedAlignmentReport:
    """Assembly report for one concatenated multi-locus alignment."""

    input_paths: list[Path]
    taxon_count: int
    locus_count: int
    alignment_length: int
    taxa: list[str]
    loci: list[ConcatenatedAlignmentLocusRow]
    partitions: tuple[LocusPartition, ...]
    occupancy_report: LocusOccupancyReport
    warnings: list[str]


class _LoadedLocus(NamedTuple):
    name: str
    path: Path
    records: list[AlignmentRecord]
    inferred_alphabet: AlignmentAlphabet
    data_type: str | None
    occupancy_by_taxon: dict[str, tuple[int, float]]

    @property
    def alignment_length(self) -> int:
        return len(self.records[0].sequence)


def _load_alignment_records(path: Path) -> list[AlignmentRecord]:
    from bijux_phylogenetics.io.fasta.core import load_fasta_alignment

    return load_fasta_alignment(path)


def _infer_loaded_alignment_alphabet(
    records: list[AlignmentRecord],
) -> AlignmentAlphabet:
    from bijux_phylogenetics.io.fasta.core import infer_alignment_alphabet

    return infer_alignment_alphabet(records)


def _validated_locus_names(
    paths: list[Path],
    locus_names: tuple[str, ...] | None,
) -> list[str]:
    if locus_names is None:
        names = [path.stem for path in paths]
    else:
        if len(locus_names) != len(paths):
            raise ValueError("locus names must match the number of input alignments")
        names = [name.strip() for name in locus_names]
    if not all(names):
        raise InvalidPartitionError(
            "every concatenated locus must have a non-empty name"
        )
    duplicate_names = sorted({name for name in names if names.count(name) > 1})
    if duplicate_names:
        raise InvalidPartitionError(
            f"concatenated locus names must be unique, duplicated: {', '.join(duplicate_names)}"
        )
    return names


def _validated_data_types(
    paths: list[Path],
    data_types: tuple[str | None, ...] | None,
) -> list[str | None]:
    if data_types is None:
        return [None] * len(paths)
    if len(data_types) != len(paths):
        raise ValueError("data types must match the number of input alignments")
    return [
        None if data_type is None else normalize_partition_data_type(data_type)
        for data_type in data_types
    ]


def _validate_unique_taxa(path: Path, records: list[AlignmentRecord]) -> None:
    seen: set[str] = set()
    for record in records:
        if record.identifier in seen:
            raise InvalidAlignmentError(
                f"alignment contains duplicate taxon '{record.identifier}': {path}"
            )
        seen.add(record.identifier)


def _partition_data_type(
    *,
    inferred_alphabet: AlignmentAlphabet,
    declared_data_type: str | None,
) -> str | None:
    if declared_data_type is not None:
        return declared_data_type
    if inferred_alphabet == "unknown":
        return None
    return normalize_partition_data_type(inferred_alphabet)


def _load_locus(
    *,
    path: Path,
    locus_name: str,
    declared_data_type: str | None,
) -> _LoadedLocus:
    records = _load_alignment_records(path)
    if not records:
        raise InvalidAlignmentError(f"alignment does not contain any records: {path}")
    _validate_unique_taxa(path, records)

    inferred_alphabet = _infer_loaded_alignment_alphabet(records)
    occupancy_by_taxon = {
        record.identifier: (len(record.sequence), 1.0) for record in records
    }
    return _LoadedLocus(
        name=locus_name,
        path=path,
        records=records,
        inferred_alphabet=inferred_alphabet,
        data_type=_partition_data_type(
            inferred_alphabet=inferred_alphabet,
            declared_data_type=declared_data_type,
        ),
        occupancy_by_taxon=occupancy_by_taxon,
    )


def concatenate_locus_alignments(
    alignment_paths: list[Path],
    *,
    locus_names: tuple[str, ...] | None = None,
    data_types: tuple[str | None, ...] | None = None,
    concatenated_alignment_path: Path | None = None,
    concatenated_partition_path: Path | None = None,
) -> tuple[
    list[AlignmentRecord], tuple[LocusPartition, ...], ConcatenatedAlignmentReport
]:
    """Concatenate aligned per-locus FASTA inputs into one partitioned supermatrix."""
    if not alignment_paths:
        raise ValueError("at least one aligned locus is required for concatenation")

    normalized_paths = [Path(path) for path in alignment_paths]
    resolved_locus_names = _validated_locus_names(normalized_paths, locus_names)
    resolved_data_types = _validated_data_types(normalized_paths, data_types)

    loaded_loci = [
        _load_locus(
            path=path,
            locus_name=locus_name,
            declared_data_type=declared_data_type,
        )
        for path, locus_name, declared_data_type in zip(
            normalized_paths,
            resolved_locus_names,
            resolved_data_types,
            strict=True,
        )
    ]

    ordered_taxa: list[str] = []
    seen_taxa: set[str] = set()
    for locus in loaded_loci:
        for record in locus.records:
            if record.identifier in seen_taxa:
                continue
            seen_taxa.add(record.identifier)
            ordered_taxa.append(record.identifier)

    partitions: list[LocusPartition] = []
    concatenated_records: list[AlignmentRecord] = []
    occupancy_cells: list[LocusOccupancyCell] = []
    occupancy_by_taxon: dict[str, dict[str, float]] = {
        taxon: {} for taxon in ordered_taxa
    }
    presence_by_taxon: dict[str, dict[str, bool]] = {
        taxon: {} for taxon in ordered_taxa
    }
    taxon_site_totals: dict[str, int] = dict.fromkeys(ordered_taxa, 0)
    locus_rows: list[LocusCoverageRow] = []
    locus_summaries: list[ConcatenatedAlignmentLocusRow] = []
    warnings: list[str] = []

    cursor = 1
    sequence_by_locus_and_taxon = {
        locus.name: {record.identifier: record.sequence for record in locus.records}
        for locus in loaded_loci
    }
    for locus in loaded_loci:
        locus_end = cursor + locus.alignment_length - 1
        partitions.append(
            LocusPartition(
                name=locus.name,
                segments=(LocusSegment(start=cursor, end=locus_end),),
                total_sites=locus.alignment_length,
                data_type=locus.data_type,
            )
        )
        cursor = locus_end + 1

        covered_taxon_count = 0
        locus_observed_site_count = 0
        missing_taxa: list[str] = []
        for taxon in ordered_taxa:
            observed_site_count, occupancy_fraction = locus.occupancy_by_taxon.get(
                taxon,
                (0, 0.0),
            )
            present = observed_site_count > 0
            if present:
                covered_taxon_count += 1
            else:
                missing_taxa.append(taxon)
            locus_observed_site_count += observed_site_count
            taxon_site_totals[taxon] += observed_site_count
            occupancy_by_taxon[taxon][locus.name] = occupancy_fraction
            presence_by_taxon[taxon][locus.name] = present
            occupancy_cells.append(
                LocusOccupancyCell(
                    taxon=taxon,
                    locus_name=locus.name,
                    observed_site_count=observed_site_count,
                    total_site_count=locus.alignment_length,
                    occupancy_fraction=occupancy_fraction,
                    present=present,
                )
            )

        locus_rows.append(
            LocusCoverageRow(
                locus_name=locus.name,
                covered_taxon_count=covered_taxon_count,
                total_taxa=len(ordered_taxa),
                taxon_coverage_fraction=covered_taxon_count / len(ordered_taxa),
                observed_site_count=locus_observed_site_count,
                total_site_count=locus.alignment_length * len(ordered_taxa),
                site_coverage_fraction=(
                    locus_observed_site_count
                    / (locus.alignment_length * len(ordered_taxa))
                ),
                low_coverage=False,
            )
        )
        locus_summaries.append(
            ConcatenatedAlignmentLocusRow(
                locus_name=locus.name,
                alignment_path=locus.path,
                inferred_alphabet=locus.inferred_alphabet,
                data_type=locus.data_type,
                taxon_count=len(locus.records),
                alignment_length=locus.alignment_length,
                missing_taxa=missing_taxa,
            )
        )
        if locus.data_type is None:
            warnings.append(
                f"locus '{locus.name}' uses an unknown alphabet and omits a declared partition data type"
            )

    total_alignment_length = sum(partition.total_sites for partition in partitions)
    for taxon in ordered_taxa:
        chunks: list[str] = []
        for locus in loaded_loci:
            chunks.append(
                sequence_by_locus_and_taxon[locus.name].get(
                    taxon,
                    _DEFAULT_MISSING_SYMBOL * locus.alignment_length,
                )
            )
        concatenated_records.append(
            AlignmentRecord(identifier=taxon, sequence="".join(chunks))
        )

    taxon_rows = [
        TaxonCoverageRow(
            taxon=taxon,
            covered_locus_count=sum(presence_by_taxon[taxon].values()),
            total_locus_count=len(partitions),
            locus_coverage_fraction=sum(presence_by_taxon[taxon].values())
            / len(partitions),
            observed_site_count=taxon_site_totals[taxon],
            total_site_count=total_alignment_length,
            site_coverage_fraction=taxon_site_totals[taxon] / total_alignment_length,
            low_coverage=False,
            occupancies=occupancy_by_taxon[taxon],
        )
        for taxon in ordered_taxa
    ]

    occupancy_report = LocusOccupancyReport(
        alignment_path=(
            Path("<concatenated-alignment>")
            if concatenated_alignment_path is None
            else concatenated_alignment_path
        ),
        partition_path=(
            Path("<concatenated-partitions>")
            if concatenated_partition_path is None
            else concatenated_partition_path
        ),
        taxon_count=len(ordered_taxa),
        locus_count=len(partitions),
        alignment_length=total_alignment_length,
        assigned_site_count=total_alignment_length,
        unassigned_site_count=0,
        partitions=tuple(partitions),
        cells=occupancy_cells,
        taxa=taxon_rows,
        loci=locus_rows,
        low_coverage_taxa=[],
        low_coverage_loci=[],
        taxon_coverage_threshold=None,
        locus_coverage_threshold=None,
        minimum_locus_occupancy=0.0,
        warnings=list(dict.fromkeys(warnings)),
    )

    return (
        concatenated_records,
        tuple(partitions),
        ConcatenatedAlignmentReport(
            input_paths=normalized_paths,
            taxon_count=len(ordered_taxa),
            locus_count=len(partitions),
            alignment_length=total_alignment_length,
            taxa=ordered_taxa,
            loci=locus_summaries,
            partitions=tuple(partitions),
            occupancy_report=occupancy_report,
            warnings=list(dict.fromkeys(warnings)),
        ),
    )
