from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re

from bijux_phylogenetics.core.alignment import AlignmentRecord
from bijux_phylogenetics.errors import InvalidAlignmentError, InvalidPartitionError
from bijux_phylogenetics.io.fasta import (
    load_fasta_alignment,
    summarise_records_as_alignment_summary,
)

__all__ = [
    "LocusCoverageRow",
    "LocusOccupancyCell",
    "LocusOccupancyReport",
    "LocusPartition",
    "LocusSegment",
    "TaxonCoverageRow",
    "build_locus_occupancy_report",
    "parse_locus_partitions",
]

_PARTITION_LINE = re.compile(
    r"^(?:(?P<data_type>[A-Za-z0-9_-]+)\s*,\s*)?"
    r"(?P<name>[^=]+?)\s*=\s*(?P<ranges>[^;]+?)\s*;?$"
)


@dataclass(frozen=True, slots=True)
class LocusSegment:
    """One inclusive 1-based coordinate block inside a locus partition."""

    start: int
    end: int

    @property
    def length(self) -> int:
        """Return the number of sites covered by this segment."""
        return self.end - self.start + 1


@dataclass(frozen=True, slots=True)
class LocusPartition:
    """One named locus assembled from one or more concatenated segments."""

    name: str
    segments: tuple[LocusSegment, ...]
    total_sites: int
    data_type: str | None = None


@dataclass(frozen=True, slots=True)
class LocusOccupancyCell:
    """Coverage state for one taxon within one named locus."""

    taxon: str
    locus_name: str
    observed_site_count: int
    total_site_count: int
    occupancy_fraction: float
    present: bool


@dataclass(frozen=True, slots=True)
class TaxonCoverageRow:
    """Per-taxon coverage across the locus set."""

    taxon: str
    covered_locus_count: int
    total_locus_count: int
    locus_coverage_fraction: float
    observed_site_count: int
    total_site_count: int
    low_coverage: bool
    occupancies: dict[str, float]


@dataclass(frozen=True, slots=True)
class LocusCoverageRow:
    """Per-locus taxon coverage across the taxon set."""

    locus_name: str
    covered_taxon_count: int
    total_taxa: int
    taxon_coverage_fraction: float
    observed_site_count: int
    total_site_count: int
    low_coverage: bool


@dataclass(slots=True)
class LocusOccupancyReport:
    """Coverage report over a concatenated multi-locus alignment."""

    alignment_path: Path
    partition_path: Path
    taxon_count: int
    locus_count: int
    alignment_length: int
    assigned_site_count: int
    unassigned_site_count: int
    partitions: tuple[LocusPartition, ...]
    cells: list[LocusOccupancyCell]
    taxa: list[TaxonCoverageRow]
    loci: list[LocusCoverageRow]
    low_coverage_taxa: list[str]
    low_coverage_loci: list[str]
    taxon_coverage_threshold: float | None
    locus_coverage_threshold: float | None
    warnings: list[str]


def _strip_partition_comments(raw_line: str) -> str:
    line = raw_line.split("#", 1)[0].strip()
    if line.lower().startswith("charset "):
        return line[8:].strip()
    return line


def _parse_segment(raw: str) -> LocusSegment:
    text = raw.strip()
    if not text:
        raise InvalidPartitionError("partition file contains an empty coordinate block")
    if "-" in text:
        left, right = text.split("-", 1)
    else:
        left = right = text
    try:
        start = int(left.strip())
        end = int(right.strip())
    except ValueError as error:
        raise InvalidPartitionError(
            f"partition coordinate '{text}' is not an integer range"
        ) from error
    if start < 1 or end < 1:
        raise InvalidPartitionError(
            f"partition coordinate '{text}' must use positive 1-based positions"
        )
    if end < start:
        raise InvalidPartitionError(
            f"partition coordinate '{text}' ends before it starts"
        )
    return LocusSegment(start=start, end=end)


def parse_locus_partitions(path: Path) -> tuple[LocusPartition, ...]:
    """Parse a simple IQ-TREE, RAxML, or NEXUS-style partition file."""
    partitions: list[LocusPartition] = []
    for line_number, raw_line in enumerate(
        path.read_text(encoding="utf-8").splitlines(), start=1
    ):
        line = _strip_partition_comments(raw_line)
        if not line or line.lower() in {"begin sets;", "end;"}:
            continue
        match = _PARTITION_LINE.match(line)
        if match is None:
            raise InvalidPartitionError(
                f"partition line {line_number} is not in NAME=START-END form"
            )
        name = match.group("name").strip()
        if not name:
            raise InvalidPartitionError(
                f"partition line {line_number} does not declare a locus name"
            )
        segments = tuple(
            _parse_segment(block) for block in match.group("ranges").split(",")
        )
        partitions.append(
            LocusPartition(
                name=name,
                segments=segments,
                total_sites=sum(segment.length for segment in segments),
                data_type=match.group("data_type"),
            )
        )
    if not partitions:
        raise InvalidPartitionError("partition file does not define any loci")
    return tuple(partitions)


def _validate_threshold(value: float | None, label: str) -> None:
    if value is None:
        return
    if not 0.0 <= value <= 1.0:
        raise ValueError(f"{label} must be between 0 and 1 inclusive")


def _validate_partitions(
    partitions: tuple[LocusPartition, ...],
    *,
    alignment_length: int,
) -> tuple[int, int]:
    seen_names: set[str] = set()
    occupied_sites: set[int] = set()
    for partition in partitions:
        if partition.name in seen_names:
            raise InvalidPartitionError(
                f"partition name '{partition.name}' appears more than once"
            )
        seen_names.add(partition.name)
        for segment in partition.segments:
            if segment.end > alignment_length:
                raise InvalidPartitionError(
                    f"partition '{partition.name}' extends beyond the alignment length"
                )
            for site in range(segment.start, segment.end + 1):
                if site in occupied_sites:
                    raise InvalidPartitionError(
                        f"partition '{partition.name}' overlaps another locus at site {site}"
                    )
                occupied_sites.add(site)
    assigned_site_count = len(occupied_sites)
    return assigned_site_count, alignment_length - assigned_site_count


def _slice_partition_sequence(sequence: str, partition: LocusPartition) -> str:
    return "".join(sequence[segment.start - 1 : segment.end] for segment in partition.segments)


def _partition_records(
    records: list[AlignmentRecord], partition: LocusPartition
) -> list[AlignmentRecord]:
    return [
        AlignmentRecord(
            identifier=record.identifier,
            sequence=_slice_partition_sequence(record.sequence, partition),
        )
        for record in records
    ]


def _observed_site_count(
    *,
    alignment_length: int,
    gap_fraction: float,
    missing_fraction: float,
) -> int:
    absent = round(alignment_length * (gap_fraction + missing_fraction))
    observed = alignment_length - absent
    return max(0, observed)


def build_locus_occupancy_report(
    alignment_path: Path,
    partition_path: Path,
    *,
    taxon_coverage_threshold: float | None = None,
    locus_coverage_threshold: float | None = None,
) -> LocusOccupancyReport:
    """Quantify locus occupancy across taxa for a concatenated alignment."""
    _validate_threshold(taxon_coverage_threshold, "taxon coverage threshold")
    _validate_threshold(locus_coverage_threshold, "locus coverage threshold")

    records = load_fasta_alignment(alignment_path)
    if not records:
        raise InvalidAlignmentError("alignment does not contain any records")

    partitions = parse_locus_partitions(partition_path)
    alignment_length = len(records[0].sequence)
    assigned_site_count, unassigned_site_count = _validate_partitions(
        partitions, alignment_length=alignment_length
    )

    warnings: list[str] = []
    if unassigned_site_count > 0:
        warnings.append(
            f"{unassigned_site_count} alignment sites are not assigned to any locus"
        )

    cells: list[LocusOccupancyCell] = []
    occupancies_by_taxon: dict[str, dict[str, float]] = {
        record.identifier: {} for record in records
    }
    taxon_site_totals: dict[str, int] = {record.identifier: 0 for record in records}
    locus_rows: list[LocusCoverageRow] = []

    for partition in partitions:
        partition_records = _partition_records(records, partition)
        summary = summarise_records_as_alignment_summary(
            path=alignment_path, records=partition_records
        )
        uncertainty_by_taxon = {
            row.identifier: row for row in summary.per_sequence_uncertainty
        }
        covered_taxon_count = 0
        locus_observed_site_count = 0
        for record in partition_records:
            profile = uncertainty_by_taxon[record.identifier]
            observed_site_count = _observed_site_count(
                alignment_length=partition.total_sites,
                gap_fraction=profile.gap_fraction,
                missing_fraction=profile.missing_fraction,
            )
            occupancy_fraction = observed_site_count / partition.total_sites
            present = observed_site_count > 0
            if present:
                covered_taxon_count += 1
            locus_observed_site_count += observed_site_count
            taxon_site_totals[record.identifier] += observed_site_count
            occupancies_by_taxon[record.identifier][partition.name] = occupancy_fraction
            cells.append(
                LocusOccupancyCell(
                    taxon=record.identifier,
                    locus_name=partition.name,
                    observed_site_count=observed_site_count,
                    total_site_count=partition.total_sites,
                    occupancy_fraction=occupancy_fraction,
                    present=present,
                )
            )
        taxon_coverage_fraction = covered_taxon_count / len(records)
        low_coverage = (
            locus_coverage_threshold is not None
            and taxon_coverage_fraction < locus_coverage_threshold
        )
        locus_rows.append(
            LocusCoverageRow(
                locus_name=partition.name,
                covered_taxon_count=covered_taxon_count,
                total_taxa=len(records),
                taxon_coverage_fraction=taxon_coverage_fraction,
                observed_site_count=locus_observed_site_count,
                total_site_count=partition.total_sites * len(records),
                low_coverage=low_coverage,
            )
        )

    taxon_rows: list[TaxonCoverageRow] = []
    for record in records:
        occupancies = occupancies_by_taxon[record.identifier]
        covered_locus_count = sum(value > 0.0 for value in occupancies.values())
        total_locus_count = len(partitions)
        locus_coverage_fraction = covered_locus_count / total_locus_count
        low_coverage = (
            taxon_coverage_threshold is not None
            and locus_coverage_fraction < taxon_coverage_threshold
        )
        taxon_rows.append(
            TaxonCoverageRow(
                taxon=record.identifier,
                covered_locus_count=covered_locus_count,
                total_locus_count=total_locus_count,
                locus_coverage_fraction=locus_coverage_fraction,
                observed_site_count=taxon_site_totals[record.identifier],
                total_site_count=assigned_site_count,
                low_coverage=low_coverage,
                occupancies=occupancies,
            )
        )

    low_coverage_taxa = [row.taxon for row in taxon_rows if row.low_coverage]
    low_coverage_loci = [row.locus_name for row in locus_rows if row.low_coverage]
    return LocusOccupancyReport(
        alignment_path=alignment_path,
        partition_path=partition_path,
        taxon_count=len(records),
        locus_count=len(partitions),
        alignment_length=alignment_length,
        assigned_site_count=assigned_site_count,
        unassigned_site_count=unassigned_site_count,
        partitions=partitions,
        cells=cells,
        taxa=taxon_rows,
        loci=locus_rows,
        low_coverage_taxa=low_coverage_taxa,
        low_coverage_loci=low_coverage_loci,
        taxon_coverage_threshold=taxon_coverage_threshold,
        locus_coverage_threshold=locus_coverage_threshold,
        warnings=warnings,
    )
