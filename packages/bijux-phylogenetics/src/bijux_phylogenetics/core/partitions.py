from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re

from bijux_phylogenetics.runtime.errors import InvalidPartitionError

__all__ = [
    "LocusPartition",
    "LocusSegment",
    "PartitionSummaryReport",
    "PartitionSummaryRow",
    "build_partition_summary_report",
    "iterate_partition_sites",
    "normalize_partition_data_type",
    "partition_coordinate_text",
    "parse_locus_partitions",
    "slice_partition_sequence",
    "validate_locus_partitions",
    "write_partition_summary_table",
    "write_locus_partitions",
]

_PARTITION_LINE = re.compile(
    r"^(?:(?P<data_type>[A-Za-z0-9_-]+)\s*,\s*)?"
    r"(?P<name>[^=]+?)\s*=\s*(?P<ranges>[^;]+?)\s*;?$"
)
_STEP_SEGMENT = re.compile(r"^(?P<start>\d+)-(?P<end>\d+)\\(?P<step>\d+)$")
_RANGE_SEGMENT = re.compile(r"^(?P<start>\d+)-(?P<end>\d+)$")
_POINT_SEGMENT = re.compile(r"^(?P<site>\d+)$")


@dataclass(frozen=True, slots=True)
class LocusSegment:
    """One inclusive 1-based coordinate block inside a locus partition."""

    start: int
    end: int
    step: int = 1

    @property
    def length(self) -> int:
        """Return the number of covered sites after applying the stride."""
        return ((self.end - self.start) // self.step) + 1


@dataclass(frozen=True, slots=True)
class LocusPartition:
    """One named locus assembled from one or more concatenated segments."""

    name: str
    segments: tuple[LocusSegment, ...]
    total_sites: int
    data_type: str | None = None


@dataclass(frozen=True, slots=True)
class PartitionSummaryRow:
    """One row in a validated partition summary table."""

    locus_name: str
    data_type: str | None
    segment_count: int
    total_sites: int
    start_site: int
    end_site: int
    coordinate_text: str


@dataclass(frozen=True, slots=True)
class PartitionSummaryReport:
    """Validated coordinate and datatype summary for one partition file."""

    alignment_length: int
    partition_count: int
    assigned_site_count: int
    unassigned_site_count: int
    mixed_data_types: bool
    declared_data_types: list[str]
    rows: list[PartitionSummaryRow]
    warnings: list[str]


def normalize_partition_data_type(data_type: str | None) -> str | None:
    """Map common partition datatype aliases to stable uppercase names."""
    if data_type is None:
        return None
    normalized = data_type.strip().upper()
    aliases = {
        "AA": "PROTEIN",
        "AMINOACID": "PROTEIN",
        "AMINO-ACID": "PROTEIN",
        "AMINO_ACID": "PROTEIN",
        "PROT": "PROTEIN",
    }
    return aliases.get(normalized, normalized)


def _strip_partition_comments(raw_line: str) -> str:
    line = raw_line.split("#", 1)[0].strip()
    if line.lower().startswith("charset "):
        return line[8:].strip()
    return line


def _parse_segment(raw: str) -> LocusSegment:
    text = raw.strip()
    if not text:
        raise InvalidPartitionError("partition file contains an empty coordinate block")
    match = _STEP_SEGMENT.match(text)
    if match is not None:
        start = int(match.group("start"))
        end = int(match.group("end"))
        step = int(match.group("step"))
    else:
        match = _RANGE_SEGMENT.match(text)
        if match is not None:
            start = int(match.group("start"))
            end = int(match.group("end"))
            step = 1
        else:
            match = _POINT_SEGMENT.match(text)
            if match is None:
                raise InvalidPartitionError(
                    f"partition coordinate '{text}' is not an integer range"
                )
            start = int(match.group("site"))
            end = start
            step = 1
    if start < 1 or end < 1:
        raise InvalidPartitionError(
            f"partition coordinate '{text}' must use positive 1-based positions"
        )
    if end < start:
        raise InvalidPartitionError(
            f"partition coordinate '{text}' ends before it starts"
        )
    if step < 1:
        raise InvalidPartitionError(
            f"partition coordinate '{text}' must use a positive stride"
        )
    return LocusSegment(start=start, end=end, step=step)


def _split_coordinate_blocks(raw: str) -> list[str]:
    blocks: list[str] = []
    for group in raw.split(","):
        blocks.extend(token for token in group.split() if token)
    return blocks


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
            _parse_segment(block)
            for block in _split_coordinate_blocks(match.group("ranges"))
        )
        partitions.append(
            LocusPartition(
                name=name,
                segments=segments,
                total_sites=sum(segment.length for segment in segments),
                data_type=normalize_partition_data_type(match.group("data_type")),
            )
        )
    if not partitions:
        raise InvalidPartitionError("partition file does not define any loci")
    return tuple(partitions)


def iterate_partition_sites(segment: LocusSegment) -> range:
    """Return the covered 1-based site coordinates for one partition segment."""
    return range(segment.start, segment.end + 1, segment.step)


def validate_locus_partitions(
    partitions: tuple[LocusPartition, ...],
    *,
    alignment_length: int,
) -> tuple[int, int]:
    """Reject duplicate names, overlaps, and coordinates beyond the alignment."""
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
            for site in iterate_partition_sites(segment):
                if site in occupied_sites:
                    raise InvalidPartitionError(
                        f"partition '{partition.name}' overlaps another locus at site {site}"
                    )
                occupied_sites.add(site)
    assigned_site_count = len(occupied_sites)
    return assigned_site_count, alignment_length - assigned_site_count


def slice_partition_sequence(sequence: str, partition: LocusPartition) -> str:
    """Extract one concatenated partition sequence from an aligned record."""
    return "".join(
        "".join(sequence[site - 1] for site in iterate_partition_sites(segment))
        for segment in partition.segments
    )


def _serialize_segment(segment: LocusSegment) -> str:
    if segment.start == segment.end:
        return str(segment.start)
    if segment.step == 1:
        return f"{segment.start}-{segment.end}"
    return f"{segment.start}-{segment.end}\\{segment.step}"


def partition_coordinate_text(partition: LocusPartition) -> str:
    """Serialize one partition into stable coordinate text."""
    return ",".join(_serialize_segment(segment) for segment in partition.segments)


def build_partition_summary_report(
    partitions: tuple[LocusPartition, ...],
    *,
    alignment_length: int,
) -> PartitionSummaryReport:
    """Validate one partition set and summarize its locus coverage."""
    assigned_site_count, unassigned_site_count = validate_locus_partitions(
        partitions,
        alignment_length=alignment_length,
    )
    declared_data_types = sorted(
        {
            partition.data_type
            for partition in partitions
            if partition.data_type is not None
        }
    )
    warnings: list[str] = []
    if unassigned_site_count > 0:
        warnings.append(
            f"{unassigned_site_count} alignment sites are not assigned to any partition"
        )
    rows = [
        PartitionSummaryRow(
            locus_name=partition.name,
            data_type=partition.data_type,
            segment_count=len(partition.segments),
            total_sites=partition.total_sites,
            start_site=min(segment.start for segment in partition.segments),
            end_site=max(segment.end for segment in partition.segments),
            coordinate_text=partition_coordinate_text(partition),
        )
        for partition in partitions
    ]
    return PartitionSummaryReport(
        alignment_length=alignment_length,
        partition_count=len(partitions),
        assigned_site_count=assigned_site_count,
        unassigned_site_count=unassigned_site_count,
        mixed_data_types=len(declared_data_types) > 1,
        declared_data_types=declared_data_types,
        rows=rows,
        warnings=warnings,
    )


def write_locus_partitions(path: Path, partitions: tuple[LocusPartition, ...]) -> Path:
    """Write a partition file for the retained loci in concatenated alignment order."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    for partition in partitions:
        prefix = f"{partition.data_type}," if partition.data_type else ""
        lines.append(
            f"{prefix}{partition.name} = {partition_coordinate_text(partition)}"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_partition_summary_table(path: Path, report: PartitionSummaryReport) -> Path:
    """Persist one stable TSV summary of the validated partition set."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "\t".join(
            [
                "locus_name",
                "data_type",
                "segment_count",
                "total_sites",
                "start_site",
                "end_site",
                "coordinate_text",
            ]
        )
    ]
    lines.extend(
        "\t".join(
            [
                row.locus_name,
                row.data_type or "",
                str(row.segment_count),
                str(row.total_sites),
                str(row.start_site),
                str(row.end_site),
                row.coordinate_text,
            ]
        )
        for row in report.rows
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
