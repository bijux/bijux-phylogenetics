from __future__ import annotations

from collections import OrderedDict
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError


@dataclass(frozen=True, slots=True)
class AlignmentSitePattern:
    """One unique aligned site pattern with its stable integer multiplicity."""

    pattern_id: str
    states: tuple[str, ...]
    weight: int
    site_positions: list[int]


@dataclass(slots=True)
class CompressedAlignmentSitePatterns:
    """One alignment summarized into unique site patterns and integer weights."""

    source_path: Path | None
    taxon_order: list[str]
    alignment_length: int
    pattern_count: int
    patterns: list[AlignmentSitePattern]


def alignment_site_columns(records: list[AlignmentRecord]) -> list[tuple[str, ...]]:
    """Return aligned columns in stable taxon order."""
    alignment_length = _validated_alignment_length(records)
    return [
        tuple(record.sequence[position] for record in records)
        for position in range(alignment_length)
    ]


def compress_alignment_site_patterns(path: Path) -> CompressedAlignmentSitePatterns:
    """Compress identical aligned columns from one FASTA alignment file."""
    from bijux_phylogenetics.io.fasta.core import load_fasta_alignment

    records = load_fasta_alignment(path)
    return compress_alignment_site_patterns_from_records(records, source_path=path)


def compress_alignment_site_patterns_from_records(
    records: list[AlignmentRecord],
    *,
    source_path: Path | None = None,
) -> CompressedAlignmentSitePatterns:
    """Compress identical aligned columns while preserving stable taxon order."""
    alignment_length = _validated_alignment_length(records)
    grouped_positions: OrderedDict[tuple[str, ...], list[int]] = OrderedDict()
    for position, states in enumerate(alignment_site_columns(records), start=1):
        grouped_positions.setdefault(states, []).append(position)

    patterns = [
        AlignmentSitePattern(
            pattern_id=f"pattern-{index}",
            states=states,
            weight=len(site_positions),
            site_positions=site_positions,
        )
        for index, (states, site_positions) in enumerate(
            grouped_positions.items(),
            start=1,
        )
    ]
    return CompressedAlignmentSitePatterns(
        source_path=source_path,
        taxon_order=[record.identifier for record in records],
        alignment_length=alignment_length,
        pattern_count=len(patterns),
        patterns=patterns,
    )


def iter_uncompressed_alignment_sites(
    records: list[AlignmentRecord],
) -> Iterable[tuple[str, ...]]:
    """Yield aligned columns in stable taxon order."""
    yield from alignment_site_columns(records)


def _validated_alignment_length(records: list[AlignmentRecord]) -> int:
    if not records:
        raise InvalidAlignmentError(
            "site-pattern compression requires at least one alignment record"
        )
    alignment_lengths = {len(record.sequence) for record in records}
    if len(alignment_lengths) != 1:
        raise InvalidAlignmentError(
            "site-pattern compression requires equal-length aligned sequences"
        )
    return next(iter(alignment_lengths))
