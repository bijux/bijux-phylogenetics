from __future__ import annotations

from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood import (
    compress_alignment_site_patterns_from_records,
)
from bijux_phylogenetics.phylo.likelihood.codon import (
    compress_codon_site_patterns_from_records,
)


def test_shared_tokenized_compressor_matches_codon_pattern_grouping() -> None:
    records = [
        AlignmentRecord(identifier="A", sequence="AAACCCAAA"),
        AlignmentRecord(identifier="B", sequence="GGGTTTGGG"),
    ]

    compressed_patterns = compress_alignment_site_patterns_from_records(
        records,
        site_tokenizer=_codon_tokens,
    )
    codon_patterns = compress_codon_site_patterns_from_records(records)

    assert compressed_patterns.taxon_order == ["A", "B"]
    assert compressed_patterns.alignment_length == 3
    assert compressed_patterns.pattern_count == 2
    assert [pattern.states for pattern in compressed_patterns.patterns] == [
        ("AAA", "GGG"),
        ("CCC", "TTT"),
    ]
    assert [pattern.weight for pattern in compressed_patterns.patterns] == [2, 1]
    assert [pattern.site_positions for pattern in compressed_patterns.patterns] == [
        [1, 3],
        [2],
    ]
    assert compressed_patterns == codon_patterns


def test_shared_tokenized_compressor_preserves_multicharacter_mk_states() -> None:
    records = [
        AlignmentRecord(identifier="A", sequence="north"),
        AlignmentRecord(identifier="B", sequence="south"),
        AlignmentRecord(identifier="C", sequence="north"),
    ]

    compressed_patterns = compress_alignment_site_patterns_from_records(
        records,
        site_tokenizer=_whole_trait_token,
    )

    assert compressed_patterns.taxon_order == ["A", "B", "C"]
    assert compressed_patterns.alignment_length == 1
    assert compressed_patterns.pattern_count == 1
    assert [pattern.states for pattern in compressed_patterns.patterns] == [
        ("north", "south", "north")
    ]
    assert [pattern.weight for pattern in compressed_patterns.patterns] == [1]
    assert [pattern.site_positions for pattern in compressed_patterns.patterns] == [[1]]


def _codon_tokens(sequence: str) -> tuple[str, ...]:
    return tuple(sequence[index : index + 3] for index in range(0, len(sequence), 3))


def _whole_trait_token(sequence: str) -> tuple[str, ...]:
    return (sequence,)
