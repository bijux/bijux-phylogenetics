from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.phylo.likelihood import compress_alignment_site_patterns

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_compress_alignment_site_patterns_groups_identical_columns_with_integer_weights() -> (
    None
):
    compressed = compress_alignment_site_patterns(
        fixture("alignments", "jc69_site_pattern_alignment.fasta")
    )

    assert compressed.taxon_order == ["A", "B", "C", "D"]
    assert compressed.alignment_length == 10
    assert compressed.pattern_count == 6
    assert [pattern.states for pattern in compressed.patterns] == [
        ("A", "A", "G", "T"),
        ("C", "C", "C", "C"),
        ("G", "T", "G", "A"),
        ("T", "A", "G", "A"),
        ("G", "A", "G", "T"),
        ("A", "C", "G", "T"),
    ]
    assert [pattern.weight for pattern in compressed.patterns] == [2, 2, 1, 1, 1, 3]
    assert [pattern.site_positions for pattern in compressed.patterns] == [
        [1, 2],
        [3, 4],
        [5],
        [6],
        [7],
        [8, 9, 10],
    ]
