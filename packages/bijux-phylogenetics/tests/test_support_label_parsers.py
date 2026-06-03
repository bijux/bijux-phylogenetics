from __future__ import annotations

from bijux_phylogenetics.io.fasttree_support import parse_fasttree_branch_support_label
from bijux_phylogenetics.io.iqtree_support import (
    parse_iqtree_branch_support_label,
    support_fraction,
)


def test_parse_fasttree_branch_support_label_accepts_scientific_notation() -> None:
    label = parse_fasttree_branch_support_label(" 9.5e-01 ")

    assert label is not None
    assert label.raw_label == "9.5e-01"
    assert label.local_support == 0.95


def test_parse_iqtree_branch_support_label_accepts_spaced_composite_values() -> None:
    label = parse_iqtree_branch_support_label(" 82 / 97 ")

    assert label is not None
    assert label.raw_label == "82 / 97"
    assert label.sh_alrt_support == 82.0
    assert label.ufboot_support == 97.0
    assert support_fraction(label.sh_alrt_support) == 0.82
    assert support_fraction(label.ufboot_support) == 0.97


def test_support_label_parsers_reject_nonfinite_values() -> None:
    assert parse_fasttree_branch_support_label("nan") is None
    assert parse_fasttree_branch_support_label("inf") is None
    assert parse_iqtree_branch_support_label("nan") is None
    assert parse_iqtree_branch_support_label("82/inf") is None
    assert support_fraction(float("nan")) is None
