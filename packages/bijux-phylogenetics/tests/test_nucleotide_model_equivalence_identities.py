from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    evaluate_jc69_tree_likelihood,
    evaluate_k80_tree_likelihood,
)

FIXTURES = Path(__file__).parent / "fixtures"
_UNIFORM_BASE_FREQUENCIES = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def _load_reference_fixture():
    tree = load_tree(fixture("trees", "likelihood_site_pattern_tree_4_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "jc69_site_pattern_alignment.fasta")
    )
    return tree, records


def _assert_likelihood_match(left: float, right: float) -> None:
    assert math.isclose(left, right, rel_tol=0.0, abs_tol=1e-12)


def _assert_likelihood_mismatch(left: float, right: float) -> None:
    assert not math.isclose(left, right, rel_tol=0.0, abs_tol=1e-9)


def test_k80_with_unit_kappa_matches_jc69() -> None:
    tree, records = _load_reference_fixture()

    jc69 = evaluate_jc69_tree_likelihood(tree, records)
    k80 = evaluate_k80_tree_likelihood(tree, records, kappa=1.0)

    _assert_likelihood_match(k80.log_likelihood, jc69.log_likelihood)


def test_k80_with_nonunit_kappa_diverges_from_jc69() -> None:
    tree, records = _load_reference_fixture()

    jc69 = evaluate_jc69_tree_likelihood(tree, records)
    k80 = evaluate_k80_tree_likelihood(tree, records, kappa=4.0)

    _assert_likelihood_mismatch(k80.log_likelihood, jc69.log_likelihood)
