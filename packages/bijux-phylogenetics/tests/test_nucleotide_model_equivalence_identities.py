from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    evaluate_gtr_tree_likelihood,
    evaluate_hky85_tree_likelihood,
    evaluate_jc69_tree_likelihood,
    evaluate_k80_tree_likelihood,
    evaluate_nucleotide_site_log_likelihoods_from_alignment,
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


def _site_log_likelihood_vector(report) -> list[float]:
    return [row.log_likelihood for row in report.site_log_likelihoods]


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


def test_hky85_with_uniform_frequencies_and_unit_kappa_matches_jc69() -> None:
    tree, records = _load_reference_fixture()

    jc69 = evaluate_jc69_tree_likelihood(tree, records)
    hky85 = evaluate_hky85_tree_likelihood(
        tree,
        records,
        kappa=1.0,
        base_frequencies=_UNIFORM_BASE_FREQUENCIES,
    )

    _assert_likelihood_match(hky85.log_likelihood, jc69.log_likelihood)


def test_hky85_with_nonunit_parameters_diverges_from_jc69() -> None:
    tree, records = _load_reference_fixture()

    jc69 = evaluate_jc69_tree_likelihood(tree, records)
    hky85 = evaluate_hky85_tree_likelihood(
        tree,
        records,
        kappa=4.0,
        base_frequencies={"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3},
    )

    _assert_likelihood_mismatch(hky85.log_likelihood, jc69.log_likelihood)


def test_gtr_with_uniform_exchangeabilities_and_frequencies_matches_jc69() -> None:
    tree, records = _load_reference_fixture()

    jc69 = evaluate_jc69_tree_likelihood(tree, records)
    gtr = evaluate_gtr_tree_likelihood(
        tree,
        records,
        exchangeabilities={
            "AC": 1.0,
            "AG": 1.0,
            "AT": 1.0,
            "CG": 1.0,
            "CT": 1.0,
            "GT": 1.0,
        },
        base_frequencies=_UNIFORM_BASE_FREQUENCIES,
    )

    _assert_likelihood_match(gtr.log_likelihood, jc69.log_likelihood)


def test_gtr_with_nonunit_parameters_diverges_from_jc69() -> None:
    tree, records = _load_reference_fixture()

    jc69 = evaluate_jc69_tree_likelihood(tree, records)
    gtr = evaluate_gtr_tree_likelihood(
        tree,
        records,
        exchangeabilities={
            "AC": 1.0,
            "AG": 4.5,
            "AT": 0.8,
            "CG": 1.6,
            "CT": 2.4,
            "GT": 3.1,
        },
        base_frequencies={"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3},
    )

    _assert_likelihood_mismatch(gtr.log_likelihood, jc69.log_likelihood)


def test_k80_site_rows_match_jc69_when_kappa_is_one() -> None:
    jc69 = evaluate_nucleotide_site_log_likelihoods_from_alignment(
        fixture("trees", "likelihood_site_pattern_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_site_pattern_alignment.fasta"),
        model_name="jc69",
    )
    k80 = evaluate_nucleotide_site_log_likelihoods_from_alignment(
        fixture("trees", "likelihood_site_pattern_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_site_pattern_alignment.fasta"),
        model_name="k80",
        kappa=1.0,
    )

    assert len(jc69.site_log_likelihoods) == len(k80.site_log_likelihoods)
    for left, right in zip(_site_log_likelihood_vector(jc69), _site_log_likelihood_vector(k80), strict=True):
        _assert_likelihood_match(left, right)


def test_hky85_site_rows_match_jc69_under_uniform_unit_parameters() -> None:
    jc69 = evaluate_nucleotide_site_log_likelihoods_from_alignment(
        fixture("trees", "likelihood_site_pattern_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_site_pattern_alignment.fasta"),
        model_name="jc69",
    )
    hky85 = evaluate_nucleotide_site_log_likelihoods_from_alignment(
        fixture("trees", "likelihood_site_pattern_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_site_pattern_alignment.fasta"),
        model_name="hky85",
        kappa=1.0,
        base_frequencies=_UNIFORM_BASE_FREQUENCIES,
    )

    assert len(jc69.site_log_likelihoods) == len(hky85.site_log_likelihoods)
    for left, right in zip(_site_log_likelihood_vector(jc69), _site_log_likelihood_vector(hky85), strict=True):
        _assert_likelihood_match(left, right)
