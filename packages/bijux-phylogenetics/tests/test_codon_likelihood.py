from __future__ import annotations

import math
from pathlib import Path

import pytest

from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.codon import (
    evaluate_codon_ctmc_tree_likelihood,
    evaluate_codon_ctmc_tree_likelihood_from_alignment,
)
from bijux_phylogenetics.phylo.likelihood.codon_states import resolve_codon_state_space
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_codon_ctmc_likelihood_returns_finite_report() -> None:
    report = evaluate_codon_ctmc_tree_likelihood_from_alignment(
        fixture("trees", "codon_likelihood_tree_2_taxa.nwk"),
        fixture("alignments", "codon_likelihood_alignment_2_taxa.fasta"),
    )

    assert report.taxa == ["A", "B"]
    assert report.site_count == 2
    assert report.pattern_count == 2
    assert report.compression_used is False
    assert report.state_count == 61
    assert report.genetic_code_id == 1
    assert report.genetic_code_name == "Standard"
    assert report.codon_frequency_source == "uniform"
    assert report.observation_policy == "reject"
    assert report.tree_newick == "(A:0.1,B:0.2);"
    assert math.isfinite(report.log_likelihood)


def test_codon_ctmc_likelihood_rejects_stop_codon_states() -> None:
    try:
        evaluate_codon_ctmc_tree_likelihood_from_alignment(
            fixture("trees", "codon_likelihood_tree_2_taxa.nwk"),
            fixture("alignments", "codon_likelihood_alignment_stop_2_taxa.fasta"),
        )
    except InvalidAlignmentError as error:
        assert "excludes stop codon states" in str(error)
        assert "TAA" in str(error)
    else:
        raise AssertionError("stop-codon codon site should be rejected")


def test_codon_ctmc_likelihood_changes_with_codon_frequencies() -> None:
    uniform_report = evaluate_codon_ctmc_tree_likelihood_from_alignment(
        fixture("trees", "codon_likelihood_tree_2_taxa.nwk"),
        fixture("alignments", "codon_likelihood_alignment_2_taxa.fasta"),
    )
    state_space = resolve_codon_state_space()
    biased_frequencies = dict.fromkeys(state_space.state_order, 1.0)
    biased_frequencies["AAA"] = 20.0
    biased_frequencies["AAG"] = 15.0
    biased_frequencies["GCT"] = 0.5
    biased_frequencies["GCC"] = 0.5

    biased_report = evaluate_codon_ctmc_tree_likelihood_from_alignment(
        fixture("trees", "codon_likelihood_tree_2_taxa.nwk"),
        fixture("alignments", "codon_likelihood_alignment_2_taxa.fasta"),
        codon_frequencies=biased_frequencies,
    )

    assert biased_report.codon_frequency_source == "provided"
    assert not math.isclose(
        biased_report.log_likelihood,
        uniform_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_codon_ctmc_likelihood_rejects_ambiguous_codons_under_reject_policy() -> None:
    with pytest.raises(
        InvalidAlignmentError,
        match="observation policy 'reject' requires resolved sense codons only",
    ):
        evaluate_codon_ctmc_tree_likelihood_from_alignment(
            fixture("trees", "codon_likelihood_tree_2_taxa.nwk"),
            fixture("alignments", "codon_likelihood_alignment_ambiguity_2_taxa.fasta"),
            observation_policy="reject",
        )


def test_codon_ctmc_observation_policies_change_likelihood_on_ambiguity_fixture() -> None:
    missing_report = evaluate_codon_ctmc_tree_likelihood_from_alignment(
        fixture("trees", "codon_likelihood_tree_2_taxa.nwk"),
        fixture("alignments", "codon_likelihood_alignment_ambiguity_2_taxa.fasta"),
        observation_policy="treat-as-missing",
    )
    ambiguity_report = evaluate_codon_ctmc_tree_likelihood_from_alignment(
        fixture("trees", "codon_likelihood_tree_2_taxa.nwk"),
        fixture("alignments", "codon_likelihood_alignment_ambiguity_2_taxa.fasta"),
        observation_policy="ambiguity-vector",
    )

    assert missing_report.observation_policy == "treat-as-missing"
    assert ambiguity_report.observation_policy == "ambiguity-vector"
    assert math.isfinite(missing_report.log_likelihood)
    assert math.isfinite(ambiguity_report.log_likelihood)
    assert not math.isclose(
        missing_report.log_likelihood,
        ambiguity_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_codon_ctmc_observation_policy_treats_missing_and_gap_codons_explicitly() -> None:
    records = [
        AlignmentRecord(identifier="A", sequence="ATG???GGG"),
        AlignmentRecord(identifier="B", sequence="ATG---GGG"),
    ]

    with pytest.raises(
        InvalidAlignmentError,
        match="does not allow missing or gap codon state",
    ):
        evaluate_codon_ctmc_tree_likelihood(
            load_tree(fixture("trees", "codon_likelihood_tree_2_taxa.nwk")),
            records,
            observation_policy="reject",
        )

    report = evaluate_codon_ctmc_tree_likelihood(
        load_tree(fixture("trees", "codon_likelihood_tree_2_taxa.nwk")),
        records,
        observation_policy="treat-as-missing",
    )

    assert report.observation_policy == "treat-as-missing"
    assert math.isfinite(report.log_likelihood)
