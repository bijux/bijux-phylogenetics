from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    evaluate_jc69_tree_likelihood,
    evaluate_nucleotide_site_log_likelihoods_from_alignment,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_site_log_likelihood_export_supports_fifth_state_observation_policy() -> None:
    direct_report = evaluate_jc69_tree_likelihood(
        load_tree(fixture("trees", "jc69_joint_ancestral_difference_tree_3_taxa.nwk")),
        load_fasta_alignment(
            fixture("alignments", "example_alignment_ambiguity.fasta")
        ),
        observation_policy="fifth-state",
    )
    site_report = evaluate_nucleotide_site_log_likelihoods_from_alignment(
        fixture("trees", "jc69_joint_ancestral_difference_tree_3_taxa.nwk"),
        fixture("alignments", "example_alignment_ambiguity.fasta"),
        model_name="jc69",
        observation_policy="fifth-state",
    )

    assert site_report.state_count == 5
    assert site_report.observation_policy == "fifth-state"
    assert len(site_report.site_log_likelihoods) == 6
    assert math.isclose(
        site_report.log_likelihood,
        direct_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
