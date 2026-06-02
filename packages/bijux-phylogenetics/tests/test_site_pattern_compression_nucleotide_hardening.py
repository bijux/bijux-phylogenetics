from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    compress_alignment_site_patterns_from_records,
    evaluate_nucleotide_site_log_likelihoods_from_alignment,
    sum_alignment_site_log_likelihoods,
    sum_compressed_site_pattern_log_likelihoods,
)
from bijux_phylogenetics.phylo.likelihood.dna import (
    UNIFORM_DNA_ROOT_PRIOR,
    evaluate_fixed_topology_dna_site_log_likelihood,
)
from bijux_phylogenetics.phylo.likelihood.dna_observation_policies import (
    normalize_dna_likelihood_records,
)
from bijux_phylogenetics.phylo.likelihood.jc69 import jc69_transition_probability_matrix

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_compressed_nucleotide_patterns_match_uncompressed_likelihood_with_ambiguity() -> (
    None
):
    tree = load_tree(fixture("trees", "likelihood_site_pattern_tree_4_taxa.nwk"))
    records = normalize_dna_likelihood_records(
        load_fasta_alignment(
            fixture("alignments", "jc69_site_pattern_hardening_alignment_4_taxa.fasta")
        ),
        model_name="JC69",
        observation_policy="ambiguity-vector",
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(records)

    def site_log_likelihood(states: tuple[str, ...]) -> float:
        return evaluate_fixed_topology_dna_site_log_likelihood(
            tree,
            states,
            taxon_order=compressed_patterns.taxon_order,
            model_name="JC69",
            observation_policy="ambiguity-vector",
            root_prior=UNIFORM_DNA_ROOT_PRIOR,
            transition_matrix_for_child=lambda child: (
                jc69_transition_probability_matrix(
                    max(float(child.branch_length or 0.0), 0.0)
                )
            ),
        )

    uncompressed_total = sum_alignment_site_log_likelihoods(
        records,
        site_log_likelihood=site_log_likelihood,
    )
    compressed_total = sum_compressed_site_pattern_log_likelihoods(
        compressed_patterns,
        site_log_likelihood=site_log_likelihood,
    )

    assert compressed_patterns.pattern_count == 4
    assert [pattern.weight for pattern in compressed_patterns.patterns] == [2, 1, 2, 1]
    assert math.isclose(
        compressed_total,
        uncompressed_total,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_nucleotide_site_log_likelihood_report_preserves_original_site_order() -> None:
    report = evaluate_nucleotide_site_log_likelihoods_from_alignment(
        fixture("trees", "likelihood_site_pattern_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_site_pattern_hardening_alignment_4_taxa.fasta"),
        model_name="jc69",
        observation_policy="ambiguity-vector",
    )

    assert [row.site_position for row in report.site_log_likelihoods] == [
        1,
        2,
        3,
        4,
        5,
        6,
    ]
    assert [row.pattern_id for row in report.site_log_likelihoods] == [
        "pattern-1",
        "pattern-2",
        "pattern-1",
        "pattern-3",
        "pattern-4",
        "pattern-3",
    ]
