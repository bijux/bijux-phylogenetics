from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    build_equal_rate_codon_ctmc_rate_matrix,
    sum_compressed_site_pattern_log_likelihoods,
)
from bijux_phylogenetics.phylo.likelihood.codon import (
    _evaluate_codon_site_log_likelihood,
    compress_codon_site_patterns_from_records,
    normalize_codon_likelihood_records,
)
from bijux_phylogenetics.phylo.likelihood.pruning import (
    build_transition_matrix_evaluator,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_compressed_codon_patterns_match_uncompressed_likelihood() -> None:
    tree = load_tree(fixture("trees", "codon_likelihood_tree_2_taxa.nwk"))
    records = normalize_codon_likelihood_records(
        load_fasta_alignment(
            fixture("alignments", "codon_site_pattern_alignment_2_taxa.fasta")
        )
    )
    compressed_patterns = compress_codon_site_patterns_from_records(records)
    state_space, rate_matrix, root_prior, _source = (
        build_equal_rate_codon_ctmc_rate_matrix()
    )
    transition_evaluator = build_transition_matrix_evaluator(rate_matrix)
    codon_sequences = [
        [
            record.sequence[index : index + 3]
            for index in range(0, len(record.sequence), 3)
        ]
        for record in records
    ]

    def site_log_likelihood(states: tuple[str, ...]) -> float:
        return _evaluate_codon_site_log_likelihood(
            tree,
            states,
            taxon_order=compressed_patterns.taxon_order,
            state_space=state_space,
            root_prior=root_prior,
            transition_evaluator=transition_evaluator,
        )

    uncompressed_total = sum(
        site_log_likelihood(states)
        for states in zip(*codon_sequences, strict=True)
    )
    compressed_total = sum_compressed_site_pattern_log_likelihoods(
        compressed_patterns,
        site_log_likelihood=site_log_likelihood,
    )

    assert compressed_patterns.pattern_count == 3
    assert [pattern.weight for pattern in compressed_patterns.patterns] == [2, 2, 2]
    assert math.isclose(
        compressed_total,
        uncompressed_total,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
