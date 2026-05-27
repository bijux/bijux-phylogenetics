from __future__ import annotations

import math
from pathlib import Path

import bijux_phylogenetics.distance as distance_module
from bijux_phylogenetics.distance.genetic_distance_matrix import (
    builder as distance_builder_module,
)
from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    evaluate_jc69_tree_likelihood,
    evaluate_jc69_tree_likelihood_from_alignment,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_jc69_fixed_tree_likelihood_matches_two_tip_analytical_fixture() -> None:
    tree_path = fixture("trees", "jc69_likelihood_tree_2_taxa.nwk")
    alignment_path = fixture("alignments", "jc69_likelihood_alignment_2_taxa.fasta")

    report = evaluate_jc69_tree_likelihood_from_alignment(tree_path, alignment_path)
    expected_probability = _expected_two_tip_fixture_probability()

    assert report.taxa == ["A", "B"]
    assert report.site_count == 4
    assert report.pattern_count == 4
    assert report.compression_used is True
    assert report.tree_newick == "(A:0.1,B:0.2);"
    assert math.isclose(
        report.log_likelihood,
        math.log(expected_probability),
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_jc69_likelihood_does_not_reuse_distance_jc69_surface(
    monkeypatch,
) -> None:
    def fail_if_called(*_args, **_kwargs):
        raise AssertionError("JC69 likelihood must not call the distance JC69 surface")

    monkeypatch.setattr(
        distance_module,
        "compute_pairwise_genetic_distance_matrix",
        fail_if_called,
    )
    monkeypatch.setattr(
        distance_builder_module,
        "compute_pairwise_genetic_distance_matrix",
        fail_if_called,
    )

    report = evaluate_jc69_tree_likelihood(
        load_tree(fixture("trees", "jc69_likelihood_tree_2_taxa.nwk")),
        load_fasta_alignment(
            fixture("alignments", "jc69_likelihood_alignment_2_taxa.fasta")
        ),
    )

    assert math.isfinite(report.log_likelihood)


def _expected_two_tip_fixture_probability() -> float:
    total_branch_length = 0.3
    decay = math.exp((-4.0 * total_branch_length) / 3.0)
    same_site_probability = 0.25 * (0.25 + (0.75 * decay))
    different_site_probability = 0.25 * (0.25 - (0.25 * decay))
    return (same_site_probability**2) * (different_site_probability**2)
