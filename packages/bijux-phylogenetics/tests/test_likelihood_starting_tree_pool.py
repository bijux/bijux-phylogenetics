from __future__ import annotations

import math
from pathlib import Path

import pytest

from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    build_nucleotide_likelihood_starting_tree_pool,
    build_nucleotide_likelihood_starting_tree_pool_from_alignment,
    validate_nucleotide_likelihood_random_start_tree_count,
    validate_nucleotide_likelihood_starting_tree_pool_model,
)
from bijux_phylogenetics.phylo.likelihood.stepwise_addition import (
    build_likelihood_stepwise_addition_tree_from_alignment,
)

pytestmark = pytest.mark.slow

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_likelihood_starting_tree_pool_reports_distinct_scored_strategies() -> None:
    report = build_nucleotide_likelihood_starting_tree_pool_from_alignment(
        fixture("trees", "jc69_likelihood_nni_start_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        random_start_tree_count=2,
        random_start_tree_seed=17,
    )

    assert report.algorithm == "nucleotide-likelihood-starting-tree-pool"
    assert report.model_name == "JC69"
    assert report.taxon_count == 4
    assert report.site_count == 12
    assert report.pattern_count == 2
    assert report.random_start_tree_count == 2
    assert report.random_start_tree_seed == 17
    assert [row.tree_id for row in report.starting_tree_summaries] == [
        "input-tree",
        "likelihood-stepwise-addition-tree",
        "random-tree-seed-17",
        "random-tree-seed-18",
    ]
    assert [row.source_strategy for row in report.starting_tree_summaries] == [
        "input-tree",
        "likelihood-stepwise-addition-tree",
        "random-tree",
        "random-tree",
    ]
    assert [row.generation_seed for row in report.starting_tree_summaries] == [
        None,
        None,
        17,
        18,
    ]
    assert len({row.topology_hash for row in report.starting_tree_summaries}) == 4
    assert all(
        math.isfinite(row.starting_log_likelihood)
        for row in report.starting_tree_summaries
    )
    assert {
        row.substitution_parameter_policy for row in report.starting_tree_summaries
    } == {"fixed-from-model"}
    assert all(
        row.substitution_parameter_values == {}
        for row in report.starting_tree_summaries
    )
    assert all(
        row.substitution_parameter_warnings == []
        for row in report.starting_tree_summaries
    )


def test_likelihood_starting_tree_pool_rejects_duplicate_strategy_topology_hashes() -> (
    None
):
    stepwise_tree, _report = build_likelihood_stepwise_addition_tree_from_alignment(
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
    )

    try:
        build_nucleotide_likelihood_starting_tree_pool(
            stepwise_tree,
            fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
            model_name="jc69",
            random_start_tree_count=1,
            random_start_tree_seed=17,
        )
    except ValueError as error:
        assert "duplicate topology hash" in str(error)
        assert "input-tree" in str(error)
        assert "likelihood-stepwise-addition-tree" in str(error)
    else:
        raise AssertionError("duplicate starting-tree topologies must fail clearly")


def test_likelihood_starting_tree_pool_validators_reject_unsupported_inputs() -> None:
    try:
        validate_nucleotide_likelihood_random_start_tree_count(0)
    except ValueError as error:
        assert str(error) == "random_start_tree_count must be at least one"
    else:
        raise AssertionError("random_start_tree_count validator must reject zero")

    try:
        validate_nucleotide_likelihood_starting_tree_pool_model("k80")
    except ValueError as error:
        assert (
            str(error)
            == "nucleotide likelihood starting-tree pool model_name must be 'jc69'"
        )
    else:
        raise AssertionError("starting-tree pool model validator must reject k80")


def test_likelihood_starting_tree_pool_preserves_declared_tree_input_path() -> None:
    report = build_nucleotide_likelihood_starting_tree_pool_from_alignment(
        fixture("trees", "jc69_likelihood_spr_start_tree_5_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_spr_alignment_5_taxa.fasta"),
        model_name="jc69",
        random_start_tree_count=1,
        random_start_tree_seed=23,
    )

    assert report.tree_path == str(
        fixture("trees", "jc69_likelihood_spr_start_tree_5_taxa.nwk")
    )
    assert report.alignment_path == str(
        fixture("alignments", "jc69_likelihood_spr_alignment_5_taxa.fasta")
    )
    assert report.taxon_count == len(
        load_tree(
            fixture("trees", "jc69_likelihood_spr_start_tree_5_taxa.nwk")
        ).tip_names
    )
