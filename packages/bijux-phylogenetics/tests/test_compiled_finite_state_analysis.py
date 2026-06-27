from __future__ import annotations

import math
from pathlib import Path

import numpy

from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.alignment.partitions import (
    LocusPartition,
    LocusSegment,
)
from bijux_phylogenetics.phylo.likelihood.finite_state import (
    FiniteStateLikelihoodSurface,
    compile_finite_state_likelihood_analysis,
    evaluate_compiled_finite_state_likelihood_analysis,
    evaluate_compiled_finite_state_partition_log_likelihood_rows,
    evaluate_compiled_finite_state_site_log_likelihood_rows,
    validate_compiled_finite_state_partition_likelihood_reconstruction,
    validate_compiled_finite_state_site_log_likelihood_reconstruction,
)
from bijux_phylogenetics.phylo.likelihood.patterns import (
    CompressedAlignmentSitePatterns,
    compress_alignment_site_patterns_from_records,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(name: str) -> Path:
    return FIXTURES / "trees" / name


def test_compiled_analysis_keeps_partition_weights_and_missing_data_policy() -> None:
    tree = load_tree(fixture("felsenstein_two_tip_tree.nwk"))
    records = [
        AlignmentRecord(identifier="A", sequence="0101"),
        AlignmentRecord(identifier="B", sequence="0100"),
    ]
    analysis = compile_finite_state_likelihood_analysis(
        tree,
        _binary_surface(),
        records=records,
        pattern_compiler=_pattern_compiler,
        missing_data_policy="gap=reject;missing=treat-as-missing",
        partition_map=(
            LocusPartition(
                name="left-half",
                segments=(LocusSegment(start=1, end=2),),
                total_sites=2,
            ),
            LocusPartition(
                name="right-half",
                segments=(LocusSegment(start=3, end=4),),
                total_sites=2,
            ),
        ),
        site_weight_by_position={2: 3.0, 4: 0.0},
    )

    assert analysis.alignment_length == 4
    assert analysis.partition_count == 2
    assert analysis.pattern_count == 4
    assert analysis.missing_data_policy == "gap=reject;missing=treat-as-missing"
    assert math.isclose(analysis.site_weight_total, 5.0, abs_tol=1e-12)

    left_partition, right_partition = analysis.partition_contexts
    assert left_partition.partition_name == "left-half"
    assert right_partition.partition_name == "right-half"
    assert left_partition.pattern_count == 2
    assert right_partition.pattern_count == 2

    left_terms = {
        term.source_site_positions: (
            term.source_site_weights,
            term.effective_weight,
        )
        for term in left_partition.compiled_pattern_terms
    }
    right_terms = {
        term.source_site_positions: (
            term.source_site_weights,
            term.effective_weight,
        )
        for term in right_partition.compiled_pattern_terms
    }

    assert left_terms[(1,)] == ((1.0,), 1.0)
    assert left_terms[(2,)] == ((3.0,), 3.0)
    assert right_terms[(3,)] == ((1.0,), 1.0)
    assert right_terms[(4,)] == ((0.0,), 0.0)


def test_compiled_analysis_expands_weighted_source_site_rows_in_original_order() -> None:
    tree = load_tree(fixture("felsenstein_two_tip_tree.nwk"))
    records = [
        AlignmentRecord(identifier="A", sequence="0101"),
        AlignmentRecord(identifier="B", sequence="0100"),
    ]
    analysis = compile_finite_state_likelihood_analysis(
        tree,
        _binary_surface(),
        records=records,
        pattern_compiler=_pattern_compiler,
        missing_data_policy="reject",
        partition_map=(
            LocusPartition(
                name="left-half",
                segments=(LocusSegment(start=1, end=2),),
                total_sites=2,
                data_type="binary",
            ),
            LocusPartition(
                name="right-half",
                segments=(LocusSegment(start=3, end=4),),
                total_sites=2,
                data_type="binary",
            ),
        ),
        site_weight_by_position={2: 3.0, 4: 0.0},
    )

    rows, total_log_likelihood = evaluate_compiled_finite_state_site_log_likelihood_rows(
        analysis
    )

    assert [row.source_site_position for row in rows] == [1, 2, 3, 4]
    assert [row.partition_name for row in rows] == [
        "left-half",
        "left-half",
        "right-half",
        "right-half",
    ]
    assert [row.pattern_site_count for row in rows] == [1, 1, 1, 1]
    assert [row.site_weight for row in rows] == [1.0, 3.0, 1.0, 0.0]
    validate_compiled_finite_state_site_log_likelihood_reconstruction(
        rows,
        expected_total_log_likelihood=total_log_likelihood,
        expected_site_count=analysis.alignment_length,
        expected_pattern_count=analysis.pattern_count,
        owner_name="compiled finite-state analysis test",
    )
    assert math.isclose(
        sum(row.site_weight * row.log_likelihood for row in rows),
        total_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_compiled_analysis_summarizes_partition_likelihood_rows() -> None:
    tree = load_tree(fixture("felsenstein_two_tip_tree.nwk"))
    records = [
        AlignmentRecord(identifier="A", sequence="0101"),
        AlignmentRecord(identifier="B", sequence="0100"),
    ]
    analysis = compile_finite_state_likelihood_analysis(
        tree,
        _binary_surface(),
        records=records,
        pattern_compiler=_pattern_compiler,
        missing_data_policy="reject",
        partition_map=(
            LocusPartition(
                name="left-half",
                segments=(LocusSegment(start=1, end=2),),
                total_sites=2,
                data_type="binary",
            ),
            LocusPartition(
                name="right-half",
                segments=(LocusSegment(start=3, end=4),),
                total_sites=2,
                data_type="binary",
            ),
        ),
        site_weight_by_position={2: 3.0, 4: 0.0},
    )

    rows, total_log_likelihood = evaluate_compiled_finite_state_partition_log_likelihood_rows(
        analysis
    )

    assert [row.partition_name for row in rows] == ["left-half", "right-half"]
    assert [row.model_name for row in rows] == [
        "binary synthetic",
        "binary synthetic",
    ]
    assert [row.context_index for row in rows] == [0, 0]
    assert [row.context_linkage_key for row in rows] == ["shared", "shared"]
    assert [row.site_count for row in rows] == [2, 2]
    assert [row.pattern_count for row in rows] == [2, 2]
    assert [row.site_weight_total for row in rows] == [4.0, 1.0]
    assert rows[0].parameter_values == {}
    assert rows[1].parameter_values == {}
    validate_compiled_finite_state_partition_likelihood_reconstruction(
        rows,
        expected_total_log_likelihood=total_log_likelihood,
        expected_partition_count=analysis.partition_count,
        owner_name="compiled finite-state partition summary test",
    )
    assert math.isclose(
        sum(row.log_likelihood for row in rows),
        total_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_compiled_analysis_reuses_cached_transition_and_pattern_state() -> None:
    tree = load_tree(fixture("felsenstein_two_tip_tree.nwk"))
    records = [
        AlignmentRecord(identifier="A", sequence="01"),
        AlignmentRecord(identifier="B", sequence="01"),
    ]
    transition_evaluation_count = 0

    def transition(branch_length: float) -> numpy.ndarray:
        nonlocal transition_evaluation_count
        transition_evaluation_count += 1
        return numpy.array(
            [
                [1.0 - branch_length, branch_length],
                [branch_length / 2.0, 1.0 - (branch_length / 2.0)],
            ],
            dtype=float,
        )

    analysis = compile_finite_state_likelihood_analysis(
        tree,
        FiniteStateLikelihoodSurface(
            model_name="binary synthetic",
            state_labels=("0", "1"),
            root_prior=numpy.array([0.5, 0.5], dtype=float),
            root_prior_source="provided",
            observation_policy="reject",
            transition_matrix_for_branch_length=transition,
            leaf_likelihood_vector_for_state=_binary_leaf_vector,
        ),
        records=records,
        pattern_compiler=_pattern_compiler,
        missing_data_policy="reject",
    )

    first_log_likelihood = evaluate_compiled_finite_state_likelihood_analysis(analysis)
    second_log_likelihood = evaluate_compiled_finite_state_likelihood_analysis(
        analysis
    )

    assert transition_evaluation_count == 2
    assert math.isclose(
        first_log_likelihood,
        second_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def _binary_surface() -> FiniteStateLikelihoodSurface:
    return FiniteStateLikelihoodSurface(
        model_name="binary synthetic",
        state_labels=("0", "1"),
        root_prior=numpy.array([0.5, 0.5], dtype=float),
        root_prior_source="provided",
        observation_policy="reject",
        transition_matrix_for_branch_length=lambda _branch_length: numpy.array(
            [[0.9, 0.1], [0.2, 0.8]],
            dtype=float,
        ),
        leaf_likelihood_vector_for_state=_binary_leaf_vector,
    )


def _pattern_compiler(
    records: list[AlignmentRecord],
    source_path: Path | None,
) -> CompressedAlignmentSitePatterns:
    return compress_alignment_site_patterns_from_records(
        records,
        source_path=source_path,
    )


def _binary_leaf_vector(taxon_name: str, observed_state: str) -> numpy.ndarray:
    if observed_state == "0":
        return numpy.array([1.0, 0.0], dtype=float)
    if observed_state == "1":
        return numpy.array([0.0, 1.0], dtype=float)
    raise ValueError(f"unexpected synthetic state '{observed_state}' at {taxon_name}")
