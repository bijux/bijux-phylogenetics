from __future__ import annotations

import csv
import json
import os
from pathlib import Path
import shutil
import subprocess

import pytest

from bijux_phylogenetics.parity import (
    list_ape_parity_cases,
    run_ape_parity_cases,
    write_ape_parity_observation_table,
    write_ape_parity_summary_table,
)
from tests.support.fake_reference_parity import fake_ape_rscript


def _r_package_available(rscript: str, package_name: str) -> bool:
    repository_root = Path(__file__).resolve().parents[3]
    environment = dict(os.environ)
    r_library = repository_root / "artifacts" / "r-lib"
    if r_library.is_dir():
        environment["R_LIBS_USER"] = str(r_library)
    result = subprocess.run(
        [
            rscript,
            "-e",
            f"cat(requireNamespace('{package_name}', quietly=TRUE), '\\n')",
        ],
        capture_output=True,
        check=False,
        cwd=repository_root,
        env=environment,
        text=True,
    )
    return result.returncode == 0 and result.stdout.strip() == "TRUE"


def test_list_ape_parity_cases_returns_governed_read_tree_registry() -> None:
    cases = list_ape_parity_cases()

    assert [case.case_id for case in cases] == [
        "read-tree-balanced-rooted-ultrametric",
        "read-tree-unrooted-branch-length",
        "read-tree-internal-node-labels",
        "read-tree-support-labels",
        "read-tree-quoted-taxon-labels",
        "read-tree-multiple-trees",
        "read-tree-malformed-newick",
        "write-tree-balanced-rooted-ultrametric",
        "write-tree-unrooted-branch-length",
        "write-tree-internal-node-labels",
        "write-tree-support-labels",
        "write-tree-quoted-taxon-labels",
        "write-tree-multiple-trees",
        "consensus-majority-conflicting-four-taxon",
        "consensus-strict-conflicting-four-taxon",
        "consensus-majority-posterior-six-taxon",
        "consensus-mismatched-taxon-set",
        "prop-clades-duplicate-conflict-four-taxon",
        "prop-clades-absent-cross-pairing-clades",
        "prop-clades-child-order-insensitive",
        "prop-clades-posterior-six-taxon",
        "prop-clades-mismatched-taxon-set",
        "root-tree-single-outgroup-tip",
        "root-tree-multiple-outgroup-tips",
        "root-tree-already-rooted",
        "root-tree-missing-outgroup",
        "root-tree-non-monophyletic-outgroup",
        "unroot-tree-balanced-rooted",
        "unroot-tree-rootable",
        "unroot-tree-after-outgroup-rooting",
        "unroot-tree-already-unrooted",
        "unroot-tree-invalid-newick",
        "drop-tip-rooted-single",
        "drop-tip-rooted-multiple",
        "drop-tip-root-change-after-outgroup-rooting",
        "drop-tip-unrooted-three-tip",
        "drop-tip-unrooted-two-tip",
        "drop-tip-unknown-tip-name",
        "keep-tip-rooted-selected-two",
        "keep-tip-rooted-order-insensitive",
        "keep-tip-root-change-after-outgroup-rooting",
        "keep-tip-unrooted-three-tip",
        "keep-tip-unrooted-two-tip",
        "extract-clade-root",
        "extract-clade-mammals",
        "extract-clade-birds",
        "extract-clade-tip-node-invalid",
        "extract-clade-node-out-of-bounds",
        "get-mrca-balanced-two-tip",
        "get-mrca-balanced-full-tip-set",
        "get-mrca-balanced-duplicate-request",
        "get-mrca-pectinate-many-tip",
        "get-mrca-rooted-polytomy",
        "get-mrca-after-outgroup-rooting",
        "get-mrca-missing-tip",
        "is-monophyletic-rooted-two-tip",
        "is-monophyletic-rooted-three-tip-reroot-false",
        "is-monophyletic-rooted-three-tip-reroot-true",
        "is-monophyletic-rooted-full-tip-set",
        "is-monophyletic-rooted-mixed-missing",
        "is-monophyletic-unrooted-two-tip",
        "is-monophyletic-unrooted-three-tip",
        "is-monophyletic-after-outgroup-rooting",
        "is-monophyletic-rooted-polytomy",
        "is-monophyletic-all-missing-rerooted",
        "cophenetic-rooted-ultrametric",
        "cophenetic-unrooted-branch-length",
        "dist-topo-rooted-identical",
        "dist-topo-rooted-child-order",
        "dist-topo-rooted-conflict",
        "dist-topo-rooted-polytomy",
        "dist-topo-unrooted-conflict",
        "dist-topo-rooted-large",
        "vcv-rooted-ultrametric",
        "vcv-rooted-non-ultrametric",
        "vcv-unrooted-branch-length",
        "vcv-zero-branch-singular",
        "ace-continuous-balanced-rooted-ultrametric",
        "ace-continuous-pectinate-non-ultrametric",
        "ace-continuous-balanced-six-taxon",
        "ace-continuous-missing-values-pruned",
        "ace-discrete-binary-balanced-rooted-ultrametric",
        "ace-discrete-multistate-balanced-rooted-ultrametric",
        "ace-discrete-multistate-pectinate-non-ultrametric",
        "ace-discrete-missing-values-pruned",
        "ace-discrete-sym-balanced-rooted-ultrametric",
        "ace-discrete-sym-pectinate-non-ultrametric",
        "ace-discrete-sym-balanced-six-taxon",
        "ace-discrete-sym-missing-values-pruned",
        "ace-discrete-ard-binary-balanced-rooted-ultrametric",
        "ace-discrete-ard-multistate-pectinate-non-ultrametric",
        "ace-discrete-ard-balanced-six-taxon",
        "ace-discrete-ard-missing-values-pruned",
        "pic-balanced-rooted-ultrametric",
        "pic-pectinate-non-ultrametric",
        "pic-balanced-six-taxon",
        "node-depth-rooted-ultrametric",
        "node-depth-rooted-non-ultrametric",
        "node-depth-zero-branch-lengths",
        "node-depth-after-outgroup-rooting",
        "branching-times-rooted-ultrametric",
        "branching-times-internal-node-labels",
        "branching-times-medium-ultrametric",
        "branching-times-zero-internal-branch",
        "gamma-stat-rooted-ultrametric",
        "gamma-stat-internal-node-labels",
        "gamma-stat-medium-ultrametric",
        "gamma-stat-zero-internal-branch",
        "rtree-rooted-six-taxon-uniform",
        "rtree-rooted-twelve-taxon-uniform",
        "rcoal-rooted-six-taxon",
        "rcoal-rooted-twelve-taxon",
        "is-ultrametric-rooted-ultrametric",
        "is-ultrametric-near-ultrametric-default",
        "is-ultrametric-near-ultrametric-tight",
        "is-ultrametric-non-ultrametric",
        "neighbor-joining-analytical-three-taxon",
        "neighbor-joining-ultrametric-four-taxon",
        "neighbor-joining-nonultrametric-four-taxon",
        "dna-dnabin-structure-clean",
        "dna-dnabin-structure-lowercase",
        "dna-dnabin-structure-gaps",
        "dna-dnabin-structure-ambiguity",
        "dna-base-frequency-lowercase",
        "dna-base-frequency-ambiguity",
        "dna-base-frequency-missing-data",
        "dna-base-frequency-all-gap-missing",
        "dna-segregating-sites-lowercase",
        "dna-segregating-sites-invariant",
        "dna-segregating-sites-one-variable",
        "dna-segregating-sites-gaps",
        "dna-segregating-sites-ambiguity",
        "dna-segregating-sites-missing-data",
        "dna-segregating-sites-all-gap-missing",
        "dna-raw-distance-clean",
        "dna-raw-distance-gaps",
        "dna-raw-distance-gaps-complete-deletion",
        "dna-raw-distance-ambiguity",
        "dna-raw-distance-identical",
        "dna-raw-distance-high-divergence",
        "dna-raw-distance-missing-data",
        "dna-raw-distance-unequal-length-invalid",
        "dna-jc69-distance-clean",
        "dna-jc69-distance-gaps",
        "dna-jc69-distance-gaps-complete-deletion",
        "dna-jc69-distance-ambiguity",
        "dna-jc69-distance-identical",
        "dna-jc69-distance-high-divergence",
        "dna-jc69-distance-missing-data",
        "dna-jc69-distance-unequal-length-invalid",
        "dna-k80-distance-clean",
        "dna-k80-distance-gaps",
        "dna-k80-distance-gaps-complete-deletion",
        "dna-k80-distance-ambiguity",
        "dna-k80-distance-identical",
        "dna-k80-distance-high-divergence",
        "dna-k80-distance-missing-data",
        "dna-k80-distance-unequal-length-invalid",
        "dna-f81-distance-clean",
        "dna-f81-distance-gaps",
        "dna-f81-distance-gaps-complete-deletion",
        "dna-f81-distance-ambiguity",
        "dna-f81-distance-identical",
        "dna-f81-distance-high-divergence",
        "dna-f81-distance-missing-data",
        "dna-f81-distance-unequal-length-invalid",
        "dna-tn93-distance-clean",
        "dna-tn93-distance-gaps",
        "dna-tn93-distance-gaps-complete-deletion",
        "dna-tn93-distance-ambiguity",
        "dna-tn93-distance-identical",
        "dna-tn93-distance-high-divergence",
        "dna-tn93-distance-missing-data",
        "dna-tn93-distance-unequal-length-invalid",
        "dna-translation-valid-frame",
        "dna-translation-ambiguous-codon",
        "dna-translation-internal-stop",
        "dna-translation-terminal-stop",
        "dna-translation-frame-error-truncation",
        "dna-translation-vertebrate-mitochondrial",
    ]
    case_map = {case.case_id: case for case in cases}
    assert (
        case_map[
            "ace-discrete-sym-balanced-rooted-ultrametric"
        ].transition_rate_tolerance
        == 1.0
    )
    assert (
        case_map["ace-discrete-sym-pectinate-non-ultrametric"].transition_rate_tolerance
        == 1.0
    )
    assert (
        case_map[
            "ace-discrete-ard-multistate-pectinate-non-ultrametric"
        ].transition_rate_tolerance
        == 0.1
    )
    assert [case.fixture_id for case in cases] == [
        "balanced_rooted_ultrametric",
        "unrooted_branch_length_tree",
        "internal_node_labels",
        "branch_support_labels",
        "quoted_taxon_labels",
        "basic_newick_tree_set",
        "malformed_unbalanced_parentheses",
        "balanced_rooted_ultrametric",
        "unrooted_branch_length_tree",
        "internal_node_labels",
        "branch_support_labels",
        "quoted_taxon_labels",
        "basic_newick_tree_set",
        "consensus_conflicting_four_taxon_tree_set",
        "consensus_conflicting_four_taxon_tree_set",
        "consensus_posterior_six_taxon_tree_set",
        "consensus_mismatched_taxon_tree_set",
        "prop_clades_duplicate_conflict_tree_set",
        "prop_clades_absent_clade_tree_set",
        "prop_clades_child_order_tree_set",
        "prop_clades_posterior_six_taxon_tree_set",
        "prop_clades_mismatched_taxon_tree_set",
        "outgroup_rootable_unrooted",
        "outgroup_rootable_unrooted",
        "outgroup_rooted_on_d",
        "outgroup_rootable_unrooted",
        "outgroup_rootable_unrooted",
        "balanced_rooted_ultrametric",
        "outgroup_rootable_unrooted",
        "outgroup_rooted_on_d",
        "unrooted_branch_length_tree",
        "malformed_unbalanced_parentheses",
        "balanced_rooted_ultrametric",
        "balanced_rooted_ultrametric",
        "outgroup_rooted_on_d",
        "unrooted_branch_length_tree",
        "unrooted_branch_length_tree",
        "balanced_rooted_ultrametric",
        "balanced_rooted_ultrametric",
        "balanced_rooted_ultrametric",
        "outgroup_rooted_on_d",
        "unrooted_branch_length_tree",
        "unrooted_branch_length_tree",
        "internal_node_labels",
        "internal_node_labels",
        "internal_node_labels",
        "internal_node_labels",
        "internal_node_labels",
        "balanced_rooted_ultrametric",
        "balanced_rooted_ultrametric",
        "balanced_rooted_ultrametric",
        "pectinate_rooted_non_ultrametric",
        "rooted_polytomy",
        "outgroup_rooted_on_d",
        "balanced_rooted_ultrametric",
        "balanced_rooted_ultrametric",
        "balanced_rooted_ultrametric",
        "balanced_rooted_ultrametric",
        "balanced_rooted_ultrametric",
        "balanced_rooted_ultrametric",
        "unrooted_branch_length_tree",
        "unrooted_branch_length_tree",
        "outgroup_rooted_on_d",
        "rooted_polytomy",
        "unrooted_branch_length_tree",
        "balanced_rooted_ultrametric",
        "unrooted_branch_length_tree",
        "topology_distance_identical_rooted_pair",
        "topology_distance_rooted_child_order_pair",
        "topology_distance_rooted_conflict_pair",
        "topology_distance_rooted_polytomy_pair",
        "topology_distance_unrooted_conflict_pair",
        "topology_distance_large_rooted_pair",
        "balanced_rooted_ultrametric",
        "pectinate_rooted_non_ultrametric",
        "unrooted_branch_length_tree",
        "zero_branch_lengths",
        "balanced_rooted_ultrametric",
        "pectinate_rooted_non_ultrametric",
        "balanced_rooted_six_taxon",
        "balanced_rooted_six_taxon",
        "balanced_rooted_ultrametric",
        "balanced_rooted_ultrametric",
        "pectinate_rooted_non_ultrametric",
        "balanced_rooted_ultrametric",
        "balanced_rooted_ultrametric",
        "pectinate_rooted_non_ultrametric",
        "balanced_rooted_six_taxon",
        "balanced_rooted_six_taxon",
        "balanced_rooted_ultrametric",
        "pectinate_rooted_non_ultrametric",
        "balanced_rooted_six_taxon",
        "balanced_rooted_six_taxon",
        "balanced_rooted_ultrametric",
        "pectinate_rooted_non_ultrametric",
        "balanced_rooted_six_taxon",
        "balanced_rooted_ultrametric",
        "pectinate_rooted_non_ultrametric",
        "zero_branch_lengths",
        "outgroup_rooted_on_d",
        "balanced_rooted_ultrametric",
        "internal_node_labels",
        "larger_binary_tree",
        "ultrametric_zero_internal_branch",
        "balanced_rooted_ultrametric",
        "internal_node_labels",
        "larger_binary_tree",
        "ultrametric_zero_internal_branch",
        "rtree_rooted_six_taxon_uniform_64",
        "rtree_rooted_twelve_taxon_uniform_128",
        "rcoal_rooted_six_taxon_64",
        "rcoal_rooted_twelve_taxon_128",
        "balanced_rooted_ultrametric",
        "near_ultrametric_branch_jitter",
        "near_ultrametric_branch_jitter",
        "pectinate_rooted_non_ultrametric",
        "analytical_three_taxon",
        "ultrametric_four_taxon",
        "nonultrametric_four_taxon",
        "clean_aligned_dna",
        "lowercase_aligned_dna",
        "dna_with_gaps",
        "dna_with_ambiguity",
        "lowercase_aligned_dna",
        "dna_with_ambiguity",
        "dna_with_missing_data",
        "all_gap_missing_alignment",
        "lowercase_aligned_dna",
        "invariant_aligned_dna",
        "one_variable_site_alignment",
        "dna_with_gaps",
        "dna_with_ambiguity",
        "dna_with_missing_data",
        "all_gap_missing_alignment",
        "clean_aligned_dna",
        "dna_with_gaps",
        "dna_with_gaps",
        "dna_with_ambiguity",
        "identical_sequences",
        "high_divergence_sequences",
        "dna_with_missing_data",
        "unequal_length_invalid_input",
        "clean_aligned_dna",
        "dna_with_gaps",
        "dna_with_gaps",
        "dna_with_ambiguity",
        "identical_sequences",
        "high_divergence_sequences",
        "dna_with_missing_data",
        "unequal_length_invalid_input",
        "clean_aligned_dna",
        "dna_with_gaps",
        "dna_with_gaps",
        "dna_with_ambiguity",
        "identical_sequences",
        "high_divergence_sequences",
        "dna_with_missing_data",
        "unequal_length_invalid_input",
        "clean_aligned_dna",
        "dna_with_gaps",
        "dna_with_gaps",
        "dna_with_ambiguity",
        "identical_sequences",
        "high_divergence_sequences",
        "dna_with_missing_data",
        "unequal_length_invalid_input",
        "clean_aligned_dna",
        "dna_with_gaps",
        "dna_with_gaps",
        "dna_with_ambiguity",
        "identical_sequences",
        "high_divergence_sequences",
        "dna_with_missing_data",
        "unequal_length_invalid_input",
        "coding_valid_reading_frame",
        "coding_ambiguous_codon",
        "coding_internal_stop",
        "coding_terminal_stop",
        "coding_frame_error",
        "coding_mitochondrial_triplet",
    ]
    assert {case.function_name for case in cases} == {
        "ape::read.tree",
        "ape::write.tree",
        "ape::root",
        "ape::unroot",
        "ape::drop.tip",
        "ape::consensus",
        "ape::prop.clades",
        "ape::dist.topo",
        "ape::keep.tip",
        "ape::getMRCA",
        "ape::is.ultrametric",
        "ape::is.monophyletic",
        "ape::cophenetic.phylo",
        "ape::node.depth.edgelength",
        "ape::branching.times",
        "ape::gammaStat",
        "ape::rcoal",
        "ape::rtree",
        "ape::vcv.phylo",
        "ape::nj",
        "ape::pic",
        "ape::as.DNAbin",
        "ape::ace",
        "ape::base.freq",
        "ape::seg.sites",
        "ape::dist.dna",
        "ape::trans",
        "ape::extract.clade",
    }
    assert {case.operation for case in cases} == {
        "read-tree-structure",
        "read-tree-set-structure",
        "write-tree-structure",
        "write-tree-set-structure",
        "tree-consensus",
        "tree-clade-support",
        "root-tree-outgroup",
        "unroot-tree",
        "drop-tree-taxa",
        "keep-tree-taxa",
        "extract-tree-clade",
        "get-tree-mrca",
        "assess-tree-monophyly",
        "tree-tip-distance",
        "tree-topology-distance",
        "tree-brownian-covariance",
        "tree-continuous-ancestral-states",
        "tree-discrete-ancestral-states",
        "tree-independent-contrasts",
        "tree-node-depth",
        "tree-branching-times",
        "tree-diversification-gamma-statistic",
        "tree-simulation-envelope",
        "tree-ultrametricity",
        "distance-matrix-neighbor-joining",
        "dna-dnabin-structure",
        "dna-base-frequency",
        "dna-segregating-sites",
        "dna-distance",
        "dna-translation",
    }


@pytest.mark.slow
def test_run_ape_parity_cases_passes_against_fake_reference_runner(
    tmp_path: Path,
) -> None:
    rscript = fake_ape_rscript(tmp_path / "fake-ape-rscript")

    report = run_ape_parity_cases(
        rscript_executable=str(rscript),
        failure_root=tmp_path / "ape-parity-failures",
    )

    assert report.all_passed is True
    assert report.case_count == 180
    assert report.passed_case_count == 180
    assert report.failed_case_count == 0
    assert report.skipped_case_count == 0
    assert [row.function_name for row in report.summary_rows] == [
        "ape::ace",
        "ape::as.DNAbin",
        "ape::base.freq",
        "ape::branching.times",
        "ape::consensus",
        "ape::cophenetic.phylo",
        "ape::dist.dna",
        "ape::dist.topo",
        "ape::drop.tip",
        "ape::extract.clade",
        "ape::gammaStat",
        "ape::getMRCA",
        "ape::is.monophyletic",
        "ape::is.ultrametric",
        "ape::keep.tip",
        "ape::nj",
        "ape::node.depth.edgelength",
        "ape::pic",
        "ape::prop.clades",
        "ape::rcoal",
        "ape::read.tree",
        "ape::root",
        "ape::rtree",
        "ape::seg.sites",
        "ape::trans",
        "ape::unroot",
        "ape::vcv.phylo",
        "ape::write.tree",
    ]
    assert all(observation.r_version == "4.6.0" for observation in report.observations)
    assert all(
        observation.ape_version == "5.0.0" for observation in report.observations
    )
    assert all(
        observation.reproducible_artifact_root is None
        for observation in report.observations
    )
    continuous_missing_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "ace-continuous-missing-values-pruned"
    )
    assert continuous_missing_case.reference_summary is not None
    assert continuous_missing_case.reference_summary["dropped_missing_taxa"] == ["B"]
    assert continuous_missing_case.reference_summary["dropped_non_numeric_taxa"] == [
        "C"
    ]
    discrete_missing_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "ace-discrete-missing-values-pruned"
    )
    assert discrete_missing_case.reference_summary is not None
    assert discrete_missing_case.reference_summary["dropped_missing_taxa"] == ["D"]
    assert discrete_missing_case.reference_summary["transition_rate_rows"]
    discrete_sym_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "ace-discrete-sym-balanced-six-taxon"
    )
    assert discrete_sym_case.reference_summary is not None
    assert discrete_sym_case.reference_summary["baseline_model"] == "equal-rates"
    assert (
        discrete_sym_case.reference_summary["preferred_model_by_aic"] == "equal-rates"
    )
    discrete_ard_case = next(
        observation
        for observation in report.observations
        if observation.case_id
        == "ace-discrete-ard-multistate-pectinate-non-ultrametric"
    )
    assert discrete_ard_case.reference_summary is not None
    assert discrete_ard_case.reference_summary["model"] == "all-rates-different"
    assert discrete_ard_case.reference_summary["overparameterized"] is True
    assert discrete_ard_case.reference_summary["baseline_model"] == "equal-rates"
    internal_label_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "read-tree-internal-node-labels"
    )
    assert internal_label_case.reference_summary is not None
    assert internal_label_case.reference_summary["tree_count"] == 1
    support_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "read-tree-support-labels"
    )
    assert support_case.reference_summary is not None
    multiple_tree_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "read-tree-multiple-trees"
    )
    assert multiple_tree_case.reference_summary is not None
    assert multiple_tree_case.reference_summary["tree_count"] == 3
    malformed_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "read-tree-malformed-newick"
    )
    assert malformed_case.status == "passed"
    assert malformed_case.reference_error is not None
    assert malformed_case.bijux_error is not None
    write_multiple_tree_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "write-tree-multiple-trees"
    )
    assert write_multiple_tree_case.reference_summary is not None
    assert write_multiple_tree_case.reference_summary["tree_count"] == 3
    consensus_majority_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "consensus-majority-conflicting-four-taxon"
    )
    assert consensus_majority_case.reference_summary is not None
    assert (
        consensus_majority_case.reference_summary["consensus_method"] == "majority-rule"
    )
    assert consensus_majority_case.reference_summary[
        "consensus_threshold"
    ] == pytest.approx(0.5)
    assert consensus_majority_case.reference_summary["included_clade_count"] == 1
    consensus_error_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "consensus-mismatched-taxon-set"
    )
    assert consensus_error_case.status == "passed"
    assert consensus_error_case.reference_error is not None
    assert consensus_error_case.bijux_error is not None
    prop_clades_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "prop-clades-posterior-six-taxon"
    )
    assert prop_clades_case.reference_summary is not None
    assert prop_clades_case.reference_summary["supported_clade_count"] == 4
    assert prop_clades_case.reference_summary["absent_clade_count"] == 0
    assert prop_clades_case.reference_summary["unscored_clade_count"] == 0
    prop_clades_error_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "prop-clades-mismatched-taxon-set"
    )
    assert prop_clades_error_case.status == "passed"
    assert prop_clades_error_case.reference_error is not None
    assert prop_clades_error_case.bijux_error is not None
    quoted_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "write-tree-quoted-taxon-labels"
    )
    assert quoted_case.reference_summary is not None
    assert quoted_case.reference_summary["tip_labels"] == [
        "A.B-1",
        "Homo sapiens",
        "Mus musculus",
    ]
    root_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "root-tree-multiple-outgroup-tips"
    )
    assert root_case.reference_summary is not None
    assert root_case.reference_summary["rooted"] is True
    unroot_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "unroot-tree-already-unrooted"
    )
    assert unroot_case.reference_summary is not None
    assert unroot_case.reference_summary["rooted"] is False
    drop_tip_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "drop-tip-unrooted-two-tip"
    )
    assert drop_tip_case.reference_summary is not None
    assert drop_tip_case.reference_summary["rooted"] is True
    assert drop_tip_case.reference_summary["dropped_taxa"] == ["C", "D"]
    assert drop_tip_case.reference_summary["absent_requested_taxa"] == []
    keep_tip_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "keep-tip-rooted-order-insensitive"
    )
    assert keep_tip_case.reference_summary is not None
    assert keep_tip_case.reference_summary["requested_taxa"] == ["A", "C"]
    assert keep_tip_case.reference_summary["dropped_taxa"] == ["B", "D"]
    extract_clade_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "extract-clade-mammals"
    )
    assert extract_clade_case.reference_summary is not None
    assert extract_clade_case.reference_summary["requested_node_id"] == 6
    assert extract_clade_case.reference_summary["matched_node_id"] == 6
    assert extract_clade_case.reference_summary["matched_node_name"] == "Mammals"
    get_mrca_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "get-mrca-balanced-duplicate-request"
    )
    assert get_mrca_case.reference_summary is not None
    assert get_mrca_case.reference_summary["duplicate_requested_taxa"] == ["A"]
    assert get_mrca_case.reference_summary["matched_node_id"] == 6
    monophyly_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "is-monophyletic-rooted-three-tip-reroot-true"
    )
    assert monophyly_case.reference_summary is not None
    assert monophyly_case.reference_summary["monophyletic"] is True
    assert monophyly_case.reference_summary["complementary_clade_used"] is True
    assert monophyly_case.reference_summary["matched_extra_taxa"] == ["D"]
    cophenetic_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "cophenetic-rooted-ultrametric"
    )
    assert cophenetic_case.reference_summary is not None
    assert cophenetic_case.reference_summary["pair_count"] == 16
    assert cophenetic_case.reference_summary["symmetric"] is True
    topology_conflict_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "dist-topo-rooted-conflict"
    )
    assert topology_conflict_case.reference_summary is not None
    assert topology_conflict_case.reference_summary["robinson_foulds_distance"] == 2
    assert topology_conflict_case.reference_summary["shared_split_count"] == 1
    topology_large_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "dist-topo-rooted-large"
    )
    assert topology_large_case.reference_summary is not None
    assert topology_large_case.reference_summary["tip_count"] == 128
    assert topology_large_case.reference_summary["left_split_count"] == 126
    assert topology_large_case.reference_summary["robinson_foulds_distance"] == 24
    vcv_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "vcv-zero-branch-singular"
    )
    assert vcv_case.reference_summary is not None
    assert vcv_case.reference_summary["singular"] is True
    assert vcv_case.reference_summary["near_singular"] is True
    assert vcv_case.reference_summary["matrix_dimension"] == 4
    node_depth_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "node-depth-after-outgroup-rooting"
    )
    assert node_depth_case.reference_summary is not None
    assert node_depth_case.reference_summary["node_count"] == 7
    assert node_depth_case.reference_summary["rooted"] is True
    assert node_depth_case.reference_summary["zero_branch_length_count"] == 1
    branching_times_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "branching-times-zero-internal-branch"
    )
    assert branching_times_case.reference_summary is not None
    assert branching_times_case.reference_summary["internal_node_count"] == 3
    assert branching_times_case.reference_summary["root_age"] == pytest.approx(0.3)
    assert branching_times_case.reference_summary["zero_branch_length_count"] == 1
    gamma_stat_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "gamma-stat-zero-internal-branch"
    )
    assert gamma_stat_case.reference_summary is not None
    assert gamma_stat_case.reference_summary["tip_count"] == 4
    assert gamma_stat_case.reference_summary["bifurcating"] is True
    assert gamma_stat_case.reference_summary["gamma_statistic"] == pytest.approx(
        -0.979795897113271,
        abs=1e-12,
    )
    rtree_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "rtree-rooted-six-taxon-uniform"
    )
    assert rtree_case.reference_summary is not None
    assert rtree_case.reference_summary["branch_length_model"] == "uniform"
    assert rtree_case.reference_summary["pooled_branch_count"] == 640
    rcoal_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "rcoal-rooted-twelve-taxon"
    )
    assert rcoal_case.reference_summary is not None
    assert (
        rcoal_case.reference_summary["branch_length_model"]
        == "coalescent-waiting-times"
    )
    assert rcoal_case.reference_summary["pooled_branch_count"] == 2816
    ultrametric_default_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "is-ultrametric-near-ultrametric-default"
    )
    assert ultrametric_default_case.reference_summary is not None
    assert ultrametric_default_case.reference_summary["ultrametric"] is True
    assert ultrametric_default_case.reference_summary["offending_taxa"] == [
        "A",
        "B",
        "C",
        "D",
    ]
    ultrametric_tight_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "is-ultrametric-near-ultrametric-tight"
    )
    assert ultrametric_tight_case.reference_summary is not None
    assert ultrametric_tight_case.reference_summary["ultrametric"] is False
    assert ultrametric_tight_case.reference_summary[
        "max_tip_depth_deviation"
    ] == pytest.approx(
        1e-9,
        abs=1e-15,
    )
    ultrametric_non_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "is-ultrametric-non-ultrametric"
    )
    assert ultrametric_non_case.reference_summary is not None
    assert ultrametric_non_case.reference_summary["offending_taxa"] == ["A", "B", "D"]
    translation_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "dna-translation-terminal-stop"
    )
    assert translation_case.reference_summary is not None
    assert translation_case.reference_summary["stop_codon_count"] == 1
    frame_error_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "dna-translation-frame-error-truncation"
    )
    assert frame_error_case.reference_summary is not None
    assert frame_error_case.reference_summary["dropped_trailing_nucleotide_count"] == 2
    assert frame_error_case.reference_summary["warning_count"] == 1


def test_run_ape_parity_cases_records_branch_length_failure_for_tree_structure_mismatch(
    tmp_path: Path,
) -> None:
    rscript = fake_ape_rscript(
        tmp_path / "fake-ape-rscript",
        normalized_tree_overrides={
            "read-tree-balanced-rooted-ultrametric": "(A:0.2,B:0.2,(C:0.1,D:0.1):0.1);\n"
        },
    )

    report = run_ape_parity_cases(
        case_ids=["read-tree-balanced-rooted-ultrametric"],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "ape-parity-failures",
    )

    assert report.all_passed is False
    observation = report.observations[0]
    assert observation.status == "failed"
    assert observation.mismatch_reason is not None
    assert "branch lengths differ" in observation.mismatch_reason
    assert observation.reproducible_artifact_root is not None
    artifact_root = observation.reproducible_artifact_root
    assert artifact_root.exists()
    comparison_payload = json.loads(
        (artifact_root / "comparison.json").read_text(encoding="utf-8")
    )
    assert "branch lengths differ" in comparison_payload["mismatch_reason"]
    observed_summary = json.loads(
        (artifact_root / "reference-summary.observed.json").read_text(encoding="utf-8")
    )
    assert observed_summary["rooted"] is True
    assert (artifact_root / "bijux-normalized.txt").exists()


def test_run_ape_parity_cases_records_covariance_tables_for_vcv_summary_mismatch(
    tmp_path: Path,
) -> None:
    rscript = fake_ape_rscript(
        tmp_path / "fake-ape-rscript",
        summary_overrides={"matrix_rank": 99},
    )

    report = run_ape_parity_cases(
        case_ids=["vcv-rooted-ultrametric"],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "ape-parity-failures",
    )

    assert report.all_passed is False
    observation = report.observations[0]
    assert observation.status == "failed"
    assert observation.mismatch_reason == "summary_mismatch"
    assert observation.reproducible_artifact_root is not None
    artifact_root = observation.reproducible_artifact_root
    assert (artifact_root / "reference-covariance-matrix.tsv").exists()
    assert (artifact_root / "reference-covariance-long.tsv").exists()
    assert (artifact_root / "reference-rows.observed.tsv").exists()
    assert (artifact_root / "bijux-rows.tsv").exists()


def test_run_ape_parity_cases_passes_expected_rooting_errors_against_fake_reference_runner(
    tmp_path: Path,
) -> None:
    rscript = fake_ape_rscript(tmp_path / "fake-ape-rscript")

    report = run_ape_parity_cases(
        case_ids=[
            "root-tree-missing-outgroup",
            "root-tree-non-monophyletic-outgroup",
        ],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "ape-parity-failures",
    )

    assert report.all_passed is True
    assert [observation.status for observation in report.observations] == [
        "passed",
        "passed",
    ]


def test_run_ape_parity_cases_marks_missing_rscript_as_skipped(tmp_path: Path) -> None:
    report = run_ape_parity_cases(
        case_ids=["read-tree-balanced-rooted-ultrametric"],
        rscript_executable=str(tmp_path / "missing-rscript"),
        failure_root=tmp_path / "ape-parity-failures",
    )

    observation = report.observations[0]
    assert observation.status == "skipped"
    assert observation.mismatch_reason == "rscript_unavailable"
    assert observation.reproducible_artifact_root is not None
    assert observation.reproducible_artifact_root.exists()


@pytest.mark.slow
def test_write_ape_parity_tables_writes_summary_and_observations(
    tmp_path: Path,
) -> None:
    rscript = fake_ape_rscript(tmp_path / "fake-ape-rscript")
    report = run_ape_parity_cases(
        rscript_executable=str(rscript),
        failure_root=tmp_path / "ape-parity-failures",
    )
    summary_path = tmp_path / "ape-parity-summary.tsv"
    observation_path = tmp_path / "ape-parity-observations.tsv"

    write_ape_parity_summary_table(summary_path, report)
    write_ape_parity_observation_table(observation_path, report)

    summary_rows = summary_path.read_text(encoding="utf-8").splitlines()
    assert summary_rows[0] == (
        "function_name\tcase_count\tpassed_case_count\tfailed_case_count\tskipped_case_count"
    )
    with observation_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 180
    assert rows[0]["function_name"] == "ape::read.tree"
    assert rows[0]["fixture_kind"] == "tree"
    assert rows[0]["fixture_id"]
    assert rows[0]["status"] == "passed"
    assert rows[0]["bijux_version"]


def test_run_ape_parity_cases_records_live_environment_status(tmp_path: Path) -> None:
    rscript = shutil.which("Rscript")
    if rscript is None:
        pytest.skip("Rscript is not available")
    if not _r_package_available(rscript, "jsonlite"):
        pytest.skip("jsonlite is required for live ape parity validation")

    report = run_ape_parity_cases(
        case_ids=["read-tree-balanced-rooted-ultrametric"],
        rscript_executable=rscript,
        failure_root=tmp_path / "ape-parity-failures",
    )

    observation = report.observations[0]
    if _r_package_available(rscript, "ape"):
        assert observation.status == "passed"
        assert observation.ape_version
    else:
        assert observation.status == "skipped"
        assert observation.mismatch_reason == "ape_package_unavailable"
        assert observation.reproducible_artifact_root is not None
