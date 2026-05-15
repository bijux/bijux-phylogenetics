from __future__ import annotations

import csv
import json
import os
from pathlib import Path
import shutil
import subprocess

import pytest

from bijux_phylogenetics.ape_parity import (
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
        "dna-base-frequency-lowercase",
        "dna-base-frequency-ambiguity",
        "dna-raw-distance-clean",
        "dna-raw-distance-gaps",
        "dna-raw-distance-identical",
        "dna-raw-distance-high-divergence",
        "dna-raw-distance-missing-data",
        "dna-translation-valid-frame",
        "dna-translation-internal-stop",
        "dna-translation-terminal-stop",
    ]
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
        "lowercase_aligned_dna",
        "dna_with_ambiguity",
        "clean_aligned_dna",
        "dna_with_gaps",
        "identical_sequences",
        "high_divergence_sequences",
        "dna_with_missing_data",
        "coding_valid_reading_frame",
        "coding_internal_stop",
        "coding_terminal_stop",
    ]
    assert {case.function_name for case in cases} == {
        "ape::read.tree",
        "ape::write.tree",
        "ape::root",
        "ape::unroot",
        "ape::drop.tip",
        "ape::keep.tip",
        "ape::getMRCA",
        "ape::base.freq",
        "ape::dist.dna",
        "ape::trans",
        "ape::extract.clade",
    }
    assert {case.operation for case in cases} == {
        "read-tree-structure",
        "read-tree-set-structure",
        "write-tree-structure",
        "write-tree-set-structure",
        "root-tree-outgroup",
        "unroot-tree",
        "drop-tree-taxa",
        "keep-tree-taxa",
        "extract-tree-clade",
        "get-tree-mrca",
        "dna-base-frequency",
        "dna-raw-distance",
        "dna-translation",
    }


def test_run_ape_parity_cases_passes_against_fake_reference_runner(
    tmp_path: Path,
) -> None:
    rscript = fake_ape_rscript(tmp_path / "fake-ape-rscript")

    report = run_ape_parity_cases(
        rscript_executable=str(rscript),
        failure_root=tmp_path / "ape-parity-failures",
    )

    assert report.all_passed is True
    assert report.case_count == 56
    assert report.passed_case_count == 56
    assert report.failed_case_count == 0
    assert report.skipped_case_count == 0
    assert [row.function_name for row in report.summary_rows] == [
        "ape::base.freq",
        "ape::dist.dna",
        "ape::drop.tip",
        "ape::extract.clade",
        "ape::getMRCA",
        "ape::keep.tip",
        "ape::read.tree",
        "ape::root",
        "ape::trans",
        "ape::unroot",
        "ape::write.tree",
    ]
    assert all(observation.r_version == "4.6.0" for observation in report.observations)
    assert all(observation.ape_version == "5.0.0" for observation in report.observations)
    assert all(observation.reproducible_artifact_root is None for observation in report.observations)
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
    translation_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "dna-translation-terminal-stop"
    )
    assert translation_case.reference_summary is not None
    assert translation_case.reference_summary["stop_codon_count"] == 1


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


def test_write_ape_parity_tables_writes_summary_and_observations(tmp_path: Path) -> None:
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
    assert len(rows) == 56
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
