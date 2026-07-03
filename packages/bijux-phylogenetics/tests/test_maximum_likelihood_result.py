from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.io.newick import loads_newick
from bijux_phylogenetics.phylo.likelihood.bootstrap_inference import (
    bootstrap_nucleotide_likelihood_tree_inference_from_alignment,
)
from bijux_phylogenetics.phylo.likelihood.maximum_likelihood_result import (
    NucleotideMaximumLikelihoodResult,
    build_nucleotide_maximum_likelihood_benchmark_metadata,
    build_nucleotide_maximum_likelihood_result,
    build_nucleotide_maximum_likelihood_support_reports,
    infer_nucleotide_maximum_likelihood_result_from_alignment,
    load_nucleotide_maximum_likelihood_result_json,
    write_nucleotide_maximum_likelihood_result_json,
)
from bijux_phylogenetics.phylo.likelihood.models import (
    LikelihoodWrapperCorrespondenceObservation,
    LikelihoodWrapperCorrespondenceReport,
    LikelihoodWrapperCorrespondenceSummaryRow,
)
from bijux_phylogenetics.phylo.likelihood.sh_like_branch_support import (
    evaluate_nucleotide_sh_like_branch_support,
)
from bijux_phylogenetics.phylo.likelihood.tree_inference import (
    infer_nucleotide_likelihood_tree_from_alignment,
)

pytestmark = pytest.mark.slow

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_maximum_likelihood_result_round_trips_inference_payload(
    tmp_path: Path,
) -> None:
    inference_report = infer_nucleotide_likelihood_tree_from_alignment(
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        search_method="nni",
        start_tree_count=3,
        start_tree_seed=17,
        upper_branch_length_bound=1.0,
    )

    result = build_nucleotide_maximum_likelihood_result(inference_report)
    written_path = write_nucleotide_maximum_likelihood_result_json(
        tmp_path / "maximum-likelihood-result.json",
        result,
    )
    reloaded_result = load_nucleotide_maximum_likelihood_result_json(written_path)

    assert isinstance(reloaded_result, NucleotideMaximumLikelihoodResult)
    assert reloaded_result.to_dict() == result.to_dict()
    assert reloaded_result.model_name == "JC69"
    assert reloaded_result.final_tree_newick == inference_report.best_final_tree_newick
    assert (
        reloaded_result.final_log_likelihood
        == inference_report.best_final_log_likelihood
    )
    assert reloaded_result.parameter_values == (
        inference_report.model_selection_report.rows[0].parameter_values
    )
    assert reloaded_result.search_trace_rows
    assert reloaded_result.support_reports.bootstrap is None
    assert reloaded_result.support_reports.sh_like is None


def test_maximum_likelihood_result_accepts_matching_support_reports() -> None:
    alignment_path = fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta")
    inference_report = infer_nucleotide_likelihood_tree_from_alignment(
        alignment_path,
        model_name="jc69",
        search_method="nni",
        start_tree_count=3,
        start_tree_seed=17,
        upper_branch_length_bound=1.0,
    )

    bootstrap_report = bootstrap_nucleotide_likelihood_tree_inference_from_alignment(
        alignment_path,
        model_name="jc69",
        search_method="nni",
        start_tree_count=3,
        start_tree_seed=17,
        replicate_count=3,
        bootstrap_seed=5,
        upper_branch_length_bound=1.0,
    )
    sh_like_report = evaluate_nucleotide_sh_like_branch_support(
        loads_newick(inference_report.best_final_tree_newick),
        alignment_path,
        model_name="jc69",
        resampling_replicate_count=3,
        resampling_seed=11,
    )

    support_reports = build_nucleotide_maximum_likelihood_support_reports(
        inference_report,
        bootstrap_report=bootstrap_report,
        sh_like_report=sh_like_report,
    )
    result = build_nucleotide_maximum_likelihood_result(
        inference_report,
        support_reports=support_reports,
    )

    assert result.support_reports.bootstrap == bootstrap_report
    assert result.support_reports.sh_like == sh_like_report
    assert sh_like_report.caution_label in result.warning_messages


def test_maximum_likelihood_result_builds_benchmark_metadata_summary() -> None:
    correspondence_report = LikelihoodWrapperCorrespondenceReport(
        observations=[
            LikelihoodWrapperCorrespondenceObservation(
                case_id="iqtree-jc69-reference",
                wrapper_engine="iqtree",
                native_surface="native-ml-tree-inference",
                wrapper_surface="iqtree-ml-tree-inference",
                comparison_policy="fixture-reference",
                status="exact-match",
                supported=True,
                blocking=False,
                tolerance=None,
                rationale="native and wrapper topologies agree on the governed fixture",
                input_fixtures=["alignments/jc69_reference.fasta"],
                expected_output={"topology": "match"},
                observed_output={"topology": "match"},
            )
        ],
        summary_rows=[
            LikelihoodWrapperCorrespondenceSummaryRow(
                status="exact-match",
                case_count=1,
                blocking_case_count=0,
                case_ids=["iqtree-jc69-reference"],
            )
        ],
        case_count=1,
        supported_case_count=1,
        exact_match_case_count=1,
        tolerance_match_case_count=0,
        expected_model_assumption_difference_case_count=0,
        unsupported_case_count=0,
        native_bug_case_count=0,
        blocking_case_count=0,
        all_supported_cases_clear=True,
    )

    metadata = build_nucleotide_maximum_likelihood_benchmark_metadata(
        correspondence_report
    )

    assert metadata.wrapper_correspondence_case_count == 1
    assert metadata.wrapper_correspondence_supported_case_count == 1
    assert metadata.wrapper_correspondence_blocking_case_count == 0
    assert metadata.wrapper_correspondence_all_supported_cases_clear is True
    assert (
        metadata.wrapper_correspondence_summary_rows
        == correspondence_report.summary_rows
    )


def test_infer_nucleotide_maximum_likelihood_result_from_alignment_can_attach_support() -> (
    None
):
    result = infer_nucleotide_maximum_likelihood_result_from_alignment(
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        search_method="nni",
        start_tree_count=3,
        start_tree_seed=17,
        bootstrap_replicate_count=3,
        bootstrap_seed=5,
        sh_like_resampling_replicate_count=3,
        sh_like_resampling_seed=11,
        upper_branch_length_bound=1.0,
    )

    assert result.model_name == "JC69"
    assert result.support_reports.bootstrap is not None
    assert result.support_reports.bootstrap.replicate_count == 3
    assert result.support_reports.sh_like is not None
    assert result.support_reports.sh_like.resampling_replicate_count == 3


def test_infer_nucleotide_maximum_likelihood_result_can_attach_benchmark_metadata() -> (
    None
):
    result = infer_nucleotide_maximum_likelihood_result_from_alignment(
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        search_method="nni",
        start_tree_count=3,
        start_tree_seed=17,
        include_wrapper_correspondence_benchmark=True,
        upper_branch_length_bound=1.0,
    )

    assert result.benchmark_metadata is not None
    assert result.benchmark_metadata.wrapper_correspondence_case_count >= 1
