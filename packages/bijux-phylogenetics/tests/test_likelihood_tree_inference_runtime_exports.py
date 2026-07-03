from __future__ import annotations

import bijux_phylogenetics.phylo.likelihood as likelihood_api
from bijux_phylogenetics.phylo.likelihood import (
    NucleotideLikelihoodTreeInferenceReport,
    NucleotideMaximumLikelihoodBenchmarkMetadata,
    NucleotideMaximumLikelihoodResult,
    NucleotideMaximumLikelihoodSupportReports,
    build_nucleotide_maximum_likelihood_benchmark_metadata,
    build_nucleotide_maximum_likelihood_result,
    build_nucleotide_maximum_likelihood_support_reports,
    default_nucleotide_likelihood_tree_inference_model_candidates,
    infer_nucleotide_likelihood_tree,
    infer_nucleotide_likelihood_tree_from_alignment,
    infer_nucleotide_maximum_likelihood_result,
    infer_nucleotide_maximum_likelihood_result_from_alignment,
    load_nucleotide_maximum_likelihood_result_json,
    validate_nucleotide_likelihood_tree_inference_model_name,
    validate_nucleotide_likelihood_tree_inference_model_selection_criterion,
    validate_nucleotide_likelihood_tree_inference_search_method,
    validate_nucleotide_likelihood_tree_inference_start_tree_count,
    write_nucleotide_likelihood_tree_inference_artifacts,
    write_nucleotide_likelihood_tree_inference_likelihood_table,
    write_nucleotide_likelihood_tree_inference_model_table,
    write_nucleotide_likelihood_tree_inference_run_json,
    write_nucleotide_maximum_likelihood_result_json,
)


def test_public_runtime_exports_tree_inference_surface() -> None:
    assert (
        likelihood_api.NucleotideMaximumLikelihoodBenchmarkMetadata
        is NucleotideMaximumLikelihoodBenchmarkMetadata
    )
    assert (
        likelihood_api.NucleotideMaximumLikelihoodResult
        is NucleotideMaximumLikelihoodResult
    )
    assert (
        likelihood_api.NucleotideMaximumLikelihoodSupportReports
        is NucleotideMaximumLikelihoodSupportReports
    )
    assert (
        likelihood_api.NucleotideLikelihoodTreeInferenceReport
        is NucleotideLikelihoodTreeInferenceReport
    )
    assert (
        likelihood_api.build_nucleotide_maximum_likelihood_benchmark_metadata
        is build_nucleotide_maximum_likelihood_benchmark_metadata
    )
    assert (
        likelihood_api.build_nucleotide_maximum_likelihood_result
        is build_nucleotide_maximum_likelihood_result
    )
    assert (
        likelihood_api.build_nucleotide_maximum_likelihood_support_reports
        is build_nucleotide_maximum_likelihood_support_reports
    )
    assert (
        likelihood_api.default_nucleotide_likelihood_tree_inference_model_candidates
        is default_nucleotide_likelihood_tree_inference_model_candidates
    )
    assert (
        likelihood_api.infer_nucleotide_maximum_likelihood_result
        is infer_nucleotide_maximum_likelihood_result
    )
    assert (
        likelihood_api.infer_nucleotide_maximum_likelihood_result_from_alignment
        is infer_nucleotide_maximum_likelihood_result_from_alignment
    )
    assert (
        likelihood_api.infer_nucleotide_likelihood_tree
        is infer_nucleotide_likelihood_tree
    )
    assert (
        likelihood_api.infer_nucleotide_likelihood_tree_from_alignment
        is infer_nucleotide_likelihood_tree_from_alignment
    )
    assert (
        likelihood_api.validate_nucleotide_likelihood_tree_inference_model_name
        is validate_nucleotide_likelihood_tree_inference_model_name
    )
    assert (
        likelihood_api.validate_nucleotide_likelihood_tree_inference_model_selection_criterion
        is validate_nucleotide_likelihood_tree_inference_model_selection_criterion
    )
    assert (
        likelihood_api.validate_nucleotide_likelihood_tree_inference_search_method
        is validate_nucleotide_likelihood_tree_inference_search_method
    )
    assert (
        likelihood_api.validate_nucleotide_likelihood_tree_inference_start_tree_count
        is validate_nucleotide_likelihood_tree_inference_start_tree_count
    )
    assert (
        likelihood_api.load_nucleotide_maximum_likelihood_result_json
        is load_nucleotide_maximum_likelihood_result_json
    )
    assert (
        likelihood_api.write_nucleotide_likelihood_tree_inference_artifacts
        is write_nucleotide_likelihood_tree_inference_artifacts
    )
    assert (
        likelihood_api.write_nucleotide_likelihood_tree_inference_likelihood_table
        is write_nucleotide_likelihood_tree_inference_likelihood_table
    )
    assert (
        likelihood_api.write_nucleotide_likelihood_tree_inference_model_table
        is write_nucleotide_likelihood_tree_inference_model_table
    )
    assert (
        likelihood_api.write_nucleotide_likelihood_tree_inference_run_json
        is write_nucleotide_likelihood_tree_inference_run_json
    )
    assert (
        likelihood_api.write_nucleotide_maximum_likelihood_result_json
        is write_nucleotide_maximum_likelihood_result_json
    )
