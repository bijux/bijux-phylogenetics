from __future__ import annotations

import bijux_phylogenetics.phylo.likelihood as likelihood_api
from bijux_phylogenetics.phylo.likelihood import (
    LikelihoodWrapperCorrespondenceObservation,
    LikelihoodWrapperCorrespondenceReport,
    LikelihoodWrapperCorrespondenceSummaryRow,
    summarize_likelihood_wrapper_correspondence,
    write_likelihood_wrapper_correspondence_artifacts,
    write_likelihood_wrapper_correspondence_observation_table,
    write_likelihood_wrapper_correspondence_run_json,
    write_likelihood_wrapper_correspondence_summary_table,
)


def test_public_likelihood_exports_wrapper_correspondence_surface() -> None:
    assert (
        likelihood_api.LikelihoodWrapperCorrespondenceObservation
        is LikelihoodWrapperCorrespondenceObservation
    )
    assert (
        likelihood_api.LikelihoodWrapperCorrespondenceSummaryRow
        is LikelihoodWrapperCorrespondenceSummaryRow
    )
    assert (
        likelihood_api.LikelihoodWrapperCorrespondenceReport
        is LikelihoodWrapperCorrespondenceReport
    )
    assert (
        likelihood_api.summarize_likelihood_wrapper_correspondence
        is summarize_likelihood_wrapper_correspondence
    )
    assert (
        likelihood_api.write_likelihood_wrapper_correspondence_summary_table
        is write_likelihood_wrapper_correspondence_summary_table
    )
    assert (
        likelihood_api.write_likelihood_wrapper_correspondence_observation_table
        is write_likelihood_wrapper_correspondence_observation_table
    )
    assert (
        likelihood_api.write_likelihood_wrapper_correspondence_run_json
        is write_likelihood_wrapper_correspondence_run_json
    )
    assert (
        likelihood_api.write_likelihood_wrapper_correspondence_artifacts
        is write_likelihood_wrapper_correspondence_artifacts
    )
