from __future__ import annotations

import bijux_phylogenetics.trees as trees_api
from bijux_phylogenetics.trees import (
    PosteriorAgreementSubtreeCandidateRow,
    PosteriorAgreementSubtreeReport,
    summarize_posterior_agreement_subtree,
    write_posterior_agreement_subtree_artifacts,
    write_posterior_agreement_subtree_removed_taxa_table,
    write_posterior_agreement_subtree_search_table,
    write_posterior_agreement_subtree_summary_table,
)


def test_package_tree_gateway_exports_posterior_agreement_subtree_surface() -> None:
    assert (
        trees_api.PosteriorAgreementSubtreeCandidateRow
        is PosteriorAgreementSubtreeCandidateRow
    )
    assert trees_api.PosteriorAgreementSubtreeReport is PosteriorAgreementSubtreeReport
    assert (
        trees_api.summarize_posterior_agreement_subtree
        is summarize_posterior_agreement_subtree
    )
    assert (
        trees_api.write_posterior_agreement_subtree_artifacts
        is write_posterior_agreement_subtree_artifacts
    )
    assert (
        trees_api.write_posterior_agreement_subtree_removed_taxa_table
        is write_posterior_agreement_subtree_removed_taxa_table
    )
    assert (
        trees_api.write_posterior_agreement_subtree_search_table
        is write_posterior_agreement_subtree_search_table
    )
    assert (
        trees_api.write_posterior_agreement_subtree_summary_table
        is write_posterior_agreement_subtree_summary_table
    )
