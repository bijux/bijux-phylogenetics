from __future__ import annotations

import bijux_phylogenetics.trees as trees_api
from bijux_phylogenetics.trees import (
    TreeSetMaximumCladeCredibilityCandidateRow,
    TreeSetMaximumCladeCredibilityReport,
    compute_maximum_clade_credibility_tree,
    write_maximum_clade_credibility_artifacts,
    write_maximum_clade_credibility_score_table,
)


def test_package_tree_gateway_exports_maximum_clade_credibility_surface() -> None:
    assert (
        trees_api.TreeSetMaximumCladeCredibilityCandidateRow
        is TreeSetMaximumCladeCredibilityCandidateRow
    )
    assert (
        trees_api.TreeSetMaximumCladeCredibilityReport
        is TreeSetMaximumCladeCredibilityReport
    )
    assert (
        trees_api.compute_maximum_clade_credibility_tree
        is compute_maximum_clade_credibility_tree
    )
    assert (
        trees_api.write_maximum_clade_credibility_artifacts
        is write_maximum_clade_credibility_artifacts
    )
    assert (
        trees_api.write_maximum_clade_credibility_score_table
        is write_maximum_clade_credibility_score_table
    )
