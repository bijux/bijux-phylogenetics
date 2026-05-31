from __future__ import annotations

import bijux_phylogenetics.trees as trees_api
from bijux_phylogenetics.trees import (
    TreeSetCredibleCladeRow,
    TreeSetCredibleCladeSetReport,
    compute_credible_clade_set,
    write_credible_clade_set_artifacts,
    write_credible_clade_set_excluded_table,
    write_credible_clade_set_included_table,
)


def test_package_tree_gateway_exports_credible_clade_set_surface() -> None:
    assert trees_api.TreeSetCredibleCladeRow is TreeSetCredibleCladeRow
    assert trees_api.TreeSetCredibleCladeSetReport is TreeSetCredibleCladeSetReport
    assert trees_api.compute_credible_clade_set is compute_credible_clade_set
    assert (
        trees_api.write_credible_clade_set_artifacts
        is write_credible_clade_set_artifacts
    )
    assert (
        trees_api.write_credible_clade_set_excluded_table
        is write_credible_clade_set_excluded_table
    )
    assert (
        trees_api.write_credible_clade_set_included_table
        is write_credible_clade_set_included_table
    )
