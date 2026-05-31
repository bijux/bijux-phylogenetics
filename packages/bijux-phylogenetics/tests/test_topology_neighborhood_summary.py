from __future__ import annotations

from pathlib import Path

import bijux_phylogenetics.phylo.topology as topology_api
from bijux_phylogenetics.io.newick import loads_newick
from bijux_phylogenetics.phylo.topology import (
    RootedSprEnumerationBudget,
    TopologyNeighborhoodSummaryReport,
    enumerate_rooted_nni_neighbors,
    enumerate_rooted_spr_neighbors,
    enumerate_rooted_tbr_neighbors,
    summarize_topology_neighborhood,
    write_topology_neighborhood_summary_table,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_topology_gateway_exports_neighborhood_summary_surface() -> None:
    assert (
        topology_api.TopologyNeighborhoodSummaryReport
        is TopologyNeighborhoodSummaryReport
    )
    assert topology_api.summarize_topology_neighborhood is summarize_topology_neighborhood
    assert (
        topology_api.write_topology_neighborhood_summary_table
        is write_topology_neighborhood_summary_table
    )


def test_rooted_nni_neighborhood_summary_reports_comparable_counts() -> None:
    summary = summarize_topology_neighborhood(
        enumerate_rooted_nni_neighbors(loads_newick("(((A,C),B),D);"))
    )

    assert summary.neighborhood_family == "rooted-nni"
    assert summary.algorithm == "rooted-nni-neighbor-enumeration"
    assert summary.candidate_count == 4
    assert summary.valid_count == 4
    assert summary.duplicate_count == 0
    assert summary.skipped_count == 0
    assert summary.skipped_reason == "none"
    assert summary.budget_reason == "unbounded"


def test_rooted_spr_neighborhood_summary_reports_duplicate_and_skip_counts() -> None:
    summary = summarize_topology_neighborhood(
        enumerate_rooted_spr_neighbors(fixture("trees", "example_tree.nwk"))
    )

    assert summary.neighborhood_family == "rooted-spr"
    assert summary.algorithm == "rooted-spr-neighbor-enumeration"
    assert summary.candidate_count == 32
    assert summary.valid_count == 12
    assert summary.duplicate_count == 12
    assert summary.skipped_count == 8
    assert summary.skipped_reason == "identity-candidate:8"
    assert summary.budget_reason == "unbounded"


def test_rooted_spr_budgeted_summary_reports_budget_reason() -> None:
    summary = summarize_topology_neighborhood(
        enumerate_rooted_spr_neighbors(
            loads_newick("(((A,C),B),D);"),
            budget=RootedSprEnumerationBudget(
                max_pruned_clade_count=1,
                max_regraft_target_count_per_pruned_clade=3,
            ),
        )
    )

    assert summary.neighborhood_family == "rooted-spr"
    assert summary.candidate_count == 30
    assert summary.valid_count == 2
    assert summary.duplicate_count == 0
    assert summary.skipped_count == 28
    assert summary.skipped_reason == "identity-candidate:1;budget-candidate:27"
    assert (
        summary.budget_reason
        == "max-pruned-clades=1;max-regraft-targets-per-pruned-clade=3"
    )


def test_rooted_tbr_neighborhood_summary_reports_identity_skips() -> None:
    summary = summarize_topology_neighborhood(
        enumerate_rooted_tbr_neighbors(fixture("trees", "example_tree.nwk"))
    )

    assert summary.neighborhood_family == "rooted-tbr"
    assert summary.algorithm == "rooted-tbr-neighbor-enumeration"
    assert summary.candidate_count == 24
    assert summary.valid_count == 0
    assert summary.duplicate_count == 0
    assert summary.skipped_count == 24
    assert summary.skipped_reason == "identity-candidate:24"
    assert summary.budget_reason == "unbounded"


def test_write_topology_neighborhood_summary_table_materializes_comparable_rows(
    tmp_path: Path,
) -> None:
    reports = [
        summarize_topology_neighborhood(
            enumerate_rooted_nni_neighbors(loads_newick("(((A,C),B),D);"))
        ),
        summarize_topology_neighborhood(
            enumerate_rooted_spr_neighbors(fixture("trees", "example_tree.nwk"))
        ),
        summarize_topology_neighborhood(
            enumerate_rooted_tbr_neighbors(
                fixture("parsimony", "spr_search_start_tree_5_taxa.nwk")
            )
        ),
    ]

    output_path = write_topology_neighborhood_summary_table(
        tmp_path / "summary.tsv",
        reports,
    )

    lines = output_path.read_text(encoding="utf-8").splitlines()
    assert lines[0] == (
        "neighborhood_family\talgorithm\tcandidate_count\tvalid_count\t"
        "duplicate_count\tskipped_count\tskipped_reason\tbudget_reason"
    )
    assert lines[1].startswith("rooted-nni\trooted-nni-neighbor-enumeration\t4\t4\t0\t0\t")
    assert lines[2].startswith("rooted-spr\trooted-spr-neighbor-enumeration\t32\t12\t12\t8\t")
    assert lines[3].startswith("rooted-tbr\trooted-tbr-neighbor-enumeration\t52\t10\t36\t6\t")
