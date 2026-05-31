from __future__ import annotations

from bijux_phylogenetics.io.newick import loads_newick
import bijux_phylogenetics.phylo.topology as topology_api
from bijux_phylogenetics.phylo.topology import (
    AffectedSubtreeReport,
    summarize_affected_subtrees,
)
from bijux_phylogenetics.phylo.topology.rooted_nni import (
    apply_rooted_nni_move,
    resolve_rooted_nni_move_candidate,
)


def test_topology_gateway_exports_affected_subtree_surface() -> None:
    assert topology_api.AffectedSubtreeReport is AffectedSubtreeReport
    assert topology_api.summarize_affected_subtrees is summarize_affected_subtrees


def test_affected_subtree_summary_reports_exact_changed_and_unchanged_nni_branches() -> (
    None
):
    original_tree = loads_newick("(((A,C),B),D);")
    candidate, _available_move_count = resolve_rooted_nni_move_candidate(original_tree, 1)
    moved_tree = apply_rooted_nni_move(original_tree, candidate)

    report = summarize_affected_subtrees(original_tree, moved_tree)

    assert report.original_branch_clade_ids == ["A", "B", "C", "D", "A|C", "A|B|C"]
    assert report.moved_branch_clade_ids == ["A", "B", "C", "D", "A|C", "A|C|D"]
    assert report.retired_branch_clade_ids == ["A|B|C"]
    assert report.introduced_branch_clade_ids == ["A|C|D"]
    assert report.affected_branch_clade_ids == ["A|B|C", "A|C|D"]
    assert report.unaffected_branch_clade_ids == ["A", "B", "C", "D", "A|C"]
