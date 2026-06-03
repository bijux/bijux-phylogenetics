from __future__ import annotations

from bijux_phylogenetics.ancestral.discrete import DiscreteAncestralReport
from bijux_phylogenetics.ancestral.discrete.state_resolution import (
    resolve_clade_consensus_state,
)
from bijux_phylogenetics.io.newick import loads_newick
from bijux_phylogenetics.phylo.topology.tree import TreeNode

from .contracts import AncestralTransitionBranchRow


def _build_transition_branch_rows(
    report: DiscreteAncestralReport,
    *,
    observed_states_by_taxon: dict[str, str] | None = None,
) -> list[AncestralTransitionBranchRow]:
    tree = loads_newick(report.analysis_tree_newick)
    estimate_by_node = {estimate.node: estimate for estimate in report.estimates}
    branch_rows: list[AncestralTransitionBranchRow] = []

    def visit(node: TreeNode) -> None:
        parent_estimate = estimate_by_node[_node_signature(node)]
        for child in node.children:
            child_estimate = estimate_by_node[_node_signature(child)]
            parent_state = _resolve_transition_state(
                estimate=parent_estimate,
                observed_states_by_taxon=observed_states_by_taxon,
            )
            child_state = _resolve_transition_state(
                estimate=child_estimate,
                observed_states_by_taxon=observed_states_by_taxon,
            )
            overlapping_states = sorted(
                set(parent_estimate.state_set) & set(child_estimate.state_set)
            )
            changed = parent_state != child_state
            branch_rows.append(
                AncestralTransitionBranchRow(
                    parent_node=parent_estimate.node,
                    child_node=child_estimate.node,
                    child_descendant_taxa=child_estimate.descendant_taxa,
                    branch_length=child.branch_length,
                    parent_most_likely_state=parent_state,
                    child_most_likely_state=child_state,
                    parent_state_set=parent_estimate.state_set,
                    child_state_set=child_estimate.state_set,
                    overlapping_states=overlapping_states,
                    changed=changed,
                    transition=(f"{parent_state}->{child_state}" if changed else ""),
                    certainty_class=_transition_certainty_class(
                        changed=changed,
                        overlapping_states=overlapping_states,
                        parent_state_set=parent_estimate.state_set,
                        child_state_set=child_estimate.state_set,
                    ),
                )
            )
            visit(child)

    visit(tree.root)
    return branch_rows


def _resolve_transition_state(
    *,
    estimate,
    observed_states_by_taxon: dict[str, str] | None,
) -> str:
    if observed_states_by_taxon is None:
        return estimate.most_likely_state
    return resolve_clade_consensus_state(
        clade_taxa=estimate.descendant_taxa,
        candidate_states=estimate.state_set,
        observed_states_by_taxon=observed_states_by_taxon,
        fallback_state=estimate.most_likely_state,
    )


def _transition_certainty_class(
    *,
    changed: bool,
    overlapping_states: list[str],
    parent_state_set: list[str],
    child_state_set: list[str],
) -> str:
    if changed:
        if not overlapping_states:
            return "certain_change"
        return "uncertain_change"
    if len(parent_state_set) == 1 and len(child_state_set) == 1:
        return "certain_no_change"
    return "uncertain_no_change"


def _node_signature(node: TreeNode) -> str:
    if node.is_leaf():
        return node.name or "<unnamed>"
    descendant_taxa: list[str] = []
    for child in node.children:
        descendant_taxa.extend(_node_signature_taxa(child))
    return "|".join(sorted(descendant_taxa))


def _node_signature_taxa(node: TreeNode) -> list[str]:
    if node.is_leaf():
        return [node.name] if node.name is not None else []
    descendant_taxa: list[str] = []
    for child in node.children:
        descendant_taxa.extend(_node_signature_taxa(child))
    return descendant_taxa
