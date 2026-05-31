from __future__ import annotations

from dataclasses import dataclass
from functools import cache
from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import load_taxon_table
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode

from .models import (
    DuplicationLossTransferAssociationRow,
    DuplicationLossTransferEventRow,
    DuplicationLossTransferReport,
)

_UNREACHABLE_COST = float("inf")


@dataclass(frozen=True, slots=True)
class _ReconciliationChoice:
    event_type: str
    mapped_species_branch: str
    left_child_species_branch: str | None
    right_child_species_branch: str | None
    transferred_child_side: str | None
    transfer_recipient_branch: str | None
    loss_branches: tuple[str, ...]
    event_cost: float
    total_cost: float


@dataclass(frozen=True, slots=True)
class _SpeciesIndex:
    tree: PhyloTree
    nodes_by_branch: dict[str, TreeNode]
    descendants_by_branch: dict[str, set[str]]
    descendants_including_self_by_branch: dict[str, tuple[str, ...]]
    descendant_tips_by_branch: dict[str, list[str]]
    ancestry_by_branch: dict[str, tuple[str, ...]]
    root_branch: str


def reconcile_duplication_loss_transfer(
    species_tree_path: Path,
    gene_tree_path: Path,
    *,
    taxon_map_path: Path | None = None,
    duplication_cost: float = 2.0,
    loss_cost: float = 1.0,
    transfer_cost: float = 3.0,
) -> DuplicationLossTransferReport:
    """Reconcile one rooted gene tree against one rooted species tree by undated DLT cost."""
    _validate_event_cost("duplication_cost", duplication_cost)
    _validate_event_cost("loss_cost", loss_cost)
    _validate_event_cost("transfer_cost", transfer_cost)

    species_tree = load_tree(species_tree_path)
    gene_tree = load_tree(gene_tree_path)
    _validate_species_tree(species_tree)
    _validate_gene_tree(gene_tree)

    species_index = _build_species_index(species_tree)
    tip_map = _resolve_gene_to_species_tip_map(
        gene_tree=gene_tree,
        species_tree=species_tree,
        taxon_map_path=taxon_map_path,
    )
    gene_nodes_by_key = {
        _gene_branch_key(node): node for node in gene_tree.iter_nodes(order="preorder")
    }
    species_branch_keys = tuple(species_index.nodes_by_branch)

    @cache
    def reconcile_gene_to_species(
        gene_branch: str,
        species_branch: str,
    ) -> _ReconciliationChoice:
        gene_node = gene_nodes_by_key[gene_branch]
        if gene_node.is_leaf():
            return _reconcile_leaf(
                gene_node,
                species_branch=species_branch,
                species_index=species_index,
                tip_map=tip_map,
                loss_cost=loss_cost,
            )
        return _reconcile_internal(
            gene_node,
            species_branch=species_branch,
            species_index=species_index,
            reconcile_gene_to_species=reconcile_gene_to_species,
            duplication_cost=duplication_cost,
            loss_cost=loss_cost,
            transfer_cost=transfer_cost,
        )

    root_gene_branch = _gene_branch_key(gene_tree.root)
    best_root_choice: _ReconciliationChoice | None = None
    for species_branch in species_branch_keys:
        node_choice = reconcile_gene_to_species(root_gene_branch, species_branch)
        if node_choice.total_cost == _UNREACHABLE_COST:
            continue
        origin_losses = _loss_branches_on_path(
            species_index,
            species_index.root_branch,
            species_branch,
        )
        total_cost = node_choice.total_cost + loss_cost * len(origin_losses)
        candidate = _ReconciliationChoice(
            event_type="origin",
            mapped_species_branch=species_index.root_branch,
            left_child_species_branch=species_branch,
            right_child_species_branch=None,
            transferred_child_side=None,
            transfer_recipient_branch=species_branch,
            loss_branches=tuple(origin_losses),
            event_cost=loss_cost * len(origin_losses),
            total_cost=total_cost,
        )
        best_root_choice = _select_better_choice(best_root_choice, candidate)

    if best_root_choice is None or best_root_choice.left_child_species_branch is None:
        raise ValueError("no valid duplication-loss-transfer reconciliation exists")

    event_rows = [
        DuplicationLossTransferEventRow(
            gene_node="origin",
            gene_node_name=None,
            descendant_gene_tips=sorted(gene_tree.tip_names),
            event_type="origin",
            mapped_species_branch=best_root_choice.mapped_species_branch,
            mapped_descendant_species=species_index.descendant_tips_by_branch[
                best_root_choice.mapped_species_branch
            ],
            left_child_gene_node=root_gene_branch,
            right_child_gene_node=None,
            left_child_species_branch=best_root_choice.left_child_species_branch,
            right_child_species_branch=None,
            transferred_child_side=None,
            transfer_recipient_branch=best_root_choice.transfer_recipient_branch,
            loss_branches=list(best_root_choice.loss_branches),
            event_cost=best_root_choice.event_cost,
        )
    ]
    event_rows.extend(
        _build_event_rows(
            gene_node=gene_tree.root,
            mapped_species_branch=best_root_choice.left_child_species_branch,
            gene_nodes_by_key=gene_nodes_by_key,
            species_index=species_index,
            reconcile_gene_to_species=reconcile_gene_to_species,
        )
    )
    observed_species_taxa = sorted(set(tip_map.values()))
    species_only_taxa = sorted(set(species_tree.tip_names) - set(observed_species_taxa))
    return DuplicationLossTransferReport(
        species_tree_path=species_tree_path,
        gene_tree_path=gene_tree_path,
        taxon_map_path=taxon_map_path,
        duplication_cost=duplication_cost,
        loss_cost=loss_cost,
        transfer_cost=transfer_cost,
        observed_species_taxa=observed_species_taxa,
        species_only_taxa=species_only_taxa,
        gene_tip_count=gene_tree.tip_count,
        root_mapping_branch=best_root_choice.left_child_species_branch,
        reconciliation_score=best_root_choice.total_cost,
        duplication_event_count=sum(
            1 for row in event_rows if row.event_type == "duplication"
        ),
        loss_event_count=sum(len(row.loss_branches) for row in event_rows),
        transfer_event_count=sum(
            1 for row in event_rows if row.event_type == "transfer"
        ),
        speciation_event_count=sum(
            1 for row in event_rows if row.event_type == "speciation"
        ),
        mapping_rows=[
            DuplicationLossTransferAssociationRow(
                gene_taxon=gene_taxon,
                species_taxon=species_taxon,
            )
            for gene_taxon, species_taxon in sorted(tip_map.items())
        ],
        event_rows=event_rows,
    )


def _reconcile_leaf(
    gene_node: TreeNode,
    *,
    species_branch: str,
    species_index: _SpeciesIndex,
    tip_map: dict[str, str],
    loss_cost: float,
) -> _ReconciliationChoice:
    associated_species = tip_map[gene_node.name or ""]
    if associated_species not in species_index.descendants_by_branch[species_branch]:
        return _unreachable_choice(species_branch)
    loss_branches = _loss_branches_on_path(
        species_index,
        species_branch,
        associated_species,
    )
    event_cost = loss_cost * len(loss_branches)
    return _ReconciliationChoice(
        event_type="leaf",
        mapped_species_branch=species_branch,
        left_child_species_branch=None,
        right_child_species_branch=None,
        transferred_child_side=None,
        transfer_recipient_branch=None,
        loss_branches=tuple(loss_branches),
        event_cost=event_cost,
        total_cost=event_cost,
    )


def _reconcile_internal(
    gene_node: TreeNode,
    *,
    species_branch: str,
    species_index: _SpeciesIndex,
    reconcile_gene_to_species,
    duplication_cost: float,
    loss_cost: float,
    transfer_cost: float,
) -> _ReconciliationChoice:
    left_gene = gene_node.children[0]
    right_gene = gene_node.children[1]
    left_gene_key = _gene_branch_key(left_gene)
    right_gene_key = _gene_branch_key(right_gene)
    best_choice: _ReconciliationChoice | None = None
    species_node = species_index.nodes_by_branch[species_branch]

    if not species_node.is_leaf():
        left_species_root = _species_branch_key(species_node.children[0])
        right_species_root = _species_branch_key(species_node.children[1])
        for left_target in species_index.descendants_including_self_by_branch[
            left_species_root
        ]:
            left_choice = reconcile_gene_to_species(left_gene_key, left_target)
            if left_choice.total_cost == _UNREACHABLE_COST:
                continue
            left_losses = _loss_branches_on_path(
                species_index,
                left_species_root,
                left_target,
            )
            for right_target in species_index.descendants_including_self_by_branch[
                right_species_root
            ]:
                right_choice = reconcile_gene_to_species(right_gene_key, right_target)
                if right_choice.total_cost == _UNREACHABLE_COST:
                    continue
                right_losses = _loss_branches_on_path(
                    species_index,
                    right_species_root,
                    right_target,
                )
                best_choice = _select_better_choice(
                    best_choice,
                    _build_event_choice(
                        event_type="speciation",
                        mapped_species_branch=species_branch,
                        left_child_species_branch=left_target,
                        right_child_species_branch=right_target,
                        transferred_child_side=None,
                        transfer_recipient_branch=None,
                        loss_branches=[*left_losses, *right_losses],
                        child_total_cost=left_choice.total_cost
                        + right_choice.total_cost,
                        event_cost=loss_cost * (len(left_losses) + len(right_losses)),
                    ),
                )
                best_choice = _select_better_choice(
                    best_choice,
                    _build_event_choice(
                        event_type="speciation",
                        mapped_species_branch=species_branch,
                        left_child_species_branch=right_target,
                        right_child_species_branch=left_target,
                        transferred_child_side=None,
                        transfer_recipient_branch=None,
                        loss_branches=[*right_losses, *left_losses],
                        child_total_cost=left_choice.total_cost
                        + right_choice.total_cost,
                        event_cost=loss_cost * (len(left_losses) + len(right_losses)),
                    ),
                )

    for left_target in species_index.descendants_including_self_by_branch[
        species_branch
    ]:
        left_choice = reconcile_gene_to_species(left_gene_key, left_target)
        if left_choice.total_cost == _UNREACHABLE_COST:
            continue
        left_losses = _loss_branches_on_path(species_index, species_branch, left_target)
        for right_target in species_index.descendants_including_self_by_branch[
            species_branch
        ]:
            right_choice = reconcile_gene_to_species(right_gene_key, right_target)
            if right_choice.total_cost == _UNREACHABLE_COST:
                continue
            right_losses = _loss_branches_on_path(
                species_index,
                species_branch,
                right_target,
            )
            best_choice = _select_better_choice(
                best_choice,
                _build_event_choice(
                    event_type="duplication",
                    mapped_species_branch=species_branch,
                    left_child_species_branch=left_target,
                    right_child_species_branch=right_target,
                    transferred_child_side=None,
                    transfer_recipient_branch=None,
                    loss_branches=[*left_losses, *right_losses],
                    child_total_cost=left_choice.total_cost + right_choice.total_cost,
                    event_cost=duplication_cost
                    + loss_cost * (len(left_losses) + len(right_losses)),
                ),
            )

    transfer_targets = [
        branch
        for branch in species_index.nodes_by_branch
        if branch != species_branch
        and not _is_descendant_or_self(
            species_index,
            species_branch,
            branch,
        )
        and not _is_descendant_or_self(
            species_index,
            branch,
            species_branch,
        )
    ]
    for donor_side, donor_key, recipient_key in (
        ("left", left_gene_key, right_gene_key),
        ("right", right_gene_key, left_gene_key),
    ):
        for donor_target in species_index.descendants_including_self_by_branch[
            species_branch
        ]:
            donor_choice = reconcile_gene_to_species(donor_key, donor_target)
            if donor_choice.total_cost == _UNREACHABLE_COST:
                continue
            donor_losses = _loss_branches_on_path(
                species_index,
                species_branch,
                donor_target,
            )
            for recipient_target in transfer_targets:
                recipient_choice = reconcile_gene_to_species(
                    recipient_key,
                    recipient_target,
                )
                if recipient_choice.total_cost == _UNREACHABLE_COST:
                    continue
                best_choice = _select_better_choice(
                    best_choice,
                    _build_event_choice(
                        event_type="transfer",
                        mapped_species_branch=species_branch,
                        left_child_species_branch=(
                            donor_target if donor_side == "left" else recipient_target
                        ),
                        right_child_species_branch=(
                            recipient_target if donor_side == "left" else donor_target
                        ),
                        transferred_child_side="right"
                        if donor_side == "left"
                        else "left",
                        transfer_recipient_branch=recipient_target,
                        loss_branches=donor_losses,
                        child_total_cost=donor_choice.total_cost
                        + recipient_choice.total_cost,
                        event_cost=transfer_cost + loss_cost * len(donor_losses),
                    ),
                )

    return best_choice or _unreachable_choice(species_branch)


def _build_event_rows(
    *,
    gene_node: TreeNode,
    mapped_species_branch: str,
    gene_nodes_by_key: dict[str, TreeNode],
    species_index: _SpeciesIndex,
    reconcile_gene_to_species,
) -> list[DuplicationLossTransferEventRow]:
    choice = reconcile_gene_to_species(
        _gene_branch_key(gene_node), mapped_species_branch
    )
    left_child = None if gene_node.is_leaf() else gene_node.children[0]
    right_child = None if gene_node.is_leaf() else gene_node.children[1]
    event_rows = [
        DuplicationLossTransferEventRow(
            gene_node=_gene_branch_key(gene_node),
            gene_node_name=gene_node.name,
            descendant_gene_tips=sorted(gene_node.descendant_taxa),
            event_type=choice.event_type,
            mapped_species_branch=choice.mapped_species_branch,
            mapped_descendant_species=species_index.descendant_tips_by_branch[
                choice.mapped_species_branch
            ],
            left_child_gene_node=None
            if left_child is None
            else _gene_branch_key(left_child),
            right_child_gene_node=None
            if right_child is None
            else _gene_branch_key(right_child),
            left_child_species_branch=choice.left_child_species_branch,
            right_child_species_branch=choice.right_child_species_branch,
            transferred_child_side=choice.transferred_child_side,
            transfer_recipient_branch=choice.transfer_recipient_branch,
            loss_branches=list(choice.loss_branches),
            event_cost=choice.event_cost,
        )
    ]
    if left_child is not None and choice.left_child_species_branch is not None:
        event_rows.extend(
            _build_event_rows(
                gene_node=left_child,
                mapped_species_branch=choice.left_child_species_branch,
                gene_nodes_by_key=gene_nodes_by_key,
                species_index=species_index,
                reconcile_gene_to_species=reconcile_gene_to_species,
            )
        )
    if right_child is not None and choice.right_child_species_branch is not None:
        event_rows.extend(
            _build_event_rows(
                gene_node=right_child,
                mapped_species_branch=choice.right_child_species_branch,
                gene_nodes_by_key=gene_nodes_by_key,
                species_index=species_index,
                reconcile_gene_to_species=reconcile_gene_to_species,
            )
        )
    return event_rows


def _build_event_choice(
    *,
    event_type: str,
    mapped_species_branch: str,
    left_child_species_branch: str | None,
    right_child_species_branch: str | None,
    transferred_child_side: str | None,
    transfer_recipient_branch: str | None,
    loss_branches: list[str],
    child_total_cost: float,
    event_cost: float,
) -> _ReconciliationChoice:
    return _ReconciliationChoice(
        event_type=event_type,
        mapped_species_branch=mapped_species_branch,
        left_child_species_branch=left_child_species_branch,
        right_child_species_branch=right_child_species_branch,
        transferred_child_side=transferred_child_side,
        transfer_recipient_branch=transfer_recipient_branch,
        loss_branches=tuple(loss_branches),
        event_cost=event_cost,
        total_cost=child_total_cost + event_cost,
    )


def _select_better_choice(
    current: _ReconciliationChoice | None,
    candidate: _ReconciliationChoice,
) -> _ReconciliationChoice:
    if current is None:
        return candidate
    current_key = _choice_sort_key(current)
    candidate_key = _choice_sort_key(candidate)
    if candidate.total_cost < current.total_cost:
        return candidate
    if candidate.total_cost > current.total_cost:
        return current
    return candidate if candidate_key < current_key else current


def _choice_sort_key(choice: _ReconciliationChoice) -> tuple[object, ...]:
    return (
        choice.event_type,
        choice.mapped_species_branch,
        ""
        if choice.left_child_species_branch is None
        else choice.left_child_species_branch,
        ""
        if choice.right_child_species_branch is None
        else choice.right_child_species_branch,
        "" if choice.transferred_child_side is None else choice.transferred_child_side,
        ""
        if choice.transfer_recipient_branch is None
        else choice.transfer_recipient_branch,
        choice.loss_branches,
    )


def _unreachable_choice(species_branch: str) -> _ReconciliationChoice:
    return _ReconciliationChoice(
        event_type="unreachable",
        mapped_species_branch=species_branch,
        left_child_species_branch=None,
        right_child_species_branch=None,
        transferred_child_side=None,
        transfer_recipient_branch=None,
        loss_branches=(),
        event_cost=_UNREACHABLE_COST,
        total_cost=_UNREACHABLE_COST,
    )


def _validate_event_cost(name: str, value: float) -> None:
    if value < 0.0:
        raise ValueError(f"{name} must be nonnegative")


def _resolve_gene_to_species_tip_map(
    *,
    gene_tree: PhyloTree,
    species_tree: PhyloTree,
    taxon_map_path: Path | None,
) -> dict[str, str]:
    species_taxa = set(species_tree.tip_names)
    if taxon_map_path is None:
        if any(
            (tip_name or "") not in species_taxa for tip_name in gene_tree.tip_names
        ):
            raise ValueError(
                "duplication-loss-transfer reconciliation requires --taxon-map when gene tips do not exactly match species-tree taxa"
            )
        return {tip_name: tip_name for tip_name in gene_tree.tip_names}

    table = load_taxon_table(taxon_map_path)
    if "species_taxon" not in table.columns:
        raise ValueError("taxon map must include a 'species_taxon' column")
    tip_map = {row[table.taxon_column]: row["species_taxon"] for row in table.rows}
    for gene_taxon in gene_tree.tip_names:
        if gene_taxon not in tip_map:
            raise ValueError(f"taxon map is missing gene tip '{gene_taxon}'")
        if tip_map[gene_taxon] not in species_taxa:
            raise ValueError(
                f"taxon map assigns gene tip '{gene_taxon}' to unknown species '{tip_map[gene_taxon]}'"
            )
    extra_gene_taxa = sorted(set(tip_map) - set(gene_tree.tip_names))
    if extra_gene_taxa:
        raise ValueError(
            "taxon map contains gene taxa absent from the gene tree: "
            + ", ".join(extra_gene_taxa)
        )
    return tip_map


def _validate_species_tree(tree: PhyloTree) -> None:
    if len(tree.root.children) != 2:
        raise ValueError(
            "duplication-loss-transfer reconciliation requires a rooted species tree"
        )
    if any(node.name is None for node in tree.iter_leaves()):
        raise ValueError(
            "duplication-loss-transfer reconciliation requires named species-tree tips"
        )
    if len(set(tree.tip_names)) != tree.tip_count:
        raise ValueError(
            "duplication-loss-transfer reconciliation requires unique species-tree tip labels"
        )
    if any(len(node.children) != 2 for node in tree.iter_internal_nodes()):
        raise ValueError(
            "duplication-loss-transfer reconciliation requires a strictly binary species tree"
        )


def _validate_gene_tree(tree: PhyloTree) -> None:
    if len(tree.root.children) != 2:
        raise ValueError(
            "duplication-loss-transfer reconciliation requires a rooted gene tree"
        )
    if any(node.name is None for node in tree.iter_leaves()):
        raise ValueError(
            "duplication-loss-transfer reconciliation requires named gene-tree tips"
        )
    if len(set(tree.tip_names)) != tree.tip_count:
        raise ValueError(
            "duplication-loss-transfer reconciliation requires unique gene-tree tip labels"
        )
    if any(len(node.children) != 2 for node in tree.iter_internal_nodes()):
        raise ValueError(
            "duplication-loss-transfer reconciliation requires a strictly binary gene tree"
        )


def _build_species_index(tree: PhyloTree) -> _SpeciesIndex:
    nodes_by_branch = {
        _species_branch_key(node): node for node in tree.iter_nodes(order="preorder")
    }
    descendants_by_branch = {
        branch: set(node.descendant_taxa) for branch, node in nodes_by_branch.items()
    }
    descendants_including_self_by_branch = {
        branch: tuple(_collect_subtree_branch_keys(node))
        for branch, node in nodes_by_branch.items()
    }
    descendant_tips_by_branch = {
        branch: sorted(node.descendant_taxa) for branch, node in nodes_by_branch.items()
    }
    ancestry_by_branch = {
        branch: tuple(_ancestor_branch_keys(node))
        for branch, node in nodes_by_branch.items()
    }
    return _SpeciesIndex(
        tree=tree,
        nodes_by_branch=nodes_by_branch,
        descendants_by_branch=descendants_by_branch,
        descendants_including_self_by_branch=descendants_including_self_by_branch,
        descendant_tips_by_branch=descendant_tips_by_branch,
        ancestry_by_branch=ancestry_by_branch,
        root_branch=_species_branch_key(tree.root),
    )


def _ancestor_branch_keys(node: TreeNode) -> list[str]:
    values: list[str] = []
    current = node
    while current is not None:
        values.append(_species_branch_key(current))
        current = current.parent
    return list(reversed(values))


def _collect_subtree_branch_keys(node: TreeNode) -> list[str]:
    values = [_species_branch_key(node)]
    for child in node.children:
        values.extend(_collect_subtree_branch_keys(child))
    return values


def _loss_branches_on_path(
    species_index: _SpeciesIndex,
    start_branch: str,
    end_branch: str,
) -> list[str]:
    if start_branch == end_branch:
        return []
    if (
        end_branch
        not in species_index.descendants_including_self_by_branch[start_branch]
    ):
        return []
    current = species_index.nodes_by_branch[start_branch]
    target = species_index.nodes_by_branch[end_branch]
    loss_branches: list[str] = []
    while current is not target:
        if current.is_leaf():
            return []
        matching_child = None
        skipped_child = None
        for child in current.children:
            child_branch = _species_branch_key(child)
            if (
                end_branch
                in species_index.descendants_including_self_by_branch[child_branch]
            ):
                matching_child = child
            else:
                skipped_child = child
        if matching_child is None:
            return []
        if skipped_child is not None:
            loss_branches.append(_species_branch_key(skipped_child))
        current = matching_child
    return loss_branches


def _is_descendant_or_self(
    species_index: _SpeciesIndex,
    ancestor_branch: str,
    descendant_branch: str,
) -> bool:
    return (
        descendant_branch
        in species_index.descendants_including_self_by_branch[ancestor_branch]
    )


def _species_branch_key(node: TreeNode) -> str:
    if node.is_leaf():
        return node.name or "<unnamed-species>"
    return "|".join(node.descendant_taxa)


def _gene_branch_key(node: TreeNode) -> str:
    return node.node_id or (
        node.name if node.is_leaf() else "|".join(node.descendant_taxa)
    )
