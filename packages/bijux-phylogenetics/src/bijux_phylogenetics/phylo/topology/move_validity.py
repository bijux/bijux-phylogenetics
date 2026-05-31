from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
import json

from bijux_phylogenetics.io.newick import load_newick, loads_newick, write_newick
from bijux_phylogenetics.phylo.topology.models import TopologyMoveValidityReport

from .rooted_nni import RootedNniMoveCandidate, iter_rooted_nni_move_candidates
from .rooted_spr import (
    RootedSprMoveCandidate,
    _candidate_is_self_regraft,
    _collect_rooted_spr_move_candidates,
)
from .rooted_tbr import RootedTbrMoveCandidate, iter_rooted_tbr_move_candidates
from .tree import PhyloTree


def summarize_rooted_nni_move_validity(
    tree: PhyloTree | Path,
    candidate: RootedNniMoveCandidate,
) -> TopologyMoveValidityReport:
    """Record whether one rooted NNI move request is valid on the supplied tree."""
    resolved_tree, input_tree_path = _resolve_move_validity_tree(tree)
    candidate_payload = asdict(candidate)
    rejected_report = _reject_invalid_rooted_move_tree(
        resolved_tree,
        input_tree_path=input_tree_path,
        move_family="rooted-nni",
        candidate_payload=candidate_payload,
    )
    if rejected_report is not None:
        return rejected_report

    available_candidates = list(iter_rooted_nni_move_candidates(resolved_tree))
    if not available_candidates:
        return _build_move_validity_report(
            resolved_tree,
            input_tree_path=input_tree_path,
            move_family="rooted-nni",
            candidate_payload=candidate_payload,
            available_move_count=0,
            validity_decision="rejected",
            rejection_code="topology_move_tree_too_small",
            rejection_reason=(
                "rooted NNI move validation requires a tree with at least one legal"
                " internal rearrangement"
            ),
            evidence={
                "tip_count": resolved_tree.tip_count,
                "internal_node_count": resolved_tree.internal_node_count,
            },
        )

    missing_node_fields = _missing_nni_node_fields(resolved_tree, candidate)
    if missing_node_fields:
        return _build_move_validity_report(
            resolved_tree,
            input_tree_path=input_tree_path,
            move_family="rooted-nni",
            candidate_payload=candidate_payload,
            available_move_count=len(available_candidates),
            validity_decision="rejected",
            rejection_code="topology_move_missing_branch",
            rejection_reason=(
                "rooted NNI move validation could not resolve one or more requested"
                " move branches on the input tree"
            ),
            evidence={"missing_node_fields": missing_node_fields},
        )

    if candidate not in available_candidates:
        return _build_move_validity_report(
            resolved_tree,
            input_tree_path=input_tree_path,
            move_family="rooted-nni",
            candidate_payload=candidate_payload,
            available_move_count=len(available_candidates),
            validity_decision="rejected",
            rejection_code="topology_move_candidate_unavailable",
            rejection_reason=(
                "rooted NNI move validation rejected a requested branch exchange"
                " that is not a legal move on the input tree"
            ),
            evidence={"available_candidate_count": len(available_candidates)},
        )

    return _build_move_validity_report(
        resolved_tree,
        input_tree_path=input_tree_path,
        move_family="rooted-nni",
        candidate_payload=candidate_payload,
        available_move_count=len(available_candidates),
        validity_decision="accepted",
        rejection_code=None,
        rejection_reason=None,
        evidence={"candidate_available": True},
    )


def summarize_rooted_spr_move_validity(
    tree: PhyloTree | Path,
    candidate: RootedSprMoveCandidate,
) -> TopologyMoveValidityReport:
    """Record whether one rooted SPR move request is valid on the supplied tree."""
    resolved_tree, input_tree_path = _resolve_move_validity_tree(tree)
    candidate_payload = asdict(candidate)
    rejected_report = _reject_invalid_rooted_move_tree(
        resolved_tree,
        input_tree_path=input_tree_path,
        move_family="rooted-spr",
        candidate_payload=candidate_payload,
    )
    if rejected_report is not None:
        return rejected_report

    (
        available_candidates,
        _skipped_pruned,
        _skipped_regraft,
        _skipped_budget_move_candidates,
    ) = (
        _collect_rooted_spr_move_candidates(resolved_tree, budget=None)
    )
    if not available_candidates:
        return _build_move_validity_report(
            resolved_tree,
            input_tree_path=input_tree_path,
            move_family="rooted-spr",
            candidate_payload=candidate_payload,
            available_move_count=0,
            validity_decision="rejected",
            rejection_code="topology_move_tree_too_small",
            rejection_reason=(
                "rooted SPR move validation requires a tree with at least one legal"
                " prune-and-regraft move"
            ),
            evidence={
                "tip_count": resolved_tree.tip_count,
                "internal_node_count": resolved_tree.internal_node_count,
            },
        )

    if _candidate_is_self_regraft(candidate):
        return _build_move_validity_report(
            resolved_tree,
            input_tree_path=input_tree_path,
            move_family="rooted-spr",
            candidate_payload=candidate_payload,
            available_move_count=len(available_candidates),
            validity_decision="rejected",
            rejection_code="topology_move_self_regraft",
            rejection_reason=(
                "rooted SPR move validation rejected a self-regraft candidate that"
                " would reattach the pruned clade onto itself"
            ),
            evidence={
                "pruned_clade_id": candidate.pruned_clade_id,
                "regraft_target_branch_id": candidate.regraft_target_branch_id,
            },
        )

    missing_branch_fields = _missing_spr_branch_fields(resolved_tree, candidate)
    if missing_branch_fields:
        return _build_move_validity_report(
            resolved_tree,
            input_tree_path=input_tree_path,
            move_family="rooted-spr",
            candidate_payload=candidate_payload,
            available_move_count=len(available_candidates),
            validity_decision="rejected",
            rejection_code="topology_move_missing_branch",
            rejection_reason=(
                "rooted SPR move validation could not resolve one or more requested"
                " prune or regraft branches on the input tree"
            ),
            evidence={"missing_branch_fields": missing_branch_fields},
        )

    if candidate not in available_candidates:
        return _build_move_validity_report(
            resolved_tree,
            input_tree_path=input_tree_path,
            move_family="rooted-spr",
            candidate_payload=candidate_payload,
            available_move_count=len(available_candidates),
            validity_decision="rejected",
            rejection_code="topology_move_candidate_unavailable",
            rejection_reason=(
                "rooted SPR move validation rejected a requested prune-and-regraft"
                " combination that is not legal on the input tree"
            ),
            evidence={"available_candidate_count": len(available_candidates)},
        )

    return _build_move_validity_report(
        resolved_tree,
        input_tree_path=input_tree_path,
        move_family="rooted-spr",
        candidate_payload=candidate_payload,
        available_move_count=len(available_candidates),
        validity_decision="accepted",
        rejection_code=None,
        rejection_reason=None,
        evidence={"candidate_available": True},
    )


def summarize_rooted_tbr_move_validity(
    tree: PhyloTree | Path,
    candidate: RootedTbrMoveCandidate,
) -> TopologyMoveValidityReport:
    """Record whether one rooted TBR move request is valid on the supplied tree."""
    resolved_tree, input_tree_path = _resolve_move_validity_tree(tree)
    candidate_payload = asdict(candidate)
    rejected_report = _reject_invalid_rooted_move_tree(
        resolved_tree,
        input_tree_path=input_tree_path,
        move_family="rooted-tbr",
        candidate_payload=candidate_payload,
    )
    if rejected_report is not None:
        return rejected_report

    available_candidates = list(iter_rooted_tbr_move_candidates(resolved_tree))
    if not available_candidates:
        return _build_move_validity_report(
            resolved_tree,
            input_tree_path=input_tree_path,
            move_family="rooted-tbr",
            candidate_payload=candidate_payload,
            available_move_count=0,
            validity_decision="rejected",
            rejection_code="topology_move_tree_too_small",
            rejection_reason=(
                "rooted TBR move validation requires a tree with at least one legal"
                " bisection-and-reconnection move"
            ),
            evidence={
                "tip_count": resolved_tree.tip_count,
                "internal_node_count": resolved_tree.internal_node_count,
            },
        )

    missing_branch_fields = _missing_tbr_branch_fields(resolved_tree, candidate)
    if missing_branch_fields:
        return _build_move_validity_report(
            resolved_tree,
            input_tree_path=input_tree_path,
            move_family="rooted-tbr",
            candidate_payload=candidate_payload,
            available_move_count=len(available_candidates),
            validity_decision="rejected",
            rejection_code="topology_move_missing_branch",
            rejection_reason=(
                "rooted TBR move validation could not resolve one or more requested"
                " cut or reconnection branches on the input tree"
            ),
            evidence={"missing_branch_fields": missing_branch_fields},
        )

    if candidate not in available_candidates:
        return _build_move_validity_report(
            resolved_tree,
            input_tree_path=input_tree_path,
            move_family="rooted-tbr",
            candidate_payload=candidate_payload,
            available_move_count=len(available_candidates),
            validity_decision="rejected",
            rejection_code="topology_move_candidate_unavailable",
            rejection_reason=(
                "rooted TBR move validation rejected a requested cut-and-reconnect"
                " combination that is not legal on the input tree"
            ),
            evidence={"available_candidate_count": len(available_candidates)},
        )

    return _build_move_validity_report(
        resolved_tree,
        input_tree_path=input_tree_path,
        move_family="rooted-tbr",
        candidate_payload=candidate_payload,
        available_move_count=len(available_candidates),
        validity_decision="accepted",
        rejection_code=None,
        rejection_reason=None,
        evidence={"candidate_available": True},
    )


def write_topology_move_validity_run_json(
    path: Path,
    report: TopologyMoveValidityReport,
) -> Path:
    """Write one machine-readable topology move-validity payload."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "algorithm": report.algorithm,
        "move_family": report.move_family,
        "input_tree_path": (
            None if report.input_tree_path is None else str(report.input_tree_path)
        ),
        "input_tree_newick": report.input_tree_newick,
        "tip_count": report.tip_count,
        "internal_node_count": report.internal_node_count,
        "rooted": report.rooted,
        "available_move_count": report.available_move_count,
        "input_validation_errors": report.input_validation_errors,
        "candidate_payload": report.candidate_payload,
        "validity_decision": report.validity_decision,
        "rejection_code": report.rejection_code,
        "rejection_reason": report.rejection_reason,
        "evidence": report.evidence,
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def write_topology_move_validity_artifacts(
    out_dir: Path,
    report: TopologyMoveValidityReport,
) -> dict[str, Path]:
    """Write the governed artifact family for one topology move-validity report."""
    out_dir.mkdir(parents=True, exist_ok=True)
    input_tree_path = write_newick(
        out_dir / "input_tree.nwk",
        loads_newick(report.input_tree_newick),
    )
    run_json_path = write_topology_move_validity_run_json(
        out_dir / "run.json",
        report,
    )
    return {
        "input_tree_path": input_tree_path,
        "run_json_path": run_json_path,
    }


def _resolve_move_validity_tree(
    tree: PhyloTree | Path,
) -> tuple[PhyloTree, Path | None]:
    resolved_tree_path = tree if isinstance(tree, Path) else None
    resolved_tree = load_newick(tree) if isinstance(tree, Path) else tree.copy()
    return resolved_tree.refresh(), resolved_tree_path


def _reject_invalid_rooted_move_tree(
    tree: PhyloTree,
    *,
    input_tree_path: Path | None,
    move_family: str,
    candidate_payload: dict[str, object],
) -> TopologyMoveValidityReport | None:
    validation_errors = tree.validation_errors()
    if validation_errors:
        return _build_move_validity_report(
            tree,
            input_tree_path=input_tree_path,
            move_family=move_family,
            candidate_payload=candidate_payload,
            available_move_count=0,
            validity_decision="rejected",
            rejection_code="topology_move_invalid_tree_structure",
            rejection_reason=(
                f"{move_family} move validation requires a structurally valid tree"
            ),
            evidence={"validation_errors": validation_errors},
        )
    if tree.rooted is False and len(tree.root.children) != 2:
        return _build_move_validity_report(
            tree,
            input_tree_path=input_tree_path,
            move_family=move_family,
            candidate_payload=candidate_payload,
            available_move_count=0,
            validity_decision="rejected",
            rejection_code="topology_move_incompatible_rootedness",
            rejection_reason=(
                f"{move_family} move validation requires a rooted binary tree"
                " representation"
            ),
            evidence={
                "rooted": tree.rooted,
                "root_child_count": len(tree.root.children),
            },
        )
    if len(tree.root.children) != 2:
        return _build_move_validity_report(
            tree,
            input_tree_path=input_tree_path,
            move_family=move_family,
            candidate_payload=candidate_payload,
            available_move_count=0,
            validity_decision="rejected",
            rejection_code="topology_move_binary_root_required",
            rejection_reason=(
                f"{move_family} move validation requires a binary root"
            ),
            evidence={"root_child_count": len(tree.root.children)},
        )
    invalid_internal_node_ids = [
        node.node_id
        for node in tree.iter_internal_nodes(order="preorder")
        if len(node.children) != 2
    ]
    if invalid_internal_node_ids:
        return _build_move_validity_report(
            tree,
            input_tree_path=input_tree_path,
            move_family=move_family,
            candidate_payload=candidate_payload,
            available_move_count=0,
            validity_decision="rejected",
            rejection_code="topology_move_strict_bifurcation_required",
            rejection_reason=(
                f"{move_family} move validation requires a strictly bifurcating tree"
            ),
            evidence={"invalid_internal_node_ids": invalid_internal_node_ids},
        )
    return None


def _build_move_validity_report(
    tree: PhyloTree,
    *,
    input_tree_path: Path | None,
    move_family: str,
    candidate_payload: dict[str, object],
    available_move_count: int,
    validity_decision: str,
    rejection_code: str | None,
    rejection_reason: str | None,
    evidence: dict[str, object],
) -> TopologyMoveValidityReport:
    return TopologyMoveValidityReport(
        algorithm="topology-move-validity",
        move_family=move_family,
        input_tree_path=input_tree_path,
        input_tree_newick=tree.to_newick(),
        tip_count=tree.tip_count,
        internal_node_count=tree.internal_node_count,
        rooted=tree.rooted,
        available_move_count=available_move_count,
        input_validation_errors=tree.validation_errors(),
        candidate_payload=candidate_payload,
        validity_decision=validity_decision,
        rejection_code=rejection_code,
        rejection_reason=rejection_reason,
        evidence=evidence,
    )


def _missing_nni_node_fields(
    tree: PhyloTree,
    candidate: RootedNniMoveCandidate,
) -> list[str]:
    missing_fields: list[str] = []
    for field_name, node_id in (
        ("parent_node_id", candidate.parent_node_id),
        ("child_node_id", candidate.child_node_id),
        ("sibling_node_id", candidate.sibling_node_id),
        ("exchanged_child_node_id", candidate.exchanged_child_node_id),
    ):
        try:
            tree.node_by_id(node_id)
        except KeyError:
            missing_fields.append(field_name)
    return missing_fields


def _missing_spr_branch_fields(
    tree: PhyloTree,
    candidate: RootedSprMoveCandidate,
) -> list[str]:
    available_candidates, _skipped_pruned, _skipped_regraft = (
        _collect_rooted_spr_move_candidates(tree, budget=None)
    )
    available_pruned_node_ids = {
        available_candidate.pruned_node_id
        for available_candidate in available_candidates
    }
    available_regraft_target_branch_ids = {
        available_candidate.regraft_target_branch_id
        for available_candidate in available_candidates
    }
    missing_fields: list[str] = []
    if candidate.pruned_node_id not in available_pruned_node_ids:
        missing_fields.append("pruned_node_id")
    if candidate.regraft_target_branch_id not in available_regraft_target_branch_ids:
        missing_fields.append("regraft_target_branch_id")
    return missing_fields


def _missing_tbr_branch_fields(
    tree: PhyloTree,
    candidate: RootedTbrMoveCandidate,
) -> list[str]:
    available_candidates = list(iter_rooted_tbr_move_candidates(tree))
    available_cut_parent_node_ids = {
        available_candidate.cut_parent_node_id
        for available_candidate in available_candidates
    }
    available_cut_child_node_ids = {
        available_candidate.cut_child_node_id
        for available_candidate in available_candidates
    }
    available_cut_edge_ids = {
        available_candidate.cut_edge_id for available_candidate in available_candidates
    }
    available_attachment_branch_ids = {
        available_candidate.left_attachment_branch_id
        for available_candidate in available_candidates
    } | {
        available_candidate.right_attachment_branch_id
        for available_candidate in available_candidates
    }
    missing_fields: list[str] = []
    if candidate.cut_parent_node_id not in available_cut_parent_node_ids:
        missing_fields.append("cut_parent_node_id")
    if candidate.cut_child_node_id not in available_cut_child_node_ids:
        missing_fields.append("cut_child_node_id")
    if candidate.cut_edge_id not in available_cut_edge_ids:
        missing_fields.append("cut_edge_id")
    if candidate.left_attachment_branch_id not in available_attachment_branch_ids:
        missing_fields.append("left_attachment_branch_id")
    if candidate.right_attachment_branch_id not in available_attachment_branch_ids:
        missing_fields.append("right_attachment_branch_id")
    return missing_fields
