from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.diagnostics.validation.structure import _load_tree
from bijux_phylogenetics.io.iqtree_support import support_fraction
from bijux_phylogenetics.phylo.topology.clades import (
    canonical_clade_id,
    informative_rooted_clade_nodes,
    node_support_value,
    split_sort_key,
)

from .models import (
    _STRONG_SUPPORT_THRESHOLD,
    _SUPPORT_DISAGREEMENT_THRESHOLD,
    _WEAK_SUPPORT_THRESHOLD,
    CladeSupportPair,
    SupportComparisonReport,
    SupportConflictRow,
)


def _split_id(signature: frozenset[str]) -> str:
    return canonical_clade_id(signature)


def _support_strength(
    value: float | None,
    *,
    strong_support_threshold: float,
    weak_support_threshold: float,
) -> str:
    fraction = support_fraction(value)
    if fraction is None:
        return "unavailable"
    if fraction >= strong_support_threshold:
        return "strong"
    if fraction >= weak_support_threshold:
        return "moderate"
    return "low"


def compare_support_values(
    left_path: Path,
    right_path: Path,
    *,
    strong_support_threshold: float = _STRONG_SUPPORT_THRESHOLD,
    weak_support_threshold: float = _WEAK_SUPPORT_THRESHOLD,
    support_disagreement_threshold: float = _SUPPORT_DISAGREEMENT_THRESHOLD,
) -> SupportComparisonReport:
    """Compare clade support values and support-aware conflicts across two trees."""
    left = _load_tree(left_path)
    right = _load_tree(right_path)
    return _build_support_comparison_report(
        left_path,
        right_path,
        left,
        right,
        strong_support_threshold=strong_support_threshold,
        weak_support_threshold=weak_support_threshold,
        support_disagreement_threshold=support_disagreement_threshold,
    )


def _build_support_comparison_report(
    left_path: Path,
    right_path: Path,
    left,
    right,
    *,
    strong_support_threshold: float,
    weak_support_threshold: float,
    support_disagreement_threshold: float,
) -> SupportComparisonReport:
    left_taxa = set(left.tip_names)
    right_taxa = set(right.tip_names)
    shared_taxa = left_taxa & right_taxa
    if len(shared_taxa) < 2:
        raise ValueError("support comparison requires at least two shared taxa")

    left_clades = informative_rooted_clade_nodes(left, shared_taxa)
    right_clades = informative_rooted_clade_nodes(right, shared_taxa)
    shared_clade_ids = sorted(
        left_clades.keys() & right_clades.keys(), key=split_sort_key
    )
    all_clade_ids = sorted(
        left_clades.keys() | right_clades.keys(),
        key=split_sort_key,
    )
    shared_clades: list[CladeSupportPair] = []
    conflicting_clades: list[SupportConflictRow] = []
    for clade_id in shared_clade_ids:
        left_support = node_support_value(left_clades[clade_id])
        right_support = node_support_value(right_clades[clade_id])
        left_fraction = support_fraction(left_support)
        right_fraction = support_fraction(right_support)
        shared_clades.append(
            CladeSupportPair(
                split_id=_split_id(clade_id),
                left_support=left_support,
                right_support=right_support,
                left_support_fraction=left_fraction,
                right_support_fraction=right_fraction,
                support_fraction_delta=(
                    None
                    if left_fraction is None or right_fraction is None
                    else abs(left_fraction - right_fraction)
                ),
                support_disagreement=(
                    left_fraction is not None
                    and right_fraction is not None
                    and abs(left_fraction - right_fraction)
                    >= support_disagreement_threshold
                ),
            )
        )
    for clade_id in all_clade_ids:
        left_present = clade_id in left_clades
        right_present = clade_id in right_clades
        if left_present and right_present:
            continue
        left_support = (
            node_support_value(left_clades[clade_id]) if left_present else None
        )
        right_support = (
            node_support_value(right_clades[clade_id]) if right_present else None
        )
        left_fraction = support_fraction(left_support)
        right_fraction = support_fraction(right_support)
        strongest_support_fraction = (
            max(value for value in (left_fraction, right_fraction) if value is not None)
            if left_fraction is not None or right_fraction is not None
            else None
        )
        support_strength = _support_strength(
            left_support if left_present else right_support,
            strong_support_threshold=strong_support_threshold,
            weak_support_threshold=weak_support_threshold,
        )
        if strongest_support_fraction is None:
            conflict_classification = "support_unavailable"
            detail = "clade conflict could not be ranked because no branch support was available"
        elif strongest_support_fraction >= strong_support_threshold:
            conflict_classification = "high_support_conflict"
            detail = "conflicting clade carried strong branch support in the tree where it was present"
        elif strongest_support_fraction >= weak_support_threshold:
            conflict_classification = "moderate_support_disagreement"
            detail = "conflicting clade carried moderate branch support in the tree where it was present"
        else:
            conflict_classification = "low_support_disagreement"
            detail = "conflicting clade was only weakly supported in the tree where it was present"
        conflicting_clades.append(
            SupportConflictRow(
                split_id=_split_id(clade_id),
                comparison_status="left_only" if left_present else "right_only",
                left_present=left_present,
                right_present=right_present,
                left_support=left_support,
                right_support=right_support,
                left_support_fraction=left_fraction,
                right_support_fraction=right_fraction,
                strongest_support_fraction=strongest_support_fraction,
                support_strength=support_strength,
                conflict_classification=conflict_classification,
                detail=detail,
            )
        )
    return SupportComparisonReport(
        left_path=left_path,
        right_path=right_path,
        shared_taxa=sorted(shared_taxa),
        strong_support_threshold=strong_support_threshold,
        weak_support_threshold=weak_support_threshold,
        support_disagreement_threshold=support_disagreement_threshold,
        shared_clades=shared_clades,
        conflicting_clades=conflicting_clades,
    )
