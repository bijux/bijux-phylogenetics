from __future__ import annotations

from bijux_phylogenetics.compare.influence import TaxonInfluenceReport
from bijux_phylogenetics.compare.topology import (
    BranchLengthComparisonReport,
    SupportComparisonReport,
    TreeComparisonReport,
)
from bijux_phylogenetics.engines.validation import InferenceTreeComparisonReport
from bijux_phylogenetics.io.iqtree_support import support_fraction

from .contracts import (
    InferenceComparisonConclusionRow,
    InferenceComparisonConclusionSummary,
    InferenceComparisonConflictRow,
    InferenceComparisonSharedCladeRow,
    InferenceComparisonWeightedConflictRow,
)

_SUPPORT_DISAGREEMENT_THRESHOLD = 0.15


def _normalize_fasttree_support(value: float | None) -> float | None:
    return value


def _normalize_iqtree_support(value: float | None) -> float | None:
    return support_fraction(value)


def _conflict_weight(
    *,
    fasttree_support_fraction: float | None,
    iqtree_support_fraction: float | None,
) -> float | None:
    values = [
        value
        for value in (fasttree_support_fraction, iqtree_support_fraction)
        if value is not None
    ]
    if not values:
        return None
    return sum(values) / len(values)


def _strongest_support(
    *,
    fasttree_support_fraction: float | None,
    iqtree_support_fraction: float | None,
) -> float | None:
    values = [
        value
        for value in (fasttree_support_fraction, iqtree_support_fraction)
        if value is not None
    ]
    if not values:
        return None
    return max(values)


def _stable_evidence_class(
    fasttree_support_fraction: float | None,
    iqtree_support_fraction: float | None,
) -> str:
    if fasttree_support_fraction is None or iqtree_support_fraction is None:
        return "support_unavailable"
    minimum_support = min(fasttree_support_fraction, iqtree_support_fraction)
    if minimum_support >= 0.9:
        return "strong_stable"
    if minimum_support >= 0.7:
        return "moderate_stable"
    return "weakly_stable"


def _shared_disagreement_severity(
    fasttree_support_fraction: float | None,
    iqtree_support_fraction: float | None,
) -> str:
    strongest_support = _strongest_support(
        fasttree_support_fraction=fasttree_support_fraction,
        iqtree_support_fraction=iqtree_support_fraction,
    )
    weakest_support = None
    values = [
        value
        for value in (fasttree_support_fraction, iqtree_support_fraction)
        if value is not None
    ]
    if values:
        weakest_support = min(values)
    if strongest_support is None:
        return "support_unavailable"
    if (
        strongest_support >= 0.9
        and weakest_support is not None
        and weakest_support >= 0.7
    ):
        return "high_support_disagreement"
    if strongest_support >= 0.7:
        return "moderate_support_disagreement"
    return "low_support_disagreement"


def _shared_disagreement_detail(row: InferenceComparisonSharedCladeRow) -> str:
    if row.support_fraction_delta is None:
        return "shared clade support could not be ranked because one or both trees lacked comparable support values"
    return (
        "shared clade is present in both trees but normalized support differs by "
        f"{row.support_fraction_delta:.3f}"
    )


def build_inference_comparison_shared_clade_rows(
    comparison: InferenceTreeComparisonReport,
    *,
    support_disagreement_threshold: float = _SUPPORT_DISAGREEMENT_THRESHOLD,
) -> list[InferenceComparisonSharedCladeRow]:
    """Convert one compared engine pair into shared-clade rows."""
    rows: list[InferenceComparisonSharedCladeRow] = []
    for pair in comparison.support.shared_clades:
        fasttree_support_fraction = _normalize_fasttree_support(pair.left_support)
        iqtree_support_fraction = _normalize_iqtree_support(pair.right_support)
        support_fraction_delta = (
            None
            if fasttree_support_fraction is None or iqtree_support_fraction is None
            else abs(fasttree_support_fraction - iqtree_support_fraction)
        )
        rows.append(
            InferenceComparisonSharedCladeRow(
                split_id=pair.split_id,
                fasttree_support=pair.left_support,
                fasttree_support_fraction=fasttree_support_fraction,
                fasttree_support_label_kind="sh-like-local-support",
                iqtree_support=pair.right_support,
                iqtree_support_fraction=iqtree_support_fraction,
                iqtree_support_label_kind="ufboot-support",
                support_fraction_delta=support_fraction_delta,
                support_disagreement=(
                    support_fraction_delta is not None
                    and support_fraction_delta >= support_disagreement_threshold
                ),
            )
        )
    return rows


def build_inference_comparison_conflict_rows(
    comparison: InferenceTreeComparisonReport,
    comparison_report: object | None = None,
    *,
    support_disagreement_threshold: float = _SUPPORT_DISAGREEMENT_THRESHOLD,
) -> list[InferenceComparisonConflictRow]:
    """Build one combined topology-plus-support conflict ledger."""
    del comparison_report
    rows: list[InferenceComparisonConflictRow] = []
    shared_rows = build_inference_comparison_shared_clade_rows(
        comparison,
        support_disagreement_threshold=support_disagreement_threshold,
    )
    for support_conflict in comparison.support.conflicting_clades:
        if support_conflict.comparison_status == "left_only":
            conflict_kind = "fasttree_only"
            fasttree_present = True
            iqtree_present = False
        elif support_conflict.comparison_status == "right_only":
            conflict_kind = "iqtree_only"
            fasttree_present = False
            iqtree_present = True
        else:
            continue
        rows.append(
            InferenceComparisonConflictRow(
                split_id=support_conflict.split_id,
                conflict_kind=conflict_kind,
                fasttree_present=fasttree_present,
                iqtree_present=iqtree_present,
                fasttree_support=support_conflict.left_support,
                fasttree_support_fraction=support_conflict.left_support_fraction,
                iqtree_support=support_conflict.right_support,
                iqtree_support_fraction=support_conflict.right_support_fraction,
                detail=support_conflict.detail,
            )
        )
    for row in shared_rows:
        if not row.support_disagreement:
            continue
        rows.append(
            InferenceComparisonConflictRow(
                split_id=row.split_id,
                conflict_kind="support_disagreement",
                fasttree_present=True,
                iqtree_present=True,
                fasttree_support=row.fasttree_support,
                fasttree_support_fraction=row.fasttree_support_fraction,
                iqtree_support=row.iqtree_support,
                iqtree_support_fraction=row.iqtree_support_fraction,
                detail=(
                    "normalized support fractions differ by at least "
                    f"{support_disagreement_threshold:.2f}"
                ),
            )
        )
    return rows


def build_inference_comparison_weighted_conflict_rows(
    support_report: SupportComparisonReport,
) -> list[InferenceComparisonWeightedConflictRow]:
    """Rank disagreements by explicit support weight and severity."""
    rows: list[InferenceComparisonWeightedConflictRow] = []
    for row in support_report.shared_clades:
        if not row.support_disagreement:
            continue
        severity_class = _shared_disagreement_severity(
            row.left_support_fraction,
            row.right_support_fraction,
        )
        rows.append(
            InferenceComparisonWeightedConflictRow(
                split_id=row.split_id,
                comparison_status="shared",
                conflict_kind="support_disagreement",
                severity_class=severity_class,
                fasttree_support_fraction=row.left_support_fraction,
                iqtree_support_fraction=row.right_support_fraction,
                support_fraction_delta=row.support_fraction_delta,
                strongest_support_fraction=_strongest_support(
                    fasttree_support_fraction=row.left_support_fraction,
                    iqtree_support_fraction=row.right_support_fraction,
                ),
                support_weight=_conflict_weight(
                    fasttree_support_fraction=row.left_support_fraction,
                    iqtree_support_fraction=row.right_support_fraction,
                ),
                serious_conflict=severity_class == "high_support_disagreement",
                detail=_shared_disagreement_detail(row),
            )
        )
    for row in support_report.conflicting_clades:
        comparison_status = (
            "fasttree_only" if row.comparison_status == "left_only" else "iqtree_only"
        )
        rows.append(
            InferenceComparisonWeightedConflictRow(
                split_id=row.split_id,
                comparison_status=comparison_status,
                conflict_kind="topology_conflict",
                severity_class=row.conflict_classification,
                fasttree_support_fraction=row.left_support_fraction,
                iqtree_support_fraction=row.right_support_fraction,
                support_fraction_delta=None,
                strongest_support_fraction=row.strongest_support_fraction,
                support_weight=_conflict_weight(
                    fasttree_support_fraction=row.left_support_fraction,
                    iqtree_support_fraction=row.right_support_fraction,
                ),
                serious_conflict=row.conflict_classification == "high_support_conflict",
                detail=row.detail,
            )
        )
    return sorted(
        rows,
        key=lambda row: (
            -(row.support_weight or -1.0),
            not row.serious_conflict,
            row.comparison_status,
            row.split_id,
        ),
    )


def build_inference_comparison_conclusion_rows(
    support_report: SupportComparisonReport,
) -> list[InferenceComparisonConclusionRow]:
    """Summarize clade stability across the two compared inference engines."""
    rows: list[InferenceComparisonConclusionRow] = []
    for row in support_report.shared_clades:
        if row.support_disagreement:
            evidence_class = _shared_disagreement_severity(
                row.left_support_fraction,
                row.right_support_fraction,
            )
            rows.append(
                InferenceComparisonConclusionRow(
                    split_id=row.split_id,
                    conclusion_class="unstable_clade",
                    evidence_class=evidence_class,
                    comparison_status="shared",
                    fasttree_present=True,
                    iqtree_present=True,
                    fasttree_support_fraction=row.left_support_fraction,
                    iqtree_support_fraction=row.right_support_fraction,
                    support_fraction_delta=row.support_fraction_delta,
                    serious_conflict=evidence_class == "high_support_disagreement",
                    detail=_shared_disagreement_detail(row),
                )
            )
            continue
        rows.append(
            InferenceComparisonConclusionRow(
                split_id=row.split_id,
                conclusion_class="stable_clade",
                evidence_class=_stable_evidence_class(
                    row.left_support_fraction,
                    row.right_support_fraction,
                ),
                comparison_status="shared",
                fasttree_present=True,
                iqtree_present=True,
                fasttree_support_fraction=row.left_support_fraction,
                iqtree_support_fraction=row.right_support_fraction,
                support_fraction_delta=row.support_fraction_delta,
                serious_conflict=False,
                detail="clade is present in both trees without a support disagreement above the governed threshold",
            )
        )
    for row in support_report.conflicting_clades:
        rows.append(
            InferenceComparisonConclusionRow(
                split_id=row.split_id,
                conclusion_class="engine_specific_clade",
                evidence_class=row.conflict_classification,
                comparison_status=(
                    "fasttree_only"
                    if row.comparison_status == "left_only"
                    else "iqtree_only"
                ),
                fasttree_present=row.left_present,
                iqtree_present=row.right_present,
                fasttree_support_fraction=row.left_support_fraction,
                iqtree_support_fraction=row.right_support_fraction,
                support_fraction_delta=None,
                serious_conflict=row.conflict_classification == "high_support_conflict",
                detail=row.detail,
            )
        )
    return sorted(
        rows,
        key=lambda row: (
            row.conclusion_class != "stable_clade",
            not row.serious_conflict,
            row.split_id,
        ),
    )


def summarize_inference_comparison_conclusions(
    topology_report: TreeComparisonReport,
    branch_length_report: BranchLengthComparisonReport,
    *,
    weighted_conflict_rows: list[InferenceComparisonWeightedConflictRow],
    conclusion_rows: list[InferenceComparisonConclusionRow],
    taxon_influence_report: TaxonInfluenceReport | None,
) -> InferenceComparisonConclusionSummary:
    """Build one compact summary of stable versus fragile biological conclusions."""
    top_conflict_driver_taxa = (
        []
        if taxon_influence_report is None
        else [
            row.taxon
            for row in taxon_influence_report.rows
            if row.influence_score > 0.0
            and (
                row.conflicting_clade_delta < 0
                or row.high_support_conflict_delta < 0
                or row.support_disagreement_delta < 0
            )
        ][:3]
    )
    return InferenceComparisonConclusionSummary(
        shared_taxa_count=len(topology_report.shared_taxa),
        robinson_foulds_distance=topology_report.robinson_foulds_distance,
        normalized_robinson_foulds=topology_report.normalized_robinson_foulds,
        branch_score_distance=branch_length_report.branch_score.branch_score_distance,
        stable_clade_count=sum(
            1 for row in conclusion_rows if row.conclusion_class == "stable_clade"
        ),
        unstable_clade_count=sum(
            1 for row in conclusion_rows if row.conclusion_class == "unstable_clade"
        ),
        engine_specific_clade_count=sum(
            1
            for row in conclusion_rows
            if row.conclusion_class == "engine_specific_clade"
        ),
        support_weighted_conflict_count=len(weighted_conflict_rows),
        low_support_disagreement_count=sum(
            1
            for row in weighted_conflict_rows
            if row.severity_class == "low_support_disagreement"
        ),
        moderate_support_disagreement_count=sum(
            1
            for row in weighted_conflict_rows
            if row.severity_class == "moderate_support_disagreement"
        ),
        high_support_conflict_count=sum(
            1
            for row in weighted_conflict_rows
            if row.severity_class == "high_support_conflict"
        ),
        high_support_disagreement_count=sum(
            1
            for row in weighted_conflict_rows
            if row.severity_class == "high_support_disagreement"
        ),
        serious_conflict_count=sum(
            1 for row in weighted_conflict_rows if row.serious_conflict
        ),
        top_conflict_driver_taxa=top_conflict_driver_taxa,
    )
