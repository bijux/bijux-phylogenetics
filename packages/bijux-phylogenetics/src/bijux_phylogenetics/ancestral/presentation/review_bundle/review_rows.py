from __future__ import annotations

import json

from bijux_phylogenetics.ancestral.continuous import ContinuousAncestralReport
from bijux_phylogenetics.ancestral.discrete import DiscreteAncestralReport
from bijux_phylogenetics.ancestral.discrete.review import (
    AncestralTransitionBranchRow,
    AncestralTransitionCountRow,
)

from .contracts import (
    AncestralContinuousChangeBranchRow,
    AncestralContinuousChangeCountRow,
)


def node_table_rows(
    report: ContinuousAncestralReport | DiscreteAncestralReport,
) -> list[list[str]]:
    if isinstance(report, ContinuousAncestralReport):
        return [
            [
                estimate.node,
                "" if estimate.node_name is None else estimate.node_name,
                "true" if estimate.is_tip else "false",
                ", ".join(estimate.descendant_taxa),
                format(estimate.estimate, ".6g"),
                format(estimate.standard_error, ".6g"),
            ]
            for estimate in report.estimates
            if not estimate.is_tip
        ]
    return [
        [
            estimate.node,
            "" if estimate.node_name is None else estimate.node_name,
            "true" if estimate.is_tip else "false",
            ", ".join(estimate.descendant_taxa),
            estimate.most_likely_state,
            format(estimate.confidence, ".6g"),
        ]
        for estimate in report.estimates
        if not estimate.is_tip
    ]


def uncertainty_table_rows(
    report: ContinuousAncestralReport | DiscreteAncestralReport,
) -> list[list[str]]:
    if isinstance(report, ContinuousAncestralReport):
        return [
            [
                estimate.node,
                ", ".join(estimate.descendant_taxa),
                (
                    f"se={format(estimate.standard_error, '.6g')}; "
                    f"95%=[{format(estimate.lower_95_interval, '.6g')}, {format(estimate.upper_95_interval, '.6g')}]"
                ),
                estimate.interpretation,
            ]
            for estimate in report.estimates
            if not estimate.is_tip
        ]
    return [
        [
            estimate.node,
            ", ".join(estimate.descendant_taxa),
            (
                f"state={estimate.most_likely_state}; "
                f"confidence={format(estimate.confidence, '.6g')}; "
                f"probabilities={json.dumps(estimate.state_probabilities, sort_keys=True)}"
            ),
            estimate.interpretation,
        ]
        for estimate in report.estimates
        if not estimate.is_tip
    ]


def transition_count_table_rows(
    reconstruction_kind: str,
    *,
    count_rows: list[AncestralContinuousChangeCountRow]
    | list[AncestralTransitionCountRow],
) -> list[list[str]]:
    if reconstruction_kind == "continuous":
        continuous_rows = count_rows
        if not all(
            isinstance(row, AncestralContinuousChangeCountRow)
            for row in continuous_rows
        ):
            raise RuntimeError(
                "continuous ancestral report package received non-continuous change counts"
            )
        return [
            [
                row.direction,
                str(row.branch_count),
                format(row.branch_fraction, ".6g"),
                (
                    f"mean_delta={format(row.mean_delta, '.6g')}; "
                    f"range=[{format(row.minimum_delta, '.6g')}, {format(row.maximum_delta, '.6g')}]"
                ),
            ]
            for row in continuous_rows
        ]
    discrete_rows = count_rows
    if not all(isinstance(row, AncestralTransitionCountRow) for row in discrete_rows):
        raise RuntimeError(
            "discrete ancestral report package received non-discrete transition counts"
        )
    return [
        [
            row.transition,
            str(row.total_change_count),
            (
                "certain"
                if row.uncertain_change_count == 0
                else "mixed"
                if row.certain_change_count > 0
                else "uncertain"
            ),
            (
                f"certain={row.certain_change_count}; "
                f"uncertain={row.uncertain_change_count}"
            ),
        ]
        for row in discrete_rows
    ]


def transition_branch_table_rows(
    reconstruction_kind: str,
    *,
    branch_rows: list[AncestralContinuousChangeBranchRow]
    | list[AncestralTransitionBranchRow],
) -> list[list[str]]:
    if reconstruction_kind == "continuous":
        continuous_rows = branch_rows
        if not all(
            isinstance(row, AncestralContinuousChangeBranchRow)
            for row in continuous_rows
        ):
            raise RuntimeError(
                "continuous ancestral report package received non-continuous branch rows"
            )
        return [
            [
                row.parent_node,
                row.child_node,
                ", ".join(row.child_descendant_taxa),
                "" if row.branch_length is None else format(row.branch_length, ".6g"),
                row.direction,
                (
                    f"delta={format(row.delta, '.6g')}; "
                    f"parent={format(row.parent_estimate, '.6g')}; "
                    f"child={format(row.child_estimate, '.6g')}"
                ),
            ]
            for row in continuous_rows
        ]
    discrete_rows = branch_rows
    if not all(isinstance(row, AncestralTransitionBranchRow) for row in discrete_rows):
        raise RuntimeError(
            "discrete ancestral report package received non-discrete branch rows"
        )
    return [
        [
            row.parent_node,
            row.child_node,
            ", ".join(row.child_descendant_taxa),
            "" if row.branch_length is None else format(row.branch_length, ".6g"),
            row.transition,
            row.certainty_class,
        ]
        for row in discrete_rows
        if row.changed
    ]


def limitations(
    reconstruction_kind: str,
    report: ContinuousAncestralReport | DiscreteAncestralReport,
) -> list[str]:
    review_limitations = list(report.warnings)
    if reconstruction_kind == "continuous":
        review_limitations.append(
            "continuous branch-change counts summarize direction of reconstructed value shifts, not discrete state-transition events"
        )
    else:
        review_limitations.append(
            "discrete transition counts are reconstructed branchwise review evidence, not stochastic mapping event histories"
        )
    return review_limitations
