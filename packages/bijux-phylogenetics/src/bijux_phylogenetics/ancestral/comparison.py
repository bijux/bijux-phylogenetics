from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.ancestral.common import write_ancestral_rows
from bijux_phylogenetics.ancestral.continuous import (
    reconstruct_continuous_ancestral_states,
)
from bijux_phylogenetics.ancestral.discrete import (
    reconstruct_discrete_ancestral_states,
)
from bijux_phylogenetics.comparative.discrete_evolution import (
    run_discrete_state_transition_model,
)


@dataclass(slots=True)
class ContinuousAncestralComparisonRow:
    """Node-by-node comparison between two continuous ancestral reconstructions."""

    node: str
    descendant_taxa: list[str]
    left_estimate: float
    right_estimate: float
    estimate_delta: float
    intervals_overlap: bool


@dataclass(slots=True)
class ContinuousAncestralComparisonReport:
    """Comparison between two continuous ancestral-state models."""

    tree_path: Path
    traits_path: Path
    trait: str
    taxon_count: int
    left_model: str
    right_model: str
    rows: list[ContinuousAncestralComparisonRow]


@dataclass(slots=True)
class DiscreteAncestralModelComparisonRow:
    """Model-comparison summary for one discrete ancestral-state model."""

    model: str
    parameter_count: int
    pseudo_log_likelihood: float
    aic: float
    selected: bool


@dataclass(slots=True)
class DiscreteAncestralModelDifference:
    """Node-level state difference between the selected and another discrete model."""

    comparison_model: str
    node: str
    descendant_taxa: list[str]
    selected_state: str
    comparison_state: str
    differs: bool


@dataclass(slots=True)
class DiscreteAncestralModelComparisonReport:
    """Comparison of ER, SYM, and ARD discrete ancestral reconstructions."""

    tree_path: Path
    traits_path: Path
    trait: str
    taxon_count: int
    selected_model: str
    rows: list[DiscreteAncestralModelComparisonRow]
    node_differences: list[DiscreteAncestralModelDifference]


@dataclass(slots=True)
class DiscreteAncestralPairComparisonRow:
    """Node-wise comparison between two discrete ancestral reconstructions."""

    node: str
    descendant_taxa: list[str]
    left_state: str
    right_state: str
    left_state_set: list[str]
    right_state_set: list[str]
    left_confidence: float
    right_confidence: float
    left_ambiguous: bool
    right_ambiguous: bool
    differs: bool
    ambiguity_changed: bool


@dataclass(slots=True)
class DiscreteAncestralPairComparisonReport:
    """Direct comparison between two discrete ancestral reconstruction models."""

    tree_path: Path
    traits_path: Path
    trait: str
    taxon_count: int
    left_model: str
    right_model: str
    left_minimal_change_count: int | None
    right_minimal_change_count: int | None
    differing_node_count: int
    ambiguity_change_count: int
    rows: list[DiscreteAncestralPairComparisonRow]


@dataclass(slots=True)
class ContinuousAncestralTreeComparisonRow:
    """Node-wise comparison of continuous ancestral estimates across two trees."""

    node: str
    descendant_taxa: list[str]
    left_estimate: float
    right_estimate: float
    estimate_delta: float


@dataclass(slots=True)
class ContinuousAncestralTreeComparisonReport:
    """Continuous ancestral reconstruction differences across alternative trees."""

    left_tree_path: Path
    right_tree_path: Path
    traits_path: Path
    trait: str
    model: str
    rows: list[ContinuousAncestralTreeComparisonRow]


@dataclass(slots=True)
class DiscreteAncestralTreeComparisonRow:
    """Node-wise comparison of discrete ancestral states across two trees."""

    node: str
    descendant_taxa: list[str]
    left_state: str
    right_state: str
    differs: bool
    left_confidence: float
    right_confidence: float


@dataclass(slots=True)
class DiscreteAncestralTreeComparisonReport:
    """Discrete ancestral reconstruction differences across alternative trees."""

    left_tree_path: Path
    right_tree_path: Path
    traits_path: Path
    trait: str
    model: str
    rows: list[DiscreteAncestralTreeComparisonRow]


def compare_continuous_ancestral_models(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    left_model: str = "brownian",
    right_model: str = "ou",
    left_alpha: float = 1.0,
    right_alpha: float = 1.0,
) -> ContinuousAncestralComparisonReport:
    """Compare two continuous ancestral reconstructions node by node."""
    left = reconstruct_continuous_ancestral_states(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=left_model,
        alpha=left_alpha,
    )
    right = reconstruct_continuous_ancestral_states(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=right_model,
        alpha=right_alpha,
    )
    right_by_node = {estimate.node: estimate for estimate in right.estimates}
    rows: list[ContinuousAncestralComparisonRow] = []
    for left_estimate in left.estimates:
        if left_estimate.node not in right_by_node:
            continue
        right_estimate = right_by_node[left_estimate.node]
        rows.append(
            ContinuousAncestralComparisonRow(
                node=left_estimate.node,
                descendant_taxa=left_estimate.descendant_taxa,
                left_estimate=left_estimate.estimate,
                right_estimate=right_estimate.estimate,
                estimate_delta=right_estimate.estimate - left_estimate.estimate,
                intervals_overlap=not (
                    left_estimate.upper_95_interval < right_estimate.lower_95_interval
                    or right_estimate.upper_95_interval
                    < left_estimate.lower_95_interval
                ),
            )
        )
    return ContinuousAncestralComparisonReport(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        taxon_count=left.taxon_count,
        left_model=left_model,
        right_model=right_model,
        rows=rows,
    )


def compare_discrete_ancestral_models(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    models: tuple[str, str, str] = ("equal-rates", "symmetric", "all-rates-different"),
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
) -> DiscreteAncestralModelComparisonReport:
    """Compare supported discrete ancestral-state likelihood models and select the lowest-AIC fit."""
    reconstructions = {
        model: reconstruct_discrete_ancestral_states(
            tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
            model=model,
            state_ordering=state_ordering,
            ordered_states=ordered_states,
        )
        for model in models
    }
    transition_reports = {
        model: run_discrete_state_transition_model(
            tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
            model=model,
            state_ordering=state_ordering,
            ordered_states=ordered_states,
        )
        for model in models
    }
    rows = [
        DiscreteAncestralModelComparisonRow(
            model=model,
            parameter_count=transition_reports[model].transition_model.parameter_count,
            pseudo_log_likelihood=transition_reports[
                model
            ].transition_model.pseudo_log_likelihood,
            aic=transition_reports[model].transition_model.aic,
            selected=False,
        )
        for model in models
    ]
    selected_model = min(rows, key=lambda row: row.aic).model
    for row in rows:
        row.selected = row.model == selected_model
    selected_report = reconstructions[selected_model]
    selected_by_node = {
        estimate.node: estimate for estimate in selected_report.estimates
    }
    node_differences: list[DiscreteAncestralModelDifference] = []
    for model, report in reconstructions.items():
        if model == selected_model:
            continue
        for estimate in report.estimates:
            selected_estimate = selected_by_node[estimate.node]
            node_differences.append(
                DiscreteAncestralModelDifference(
                    comparison_model=model,
                    node=estimate.node,
                    descendant_taxa=estimate.descendant_taxa,
                    selected_state=selected_estimate.most_likely_state,
                    comparison_state=estimate.most_likely_state,
                    differs=selected_estimate.most_likely_state
                    != estimate.most_likely_state,
                )
            )
    return DiscreteAncestralModelComparisonReport(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        taxon_count=selected_report.taxon_count,
        selected_model=selected_model,
        rows=rows,
        node_differences=node_differences,
    )


def compare_discrete_ancestral_reconstructions(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    left_model: str = "fitch",
    right_model: str = "equal-rates",
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
    root_prior_mode: str = "equal",
    fixed_root_state: str | None = None,
) -> DiscreteAncestralPairComparisonReport:
    """Compare two discrete ancestral reconstructions node by node."""
    if left_model == right_model:
        raise ValueError("discrete ancestral comparison requires distinct models")
    left = reconstruct_discrete_ancestral_states(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=left_model,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
        root_prior_mode=root_prior_mode,
        fixed_root_state=fixed_root_state,
    )
    right = reconstruct_discrete_ancestral_states(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=right_model,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
        root_prior_mode=root_prior_mode,
        fixed_root_state=fixed_root_state,
    )
    right_by_node = {
        estimate.node: estimate for estimate in right.estimates if not estimate.is_tip
    }
    rows: list[DiscreteAncestralPairComparisonRow] = []
    differing_node_count = 0
    ambiguity_change_count = 0
    for left_estimate in left.estimates:
        if left_estimate.is_tip or left_estimate.node not in right_by_node:
            continue
        right_estimate = right_by_node[left_estimate.node]
        differs = left_estimate.most_likely_state != right_estimate.most_likely_state
        ambiguity_changed = left_estimate.ambiguous != right_estimate.ambiguous
        if differs:
            differing_node_count += 1
        if ambiguity_changed:
            ambiguity_change_count += 1
        rows.append(
            DiscreteAncestralPairComparisonRow(
                node=left_estimate.node,
                descendant_taxa=left_estimate.descendant_taxa,
                left_state=left_estimate.most_likely_state,
                right_state=right_estimate.most_likely_state,
                left_state_set=left_estimate.state_set,
                right_state_set=right_estimate.state_set,
                left_confidence=left_estimate.confidence,
                right_confidence=right_estimate.confidence,
                left_ambiguous=left_estimate.ambiguous,
                right_ambiguous=right_estimate.ambiguous,
                differs=differs,
                ambiguity_changed=ambiguity_changed,
            )
        )
    return DiscreteAncestralPairComparisonReport(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        taxon_count=left.taxon_count,
        left_model=left.model,
        right_model=right.model,
        left_minimal_change_count=left.minimal_change_count,
        right_minimal_change_count=right.minimal_change_count,
        differing_node_count=differing_node_count,
        ambiguity_change_count=ambiguity_change_count,
        rows=rows,
    )


def write_discrete_ancestral_comparison_table(
    path: Path,
    report: DiscreteAncestralPairComparisonReport,
) -> Path:
    """Write a node-wise discrete ancestral comparison ledger."""
    return write_ancestral_rows(
        path,
        columns=[
            "node",
            "descendant_taxa",
            "left_model",
            "right_model",
            "left_state",
            "right_state",
            "left_state_set",
            "right_state_set",
            "left_confidence",
            "right_confidence",
            "left_ambiguous",
            "right_ambiguous",
            "differs",
            "ambiguity_changed",
        ],
        rows=[
            {
                "node": row.node,
                "descendant_taxa": ",".join(row.descendant_taxa),
                "left_model": report.left_model,
                "right_model": report.right_model,
                "left_state": row.left_state,
                "right_state": row.right_state,
                "left_state_set": ",".join(row.left_state_set),
                "right_state_set": ",".join(row.right_state_set),
                "left_confidence": str(row.left_confidence),
                "right_confidence": str(row.right_confidence),
                "left_ambiguous": str(row.left_ambiguous).lower(),
                "right_ambiguous": str(row.right_ambiguous).lower(),
                "differs": str(row.differs).lower(),
                "ambiguity_changed": str(row.ambiguity_changed).lower(),
            }
            for row in report.rows
        ],
    )


def compare_continuous_ancestral_trees(
    left_tree_path: Path,
    right_tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    model: str = "brownian",
    alpha: float = 1.0,
) -> ContinuousAncestralTreeComparisonReport:
    """Compare continuous ancestral estimates across two alternative trees."""
    left = reconstruct_continuous_ancestral_states(
        left_tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=model,
        alpha=alpha,
    )
    right = reconstruct_continuous_ancestral_states(
        right_tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=model,
        alpha=alpha,
    )
    right_by_node = {
        estimate.node: estimate for estimate in right.estimates if not estimate.is_tip
    }
    rows = [
        ContinuousAncestralTreeComparisonRow(
            node=estimate.node,
            descendant_taxa=estimate.descendant_taxa,
            left_estimate=estimate.estimate,
            right_estimate=right_by_node[estimate.node].estimate,
            estimate_delta=right_by_node[estimate.node].estimate - estimate.estimate,
        )
        for estimate in left.estimates
        if not estimate.is_tip and estimate.node in right_by_node
    ]
    return ContinuousAncestralTreeComparisonReport(
        left_tree_path=left_tree_path,
        right_tree_path=right_tree_path,
        traits_path=traits_path,
        trait=trait,
        model=model,
        rows=rows,
    )


def compare_discrete_ancestral_trees(
    left_tree_path: Path,
    right_tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    model: str = "fitch",
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
) -> DiscreteAncestralTreeComparisonReport:
    """Compare discrete ancestral states across two alternative trees."""
    left = reconstruct_discrete_ancestral_states(
        left_tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=model,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
    )
    right = reconstruct_discrete_ancestral_states(
        right_tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=model,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
    )
    right_by_node = {
        estimate.node: estimate for estimate in right.estimates if not estimate.is_tip
    }
    rows = [
        DiscreteAncestralTreeComparisonRow(
            node=estimate.node,
            descendant_taxa=estimate.descendant_taxa,
            left_state=estimate.most_likely_state,
            right_state=right_by_node[estimate.node].most_likely_state,
            differs=estimate.most_likely_state
            != right_by_node[estimate.node].most_likely_state,
            left_confidence=estimate.confidence,
            right_confidence=right_by_node[estimate.node].confidence,
        )
        for estimate in left.estimates
        if not estimate.is_tip and estimate.node in right_by_node
    ]
    return DiscreteAncestralTreeComparisonReport(
        left_tree_path=left_tree_path,
        right_tree_path=right_tree_path,
        traits_path=traits_path,
        trait=trait,
        model=model,
        rows=rows,
    )
