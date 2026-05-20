"""Ancestral-state reconstruction methods and helpers."""

from __future__ import annotations

from importlib import import_module

_PUBLIC_SURFACES = (
    (
        "confidence",
        (
            "AncestralConfidenceSummary",
            "ContinuousAncestralConfidenceRow",
            "ContinuousAncestralTreeSetConfidenceRow",
            "DiscreteAncestralConfidenceRow",
            "DiscreteAncestralTreeSetConfidenceRow",
            "build_continuous_ancestral_confidence_rows",
            "build_continuous_ancestral_tree_set_confidence_rows",
            "build_discrete_ancestral_confidence_rows",
            "build_discrete_ancestral_tree_set_confidence_rows",
            "summarize_continuous_ancestral_confidence",
            "summarize_continuous_ancestral_tree_set_confidence",
            "summarize_discrete_ancestral_confidence",
            "summarize_discrete_ancestral_tree_set_confidence",
            "write_ancestral_confidence_summary_table",
            "write_continuous_ancestral_confidence_table",
            "write_continuous_ancestral_tree_set_confidence_table",
            "write_discrete_ancestral_confidence_table",
            "write_discrete_ancestral_tree_set_confidence_table",
        ),
    ),
    (
        "continuous",
        (
            "ContinuousAncestralBrownianFitDiagnostics",
            "ContinuousAncestralEstimate",
            "ContinuousAncestralExclusion",
            "ContinuousAncestralReport",
            "ContinuousAncestralSummary",
            "continuous_ancestral_exclusions",
            "reconstruct_continuous_ancestral_states",
            "reconstruct_continuous_ancestral_states_from_dataset",
            "summarize_continuous_ancestral_report",
            "write_continuous_ancestral_exclusion_table",
            "write_continuous_ancestral_summary_table",
            "write_continuous_ancestral_uncertainty_table",
        ),
    ),
    (
        "discrete",
        (
            "DiscreteAncestralEstimate",
            "DiscreteAncestralExclusion",
            "DiscreteAncestralReport",
            "DiscreteAncestralSummary",
            "DiscreteModelBaselineComparison",
            "DiscreteOptimizerDiagnostics",
            "DiscreteRerootingMethodCompatibility",
            "DiscreteTransitionRateRow",
            "discrete_ancestral_exclusions",
            "reconstruct_discrete_ancestral_states",
            "reconstruct_discrete_ancestral_states_from_dataset",
            "summarize_discrete_ancestral_report",
            "write_discrete_ancestral_exclusion_table",
            "write_discrete_ancestral_fit_table",
            "write_discrete_ancestral_probability_table",
            "write_discrete_ancestral_summary_table",
            "write_discrete_ancestral_transition_table",
        ),
    ),
    (
        "discrete_reference",
        (
            "DiscreteAncestralReferenceObservation",
            "DiscreteAncestralReferenceProbabilityRow",
            "DiscreteAncestralReferenceValidationReport",
            "validate_discrete_ancestral_reference_examples",
        ),
    ),
    (
        "irreversible_discrete",
        (
            "IrreversibleDiscreteFitRow",
            "IrreversibleDiscreteNodeComparisonRow",
            "IrreversibleDiscreteReport",
            "IrreversibleDiscreteSummary",
            "IrreversibleDiscreteTransitionComparisonRow",
            "summarize_irreversible_discrete_reconstruction",
            "summarize_irreversible_discrete_report",
            "write_irreversible_discrete_fit_table",
            "write_irreversible_discrete_node_table",
            "write_irreversible_discrete_summary_table",
            "write_irreversible_discrete_transition_table",
        ),
    ),
    (
        "ordered_discrete",
        (
            "OrderedDiscreteFitRow",
            "OrderedDiscreteNodeComparisonRow",
            "OrderedDiscreteReport",
            "OrderedDiscreteSummary",
            "OrderedDiscreteTransitionComparisonRow",
            "summarize_ordered_discrete_reconstruction",
            "summarize_ordered_discrete_report",
            "write_ordered_discrete_fit_table",
            "write_ordered_discrete_node_table",
            "write_ordered_discrete_summary_table",
            "write_ordered_discrete_transition_table",
        ),
    ),
    (
        "presentation.methods_text",
        (
            "AncestralMethodsSummaryTextResult",
            "build_ancestral_methods_summary_text",
            "write_ancestral_methods_summary_text",
        ),
    ),
    (
        "presentation.review_bundle",
        (
            "AncestralContinuousChangeBranchRow",
            "AncestralContinuousChangeCountRow",
            "AncestralReportPackageResult",
            "build_ancestral_report_package",
            "summarize_continuous_change_branches",
            "summarize_continuous_change_counts",
            "write_continuous_change_branch_table",
            "write_continuous_change_count_table",
        ),
    ),
    (
        "sensitivity",
        (
            "AncestralSensitivityReport",
            "AncestralSensitivitySummary",
            "RootSensitivityAssumptionRow",
            "RootSensitivityNodeRow",
            "RootSensitivityReport",
            "RootSensitivitySummary",
            "build_ancestral_sensitivity_report",
            "summarize_ancestral_root_sensitivity",
            "summarize_ancestral_root_sensitivity_report",
            "write_ancestral_root_assumption_table",
            "write_ancestral_root_sensitivity_node_table",
            "write_ancestral_root_sensitivity_summary_table",
        ),
    ),
    (
        "service",
        (
            "AncestralStateReportBuildResult",
            "ContinuousAncestralComparisonReport",
            "ContinuousAncestralComparisonRow",
            "ContinuousEvolutionaryModeAncestralReport",
            "DiscreteAncestralPairComparisonReport",
            "DiscreteAncestralPairComparisonRow",
            "compare_continuous_ancestral_models",
            "compare_discrete_ancestral_reconstructions",
            "reconstruct_continuous_evolutionary_mode_states",
            "render_ancestral_state_report",
            "render_ancestral_state_tree",
            "write_ancestral_state_table",
            "write_discrete_ancestral_comparison_table",
        ),
    ),
    (
        "transitions",
        (
            "AncestralTransitionBranchRow",
            "AncestralTransitionCountRow",
            "AncestralTransitionExclusion",
            "AncestralTransitionReport",
            "AncestralTransitionSummary",
            "AncestralTransitionTreeRow",
            "AncestralTransitionTreeSetBranchRow",
            "AncestralTransitionTreeSetCountRow",
            "AncestralTransitionTreeSetReport",
            "AncestralTransitionTreeSetSummary",
            "summarize_ancestral_transition_report",
            "summarize_ancestral_transition_tree_set",
            "summarize_ancestral_transition_tree_set_report",
            "summarize_ancestral_transitions",
            "write_ancestral_transition_branch_table",
            "write_ancestral_transition_count_table",
            "write_ancestral_transition_exclusion_table",
            "write_ancestral_transition_summary_table",
            "write_ancestral_transition_tree_set_branch_table",
            "write_ancestral_transition_tree_set_count_table",
            "write_ancestral_transition_tree_set_summary_table",
            "write_ancestral_transition_tree_set_tree_table",
        ),
    ),
    (
        "tree_set",
        (
            "AncestralTreeSetExclusion",
            "AncestralTreeSetTreeRow",
            "ContinuousAncestralTreeSetCladeSummaryRow",
            "ContinuousAncestralTreeSetNodeRow",
            "ContinuousAncestralTreeSetReport",
            "ContinuousAncestralTreeSetSummary",
            "DiscreteAncestralTreeSetCladeSummaryRow",
            "DiscreteAncestralTreeSetNodeRow",
            "DiscreteAncestralTreeSetReport",
            "DiscreteAncestralTreeSetSummary",
            "summarize_continuous_ancestral_tree_set",
            "summarize_continuous_ancestral_tree_set_report",
            "summarize_discrete_ancestral_tree_set",
            "summarize_discrete_ancestral_tree_set_report",
            "write_ancestral_tree_set_exclusion_table",
            "write_ancestral_tree_set_tree_table",
            "write_continuous_ancestral_tree_set_clade_table",
            "write_continuous_ancestral_tree_set_node_table",
            "write_continuous_ancestral_tree_set_summary_table",
            "write_discrete_ancestral_tree_set_clade_table",
            "write_discrete_ancestral_tree_set_node_table",
            "write_discrete_ancestral_tree_set_summary_table",
        ),
    ),
    (
        "presentation.figure_bundle",
        (
            "AncestralFigurePackageResult",
            "build_ancestral_figure_package",
        ),
    ),
    (
        "presentation.visualization",
        (
            "AncestralVisualizationResult",
            "render_ancestral_state_visualization",
        ),
    ),
)

__all__ = [name for _, names in _PUBLIC_SURFACES for name in names]

_NAME_TO_MODULE = {
    name: module_name
    for module_name, names in _PUBLIC_SURFACES
    for name in names
}


def __getattr__(name: str):
    """Resolve ancestral exports lazily from their owning submodules."""
    try:
        module_name = _NAME_TO_MODULE[name]
    except KeyError as error:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from error
    value = getattr(import_module(f"{__name__}.{module_name}"), name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
