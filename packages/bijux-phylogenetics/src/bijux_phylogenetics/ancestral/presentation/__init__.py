from __future__ import annotations

from importlib import import_module

_PUBLIC_SURFACES = (
    (
        "confidence_review",
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
        "report_rendering",
        (
            "AncestralStateReportBuildResult",
            "render_ancestral_state_report",
            "render_ancestral_state_tree",
            "write_ancestral_state_table",
        ),
    ),
    (
        "figure_bundle",
        (
            "AncestralFigurePackageResult",
            "build_ancestral_figure_package",
        ),
    ),
    (
        "visualization",
        (
            "AncestralVisualizationResult",
            "render_ancestral_state_visualization",
        ),
    ),
    (
        "methods_text",
        (
            "AncestralMethodsSummaryTextResult",
            "build_ancestral_methods_summary_text",
            "write_ancestral_methods_summary_text",
        ),
    ),
    (
        "review_bundle",
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
)

__all__ = [name for _, names in _PUBLIC_SURFACES for name in names]

_NAME_TO_MODULE = {
    name: module_name for module_name, names in _PUBLIC_SURFACES for name in names
}


def __getattr__(name: str):
    """Resolve presentation exports lazily from their owning submodules."""
    try:
        module_name = _NAME_TO_MODULE[name]
    except KeyError as error:
        raise AttributeError(
            f"module {__name__!r} has no attribute {name!r}"
        ) from error
    value = getattr(import_module(f"{__name__}.{module_name}"), name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
