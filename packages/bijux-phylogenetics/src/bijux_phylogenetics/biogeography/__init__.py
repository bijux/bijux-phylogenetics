from __future__ import annotations

from importlib import import_module

_PUBLIC_SURFACES = (
    (
        "state_models",
        (
            "ConstrainedGeographicFitRow",
            "ConstrainedGeographicReport",
            "ConstrainedGeographicSummary",
            "ConstrainedGeographicTransitionRow",
            "GeographicExcludedTaxonRow",
            "GeographicRegionProbabilityRow",
            "GeographicSamplingBiasNodeRow",
            "GeographicSamplingBiasReport",
            "GeographicSamplingBiasSummary",
            "GeographicSamplingBiasTransitionRow",
            "GeographicSamplingCountRow",
            "GeographicStateModelReport",
            "GeographicStateSummary",
            "GeographicTransitionEventRow",
            "GeographicTransitionRateRow",
            "TimeBinDefinition",
            "TimeStratifiedBranchRow",
            "TimeStratifiedTransitionMatrixRow",
            "TimeStratifiedTransitionReport",
            "TimeStratifiedTransitionSummary",
            "UnsupportedGeographicTransitionClaimRow",
            "summarize_constrained_geographic_model",
            "summarize_constrained_geographic_report",
            "summarize_geographic_sampling_bias",
            "summarize_geographic_state_model",
            "summarize_time_stratified_geographic_transitions",
            "write_constrained_geographic_exclusion_table",
            "write_constrained_geographic_fit_table",
            "write_constrained_geographic_summary_table",
            "write_constrained_geographic_transition_table",
            "write_geographic_exclusion_table",
            "write_geographic_region_probability_table",
            "write_geographic_sampling_bias_exclusion_table",
            "write_geographic_sampling_bias_node_table",
            "write_geographic_sampling_bias_summary_table",
            "write_geographic_sampling_bias_transition_table",
            "write_geographic_sampling_count_table",
            "write_geographic_state_summary_table",
            "write_geographic_transition_event_table",
            "write_geographic_transition_rate_table",
            "write_time_stratified_branch_table",
            "write_time_stratified_exclusion_table",
            "write_time_stratified_transition_matrix_table",
            "write_time_stratified_transition_summary_table",
            "write_unsupported_geographic_transition_claim_table",
        ),
    ),
    (
        "migration",
        (
            "DatedBiogeographyEventRow",
            "DatedBiogeographyNodeRow",
            "DatedBiogeographyReport",
            "DatedBiogeographySummary",
            "DatedBiogeographyTimeBinRow",
            "GeographicMigrationEventReport",
            "GeographicMigrationEventRow",
            "GeographicMigrationEventSummary",
            "GeographicMigrationTreeRow",
            "GeographicMigrationTreeSetEventRow",
            "GeographicMigrationTreeSetEventSummaryRow",
            "GeographicMigrationTreeSetReport",
            "GeographicMigrationTreeSetSummary",
            "summarize_biogeographic_transition_chronology",
            "summarize_geographic_migration_event_tree_set",
            "summarize_geographic_migration_events",
            "write_dated_biogeography_event_table",
            "write_dated_biogeography_exclusion_table",
            "write_dated_biogeography_node_table",
            "write_dated_biogeography_summary_table",
            "write_dated_biogeography_time_bin_table",
            "write_geographic_migration_event_summary_table",
            "write_geographic_migration_event_table",
            "write_geographic_migration_exclusion_table",
            "write_geographic_migration_tree_set_event_summary_table",
            "write_geographic_migration_tree_set_event_table",
            "write_geographic_migration_tree_set_exclusion_table",
            "write_geographic_migration_tree_set_summary_table",
            "write_geographic_migration_tree_set_tree_table",
        ),
    ),
    (
        "presentation",
        (
            "BiogeographyRegionCountRow",
            "BiogeographyReportExclusionRow",
            "BiogeographyReportPackageResult",
            "build_biogeography_report_package",
            "summarize_biogeography_region_counts",
            "write_biogeography_region_count_table",
            "write_biogeography_report_exclusion_table",
        ),
    ),
)

__all__ = [name for _, names in _PUBLIC_SURFACES for name in names]

_NAME_TO_MODULE = {
    name: module_name for module_name, names in _PUBLIC_SURFACES for name in names
}


def __getattr__(name: str):
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
