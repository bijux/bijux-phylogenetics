"""Governed `geiger` parity surfaces."""

from .boundary_warning_registry import (
    GeigerBoundaryWarningRow,
    write_geiger_boundary_warning_table,
)
from .generated_report import (
    GeigerBenchmarkSummaryRow,
    GeigerBoundaryWarningSummaryRow,
    GeigerExcludedModelRow,
    GeigerGoalCoverageRow,
    GeigerOptimizerMismatchCategoryRow,
    GeigerSimulationRecoveryRow,
    GeigerToleranceRuleRow,
    GeneratedGeigerParityReport,
    build_generated_geiger_parity_report,
    write_generated_geiger_parity_report_json,
    write_generated_geiger_parity_report_markdown,
)
from .likelihood_policy import (
    GeigerLikelihoodPolicyRow,
    write_geiger_likelihood_policy_table,
)
from .model_confidence import (
    GeigerModelConfidenceRow,
    write_geiger_model_confidence_table,
)
from .optimizer_triage import (
    GeigerOptimizerTriageRow,
    write_geiger_optimizer_triage_table,
)
from .parameterization_registry import (
    GeigerParameterizationRegistryRow,
    write_geiger_parameterization_registry_table,
)
from .registry import GeigerParityCase, list_geiger_parity_cases
from .runner import (
    GeigerParityObservation,
    GeigerParityReport,
    GeigerParitySummaryRow,
    run_geiger_parity_cases,
    write_geiger_parity_observation_table,
    write_geiger_parity_summary_table,
)

__all__ = [
    "GeigerParityCase",
    "GeigerParityObservation",
    "GeigerParityReport",
    "GeigerParitySummaryRow",
    "GeigerOptimizerTriageRow",
    "GeigerBoundaryWarningRow",
    "GeneratedGeigerParityReport",
    "GeigerGoalCoverageRow",
    "GeigerExcludedModelRow",
    "GeigerToleranceRuleRow",
    "GeigerOptimizerMismatchCategoryRow",
    "GeigerBoundaryWarningSummaryRow",
    "GeigerSimulationRecoveryRow",
    "GeigerBenchmarkSummaryRow",
    "GeigerLikelihoodPolicyRow",
    "GeigerModelConfidenceRow",
    "GeigerParameterizationRegistryRow",
    "list_geiger_parity_cases",
    "build_generated_geiger_parity_report",
    "run_geiger_parity_cases",
    "write_geiger_parity_observation_table",
    "write_geiger_parity_summary_table",
    "write_geiger_optimizer_triage_table",
    "write_geiger_boundary_warning_table",
    "write_geiger_likelihood_policy_table",
    "write_geiger_model_confidence_table",
    "write_geiger_parameterization_registry_table",
    "write_generated_geiger_parity_report_json",
    "write_generated_geiger_parity_report_markdown",
]
