from __future__ import annotations

"""Compatibility facade for comparative validation and assessment surfaces."""

from bijux_phylogenetics.comparative.assessment import (
    ComparativeMethodMaturityReport,
    ComparativeResidualDiagnosticSurface,
    ComparativeSensitivityReport,
    ComparativeSensitivitySummary,
    LeaveOneTaxonOutRow,
    assess_comparative_method_maturity,
    run_comparative_sensitivity_analysis,
)
from bijux_phylogenetics.comparative.validation import (
    ComparativeOUIdentifiabilityAudit,
    ComparativeOUIdentifiabilityCase,
    ComparativeParameterIntervalAuditRow,
    ComparativeParameterUncertaintyAudit,
    ComparativeReferenceObservation,
    ComparativeReferenceValidationReport,
    audit_comparative_parameter_uncertainty,
    audit_ou_identifiability_reference_examples,
    validate_comparative_reference_examples,
)

__all__ = [
    "ComparativeMethodMaturityReport",
    "ComparativeOUIdentifiabilityAudit",
    "ComparativeOUIdentifiabilityCase",
    "ComparativeParameterIntervalAuditRow",
    "ComparativeParameterUncertaintyAudit",
    "ComparativeReferenceObservation",
    "ComparativeReferenceValidationReport",
    "ComparativeResidualDiagnosticSurface",
    "ComparativeSensitivityReport",
    "ComparativeSensitivitySummary",
    "LeaveOneTaxonOutRow",
    "assess_comparative_method_maturity",
    "audit_comparative_parameter_uncertainty",
    "audit_ou_identifiability_reference_examples",
    "run_comparative_sensitivity_analysis",
    "validate_comparative_reference_examples",
]
