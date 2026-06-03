"""Evidence bundle and evidence-book helpers."""

from __future__ import annotations

from importlib import import_module

_MODULE_EXPORTS = {
    ".book": (
        "EvidenceBookValidationIssue",
        "EvidenceBookValidationReport",
        "build_evidence_book_index",
        "build_evidence_claim_map",
        "build_evidence_parity_dashboard",
        "evidence_book_root",
        "render_evidence_catalog",
        "render_evidence_parity_dashboard",
        "validate_evidence_book",
        "write_evidence_book_index",
    ),
    ".closure": (
        "build_analytical_surface_coverage",
        "build_claim_reaudit",
        "build_closure_criteria",
        "build_completion_gates",
        "build_evidence_maturity_scorecard",
        "build_evidence_review_ritual",
        "render_analytical_surface_coverage",
        "render_claim_reaudit",
        "render_closure_criteria",
        "render_completion_gates",
        "render_evidence_maturity_scorecard",
        "render_evidence_review_ritual",
    ),
    ".coverage": (
        "build_evidence_coverage_gap_report",
        "render_evidence_coverage_gap_report",
    ),
    ".freshness": (
        "build_evidence_freshness_report",
        "render_evidence_freshness_report",
    ),
    ".integrity": (
        "build_evidence_integrity_report",
        "render_evidence_integrity_report",
    ),
    ".portability": (
        "EvidencePathIssue",
        "EvidencePathValue",
        "audit_payload_path_values",
        "classify_locator_kind",
        "collect_payload_path_values",
        "is_portable_locator",
        "render_portability_rules_markdown",
    ),
    ".scaffolding": (
        "EvidenceBundleTemplateSpec",
        "build_evidence_bundle_template",
        "write_evidence_bundle_template",
    ),
    ".study_registry": (
        "EvidenceStudyBuildReport",
        "EvidenceStudyRegistration",
        "EvidenceStudyRerunReport",
        "build_registered_study",
        "get_study_registration",
        "rerun_selected_evidence",
        "study_registrations",
    ),
    ".workbench": (
        "EvidenceBookRefreshReport",
        "EvidenceBookSelectionBuildReport",
        "EvidenceBookStudyBuildReport",
        "EvidenceBookStudyRerunReport",
        "build_evidence_book_selection",
        "build_evidence_book_study",
        "list_registered_evidence_studies",
        "refresh_evidence_book",
        "rerun_evidence_book_selection",
        "write_bundle_reviewer_summaries",
        "write_docs_evidence_overview",
    ),
}

__all__ = [export for exports in _MODULE_EXPORTS.values() for export in exports]

_ATTRIBUTE_MODULES = {
    export: module_name
    for module_name, exports in _MODULE_EXPORTS.items()
    for export in exports
}


def __getattr__(name: str) -> object:
    module_name = _ATTRIBUTE_MODULES.get(name)
    if module_name is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module = import_module(module_name, __name__)
    value = getattr(module, name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted([*globals(), *__all__])
