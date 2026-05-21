from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.bayesian.beast.models import (
    FossilCalibrationValidationReport,
    TipDatingValidationReport,
)
from bijux_phylogenetics.bayesian.beast.validation import (
    validate_fossil_calibration_table,
    validate_tip_dating_metadata,
)
from bijux_phylogenetics.datasets.study_inputs import (
    load_taxon_table,
)
from bijux_phylogenetics.io.fasta.quality import build_alignment_forensic_report
from bijux_phylogenetics.io.fasta.records import link_alignment_to_tree
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment import AlignmentForensicReport

from .audit_policy import (
    _analysis_decisions,
    _append_finding,
    _build_exclusion_table,
    _categorize_findings,
    _dataset_risk_score,
    _group_imbalance_warnings,
    _minimal_fix_plan,
    _pruning_steps,
    _readiness_levels,
    _reviewer_checklist,
)
from .context import _load_dataset_context
from .crosswalk import (
    build_dataset_completeness_matrix,
    build_dataset_crosswalk,
    build_dataset_mismatch_report,
)
from .models import DatasetAuditFinding, DatasetAuditReport
from .ordering import audit_dataset_taxon_ordering
from .readiness import summarize_dataset_readiness


def audit_dataset_inputs(
    tree_path: Path,
    metadata_path: Path,
    traits_path: Path,
    *,
    alignment_path: Path | None = None,
    tip_dates_path: Path | None = None,
    calibration_path: Path | None = None,
    trait_missingness_threshold: float = 0.25,
) -> DatasetAuditReport:
    """Produce one integrated readiness decision across major phylogenetic input surfaces."""
    readiness = summarize_dataset_readiness(
        tree_path,
        metadata_path,
        traits_path,
        trait_missingness_threshold=trait_missingness_threshold,
    )
    tree = load_tree(tree_path)
    tree_taxa = set(tree.tip_names)
    metadata_table = load_taxon_table(metadata_path)
    traits_table = load_taxon_table(traits_path)
    metadata_taxa = set(metadata_table.taxa)
    trait_taxa = set(traits_table.taxa)

    findings: list[DatasetAuditFinding] = []
    for message in readiness.blockers:
        category = (
            "traits"
            if "trait" in message
            else ("metadata" if "metadata" in message else "tree")
        )
        _append_finding(
            findings,
            severity="blocker",
            category=category,
            code=message.replace(" ", "-"),
            message=message,
            affected_analyses=["comparative", "publication"],
        )
    for message in readiness.warnings:
        category = (
            "traits"
            if "trait" in message
            else ("metadata" if "metadata" in message else "dataset")
        )
        _append_finding(
            findings,
            severity="warning",
            category=category,
            code=message.replace(" ", "-"),
            message=message,
            affected_analyses=["comparative", "publication"],
        )

    alignment_forensic: AlignmentForensicReport | None = None
    alignment_taxa = set(tree_taxa)
    if alignment_path is not None:
        alignment_linkage = link_alignment_to_tree(tree_path, alignment_path)
        alignment_forensic = build_alignment_forensic_report(alignment_path)
        alignment_taxa = set(alignment_linkage.usable_taxa)
        if alignment_linkage.missing_from_alignment:
            _append_finding(
                findings,
                severity="blocker",
                category="alignment",
                code="missing-tree-taxa",
                message="alignment is missing one or more tree taxa",
                affected_analyses=[
                    "distance",
                    "maximum_likelihood",
                    "bayesian",
                    "coding",
                    "time_tree",
                    "publication",
                ],
            )
        if alignment_linkage.extra_alignment_ids:
            _append_finding(
                findings,
                severity="warning",
                category="alignment",
                code="extra-alignment-taxa",
                message="alignment contains taxa absent from the tree",
                affected_analyses=[
                    "distance",
                    "maximum_likelihood",
                    "bayesian",
                    "coding",
                    "publication",
                ],
            )
        if (
            not alignment_forensic.safe_for_distance_analysis
            and not alignment_forensic.safe_for_maximum_likelihood
        ):
            _append_finding(
                findings,
                severity="blocker",
                category="alignment",
                code="inference-unsafe-alignment",
                message="alignment is not currently safe for core inference workflows",
                affected_analyses=[
                    "distance",
                    "maximum_likelihood",
                    "bayesian",
                    "publication",
                ],
            )
        for warning in alignment_forensic.warnings:
            affected = ["publication"]
            if "coding" in warning:
                affected.append("coding")
            if "alignment" in warning or "sites" in warning or "sequence" in warning:
                affected.extend(["distance", "maximum_likelihood", "bayesian"])
            _append_finding(
                findings,
                severity="warning",
                category="alignment",
                code=warning.replace(" ", "-"),
                message=warning,
                affected_analyses=affected,
            )

    tip_dates: TipDatingValidationReport | None = None
    if tip_dates_path is not None:
        tip_dates = validate_tip_dating_metadata(
            tree_path,
            tip_dates_path,
            alignment_path=alignment_path,
        )
        if tip_dates.invalid_tip_count > 0 or tip_dates.missing_tree_taxa:
            _append_finding(
                findings,
                severity="blocker",
                category="tip_dates",
                code="invalid-tip-dates",
                message="tip-date metadata contains invalid or missing tree-tip dates",
                affected_analyses=["time_tree", "publication"],
            )
        if tip_dates.extra_tip_taxa or tip_dates.extra_alignment_taxa:
            _append_finding(
                findings,
                severity="warning",
                category="tip_dates",
                code="extra-tip-date-taxa",
                message="tip-date metadata contains taxa absent from the tree or alignment",
                affected_analyses=["time_tree", "publication"],
            )

    calibrations: FossilCalibrationValidationReport | None = None
    if calibration_path is not None:
        calibrations = validate_fossil_calibration_table(tree_path, calibration_path)
        if calibrations.invalid_calibration_count > 0:
            _append_finding(
                findings,
                severity="blocker",
                category="calibration",
                code="invalid-calibrations",
                message="calibration table contains invalid fossil calibration targets or ages",
                affected_analyses=["time_tree", "publication"],
            )

    context = _load_dataset_context(
        tree_path,
        metadata_path,
        traits_path,
        alignment_path=alignment_path,
        tip_dates_path=tip_dates_path,
        calibration_path=calibration_path,
        alignment_forensic=alignment_forensic,
        tip_dates_report=tip_dates,
        calibration_report=calibrations,
    )
    pruning_steps, analysis_taxa_set = _pruning_steps(context)
    analysis_taxa = sorted(
        tree_taxa & metadata_taxa & trait_taxa & alignment_taxa & analysis_taxa_set
    )
    if len(analysis_taxa) < 2:
        _append_finding(
            findings,
            severity="blocker",
            category="dataset",
            code="too-few-intersection-taxa",
            message="fewer than two taxa remain after intersecting all requested dataset surfaces",
            affected_analyses=[
                "comparative",
                "distance",
                "maximum_likelihood",
                "bayesian",
                "coding",
                "time_tree",
                "publication",
            ],
        )

    ordering_audit = audit_dataset_taxon_ordering(
        tree_path,
        metadata_path,
        traits_path,
        alignment_path=alignment_path,
        tip_dates_path=tip_dates_path,
    )
    if not ordering_audit.consistent:
        _append_finding(
            findings,
            severity="warning",
            category="metadata",
            code="ordering-drift",
            message="one or more dataset surfaces silently reorder shared taxa relative to the tree",
            affected_analyses=["comparative", "publication"],
        )

    group_imbalance_warnings = _group_imbalance_warnings(context, set(analysis_taxa))
    for warning in group_imbalance_warnings:
        _append_finding(
            findings,
            severity="warning",
            category=warning.surface,
            code="group-imbalance",
            message=warning.message,
            affected_analyses=["comparative", "publication"],
        )

    crosswalk = build_dataset_crosswalk(
        tree_path,
        metadata_path,
        traits_path,
        alignment_path=alignment_path,
        tip_dates_path=tip_dates_path,
        calibration_path=calibration_path,
    )
    completeness_matrix = build_dataset_completeness_matrix(
        tree_path,
        metadata_path,
        traits_path,
        alignment_path=alignment_path,
        tip_dates_path=tip_dates_path,
        calibration_path=calibration_path,
    )
    exclusion_table = _build_exclusion_table(completeness_matrix)
    mismatch_report = build_dataset_mismatch_report(
        tree_path,
        metadata_path,
        traits_path,
        alignment_path=alignment_path,
        tip_dates_path=tip_dates_path,
        calibration_path=calibration_path,
    )
    analysis_decisions = _analysis_decisions(
        findings,
        alignment_forensic=alignment_forensic,
        tip_dates=tip_dates,
        calibrations=calibrations,
    )
    readiness_levels = _readiness_levels(analysis_decisions)

    blockers = sorted(
        dict.fromkeys(
            finding.message for finding in findings if finding.severity == "blocker"
        )
    )
    warnings = sorted(
        dict.fromkeys(
            finding.message for finding in findings if finding.severity == "warning"
        )
    )
    allowed_analyses = sorted(
        decision.analysis
        for decision in analysis_decisions
        if decision.decision != "blocked"
    )
    blocked_analyses = sorted(
        decision.analysis
        for decision in analysis_decisions
        if decision.decision == "blocked"
    )
    readiness_decision = (
        "blocked" if blockers else ("ready_with_warnings" if warnings else "ready")
    )
    risk_score = _dataset_risk_score(
        findings, tip_dates=tip_dates, calibrations=calibrations
    )
    minimal_fix_plan = _minimal_fix_plan(findings, analysis_decisions)
    reviewer_checklist = _reviewer_checklist(
        readiness_decision=readiness_decision,
        blockers=blockers,
        warnings=warnings,
        mismatch_report=mismatch_report,
        risk_score=risk_score,
        decisions=analysis_decisions,
    )

    return DatasetAuditReport(
        tree_path=tree_path,
        metadata_path=metadata_path,
        traits_path=traits_path,
        alignment_path=alignment_path,
        tip_dates_path=tip_dates_path,
        calibration_path=calibration_path,
        tree_taxa=len(tree_taxa),
        analysis_taxa=analysis_taxa,
        readiness_decision=readiness_decision,
        allowed_analyses=allowed_analyses,
        blocked_analyses=blocked_analyses,
        blockers=blockers,
        warnings=warnings,
        blocker_categories=_categorize_findings(findings, "blocker"),
        warning_categories=_categorize_findings(findings, "warning"),
        findings=findings,
        analysis_decisions=analysis_decisions,
        readiness_levels=readiness_levels,
        crosswalk=crosswalk,
        completeness_matrix=completeness_matrix,
        exclusion_table=exclusion_table,
        ordering_audit=ordering_audit,
        pruning_steps=pruning_steps,
        group_imbalance_warnings=group_imbalance_warnings,
        dataset_readiness=readiness,
        alignment_forensic=alignment_forensic,
        tip_dates=tip_dates,
        calibrations=calibrations,
        mismatch_report=mismatch_report,
        risk_score=risk_score,
        minimal_fix_plan=minimal_fix_plan,
        reviewer_checklist=reviewer_checklist,
    )
