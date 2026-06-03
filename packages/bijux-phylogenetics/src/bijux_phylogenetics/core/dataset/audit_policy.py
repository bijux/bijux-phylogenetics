from __future__ import annotations

from bijux_phylogenetics.bayesian.beast.models import (
    FossilCalibrationValidationReport,
    TipDatingValidationReport,
)
from bijux_phylogenetics.phylo.alignment import AlignmentForensicReport

from .context import _DatasetContext, _infer_group_columns, _ordered_taxa
from .models import (
    DatasetAnalysisDecision,
    DatasetAuditFinding,
    DatasetCompletenessMatrix,
    DatasetExclusionRow,
    DatasetExclusionTable,
    DatasetFixRecommendation,
    DatasetGroupImbalanceWarning,
    DatasetMinimalFixPlan,
    DatasetMismatchReport,
    DatasetPruningStepSummary,
    DatasetReadinessLevel,
    DatasetReviewerChecklist,
    DatasetReviewerChecklistItem,
    DatasetRiskComponent,
    DatasetRiskScoreReport,
)


def _append_finding(
    findings: list[DatasetAuditFinding],
    *,
    severity: str,
    category: str,
    code: str,
    message: str,
    affected_analyses: list[str],
) -> None:
    finding = DatasetAuditFinding(
        severity=severity,
        category=category,
        code=code,
        message=message,
        affected_analyses=sorted(dict.fromkeys(affected_analyses)),
    )
    if finding not in findings:
        findings.append(finding)


def _categorize_findings(
    findings: list[DatasetAuditFinding], severity: str
) -> dict[str, list[str]]:
    categories: dict[str, list[str]] = {}
    for finding in findings:
        if finding.severity != severity:
            continue
        categories.setdefault(finding.category, []).append(finding.message)
    return dict(sorted(categories.items()))


def _pruning_steps(
    context: _DatasetContext,
) -> tuple[list[DatasetPruningStepSummary], set[str]]:
    current = set(context.tree_taxa)
    steps = [
        DatasetPruningStepSummary(
            step="tree",
            input_taxa=len(current),
            retained_taxa=len(current),
            excluded_taxa=0,
            reason="starting canonical taxon set from the tree tips",
        )
    ]

    for step_name, next_taxa, reason in (
        (
            "alignment",
            set(context.alignment_ids)
            if context.alignment_path is not None
            else set(context.tree_taxa),
            "retain taxa present in the alignment",
        ),
        (
            "metadata",
            set(_ordered_taxa(context.metadata_table)),
            "retain taxa present in the metadata table",
        ),
        (
            "traits",
            set(_ordered_taxa(context.traits_table)),
            "retain taxa present in the trait table",
        ),
        (
            "tip_dates",
            set(context.tip_date_taxa)
            if context.tip_dates_path is not None
            else current,
            "retain taxa with valid tip-date rows",
        ),
    ):
        if step_name == "alignment" and context.alignment_path is None:
            continue
        if step_name == "tip_dates" and context.tip_dates_path is None:
            continue
        input_count = len(current)
        current = current & set(next_taxa)
        steps.append(
            DatasetPruningStepSummary(
                step=step_name,
                input_taxa=input_count,
                retained_taxa=len(current),
                excluded_taxa=input_count - len(current),
                reason=reason,
            )
        )
    return steps, current


def _group_imbalance_warnings(
    context: _DatasetContext,
    retained_taxa: set[str],
) -> list[DatasetGroupImbalanceWarning]:
    warnings: list[DatasetGroupImbalanceWarning] = []
    for surface, table in (
        ("metadata", context.metadata_table),
        ("traits", context.traits_table),
    ):
        for column in _infer_group_columns(table):
            groups: dict[str, set[str]] = {}
            for row in table.rows:
                value = row[column].strip()
                if value:
                    groups.setdefault(value, set()).add(row[table.taxon_column])
            for group, taxa in groups.items():
                original_count = len(taxa)
                retained_count = len(taxa & retained_taxa)
                removed_count = original_count - retained_count
                removed_fraction = (
                    0.0 if original_count == 0 else removed_count / original_count
                )
                if original_count >= 1 and removed_fraction >= 0.5:
                    warnings.append(
                        DatasetGroupImbalanceWarning(
                            surface=surface,
                            group_column=column,
                            group=group,
                            original_count=original_count,
                            retained_count=retained_count,
                            removed_count=removed_count,
                            removed_fraction=removed_fraction,
                            message=f"{surface} group '{group}' in column '{column}' loses most taxa after pruning",
                        )
                    )
    return warnings


def _affected_analyses_for_missing_surface(surface: str) -> list[str]:
    return {
        "tree": [
            "inspection",
            "distance",
            "maximum_likelihood",
            "bayesian",
            "coding",
            "comparative",
            "time_tree",
            "publication",
        ],
        "alignment": [
            "distance",
            "maximum_likelihood",
            "bayesian",
            "coding",
            "time_tree",
            "publication",
        ],
        "metadata": ["comparative", "time_tree", "publication"],
        "traits": ["comparative", "publication"],
        "tip_dates": ["time_tree", "publication"],
    }.get(surface, ["publication"])


def _build_exclusion_table(
    completeness_matrix: DatasetCompletenessMatrix,
) -> DatasetExclusionTable:
    rows: list[DatasetExclusionRow] = []
    for row in completeness_matrix.rows:
        causes: list[str] = []
        if not row.in_tree:
            causes.append("absent_from_tree")
        if completeness_matrix.alignment_path is not None and not row.in_alignment:
            causes.append("absent_from_alignment")
        if not row.in_metadata:
            causes.append("absent_from_metadata")
        if not row.in_traits:
            causes.append("absent_from_traits")
        if completeness_matrix.tip_dates_path is not None and not row.in_tip_dates:
            causes.append("absent_from_tip_dates")
        if not causes:
            continue
        first_surface = causes[0].removeprefix("absent_from_")
        affected_analyses: list[str] = []
        for cause in causes:
            affected_analyses.extend(
                _affected_analyses_for_missing_surface(
                    cause.removeprefix("absent_from_")
                )
            )
        rows.append(
            DatasetExclusionRow(
                taxon=row.taxon,
                causes=causes,
                first_failed_surface=first_surface,
                affected_analyses=sorted(dict.fromkeys(affected_analyses)),
            )
        )
    return DatasetExclusionTable(rows=rows)


def _analysis_decisions(
    findings: list[DatasetAuditFinding],
    *,
    alignment_forensic: AlignmentForensicReport | None,
    tip_dates: TipDatingValidationReport | None,
    calibrations: FossilCalibrationValidationReport | None,
) -> list[DatasetAnalysisDecision]:
    analyses = [
        "inspection",
        "distance",
        "maximum_likelihood",
        "bayesian",
        "coding",
        "comparative",
        "time_tree",
        "publication",
    ]
    blockers_by_analysis = {analysis: [] for analysis in analyses}
    warnings_by_analysis = {analysis: [] for analysis in analyses}
    for finding in findings:
        target = (
            blockers_by_analysis
            if finding.severity == "blocker"
            else warnings_by_analysis
        )
        for analysis in finding.affected_analyses:
            if analysis in target:
                target[analysis].append(finding.message)

    decisions: list[DatasetAnalysisDecision] = []
    for analysis in analyses:
        reasons = sorted(
            dict.fromkeys(
                blockers_by_analysis[analysis] + warnings_by_analysis[analysis]
            )
        )
        if blockers_by_analysis[analysis]:
            decision = "blocked"
        elif warnings_by_analysis[analysis]:
            decision = "risky"
        else:
            decision = "allowed"
        if (
            analysis in {"distance", "maximum_likelihood", "bayesian", "coding"}
            and alignment_forensic is None
        ):
            decision = "blocked"
            reasons = ["analysis requires an alignment input"]
        if analysis == "time_tree" and tip_dates is None and calibrations is None:
            decision = "blocked"
            reasons = ["analysis requires tip dates or fossil calibrations"]
        decisions.append(
            DatasetAnalysisDecision(
                analysis=analysis, decision=decision, reasons=reasons
            )
        )
    return decisions


def _readiness_levels(
    decisions: list[DatasetAnalysisDecision],
) -> list[DatasetReadinessLevel]:
    decision_by_analysis = {row.analysis: row for row in decisions}

    def summarize(level: str, analyses: list[str]) -> DatasetReadinessLevel:
        members = [decision_by_analysis[analysis] for analysis in analyses]
        reasons = sorted(
            dict.fromkeys(reason for member in members for reason in member.reasons)
        )
        if any(member.decision == "blocked" for member in members):
            decision = "blocked"
        elif any(member.decision == "risky" for member in members):
            decision = "risky"
        else:
            decision = "ready"
        return DatasetReadinessLevel(level=level, decision=decision, reasons=reasons)

    return [
        summarize("inspection_ready", ["inspection"]),
        summarize("inference_ready", ["distance", "maximum_likelihood", "bayesian"]),
        summarize("comparative_ready", ["comparative"]),
        summarize("time_tree_ready", ["time_tree"]),
        summarize("publication_ready", ["publication"]),
    ]


def _dataset_risk_score(
    findings: list[DatasetAuditFinding],
    *,
    tip_dates: TipDatingValidationReport | None,
    calibrations: FossilCalibrationValidationReport | None,
) -> DatasetRiskScoreReport:
    components: list[DatasetRiskComponent] = []
    component_order = ["tree", "alignment", "taxon", "trait", "metadata", "calibration"]
    category_aliases = {
        "tree": {"tree", "dataset"},
        "alignment": {"alignment"},
        "taxon": {"taxon", "tip_dates"},
        "trait": {"traits"},
        "metadata": {"metadata"},
        "calibration": {"calibration"},
    }
    for component in component_order:
        component_findings = [
            finding
            for finding in findings
            if finding.category in category_aliases[component]
        ]
        score = 0.0
        reasons = [finding.message for finding in component_findings]
        score += 25.0 * sum(
            1 for finding in component_findings if finding.severity == "blocker"
        )
        score += 10.0 * sum(
            1 for finding in component_findings if finding.severity == "warning"
        )
        if component == "calibration" and tip_dates is None and calibrations is None:
            reasons.append("time-tree surfaces were not supplied")
        components.append(
            DatasetRiskComponent(
                component=component,
                score=min(100.0, score),
                reasons=sorted(dict.fromkeys(reasons)),
            )
        )
    total_score = round(
        sum(component.score for component in components) / len(components), 15
    )
    if total_score >= 50.0:
        risk_level = "high"
    elif total_score >= 20.0:
        risk_level = "moderate"
    else:
        risk_level = "low"
    return DatasetRiskScoreReport(
        total_score=total_score, risk_level=risk_level, components=components
    )


def _minimal_fix_plan(
    findings: list[DatasetAuditFinding],
    decisions: list[DatasetAnalysisDecision],
) -> DatasetMinimalFixPlan:
    recommendations: list[DatasetFixRecommendation] = []

    def add(
        priority: int, summary: str, affected_surfaces: list[str], unlocks: list[str]
    ) -> None:
        row = DatasetFixRecommendation(
            priority=priority,
            summary=summary,
            affected_surfaces=affected_surfaces,
            unlocks_analyses=unlocks,
        )
        if row not in recommendations:
            recommendations.append(row)

    messages = {finding.message: finding for finding in findings}
    if "tree requires complete branch lengths" in messages:
        add(
            1,
            "supply complete branch lengths for every non-root edge",
            ["tree"],
            ["comparative", "publication"],
        )
    if "metadata table is missing one or more tree taxa" in messages:
        add(
            1,
            "fill metadata rows for every tree tip or prune unmatched tips explicitly",
            ["metadata"],
            ["comparative", "publication"],
        )
    if "trait table is missing one or more tree taxa" in messages:
        add(
            1,
            "complete the trait table for every retained tree tip",
            ["traits"],
            ["comparative", "publication"],
        )
    if "alignment is missing one or more tree taxa" in messages:
        add(
            1,
            "reconcile tree and alignment taxon sets before inference",
            ["alignment", "tree"],
            [
                "distance",
                "maximum_likelihood",
                "bayesian",
                "coding",
                "time_tree",
                "publication",
            ],
        )
    if "tip-date metadata contains invalid or missing tree-tip dates" in messages:
        add(
            1,
            "correct invalid or missing tip-date rows",
            ["tip_dates"],
            ["time_tree", "publication"],
        )
    if (
        "calibration table contains invalid fossil calibration targets or ages"
        in messages
    ):
        add(
            1,
            "repair invalid fossil calibrations or remove unsupported targets",
            ["calibration"],
            ["time_tree", "publication"],
        )
    if any(
        decision.analysis in {"distance", "maximum_likelihood", "bayesian", "coding"}
        for decision in decisions
        if decision.decision == "blocked"
        and "alignment input" in " ".join(decision.reasons)
    ):
        add(
            2,
            "add an alignment input for sequence-based inference workflows",
            ["alignment"],
            ["distance", "maximum_likelihood", "bayesian", "coding"],
        )
    if any(
        decision.analysis == "time_tree"
        for decision in decisions
        if decision.decision == "blocked"
        and "tip dates or fossil calibrations" in " ".join(decision.reasons)
    ):
        add(
            2,
            "supply tip dates or fossil calibrations for time-tree analysis",
            ["tip_dates", "calibration"],
            ["time_tree"],
        )

    recommendations.sort(key=lambda row: (row.priority, row.summary))
    return DatasetMinimalFixPlan(recommendations=recommendations)


def _reviewer_checklist(
    *,
    readiness_decision: str,
    blockers: list[str],
    warnings: list[str],
    mismatch_report: DatasetMismatchReport,
    risk_score: DatasetRiskScoreReport,
    decisions: list[DatasetAnalysisDecision],
) -> DatasetReviewerChecklist:
    def status_for_messages(messages: list[str]) -> str:
        if messages:
            return "blocked"
        return "pass"

    risky_analyses = [
        decision.analysis for decision in decisions if decision.decision == "risky"
    ]
    blocked_analyses = [
        decision.analysis for decision in decisions if decision.decision == "blocked"
    ]
    items = [
        DatasetReviewerChecklistItem(
            section="overall_readiness",
            status="blocked"
            if readiness_decision == "blocked"
            else ("risk" if warnings else "pass"),
            summary=f"dataset readiness is {readiness_decision}",
            evidence=blockers or warnings,
        ),
        DatasetReviewerChecklistItem(
            section="taxon_mismatch",
            status="risk" if mismatch_report.rows else "pass",
            summary="cross-surface taxon mismatches were audited",
            evidence=[row.message for row in mismatch_report.rows[:10]],
        ),
        DatasetReviewerChecklistItem(
            section="dataset_risk",
            status="risk" if risk_score.risk_level != "low" else "pass",
            summary=f"transparent dataset risk level is {risk_score.risk_level}",
            evidence=[
                f"{row.component}: {row.score:g}" for row in risk_score.components
            ],
        ),
        DatasetReviewerChecklistItem(
            section="analysis_eligibility",
            status=status_for_messages(blocked_analyses),
            summary="downstream analysis eligibility was classified explicitly",
            evidence=blocked_analyses + risky_analyses,
        ),
    ]
    return DatasetReviewerChecklist(items=items)
