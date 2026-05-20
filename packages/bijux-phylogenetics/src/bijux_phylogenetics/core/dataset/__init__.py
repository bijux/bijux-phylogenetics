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
from bijux_phylogenetics.core.alignment import AlignmentForensicReport
from bijux_phylogenetics.core.metadata import inspect_metadata_table, load_taxon_table
from bijux_phylogenetics.core.traits import (
    check_tree_and_trait_taxon_names,
    detect_unusable_trait_columns,
    link_tree_to_traits,
)
from bijux_phylogenetics.diagnostics.validation import validate_tree_path
from bijux_phylogenetics.io.fasta.quality import build_alignment_forensic_report
from bijux_phylogenetics.io.fasta.records import link_alignment_to_tree
from bijux_phylogenetics.io.trees import load_tree
from .context import (
    _DatasetContext,
    _collect_external_ids,
    _infer_group_columns,
    _load_dataset_context,
    _ordered_taxa,
)
from .models import (
    DatasetAnalysisDecision,
    DatasetAuditFinding,
    DatasetAuditReport,
    DatasetCompletenessMatrix,
    DatasetCompletenessRow,
    DatasetCrosswalkReport,
    DatasetCrosswalkRow,
    DatasetExclusionRow,
    DatasetExclusionTable,
    DatasetFixRecommendation,
    DatasetGroupImbalanceWarning,
    DatasetMinimalFixPlan,
    DatasetMismatchReport,
    DatasetMismatchRow,
    DatasetOrderingAudit,
    DatasetOrderingConflict,
    DatasetPruningStepSummary,
    DatasetReadinessLevel,
    DatasetReadinessSummary,
    DatasetReviewerChecklist,
    DatasetReviewerChecklistItem,
    DatasetRiskComponent,
    DatasetRiskScoreReport,
)


def build_dataset_crosswalk(
    tree_path: Path,
    metadata_path: Path,
    traits_path: Path,
    *,
    alignment_path: Path | None = None,
    tip_dates_path: Path | None = None,
    calibration_path: Path | None = None,
) -> DatasetCrosswalkReport:
    """Generate an explicit taxon crosswalk across the main dataset surfaces."""
    context = _load_dataset_context(
        tree_path,
        metadata_path,
        traits_path,
        alignment_path=alignment_path,
        tip_dates_path=tip_dates_path,
        calibration_path=calibration_path,
    )
    metadata_taxa = set(_ordered_taxa(context.metadata_table))
    trait_taxa = set(_ordered_taxa(context.traits_table))
    alignment_taxa = set(context.alignment_ids)
    tip_date_taxa = set(context.tip_date_taxa)
    tree_taxa = set(context.tree_taxa)
    union_taxa = sorted(
        tree_taxa | metadata_taxa | trait_taxa | alignment_taxa | tip_date_taxa
    )
    rows = [
        DatasetCrosswalkRow(
            taxon=taxon,
            tree_tip=taxon if taxon in tree_taxa else None,
            alignment_id=taxon if taxon in alignment_taxa else None,
            metadata_id=taxon if taxon in metadata_taxa else None,
            trait_id=taxon if taxon in trait_taxa else None,
            tip_date_id=taxon if taxon in tip_date_taxa else None,
            geography_source="geography" if taxon in context.geography_taxa else None,
            calibration_targets=context.calibration_taxa_to_targets.get(taxon, []),
            external_taxonomy_ids=_collect_external_ids(
                taxon,
                context.metadata_table,
                context.traits_table,
            ),
        )
        for taxon in union_taxa
    ]
    return DatasetCrosswalkReport(
        tree_path=tree_path,
        metadata_path=metadata_path,
        traits_path=traits_path,
        alignment_path=alignment_path,
        tip_dates_path=tip_dates_path,
        calibration_path=calibration_path,
        rows=rows,
    )


def build_dataset_completeness_matrix(
    tree_path: Path,
    metadata_path: Path,
    traits_path: Path,
    *,
    alignment_path: Path | None = None,
    tip_dates_path: Path | None = None,
    calibration_path: Path | None = None,
) -> DatasetCompletenessMatrix:
    """Build a taxon-by-surface completeness matrix for one dataset."""
    crosswalk = build_dataset_crosswalk(
        tree_path,
        metadata_path,
        traits_path,
        alignment_path=alignment_path,
        tip_dates_path=tip_dates_path,
        calibration_path=calibration_path,
    )
    rows = [
        DatasetCompletenessRow(
            taxon=row.taxon,
            in_tree=row.tree_tip is not None,
            in_alignment=row.alignment_id is not None,
            in_metadata=row.metadata_id is not None,
            in_traits=row.trait_id is not None,
            in_tip_dates=row.tip_date_id is not None,
            in_geography=row.geography_source is not None,
            in_calibrations=bool(row.calibration_targets),
        )
        for row in crosswalk.rows
    ]
    surface_counts = {
        "tree": sum(1 for row in rows if row.in_tree),
        "alignment": sum(1 for row in rows if row.in_alignment),
        "metadata": sum(1 for row in rows if row.in_metadata),
        "traits": sum(1 for row in rows if row.in_traits),
        "tip_dates": sum(1 for row in rows if row.in_tip_dates),
        "geography": sum(1 for row in rows if row.in_geography),
        "calibrations": sum(1 for row in rows if row.in_calibrations),
    }
    context = _load_dataset_context(
        tree_path,
        metadata_path,
        traits_path,
        alignment_path=alignment_path,
        tip_dates_path=tip_dates_path,
        calibration_path=calibration_path,
    )
    return DatasetCompletenessMatrix(
        tree_path=tree_path,
        metadata_path=metadata_path,
        traits_path=traits_path,
        alignment_path=alignment_path,
        tip_dates_path=tip_dates_path,
        calibration_path=calibration_path,
        geography_columns=context.geography_columns,
        rows=rows,
        surface_counts=surface_counts,
    )


def build_dataset_mismatch_report(
    tree_path: Path,
    metadata_path: Path,
    traits_path: Path,
    *,
    alignment_path: Path | None = None,
    tip_dates_path: Path | None = None,
    calibration_path: Path | None = None,
) -> DatasetMismatchReport:
    """Show which taxa are missing from which requested dataset surfaces."""
    matrix = build_dataset_completeness_matrix(
        tree_path,
        metadata_path,
        traits_path,
        alignment_path=alignment_path,
        tip_dates_path=tip_dates_path,
        calibration_path=calibration_path,
    )
    requested_surfaces = ["tree", "metadata", "traits"]
    if alignment_path is not None:
        requested_surfaces.append("alignment")
    if tip_dates_path is not None:
        requested_surfaces.append("tip_dates")
    if calibration_path is not None:
        requested_surfaces.append("calibrations")

    rows: list[DatasetMismatchRow] = []
    mismatch_counts = dict.fromkeys(requested_surfaces, 0)
    for row in matrix.rows:
        present_surfaces: list[str] = []
        missing_surfaces: list[str] = []
        for surface in requested_surfaces:
            present = {
                "tree": row.in_tree,
                "metadata": row.in_metadata,
                "traits": row.in_traits,
                "alignment": row.in_alignment,
                "tip_dates": row.in_tip_dates,
                "calibrations": row.in_calibrations,
            }[surface]
            if present:
                present_surfaces.append(surface)
            else:
                missing_surfaces.append(surface)
                mismatch_counts[surface] += 1
        if missing_surfaces and present_surfaces:
            rows.append(
                DatasetMismatchRow(
                    taxon=row.taxon,
                    present_surfaces=present_surfaces,
                    missing_surfaces=missing_surfaces,
                    message=f"taxon is missing from {', '.join(missing_surfaces)} while present in {', '.join(present_surfaces)}",
                )
            )
    return DatasetMismatchReport(
        requested_surfaces=requested_surfaces,
        rows=rows,
        mismatch_counts=mismatch_counts,
    )


def _ordering_conflicts_for_surface(
    *,
    surface: str,
    canonical_order: list[str],
    observed_order: list[str],
) -> list[DatasetOrderingConflict]:
    canonical_shared = [
        taxon for taxon in canonical_order if taxon in set(observed_order)
    ]
    observed_shared = [
        taxon for taxon in observed_order if taxon in set(canonical_order)
    ]
    if canonical_shared == observed_shared:
        return []

    expected_index = {
        taxon: index for index, taxon in enumerate(canonical_shared, start=1)
    }
    return [
        DatasetOrderingConflict(
            surface=surface,
            taxon=taxon,
            expected_index=expected_index[taxon],
            observed_index=index,
        )
        for index, taxon in enumerate(observed_shared, start=1)
        if expected_index[taxon] != index
    ]


def audit_dataset_taxon_ordering(
    tree_path: Path,
    metadata_path: Path,
    traits_path: Path,
    *,
    alignment_path: Path | None = None,
    tip_dates_path: Path | None = None,
) -> DatasetOrderingAudit:
    """Detect silent taxon-order drift across dataset surfaces."""
    context = _load_dataset_context(
        tree_path,
        metadata_path,
        traits_path,
        alignment_path=alignment_path,
        tip_dates_path=tip_dates_path,
    )
    conflicts: list[DatasetOrderingConflict] = []
    drifted_surfaces: list[str] = []
    surface_orders = {
        "metadata": _ordered_taxa(context.metadata_table),
        "traits": _ordered_taxa(context.traits_table),
    }
    if alignment_path is not None:
        surface_orders["alignment"] = context.alignment_ids
    if tip_dates_path is not None:
        surface_orders["tip_dates"] = context.tip_date_taxa

    for surface, observed_order in surface_orders.items():
        surface_conflicts = _ordering_conflicts_for_surface(
            surface=surface,
            canonical_order=context.tree_taxa,
            observed_order=observed_order,
        )
        if surface_conflicts:
            drifted_surfaces.append(surface)
            conflicts.extend(surface_conflicts)

    return DatasetOrderingAudit(
        canonical_surface="tree",
        consistent=not conflicts,
        drifted_surfaces=sorted(drifted_surfaces),
        conflicts=conflicts,
    )


def summarize_dataset_readiness(
    tree_path: Path,
    metadata_path: Path,
    traits_path: Path,
    *,
    trait_missingness_threshold: float = 0.25,
) -> DatasetReadinessSummary:
    """Summarize whether a tree plus linked tables are ready for comparative analysis."""
    tree_validation = validate_tree_path(tree_path, strict=True)
    metadata = inspect_metadata_table(metadata_path)
    tree_taxa = set(load_tree(tree_path).tip_names)
    metadata_table = load_taxon_table(metadata_path)
    metadata_taxa = set(metadata_table.taxa)
    trait_name_check = check_tree_and_trait_taxon_names(tree_path, traits_path)
    traits_linkage = link_tree_to_traits(tree_path, traits_path)
    unusable_trait_columns = detect_unusable_trait_columns(
        traits_path,
        missingness_threshold=trait_missingness_threshold,
    )

    analysis_taxa = sorted(
        tree_taxa.intersection(metadata_taxa, traits_linkage.usable_taxa)
    )
    blockers: list[str] = []
    warnings: list[str] = []

    if tree_validation.branch_length_status != "complete":
        blockers.append("tree requires complete branch lengths")
    missing_metadata_taxa = sorted(tree_taxa - metadata_taxa)
    metadata_only_taxa = sorted(metadata_taxa - tree_taxa)
    if missing_metadata_taxa:
        blockers.append("metadata table is missing one or more tree taxa")
    if trait_name_check.tree_not_data:
        blockers.append("trait table is missing one or more tree taxa")
    if unusable_trait_columns:
        blockers.append("one or more trait columns exceed the missingness threshold")
    if len(analysis_taxa) < 2:
        blockers.append(
            "fewer than two taxa remain after intersecting tree, metadata, and traits"
        )
    if metadata_only_taxa:
        warnings.append("metadata table contains taxa absent from the tree")
    if trait_name_check.data_not_tree:
        warnings.append("trait table contains taxa absent from the tree")

    return DatasetReadinessSummary(
        tree_path=tree_path,
        metadata_path=metadata_path,
        traits_path=traits_path,
        tree_taxa=tree_validation.tip_count,
        analysis_taxa=analysis_taxa,
        missing_metadata_taxa=missing_metadata_taxa,
        missing_trait_taxa=trait_name_check.tree_not_data,
        metadata_only_taxa=metadata_only_taxa,
        trait_only_taxa=trait_name_check.data_not_tree,
        metadata_column_completeness=metadata.column_completeness,
        unusable_trait_columns=[column.name for column in unusable_trait_columns],
        ready_for_comparative_analysis=not blockers,
        blockers=blockers,
        warnings=warnings,
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
