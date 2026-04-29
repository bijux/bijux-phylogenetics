from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.bayesian.beast import (
    FossilCalibrationValidationReport,
    TipDatingValidationReport,
    validate_fossil_calibration_table,
    validate_tip_dating_metadata,
)
from bijux_phylogenetics.core.alignment import AlignmentForensicReport
from bijux_phylogenetics.core.metadata import MetadataColumnCompleteness, inspect_metadata_table, load_taxon_table
from bijux_phylogenetics.core.traits import detect_unusable_trait_columns, link_tree_to_traits
from bijux_phylogenetics.diagnostics.validation import validate_tree_path
from bijux_phylogenetics.io.fasta import build_alignment_forensic_report, link_alignment_to_tree
from bijux_phylogenetics.io.trees import load_tree


@dataclass(slots=True)
class DatasetReadinessSummary:
    """Readiness summary for a tree plus linked metadata and trait tables."""

    tree_path: Path
    metadata_path: Path
    traits_path: Path
    tree_taxa: int
    analysis_taxa: list[str]
    missing_metadata_taxa: list[str]
    missing_trait_taxa: list[str]
    metadata_only_taxa: list[str]
    trait_only_taxa: list[str]
    metadata_column_completeness: list[MetadataColumnCompleteness]
    unusable_trait_columns: list[str]
    ready_for_comparative_analysis: bool
    blockers: list[str]
    warnings: list[str]


@dataclass(slots=True)
class DatasetAuditReport:
    """Integrated audit across tree, alignment, metadata, traits, and optional dating surfaces."""

    tree_path: Path
    metadata_path: Path
    traits_path: Path
    alignment_path: Path | None
    tip_dates_path: Path | None
    calibration_path: Path | None
    tree_taxa: int
    analysis_taxa: list[str]
    readiness_decision: str
    allowed_analyses: list[str]
    blocked_analyses: list[str]
    blockers: list[str]
    warnings: list[str]
    dataset_readiness: DatasetReadinessSummary
    alignment_forensic: AlignmentForensicReport | None
    tip_dates: TipDatingValidationReport | None
    calibrations: FossilCalibrationValidationReport | None


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
    if traits_linkage.missing_from_traits:
        blockers.append("trait table is missing one or more tree taxa")
    if unusable_trait_columns:
        blockers.append("one or more trait columns exceed the missingness threshold")
    if len(analysis_taxa) < 2:
        blockers.append("fewer than two taxa remain after intersecting tree, metadata, and traits")
    if metadata_only_taxa:
        warnings.append("metadata table contains taxa absent from the tree")
    if traits_linkage.extra_trait_taxa:
        warnings.append("trait table contains taxa absent from the tree")

    return DatasetReadinessSummary(
        tree_path=tree_path,
        metadata_path=metadata_path,
        traits_path=traits_path,
        tree_taxa=tree_validation.tip_count,
        analysis_taxa=analysis_taxa,
        missing_metadata_taxa=missing_metadata_taxa,
        missing_trait_taxa=traits_linkage.missing_from_traits,
        metadata_only_taxa=metadata_only_taxa,
        trait_only_taxa=traits_linkage.extra_trait_taxa,
        metadata_column_completeness=metadata.column_completeness,
        unusable_trait_columns=[column.name for column in unusable_trait_columns],
        ready_for_comparative_analysis=not blockers,
        blockers=blockers,
        warnings=warnings,
    )


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
    tree_taxa = set(load_tree(tree_path).tip_names)
    metadata_taxa = set(load_taxon_table(metadata_path).taxa)
    trait_taxa = set(load_taxon_table(traits_path).taxa)
    blockers = list(readiness.blockers)
    warnings = list(readiness.warnings)

    alignment_forensic: AlignmentForensicReport | None = None
    alignment_taxa = set(tree_taxa)
    if alignment_path is not None:
        alignment_linkage = link_alignment_to_tree(tree_path, alignment_path)
        alignment_forensic = build_alignment_forensic_report(alignment_path)
        alignment_taxa = set(alignment_linkage.usable_taxa)
        if alignment_linkage.missing_from_alignment:
            blockers.append("alignment is missing one or more tree taxa")
        if alignment_linkage.extra_alignment_ids:
            warnings.append("alignment contains taxa absent from the tree")
        if not alignment_forensic.safe_for_distance_analysis and not alignment_forensic.safe_for_maximum_likelihood:
            blockers.append("alignment is not currently safe for core inference workflows")
        elif alignment_forensic.warnings:
            warnings.extend(alignment_forensic.warnings)

    tip_dates: TipDatingValidationReport | None = None
    if tip_dates_path is not None:
        tip_dates = validate_tip_dating_metadata(
            tree_path,
            tip_dates_path,
            alignment_path=alignment_path,
        )
        if tip_dates.invalid_tip_count > 0 or tip_dates.missing_tree_taxa:
            blockers.append("tip-date metadata contains invalid or missing tree-tip dates")
        if tip_dates.extra_tip_taxa or tip_dates.extra_alignment_taxa:
            warnings.append("tip-date metadata contains taxa absent from the tree or alignment")

    calibrations: FossilCalibrationValidationReport | None = None
    if calibration_path is not None:
        calibrations = validate_fossil_calibration_table(tree_path, calibration_path)
        if calibrations.invalid_calibration_count > 0:
            blockers.append("calibration table contains invalid fossil calibration targets or ages")

    analysis_taxa = sorted(tree_taxa & metadata_taxa & trait_taxa & alignment_taxa)
    if len(analysis_taxa) < 2:
        blockers.append("fewer than two taxa remain after intersecting all requested dataset surfaces")

    blocked_analyses: list[str] = []
    allowed_analyses: list[str] = []
    if blockers:
        blocked_analyses.extend(["comparative", "publication"])
    else:
        allowed_analyses.extend(["comparative", "publication"])

    if alignment_forensic is None:
        blocked_analyses.extend(["distance", "maximum_likelihood", "bayesian"])
    else:
        if alignment_forensic.safe_for_distance_analysis:
            allowed_analyses.append("distance")
        else:
            blocked_analyses.append("distance")
        if alignment_forensic.safe_for_maximum_likelihood:
            allowed_analyses.append("maximum_likelihood")
        else:
            blocked_analyses.append("maximum_likelihood")
        if alignment_forensic.safe_for_bayesian_inference:
            allowed_analyses.append("bayesian")
        else:
            blocked_analyses.append("bayesian")
        if alignment_forensic.safe_for_coding_analysis:
            allowed_analyses.append("coding")
        else:
            blocked_analyses.append("coding")

    if tip_dates is not None or calibrations is not None:
        if tip_dates is not None and calibrations is not None and tip_dates.invalid_tip_count == 0 and calibrations.invalid_calibration_count == 0:
            allowed_analyses.append("time_tree")
        else:
            blocked_analyses.append("time_tree")

    readiness_decision = "blocked" if blockers else ("ready_with_warnings" if warnings else "ready")
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
        allowed_analyses=sorted(dict.fromkeys(allowed_analyses)),
        blocked_analyses=sorted(dict.fromkeys(blocked_analyses)),
        blockers=blockers,
        warnings=sorted(dict.fromkeys(warnings)),
        dataset_readiness=readiness,
        alignment_forensic=alignment_forensic,
        tip_dates=tip_dates,
        calibrations=calibrations,
    )
