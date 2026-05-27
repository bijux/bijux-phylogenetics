from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.phylo.alignment import DnaBinAlignment
from bijux_phylogenetics.phylo.topology.bionj import build_bionj_tree as _build_bionj_tree_impl
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

from .models import (
    AmbiguityPolicy as AmbiguityPolicy,
    BalancedMinimumEvolutionNniSearchReport as BalancedMinimumEvolutionNniSearchReport,
    BalancedMinimumEvolutionNniTraceRow as BalancedMinimumEvolutionNniTraceRow,
    DistanceBootstrapReport as DistanceBootstrapReport,
    DistanceBootstrapSupportSummary as DistanceBootstrapSupportSummary,
    DistanceGapPolicySensitivityReport as DistanceGapPolicySensitivityReport,
    DistanceMatrixQualityReport as DistanceMatrixQualityReport,
    DistanceMethodAssessment as DistanceMethodAssessment,
    DistanceMethodAssumptionReport as DistanceMethodAssumptionReport,
    DistanceMethodMaturityGateReport as DistanceMethodMaturityGateReport,
    DistanceMethodReport as DistanceMethodReport,
    DistanceModel as DistanceModel,
    DistanceModelComparisonReport as DistanceModelComparisonReport,
    DistanceOutlierPair as DistanceOutlierPair,
    DistanceReferenceObservation as DistanceReferenceObservation,
    DistanceReferenceValidationReport as DistanceReferenceValidationReport,
    DistanceReproducibilityBundleReport as DistanceReproducibilityBundleReport,
    DistanceTreeBuildReport as DistanceTreeBuildReport,
    DistanceTreeMethodPolicy as DistanceTreeMethodPolicy,
    DistanceTreeReferenceComparisonReport as DistanceTreeReferenceComparisonReport,
    DistanceTreeReferenceObservation as DistanceTreeReferenceObservation,
    DistanceTreeTopologyComparison as DistanceTreeTopologyComparison,
    FitchMargoliashBranchFit as FitchMargoliashBranchFit,
    FitchMargoliashFitReport as FitchMargoliashFitReport,
    GapHandlingMode as GapHandlingMode,
    GeneticDistanceMatrix as GeneticDistanceMatrix,
    GeneticDistanceModelParameters as GeneticDistanceModelParameters,
    ImportedDistanceEntry as ImportedDistanceEntry,
    ImportedDistanceMatrixQualityReport as ImportedDistanceMatrixQualityReport,
    ImportedDistanceMatrixReport as ImportedDistanceMatrixReport,
    ImportedDistanceTreeBuildReport as ImportedDistanceTreeBuildReport,
    LowInformationPair as LowInformationPair,
    MinimumEvolutionBranchFit as MinimumEvolutionBranchFit,
    MinimumEvolutionScoreReport as MinimumEvolutionScoreReport,
    NonMetricDistanceObservation as NonMetricDistanceObservation,
    NonnegativeLeastSquaresActiveConstraint as NonnegativeLeastSquaresActiveConstraint,
    NonnegativeLeastSquaresBranchFit as NonnegativeLeastSquaresBranchFit,
    NonnegativeLeastSquaresFitReport as NonnegativeLeastSquaresFitReport,
    OrdinaryLeastSquaresBranchFit as OrdinaryLeastSquaresBranchFit,
    OrdinaryLeastSquaresFitReport as OrdinaryLeastSquaresFitReport,
    PatristicResidualDiagnosticsReport as PatristicResidualDiagnosticsReport,
    PatristicResidualRow as PatristicResidualRow,
    SaturatedDistancePair as SaturatedDistancePair,
    UPGMAUltrametricViolation as UPGMAUltrametricViolation,
)
from .shared import (
    _allowed_models_for_alphabet as _allowed_models_for_alphabet,
    _file_sha256 as _file_sha256,
    _iter_ultrametric_violations as _iter_ultrametric_violations,
    _normalize_distance_model as _normalize_distance_model,
    _pair_key as _pair_key,
    _require_supported_distance_tree_method as _require_supported_distance_tree_method,
    _unique_genetic_distance_pairs as _unique_genetic_distance_pairs,
    list_distance_tree_method_policies as list_distance_tree_method_policies,
    resolve_distance_tree_method_policy as resolve_distance_tree_method_policy,
)
from .upgma import (
    UPGMABuildReport as UPGMABuildReport,
    UPGMAClusterHeightRow as UPGMAClusterHeightRow,
    UPGMAMergeRow as UPGMAMergeRow,
)
from .single_linkage import (
    SingleLinkageBuildReport as SingleLinkageBuildReport,
    SingleLinkageClusterHeightRow as SingleLinkageClusterHeightRow,
    SingleLinkageMergeRow as SingleLinkageMergeRow,
)
from .complete_linkage import (
    CompleteLinkageBuildReport as CompleteLinkageBuildReport,
    CompleteLinkageClusterHeightRow as CompleteLinkageClusterHeightRow,
    CompleteLinkageMergeRow as CompleteLinkageMergeRow,
)
from .wpgma import (
    WPGMABuildReport as WPGMABuildReport,
    WPGMAClusterHeightRow as WPGMAClusterHeightRow,
    WPGMAMergeRow as WPGMAMergeRow,
)


def _load_alignment_for_model(path: Path, *, model: DistanceModel):
    from .matrix import _load_alignment_for_model as load_alignment_impl

    return load_alignment_impl(path, model=model)


def compute_pairwise_genetic_distance_matrix_from_dna_bin_alignment(
    alignment: DnaBinAlignment,
    *,
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
) -> GeneticDistanceMatrix:
    from .matrix import (
        compute_pairwise_genetic_distance_matrix_from_dna_bin_alignment as compute_impl,
    )

    return compute_impl(
        alignment,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )


def _build_alignment_distance_lookup(
    report: GeneticDistanceMatrix,
) -> dict[tuple[str, str], float]:
    from .matrix import _build_alignment_distance_lookup as lookup_impl

    return lookup_impl(report)


def assess_distance_method_assumptions(
    path: Path,
    *,
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
    ultrametric_tolerance: float = 1e-6,
) -> DistanceMethodAssumptionReport:
    from .quality import assess_distance_method_assumptions as assess_impl

    return assess_impl(
        path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
        ultrametric_tolerance=ultrametric_tolerance,
    )


def assess_distance_method_assumptions_from_genetic_distance_matrix(
    report: GeneticDistanceMatrix,
    *,
    ultrametric_tolerance: float = 1e-6,
) -> DistanceMethodAssumptionReport:
    from .quality import (
        assess_distance_method_assumptions_from_genetic_distance_matrix as assess_impl,
    )

    return assess_impl(report, ultrametric_tolerance=ultrametric_tolerance)


def _bio_distance_matrix(report: GeneticDistanceMatrix):
    from .matrix import _bio_distance_matrix as bio_distance_matrix_impl

    return bio_distance_matrix_impl(report)


def _distance_lookup(report: GeneticDistanceMatrix) -> dict[tuple[str, str], float]:
    from .matrix import _distance_lookup as distance_lookup_impl

    return distance_lookup_impl(report)


def build_upgma_tree(
    identifiers: list[str],
    distance_lookup: dict[tuple[str, str], float],
) -> tuple[PhyloTree, UPGMABuildReport]:
    from .upgma import build_upgma_tree as build_upgma_impl

    return build_upgma_impl(identifiers, distance_lookup)


def build_bionj_tree(
    identifiers: list[str],
    distance_lookup: dict[tuple[str, str], float],
) -> PhyloTree:
    return _build_bionj_tree_impl(identifiers, distance_lookup)


def build_wpgma_tree(
    identifiers: list[str],
    distance_lookup: dict[tuple[str, str], float],
) -> tuple[PhyloTree, WPGMABuildReport]:
    from .wpgma import build_wpgma_tree as build_wpgma_impl

    return build_wpgma_impl(identifiers, distance_lookup)


def build_single_linkage_tree(
    identifiers: list[str],
    distance_lookup: dict[tuple[str, str], float],
) -> tuple[PhyloTree, SingleLinkageBuildReport]:
    from .single_linkage import (
        build_single_linkage_tree as build_single_linkage_impl,
    )

    return build_single_linkage_impl(identifiers, distance_lookup)


def build_complete_linkage_tree(
    identifiers: list[str],
    distance_lookup: dict[tuple[str, str], float],
) -> tuple[PhyloTree, CompleteLinkageBuildReport]:
    from .complete_linkage import (
        build_complete_linkage_tree as build_complete_linkage_impl,
    )

    return build_complete_linkage_impl(identifiers, distance_lookup)


def load_imported_distance_matrix(path: Path) -> list[ImportedDistanceEntry]:
    from .imported import load_imported_distance_matrix as load_impl

    return load_impl(path)


def assess_imported_distance_method_assumptions(
    path: Path,
    *,
    ultrametric_tolerance: float = 1e-6,
) -> DistanceMethodAssumptionReport:
    from .imported import (
        assess_imported_distance_method_assumptions as assess_imported_impl,
    )

    return assess_imported_impl(path, ultrametric_tolerance=ultrametric_tolerance)


def validate_imported_distance_matrix(path: Path) -> ImportedDistanceMatrixReport:
    from .imported import validate_imported_distance_matrix as validate_impl

    return validate_impl(path)


def inspect_imported_distance_matrix_quality(
    path: Path,
) -> ImportedDistanceMatrixQualityReport:
    from .imported import (
        inspect_imported_distance_matrix_quality as inspect_quality_impl,
    )

    return inspect_quality_impl(path)


def build_tree_from_imported_distance_matrix(
    path: Path,
    *,
    method: str,
) -> tuple[PhyloTree, ImportedDistanceTreeBuildReport]:
    from .imported import (
        build_tree_from_imported_distance_matrix as build_imported_tree_impl,
    )

    return build_imported_tree_impl(path, method=method)


def fit_minimum_evolution_tree(
    tree: PhyloTree,
    identifiers: list[str],
    distance_lookup: dict[tuple[str, str], float],
) -> tuple[PhyloTree, MinimumEvolutionScoreReport]:
    from .minimum_evolution import fit_minimum_evolution_tree as fit_impl

    return fit_impl(tree, identifiers, distance_lookup)


def fit_minimum_evolution_tree_from_imported_distance_matrix(
    matrix_path: Path,
    tree_path: Path,
) -> tuple[PhyloTree, MinimumEvolutionScoreReport]:
    from .minimum_evolution import (
        fit_minimum_evolution_tree_from_imported_distance_matrix as fit_impl,
    )

    return fit_impl(matrix_path, tree_path)


def fit_fitch_margoliash_tree(
    tree: PhyloTree,
    identifiers: list[str],
    distance_lookup: dict[tuple[str, str], float],
    *,
    weighting_power: float = 2.0,
) -> tuple[PhyloTree, FitchMargoliashFitReport]:
    from .fitch_margoliash import fit_fitch_margoliash_tree as fit_impl

    return fit_impl(
        tree,
        identifiers,
        distance_lookup,
        weighting_power=weighting_power,
    )


def fit_fitch_margoliash_tree_from_imported_distance_matrix(
    matrix_path: Path,
    tree_path: Path,
    *,
    weighting_power: float = 2.0,
) -> tuple[PhyloTree, FitchMargoliashFitReport]:
    from .fitch_margoliash import (
        fit_fitch_margoliash_tree_from_imported_distance_matrix as fit_impl,
    )

    return fit_impl(
        matrix_path,
        tree_path,
        weighting_power=weighting_power,
    )


def fit_ordinary_least_squares_tree(
    tree: PhyloTree,
    identifiers: list[str],
    distance_lookup: dict[tuple[str, str], float],
) -> tuple[PhyloTree, OrdinaryLeastSquaresFitReport]:
    from .ordinary_least_squares import fit_ordinary_least_squares_tree as fit_impl

    return fit_impl(tree, identifiers, distance_lookup)


def fit_ordinary_least_squares_tree_from_imported_distance_matrix(
    matrix_path: Path,
    tree_path: Path,
) -> tuple[PhyloTree, OrdinaryLeastSquaresFitReport]:
    from .ordinary_least_squares import (
        fit_ordinary_least_squares_tree_from_imported_distance_matrix as fit_impl,
    )

    return fit_impl(matrix_path, tree_path)


def fit_nonnegative_least_squares_tree(
    tree: PhyloTree,
    identifiers: list[str],
    distance_lookup: dict[tuple[str, str], float],
) -> tuple[PhyloTree, NonnegativeLeastSquaresFitReport]:
    from .nonnegative_least_squares import (
        fit_nonnegative_least_squares_tree as fit_impl,
    )

    return fit_impl(tree, identifiers, distance_lookup)


def fit_nonnegative_least_squares_tree_from_imported_distance_matrix(
    matrix_path: Path,
    tree_path: Path,
) -> tuple[PhyloTree, NonnegativeLeastSquaresFitReport]:
    from .nonnegative_least_squares import (
        fit_nonnegative_least_squares_tree_from_imported_distance_matrix as fit_impl,
    )

    return fit_impl(matrix_path, tree_path)


def compute_patristic_residual_diagnostics(
    tree: PhyloTree,
    identifiers: list[str],
    distance_lookup: dict[tuple[str, str], float],
) -> PatristicResidualDiagnosticsReport:
    from .patristic_residuals import (
        compute_patristic_residual_diagnostics as compute_impl,
    )

    return compute_impl(tree, identifiers, distance_lookup)


def compute_patristic_residual_diagnostics_from_imported_distance_matrix(
    matrix_path: Path,
    tree_path: Path,
) -> PatristicResidualDiagnosticsReport:
    from .patristic_residuals import (
        compute_patristic_residual_diagnostics_from_imported_distance_matrix as compute_impl,
    )

    return compute_impl(matrix_path, tree_path)


def write_patristic_residual_artifacts(
    out_dir: Path,
    report: PatristicResidualDiagnosticsReport,
) -> dict[str, Path]:
    from .patristic_residuals import (
        write_patristic_residual_artifacts as write_impl,
    )

    return write_impl(out_dir, report)


def score_balanced_minimum_evolution(
    tree: PhyloTree,
    identifiers: list[str],
    distance_lookup: dict[tuple[str, str], float],
) -> float:
    from .balanced_minimum_evolution import score_balanced_minimum_evolution as score_impl

    return score_impl(tree, identifiers, distance_lookup)


def search_balanced_minimum_evolution_nni(
    identifiers: list[str],
    distance_lookup: dict[tuple[str, str], float],
    *,
    start_method: str,
) -> BalancedMinimumEvolutionNniSearchReport:
    from .balanced_minimum_evolution_nni import (
        search_balanced_minimum_evolution_nni as search_impl,
    )

    return search_impl(identifiers, distance_lookup, start_method=start_method)


def search_balanced_minimum_evolution_nni_from_imported_distance_matrix(
    matrix_path: Path,
    *,
    start_method: str,
) -> BalancedMinimumEvolutionNniSearchReport:
    from .balanced_minimum_evolution_nni import (
        search_balanced_minimum_evolution_nni_from_imported_distance_matrix as search_impl,
    )

    return search_impl(matrix_path, start_method=start_method)


def write_balanced_minimum_evolution_nni_artifacts(
    out_dir: Path,
    report: BalancedMinimumEvolutionNniSearchReport,
) -> dict[str, Path]:
    from .balanced_minimum_evolution_nni import (
        write_balanced_minimum_evolution_nni_artifacts as write_impl,
    )

    return write_impl(out_dir, report)


def build_distance_tree(
    path: Path,
    *,
    method: str,
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
) -> tuple[PhyloTree, DistanceTreeBuildReport]:
    from .tree_inference import build_distance_tree as build_tree_impl

    return build_tree_impl(
        path,
        method=method,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )


def build_distance_tree_from_genetic_distance_matrix(
    report: GeneticDistanceMatrix,
    *,
    method: str,
) -> tuple[PhyloTree, DistanceTreeBuildReport]:
    from .tree_inference import (
        build_distance_tree_from_genetic_distance_matrix as build_from_matrix_impl,
    )

    return build_from_matrix_impl(report, method=method)


def compare_distance_tree_topologies(
    path: Path,
    *,
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
) -> DistanceTreeTopologyComparison:
    from .tree_inference import compare_distance_tree_topologies as compare_trees_impl

    return compare_trees_impl(
        path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )


def bootstrap_distance_trees(
    path: Path,
    *,
    method: str,
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
    replicates: int = 100,
    seed: int = 1,
) -> tuple[list[PhyloTree], DistanceBootstrapReport]:
    from .tree_inference import bootstrap_distance_trees as bootstrap_impl

    return bootstrap_impl(
        path,
        method=method,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
        replicates=replicates,
        seed=seed,
    )


def write_distance_bootstrap_support(
    path: Path, report: DistanceBootstrapReport
) -> Path:
    from .tree_inference import (
        write_distance_bootstrap_support as write_support_impl,
    )

    return write_support_impl(path, report)


def summarize_distance_bootstrap_support(
    report: DistanceBootstrapReport,
    *,
    weak_frequency_threshold: float = 0.5,
) -> DistanceBootstrapSupportSummary:
    from .tree_inference import (
        summarize_distance_bootstrap_support as summarize_support_impl,
    )

    return summarize_support_impl(
        report, weak_frequency_threshold=weak_frequency_threshold
    )


def compare_distance_tree_to_reference_tree(
    path: Path,
    reference_tree_path: Path,
    *,
    method: str,
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
) -> DistanceTreeReferenceComparisonReport:
    from .tree_inference import (
        compare_distance_tree_to_reference_tree as compare_reference_impl,
    )

    return compare_reference_impl(
        path,
        reference_tree_path,
        method=method,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )


def compute_pairwise_genetic_distance_matrix(
    path: Path,
    *,
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
) -> GeneticDistanceMatrix:
    from .matrix import compute_pairwise_genetic_distance_matrix as compute_impl

    return compute_impl(
        path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )


def validate_distance_reference_examples() -> DistanceReferenceValidationReport:
    from .validation import validate_distance_reference_examples as validate_impl

    return validate_impl()


def inspect_distance_matrix_quality(
    path: Path,
    *,
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
) -> DistanceMatrixQualityReport:
    from .quality import inspect_distance_matrix_quality as inspect_impl

    return inspect_impl(
        path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )


def write_genetic_distance_matrix(path: Path, report: GeneticDistanceMatrix) -> Path:
    from .matrix import write_genetic_distance_matrix as write_matrix_impl

    return write_matrix_impl(path, report)


def write_genetic_distance_component_table(
    path: Path, report: GeneticDistanceMatrix
) -> Path:
    from .matrix import (
        write_genetic_distance_component_table as write_component_table_impl,
    )

    return write_component_table_impl(path, report)


def write_genetic_distance_parameter_table(
    path: Path, report: GeneticDistanceMatrix
) -> Path:
    from .matrix import (
        write_genetic_distance_parameter_table as write_parameter_table_impl,
    )

    return write_parameter_table_impl(path, report)


def compare_distance_models(
    path: Path,
    *,
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
) -> DistanceModelComparisonReport:
    from .reporting import compare_distance_models as compare_models_impl

    return compare_models_impl(
        path,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )


def compare_distance_gap_policies(
    path: Path,
    *,
    model: DistanceModel = "p-distance",
    ambiguity_policy: AmbiguityPolicy = "ignore",
) -> DistanceGapPolicySensitivityReport:
    from .reporting import (
        compare_distance_gap_policies as compare_gap_policies_impl,
    )

    return compare_gap_policies_impl(
        path,
        model=model,
        ambiguity_policy=ambiguity_policy,
    )


def assess_distance_method_maturity(
    path: Path,
    *,
    method: str = "neighbor-joining",
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
    bootstrap_replicates: int = 25,
    bootstrap_seed: int = 1,
    validate_bundle: bool = False,
) -> DistanceMethodMaturityGateReport:
    from .reporting import (
        assess_distance_method_maturity as assess_maturity_impl,
    )

    return assess_maturity_impl(
        path,
        method=method,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
        bootstrap_replicates=bootstrap_replicates,
        bootstrap_seed=bootstrap_seed,
        validate_bundle=validate_bundle,
    )


def build_distance_method_report(
    path: Path,
    *,
    method: str = "neighbor-joining",
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
    bootstrap_replicates: int = 25,
    bootstrap_seed: int = 1,
) -> DistanceMethodReport:
    from .reporting import build_distance_method_report as build_report_impl

    return build_report_impl(
        path,
        method=method,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
        bootstrap_replicates=bootstrap_replicates,
        bootstrap_seed=bootstrap_seed,
    )


def write_distance_reproducibility_bundle(
    out_dir: Path,
    *,
    alignment_path: Path,
    method: str,
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
    replicates: int = 100,
    seed: int = 1,
) -> DistanceReproducibilityBundleReport:
    from .reporting import (
        write_distance_reproducibility_bundle as write_bundle_impl,
    )

    return write_bundle_impl(
        out_dir,
        alignment_path=alignment_path,
        method=method,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
        replicates=replicates,
        seed=seed,
    )
