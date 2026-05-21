from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.comparative.common import load_comparative_dataset
from bijux_phylogenetics.comparative.evolutionary_modes.comparison import (
    compare_selected_continuous_modes as _compare_selected_continuous_modes,
)
from bijux_phylogenetics.comparative.evolutionary_modes.fitting import (
    fit_continuous_evolutionary_mode,
)
from bijux_phylogenetics.comparative.evolutionary_modes.fitting import (
    fit_evolutionary_mode_from_dataset as _fit_evolutionary_mode_from_dataset,
)
from bijux_phylogenetics.comparative.evolutionary_modes.models import (
    CONTINUOUS_GAUSSIAN_LIKELIHOOD_COMPARISON_POLICY,
    CONTINUOUS_GAUSSIAN_LIKELIHOOD_CONSTANT_POLICY,
    FITCONTINUOUS_MODEL_COMPARISON_ORDER,
    FITCONTINUOUS_MODEL_CONFIDENCE_DELTA_THRESHOLD,
    FITCONTINUOUS_MODEL_CONFIDENCE_WEIGHT_BASIS,
    FITCONTINUOUS_MODEL_RANKING_POLICY,
    ComparativeTreeRescalingReport,
    ContinuousEvolutionaryModeComparisonReport,
    ContinuousEvolutionaryModeFitReport,
    ContinuousModeBoundaryAssessment,
    ContinuousModeOptimizerDiagnostics,
    ContinuousModeOptimizerProfileRow,
    ContinuousModeSearchControls,
    EvolutionaryModeBranchLengthRow,
    EvolutionaryModeIdentifiabilityWarning,
    LikelihoodRatioTestResult,
)
from bijux_phylogenetics.comparative.evolutionary_modes.tree_transforms import (
    rescale_tree_early_burst,
    rescale_tree_ornstein_uhlenbeck,
    rescale_tree_pagel_delta,
    rescale_tree_pagel_kappa,
    rescale_tree_pagel_lambda,
    rescale_tree_white_noise,
    transform_tree_for_evolutionary_mode,
)

__all__ = [
    "CONTINUOUS_GAUSSIAN_LIKELIHOOD_COMPARISON_POLICY",
    "CONTINUOUS_GAUSSIAN_LIKELIHOOD_CONSTANT_POLICY",
    "FITCONTINUOUS_MODEL_CONFIDENCE_DELTA_THRESHOLD",
    "FITCONTINUOUS_MODEL_CONFIDENCE_WEIGHT_BASIS",
    "FITCONTINUOUS_MODEL_RANKING_POLICY",
    "ComparativeTreeRescalingReport",
    "ContinuousEvolutionaryModeComparisonReport",
    "ContinuousEvolutionaryModeFitReport",
    "ContinuousModeBoundaryAssessment",
    "ContinuousModeOptimizerDiagnostics",
    "ContinuousModeOptimizerProfileRow",
    "ContinuousModeSearchControls",
    "EvolutionaryModeBranchLengthRow",
    "EvolutionaryModeIdentifiabilityWarning",
    "LikelihoodRatioTestResult",
    "compare_continuous_evolutionary_modes",
    "compare_fitcontinuous_model_ranking",
    "fit_continuous_evolutionary_mode",
    "rescale_tree_early_burst",
    "rescale_tree_ornstein_uhlenbeck",
    "rescale_tree_pagel_delta",
    "rescale_tree_pagel_kappa",
    "rescale_tree_pagel_lambda",
    "rescale_tree_white_noise",
    "transform_tree_for_evolutionary_mode",
]


def compare_continuous_evolutionary_modes(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    standard_error_trait: str | None = None,
    lambda_bounds: tuple[float, float] = (0.0, 1.0),
    kappa_bounds: tuple[float, float] = (0.0, 3.0),
    delta_bounds: tuple[float, float] = (0.0, 3.0),
    ou_bounds: tuple[float, float] = (0.0, 10.0),
    early_burst_bounds: tuple[float, float] = (0.0, 50.0),
) -> ContinuousEvolutionaryModeComparisonReport:
    """Compare the legacy Brownian, OU, and early-burst mode trio."""
    return _compare_selected_continuous_modes(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        modes=("brownian", "ornstein-uhlenbeck", "early-burst"),
        standard_error_trait=standard_error_trait,
        lambda_bounds=lambda_bounds,
        kappa_bounds=kappa_bounds,
        delta_bounds=delta_bounds,
        ou_bounds=ou_bounds,
        early_burst_bounds=early_burst_bounds,
        dataset_loader=load_comparative_dataset,
        fit_mode_from_dataset=_fit_evolutionary_mode_from_dataset,
    )


def compare_fitcontinuous_model_ranking(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    modes: tuple[str, ...] | None = None,
    standard_error_trait: str | None = None,
    lambda_bounds: tuple[float, float] = (0.0, 1.0),
    kappa_bounds: tuple[float, float] = (0.0, 3.0),
    delta_bounds: tuple[float, float] = (0.0, 3.0),
    ou_bounds: tuple[float, float] = (0.0, 10.0),
    early_burst_bounds: tuple[float, float] = (0.0, 50.0),
) -> ContinuousEvolutionaryModeComparisonReport:
    """Compare the governed `fitContinuous` model set by AIC and AICc."""
    return _compare_selected_continuous_modes(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        modes=FITCONTINUOUS_MODEL_COMPARISON_ORDER if modes is None else modes,
        standard_error_trait=standard_error_trait,
        lambda_bounds=lambda_bounds,
        kappa_bounds=kappa_bounds,
        delta_bounds=delta_bounds,
        ou_bounds=ou_bounds,
        early_burst_bounds=early_burst_bounds,
        dataset_loader=load_comparative_dataset,
        fit_mode_from_dataset=_fit_evolutionary_mode_from_dataset,
    )
