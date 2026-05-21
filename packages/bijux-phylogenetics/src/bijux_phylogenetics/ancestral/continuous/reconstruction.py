from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.ancestral.common import (
    AncestralContinuousDataset,
    dump_pruned_tree,
    load_continuous_dataset,
    node_signature,
    stable_value,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

from .diagnostics import (
    _apply_anc_ml_profile_to_brownian_diagnostics,
    _summarize_brownian_fit_diagnostics,
)
from .estimators import (
    _build_anc_ml_estimates,
    _build_fast_anc_estimates,
    _build_local_continuous_estimates,
    _ContinuousAncMlProfileFit,
    _fit_continuous_anc_ml_profile,
    _sample_standard_deviation,
)
from .models import (
    ContinuousAncestralBrownianFitDiagnostics,
    ContinuousAncestralEstimate,
    ContinuousAncestralOptimizerDiagnostics,
    ContinuousAncestralReport,
)


def reconstruct_continuous_ancestral_states(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    model: str = "brownian",
    estimator: str | None = None,
    alpha: float = 1.0,
) -> ContinuousAncestralReport:
    """Reconstruct continuous ancestral states under a Brownian or OU-style model."""
    if model not in {"brownian", "ou"}:
        raise ValueError(f"unsupported continuous ancestral model: {model}")
    if alpha <= 0:
        raise ValueError(
            f"alpha must be positive for continuous ancestral reconstruction, got {alpha}"
        )
    dataset = load_continuous_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    return reconstruct_continuous_ancestral_states_from_dataset(
        dataset,
        model=model,
        estimator=estimator,
        alpha=alpha,
    )


def reconstruct_continuous_ancestral_states_from_dataset(
    dataset: AncestralContinuousDataset,
    *,
    model: str = "brownian",
    estimator: str | None = None,
    alpha: float = 1.0,
) -> ContinuousAncestralReport:
    """Reconstruct continuous ancestral states from one native ancestral dataset."""
    if model not in {"brownian", "ou"}:
        raise ValueError(f"unsupported continuous ancestral model: {model}")
    if alpha <= 0:
        raise ValueError(
            f"alpha must be positive for continuous ancestral reconstruction, got {alpha}"
        )
    resolved_estimator = _resolve_continuous_estimator(model, estimator)
    brownian_fit_diagnostics = (
        _summarize_brownian_fit_diagnostics(dataset) if model == "brownian" else None
    )
    optimizer_diagnostics: ContinuousAncestralOptimizerDiagnostics | None = None
    anc_ml_profile_fit: _ContinuousAncMlProfileFit | None = None
    if model == "brownian" and resolved_estimator == "anc-ml":
        anc_ml_profile_fit = _fit_continuous_anc_ml_profile(dataset)
        brownian_fit_diagnostics = _apply_anc_ml_profile_to_brownian_diagnostics(
            brownian_fit_diagnostics,
            log_likelihood=anc_ml_profile_fit.log_likelihood,
            sigma_squared=anc_ml_profile_fit.sigma_squared,
        )
        optimizer_diagnostics = anc_ml_profile_fit.optimizer_diagnostics
    return _reconstruct_continuous_from_dataset(
        dataset,
        working_tree=dataset.tree,
        model=model,
        estimator=resolved_estimator,
        alpha=alpha,
        brownian_fit_diagnostics=brownian_fit_diagnostics,
        optimizer_diagnostics=optimizer_diagnostics,
        anc_ml_profile_fit=anc_ml_profile_fit,
    )


def _reconstruct_continuous_from_dataset(
    dataset: AncestralContinuousDataset,
    *,
    working_tree: PhyloTree,
    model: str,
    estimator: str,
    alpha: float,
    brownian_fit_diagnostics: ContinuousAncestralBrownianFitDiagnostics | None,
    optimizer_diagnostics: ContinuousAncestralOptimizerDiagnostics | None,
    anc_ml_profile_fit: _ContinuousAncMlProfileFit | None,
) -> ContinuousAncestralReport:
    global_mean = sum(dataset.values_by_taxon[taxon] for taxon in dataset.taxa) / len(
        dataset.taxa
    )
    sigma = _sample_standard_deviation(
        [dataset.values_by_taxon[taxon] for taxon in dataset.taxa]
    )
    trait_range = (
        max(dataset.values_by_taxon.values()) - min(dataset.values_by_taxon.values())
        if dataset.values_by_taxon
        else 0.0
    )
    if estimator == "fast-anc":
        estimates = _build_fast_anc_estimates(
            dataset,
            trait_range=trait_range,
        )
    elif estimator == "anc-ml":
        if anc_ml_profile_fit is None:
            raise ValueError(
                "anc-ml continuous ancestral reconstruction requires one anc-ml profile fit"
            )
        estimates = _build_anc_ml_estimates(
            dataset,
            trait_range=trait_range,
            profile_fit=anc_ml_profile_fit,
        )
    else:
        estimates = _build_local_continuous_estimates(
            dataset,
            working_tree=working_tree,
            model=model,
            alpha=alpha,
            global_mean=global_mean,
            sigma=sigma,
            trait_range=trait_range,
        )
    ordered_estimates = _ordered_estimates(dataset, estimates)
    unstable_nodes = [
        estimate.node
        for estimate in ordered_estimates
        if not estimate.is_tip and estimate.unstable
    ]
    weak_support_nodes = [
        estimate.node
        for estimate in ordered_estimates
        if not estimate.is_tip and estimate.confidence < 0.75
    ]
    warnings = list(dataset.warnings)
    if unstable_nodes:
        warnings.append(
            "one or more continuous ancestral estimates have broad uncertainty intervals"
        )
    if weak_support_nodes:
        warnings.append(
            "low-confidence ancestral estimates should not be overinterpreted for evolutionary timing or trait polarity"
        )
    if (
        brownian_fit_diagnostics is not None
        and brownian_fit_diagnostics.covariance_near_singular
    ):
        warnings.append(
            "Brownian covariance diagnostics indicate a singular or ill-conditioned fit, so ancestral uncertainty should be interpreted cautiously"
        )
    if (
        brownian_fit_diagnostics is not None
        and brownian_fit_diagnostics.solver_regularized
    ):
        warnings.append(
            "Brownian covariance inversion required light diagonal regularization for numerical stability"
        )
    return ContinuousAncestralReport(
        tree_path=dataset.tree_path,
        traits_path=dataset.traits_path,
        taxon_column=dataset.taxon_column,
        trait=dataset.trait,
        model=model,
        estimator=estimator,
        alpha=stable_value(alpha),
        taxon_count=len(dataset.taxa),
        analysis_tree_newick=dump_pruned_tree(working_tree),
        missing_from_traits_taxa=dataset.missing_from_traits_taxa,
        dropped_missing_taxa=dataset.dropped_missing_taxa,
        dropped_non_numeric_taxa=dataset.dropped_non_numeric_taxa,
        warnings=warnings,
        unstable_nodes=unstable_nodes,
        weak_support_nodes=weak_support_nodes,
        brownian_fit_diagnostics=brownian_fit_diagnostics,
        optimizer_diagnostics=optimizer_diagnostics,
        estimates=ordered_estimates,
    )


def _ordered_estimates(
    dataset: AncestralContinuousDataset,
    estimates: list[ContinuousAncestralEstimate],
) -> list[ContinuousAncestralEstimate]:
    node_order = {
        node_signature(node): index
        for index, node in enumerate(dataset.tree.iter_nodes())
    }
    return sorted(estimates, key=lambda estimate: node_order[estimate.node])


def _resolve_continuous_estimator(model: str, estimator: str | None) -> str:
    default_estimators = {
        "brownian": "ace-pic",
        "ou": "generalized-least-squares",
    }
    allowed_estimators = {
        "brownian": {"ace-pic", "anc-ml", "fast-anc"},
        "ou": {"generalized-least-squares"},
    }
    resolved = default_estimators[model] if estimator is None else estimator
    if resolved not in allowed_estimators[model]:
        supported = ", ".join(sorted(allowed_estimators[model]))
        raise ValueError(
            f"unsupported continuous ancestral estimator '{resolved}' for model '{model}'; expected one of: {supported}"
        )
    return resolved
