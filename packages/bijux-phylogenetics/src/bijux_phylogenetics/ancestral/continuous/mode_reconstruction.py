from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.ancestral.common import load_continuous_dataset
from bijux_phylogenetics.ancestral.continuous.models import (
    ContinuousAncestralReport,
)
from bijux_phylogenetics.ancestral.continuous.reconstruction import (
    _reconstruct_continuous_from_dataset,
    reconstruct_continuous_ancestral_states,
)


@dataclass(slots=True)
class ContinuousEvolutionaryModeAncestralReport:
    """Continuous ancestral reconstruction over a governed transformed-tree mode."""

    tree_path: Path
    traits_path: Path
    trait: str
    mode: str
    parameter_name: str | None
    parameter_value: float | None
    transformed_tree_newick: str
    reconstruction: ContinuousAncestralReport


def reconstruct_continuous_evolutionary_mode_states(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    mode: str = "brownian",
    rate_change: float = 0.0,
) -> ContinuousEvolutionaryModeAncestralReport:
    """Reconstruct continuous ancestral states under Brownian or early-burst tree modes."""
    if mode == "brownian":
        reconstruction = reconstruct_continuous_ancestral_states(
            tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
            model="brownian",
        )
        return ContinuousEvolutionaryModeAncestralReport(
            tree_path=tree_path,
            traits_path=traits_path,
            trait=trait,
            mode="brownian",
            parameter_name=None,
            parameter_value=None,
            transformed_tree_newick=reconstruction.analysis_tree_newick,
            reconstruction=reconstruction,
        )
    if mode != "early-burst":
        raise ValueError(
            "unsupported evolutionary mode for continuous ancestral reconstruction"
        )
    dataset = load_continuous_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    from bijux_phylogenetics.comparative.evolutionary_modes.tree_transforms import (
        transform_tree_for_evolutionary_mode,
    )

    transformed_tree = transform_tree_for_evolutionary_mode(
        dataset.tree,
        mode="early-burst",
        parameter_value=rate_change,
    )
    reconstruction = _reconstruct_continuous_from_dataset(
        dataset,
        working_tree=transformed_tree,
        model="brownian",
        estimator="fast-anc",
        alpha=1.0,
        brownian_fit_diagnostics=None,
        optimizer_diagnostics=None,
        anc_ml_profile_fit=None,
    )
    return ContinuousEvolutionaryModeAncestralReport(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        mode="early-burst",
        parameter_name="rate_change",
        parameter_value=rate_change,
        transformed_tree_newick=reconstruction.analysis_tree_newick,
        reconstruction=reconstruction,
    )
