from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.ancestral.discrete import (
    reconstruct_discrete_ancestral_states,
)
from bijux_phylogenetics.runtime.errors import AncestralReconstructionError

from .analysis import (
    build_branch_rows,
    build_clade_rows,
    build_count_rows,
    build_exclusion_rows,
    build_node_rows,
    build_rate_rows,
    build_summary,
)
from .contracts import NicheTransitionReport
from .shared import resolve_internal_model


def summarize_niche_transitions(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    model: str = "er",
) -> NicheTransitionReport:
    """Fit one ecological niche transition model and summarize clade shifts."""
    internal_model = resolve_internal_model(model)
    reconstruction = reconstruct_discrete_ancestral_states(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=internal_model,
    )
    if (
        reconstruction.log_likelihood is None
        or reconstruction.parameter_count is None
        or reconstruction.aic is None
    ):
        raise AncestralReconstructionError(
            "ecological niche transition analysis requires a likelihood discrete ancestral model"
        )
    node_rows = build_node_rows(reconstruction)
    rate_rows = build_rate_rows(reconstruction)
    branch_rows = build_branch_rows(reconstruction)
    count_rows = build_count_rows(branch_rows)
    clade_rows = build_clade_rows(reconstruction, branch_rows)
    exclusion_rows = build_exclusion_rows(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    summary = build_summary(
        reconstruction=reconstruction,
        model=model,
        node_rows=node_rows,
        rate_rows=rate_rows,
        branch_rows=branch_rows,
        clade_rows=clade_rows,
        exclusion_rows=exclusion_rows,
        count_rows=count_rows,
    )
    return NicheTransitionReport(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        taxon_column=reconstruction.taxon_column,
        model=model,
        internal_model=reconstruction.model,
        summary=summary,
        node_rows=node_rows,
        rate_rows=rate_rows,
        branch_rows=branch_rows,
        count_rows=count_rows,
        clade_rows=clade_rows,
        exclusion_rows=exclusion_rows,
        warnings=list(reconstruction.warnings),
    )
