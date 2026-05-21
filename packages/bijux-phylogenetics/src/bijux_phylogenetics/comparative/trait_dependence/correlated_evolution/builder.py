from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.ancestral.discrete.policy import (
    resolve_discrete_model_name as _resolve_discrete_model_name,
)
from bijux_phylogenetics.datasets.study_inputs import load_taxon_table
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.runtime.errors import ComparativeMethodError

from .binary_coupling import _summarize_binary_trait_coupling
from .continuous_coupling import _summarize_continuous_trait_coupling
from .contracts import CorrelatedTraitEvolutionReport
from .preparation import _prepare_shared_trait_rows


def summarize_correlated_trait_evolution(
    tree_path: Path,
    traits_path: Path,
    *,
    left_trait: str,
    right_trait: str,
    taxon_column: str | None = None,
    analysis_kind: str = "auto",
    binary_model: str = "all-rates-different",
) -> CorrelatedTraitEvolutionReport:
    """Summarize coupling between two traits across one phylogenetic tree."""
    if left_trait == right_trait:
        raise ComparativeMethodError(
            "correlated trait evolution requires two distinct trait columns"
        )
    if analysis_kind not in {"auto", "continuous", "binary"}:
        raise ComparativeMethodError(
            "analysis_kind must be one of: auto, continuous, binary"
        )
    if binary_model == "meristic":
        _resolve_discrete_model_name(binary_model)
    if binary_model not in {"equal-rates", "symmetric", "all-rates-different"}:
        raise ComparativeMethodError(
            "binary correlated-trait analysis requires a supported discrete model"
        )

    tree = load_tree(tree_path)
    table = load_taxon_table(traits_path, taxon_column=taxon_column)
    for trait in (left_trait, right_trait):
        if trait not in table.columns:
            raise ComparativeMethodError(
                f"trait table does not contain required column '{trait}'"
            )

    prepared = _prepare_shared_trait_rows(
        tree=tree,
        table=table,
        left_trait=left_trait,
        right_trait=right_trait,
        analysis_kind=analysis_kind,
    )
    if prepared.analysis_kind == "continuous-brownian-contrasts":
        return _summarize_continuous_trait_coupling(
            tree_path=tree_path,
            traits_path=traits_path,
            tree_taxon_count=tree.tip_count,
            taxon_column=table.taxon_column,
            left_trait=left_trait,
            right_trait=right_trait,
            prepared=prepared,
        )
    return _summarize_binary_trait_coupling(
        tree_path=tree_path,
        traits_path=traits_path,
        tree_taxon_count=tree.tip_count,
        taxon_column=table.taxon_column,
        left_trait=left_trait,
        right_trait=right_trait,
        prepared=prepared,
        binary_model=binary_model,
    )
