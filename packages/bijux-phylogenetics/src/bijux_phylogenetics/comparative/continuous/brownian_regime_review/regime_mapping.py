from __future__ import annotations

from bijux_phylogenetics.comparative.common import ComparativeReadinessReport
from bijux_phylogenetics.comparative.traits.regime_mapping import (
    build_branch_identity_lookup,
    resolve_branch_regime_id_column,
)
from bijux_phylogenetics.datasets.study_inputs import load_taxon_table
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import ComparativeMethodError

from .contracts import BrownianRegimeBranchRow, BrownianRegimeExclusion


def load_branch_regime_rows(
    regime_map_path,
    *,
    tree: PhyloTree,
    analyzed_taxa: list[str],
    branch_id_column: str | None,
    regime_column: str,
) -> tuple[list[BrownianRegimeBranchRow], str]:
    resolved_branch_id_column = resolve_branch_regime_id_column(
        regime_map_path,
        requested=branch_id_column,
    )
    table = load_taxon_table(
        regime_map_path,
        taxon_column=resolved_branch_id_column,
    )
    if regime_column not in table.columns:
        raise ComparativeMethodError(
            f"regime map does not contain column '{regime_column}'"
        )
    branch_lookup = build_branch_identity_lookup(tree, analyzed_taxa=analyzed_taxa)
    mapped_branch_ids = {row[table.taxon_column] for row in table.rows}
    expected_branch_ids = set(branch_lookup)
    missing = sorted(expected_branch_ids - mapped_branch_ids)
    if missing:
        raise ComparativeMethodError(
            "regime map is missing one or more non-root branches: " + ", ".join(missing)
        )
    extra = sorted(mapped_branch_ids - expected_branch_ids)
    if extra:
        raise ComparativeMethodError(
            "regime map contains branches absent from the tree: " + ", ".join(extra)
        )
    rows: list[BrownianRegimeBranchRow] = []
    for row in table.rows:
        branch = branch_lookup[row[table.taxon_column]]
        regime = row[regime_column]
        if not regime:
            raise ComparativeMethodError(
                f"regime map branch '{branch.branch_id}' has an empty '{regime_column}' value"
            )
        rows.append(
            BrownianRegimeBranchRow(
                branch_id=branch.branch_id,
                regime=regime,
                branch_length=branch.branch_length,
                descendant_taxa=branch.descendant_taxa,
                analyzed_descendant_taxa=branch.analyzed_descendant_taxa,
                contributes_to_analysis=branch.contributes_to_analysis,
            )
        )
    contributing_regimes = {
        row.regime
        for row in rows
        if row.contributes_to_analysis and row.branch_length > 0.0
    }
    if len(contributing_regimes) < 2:
        raise ComparativeMethodError(
            "regime map must expose at least two regimes that contribute analyzed branch length"
        )
    return sorted(rows, key=lambda row: row.branch_id), resolved_branch_id_column


def build_regime_covariance_components(
    taxa: list[str],
    branch_rows: list[BrownianRegimeBranchRow],
) -> dict[str, list[list[float]]]:
    index = {taxon: position for position, taxon in enumerate(taxa)}
    components: dict[str, list[list[float]]] = {}
    for branch in branch_rows:
        if not branch.contributes_to_analysis or branch.branch_length <= 0.0:
            continue
        matrix = components.setdefault(
            branch.regime,
            [[0.0] * len(taxa) for _ in range(len(taxa))],
        )
        for left_taxon in branch.analyzed_descendant_taxa:
            for right_taxon in branch.analyzed_descendant_taxa:
                matrix[index[left_taxon]][index[right_taxon]] += branch.branch_length
    return components


def build_excluded_taxa(
    readiness: ComparativeReadinessReport,
) -> list[BrownianRegimeExclusion]:
    rows: list[BrownianRegimeExclusion] = []
    rows.extend(
        BrownianRegimeExclusion(taxon=taxon, reason="missing_from_trait_table")
        for taxon in readiness.missing_from_traits
    )
    rows.extend(
        BrownianRegimeExclusion(taxon=taxon, reason="missing_trait_value")
        for taxon in readiness.pruned_missing_value_taxa
    )
    rows.extend(
        BrownianRegimeExclusion(taxon=taxon, reason="non_numeric_trait_value")
        for taxon in readiness.pruned_non_numeric_taxa
    )
    rows.extend(
        BrownianRegimeExclusion(taxon=taxon, reason="absent_from_tree")
        for taxon in readiness.extra_trait_taxa
    )
    return rows
