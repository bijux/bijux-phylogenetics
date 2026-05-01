from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.core.metadata import load_taxon_table, write_taxon_rows
from bijux_phylogenetics.core.pruning import prune_tree_to_requested_taxa
from bijux_phylogenetics.core.traits import validate_traits_table
from bijux_phylogenetics.core.tree import PhyloTree, TreeNode
from bijux_phylogenetics.errors import AncestralReconstructionError
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree


@dataclass(slots=True)
class AncestralContinuousDataset:
    """Continuous-trait dataset normalized for ancestral reconstruction."""

    tree_path: Path
    traits_path: Path
    taxon_column: str
    trait: str
    tree: PhyloTree
    taxa: list[str]
    values_by_taxon: dict[str, float]
    dropped_missing_taxa: list[str]
    dropped_non_numeric_taxa: list[str]
    warnings: list[str]


@dataclass(slots=True)
class AncestralDiscreteDataset:
    """Discrete-trait dataset normalized for ancestral reconstruction."""

    tree_path: Path
    traits_path: Path
    taxon_column: str
    trait: str
    tree: PhyloTree
    taxa: list[str]
    states_by_taxon: dict[str, str]
    observed_states: list[str]
    state_counts: dict[str, int]
    sparse_states: list[str]
    dropped_missing_taxa: list[str]
    warnings: list[str]


def node_descendant_taxa(node: TreeNode) -> list[str]:
    """Return the sorted descendant-tip labels for a node."""
    if node.is_leaf():
        return [node.name] if node.name is not None else []
    taxa: list[str] = []
    for child in node.children:
        taxa.extend(node_descendant_taxa(child))
    return sorted(taxa)


def node_signature(node: TreeNode) -> str:
    """Return a deterministic internal-node signature."""
    taxa = node_descendant_taxa(node)
    if taxa:
        return "|".join(taxa)
    return node.name or "<unnamed>"


def stable_value(value: float) -> float:
    """Round a float into the repository's deterministic reporting style."""
    return float(format(round(value, 15), ".15g"))


def load_continuous_dataset(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    warning_taxon_threshold: int = 4,
) -> AncestralContinuousDataset:
    """Load a rooted-tree continuous-trait dataset over overlapping usable taxa."""
    tree = load_tree(tree_path)
    if len(tree.root.children) != 2:
        raise AncestralReconstructionError(
            "ancestral-state reconstruction requires a rooted tree"
        )
    if any(
        node.branch_length is None
        for node in tree.iter_nodes()
        if node is not tree.root
    ):
        raise AncestralReconstructionError(
            "continuous ancestral reconstruction requires complete branch lengths"
        )

    table = load_taxon_table(traits_path, taxon_column=taxon_column)
    if trait not in table.columns:
        raise AncestralReconstructionError(
            f"trait table does not contain column '{trait}'"
        )
    trait_report = validate_traits_table(traits_path, taxon_column=taxon_column)
    column_kind = next(
        (column.kind for column in trait_report.trait_columns if column.name == trait),
        None,
    )
    if column_kind != "numeric":
        raise AncestralReconstructionError(
            f"trait column '{trait}' must be numeric for continuous ancestral reconstruction"
        )

    rows_by_taxon = {row[table.taxon_column]: row for row in table.rows}
    kept_taxa: list[str] = []
    dropped_missing_taxa: list[str] = []
    dropped_non_numeric_taxa: list[str] = []
    values_by_taxon: dict[str, float] = {}
    for taxon in tree.tip_names:
        row = rows_by_taxon.get(taxon)
        if row is None or not row[trait]:
            dropped_missing_taxa.append(taxon)
            continue
        try:
            values_by_taxon[taxon] = float(row[trait])
        except ValueError:
            dropped_non_numeric_taxa.append(taxon)
            continue
        kept_taxa.append(taxon)
    if len(kept_taxa) < 2:
        raise AncestralReconstructionError(
            "continuous ancestral reconstruction requires at least two taxa with usable numeric trait values"
        )

    pruned_tree, _ = prune_tree_to_requested_taxa(tree_path, kept_taxa)
    warnings: list[str] = []
    if len(kept_taxa) < warning_taxon_threshold:
        warnings.append(
            f"continuous trait reconstruction is using only {len(kept_taxa)} taxa; results may be unstable"
        )
    if dropped_missing_taxa:
        warnings.append(
            "one or more taxa were excluded because the continuous trait value was missing"
        )
    if dropped_non_numeric_taxa:
        warnings.append(
            "one or more taxa were excluded because the continuous trait value was not numeric"
        )
    return AncestralContinuousDataset(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=table.taxon_column,
        trait=trait,
        tree=pruned_tree,
        taxa=pruned_tree.tip_names,
        values_by_taxon=values_by_taxon,
        dropped_missing_taxa=sorted(dropped_missing_taxa),
        dropped_non_numeric_taxa=sorted(dropped_non_numeric_taxa),
        warnings=warnings,
    )


def load_discrete_dataset(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
) -> AncestralDiscreteDataset:
    """Load a rooted-tree discrete-trait dataset over overlapping usable taxa."""
    tree = load_tree(tree_path)
    if len(tree.root.children) != 2:
        raise AncestralReconstructionError(
            "ancestral-state reconstruction requires a rooted tree"
        )

    table = load_taxon_table(traits_path, taxon_column=taxon_column)
    if trait not in table.columns:
        raise AncestralReconstructionError(
            f"trait table does not contain column '{trait}'"
        )

    rows_by_taxon = {row[table.taxon_column]: row for row in table.rows}
    kept_taxa: list[str] = []
    dropped_missing_taxa: list[str] = []
    states_by_taxon: dict[str, str] = {}
    for taxon in tree.tip_names:
        row = rows_by_taxon.get(taxon)
        if row is None or not row[trait]:
            dropped_missing_taxa.append(taxon)
            continue
        state = row[trait].strip()
        states_by_taxon[taxon] = state
        kept_taxa.append(taxon)
    observed_states = sorted(set(states_by_taxon.values()))
    if len(kept_taxa) < 2:
        raise AncestralReconstructionError(
            "discrete ancestral reconstruction requires at least two taxa with observed states"
        )
    if len(observed_states) < 2:
        raise AncestralReconstructionError(
            "discrete ancestral reconstruction requires at least two observed states"
        )

    pruned_tree, _ = prune_tree_to_requested_taxa(tree_path, kept_taxa)
    state_counts = {
        state: sum(1 for value in states_by_taxon.values() if value == state)
        for state in observed_states
    }
    sparse_states = sorted(state for state, count in state_counts.items() if count < 2)
    warnings: list[str] = []
    if dropped_missing_taxa:
        warnings.append(
            "one or more taxa were excluded because the discrete trait state was missing"
        )
    if sparse_states:
        warnings.append(
            "one or more discrete states are represented by fewer than two taxa and should be interpreted cautiously"
        )
    return AncestralDiscreteDataset(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=table.taxon_column,
        trait=trait,
        tree=pruned_tree,
        taxa=pruned_tree.tip_names,
        states_by_taxon=states_by_taxon,
        observed_states=observed_states,
        state_counts=state_counts,
        sparse_states=sparse_states,
        dropped_missing_taxa=sorted(dropped_missing_taxa),
        warnings=warnings,
    )


def write_ancestral_rows(
    path: Path, *, columns: list[str], rows: list[dict[str, str]]
) -> Path:
    """Write a deterministic ancestral-state table."""
    return write_taxon_rows(path, columns=columns, rows=rows)


def reconstruction_manifest(
    *,
    report_kind: str,
    title: str,
    tree_path: Path,
    traits_path: Path,
    trait: str,
    model: str,
    rendered_tree: str | None = None,
) -> dict[str, object]:
    """Build a minimal machine manifest for ancestral-state outputs."""
    payload: dict[str, object] = {
        "report_kind": report_kind,
        "title": title,
        "tree_path": str(tree_path),
        "traits_path": str(traits_path),
        "trait": trait,
        "model": model,
    }
    if rendered_tree is not None:
        payload["rendered_tree"] = rendered_tree
    return payload


def dump_pruned_tree(tree: PhyloTree) -> str:
    """Expose the pruned analysis tree in canonical Newick."""
    return dumps_newick(tree)
