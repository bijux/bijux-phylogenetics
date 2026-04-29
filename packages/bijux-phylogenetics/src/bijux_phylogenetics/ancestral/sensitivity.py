from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import tempfile

from bijux_phylogenetics.ancestral.continuous import reconstruct_continuous_ancestral_states
from bijux_phylogenetics.ancestral.discrete import reconstruct_discrete_ancestral_states
from bijux_phylogenetics.ancestral.service import (
    compare_continuous_ancestral_models,
    compare_continuous_ancestral_trees,
    compare_discrete_ancestral_models,
    compare_discrete_ancestral_trees,
)
from bijux_phylogenetics.core.metadata import load_taxon_table, write_taxon_rows
from bijux_phylogenetics.core.pruning import prune_tree_to_requested_taxa
from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.io.trees import load_tree


@dataclass(slots=True)
class AncestralSensitivitySummary:
    label: str
    changed_node_count: int
    notes: list[str]


@dataclass(slots=True)
class AncestralSensitivityReport:
    tree_path: Path
    traits_path: Path
    trait: str
    reconstruction_kind: str
    baseline_model: str
    baseline_node_count: int
    model_sensitivity: AncestralSensitivitySummary | None
    tree_sensitivity: AncestralSensitivitySummary | None
    pruning_sensitivity: AncestralSensitivitySummary | None
    trait_coding_sensitivity: AncestralSensitivitySummary | None


def build_ancestral_sensitivity_report(
    *,
    tree_path: Path,
    traits_path: Path,
    trait: str,
    reconstruction_kind: str,
    model: str,
    taxon_column: str | None = None,
    alpha: float = 1.0,
    compare_tree_path: Path | None = None,
    compare_model: str | None = None,
    drop_taxa: list[str] | None = None,
    coding_map: dict[str, str] | None = None,
) -> AncestralSensitivityReport:
    """Summarize sensitivity of ancestral reconstruction to model, tree, pruning, and coding choices."""
    baseline_node_count = _baseline_node_count(
        tree_path,
        traits_path,
        trait=trait,
        reconstruction_kind=reconstruction_kind,
        model=model,
        taxon_column=taxon_column,
        alpha=alpha,
    )
    model_sensitivity = _model_sensitivity_summary(
        tree_path,
        traits_path,
        trait=trait,
        reconstruction_kind=reconstruction_kind,
        model=model,
        taxon_column=taxon_column,
        alpha=alpha,
        compare_model=compare_model,
    )
    tree_sensitivity = _tree_sensitivity_summary(
        tree_path,
        traits_path,
        trait=trait,
        reconstruction_kind=reconstruction_kind,
        model=model,
        taxon_column=taxon_column,
        alpha=alpha,
        compare_tree_path=compare_tree_path,
    )
    pruning_sensitivity = _pruning_sensitivity_summary(
        tree_path,
        traits_path,
        trait=trait,
        reconstruction_kind=reconstruction_kind,
        model=model,
        taxon_column=taxon_column,
        alpha=alpha,
        drop_taxa=drop_taxa,
    )
    trait_coding_sensitivity = _trait_coding_sensitivity_summary(
        tree_path,
        traits_path,
        trait=trait,
        reconstruction_kind=reconstruction_kind,
        model=model,
        taxon_column=taxon_column,
        coding_map=coding_map,
    )
    return AncestralSensitivityReport(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        reconstruction_kind=reconstruction_kind,
        baseline_model=model,
        baseline_node_count=baseline_node_count,
        model_sensitivity=model_sensitivity,
        tree_sensitivity=tree_sensitivity,
        pruning_sensitivity=pruning_sensitivity,
        trait_coding_sensitivity=trait_coding_sensitivity,
    )


def _baseline_node_count(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    reconstruction_kind: str,
    model: str,
    taxon_column: str | None,
    alpha: float,
) -> int:
    if reconstruction_kind == "continuous":
        report = reconstruct_continuous_ancestral_states(
            tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
            model=model,
            alpha=alpha,
        )
    else:
        report = reconstruct_discrete_ancestral_states(
            tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
            model=model,
        )
    return sum(1 for estimate in report.estimates if not estimate.is_tip)


def _model_sensitivity_summary(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    reconstruction_kind: str,
    model: str,
    taxon_column: str | None,
    alpha: float,
    compare_model: str | None,
) -> AncestralSensitivitySummary | None:
    if compare_model is None:
        return None
    if reconstruction_kind == "continuous":
        comparison = compare_continuous_ancestral_models(
            tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
            left_model=model,
            right_model=compare_model,
            left_alpha=alpha,
            right_alpha=alpha,
        )
        changed = sum(1 for row in comparison.rows if abs(row.estimate_delta) > 1e-9)
        return AncestralSensitivitySummary(
            label="model",
            changed_node_count=changed,
            notes=[f"compared continuous models {model} and {compare_model}"],
        )
    comparison = compare_discrete_ancestral_models(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        models=(model, compare_model, "all-rates-different" if compare_model != "all-rates-different" else "equal-rates"),
    ) if model != "fitch" and compare_model != "fitch" else None
    if comparison is not None:
        changed = sum(1 for row in comparison.node_differences if row.differs)
        return AncestralSensitivitySummary(
            label="model",
            changed_node_count=changed,
            notes=[f"selected discrete model {comparison.selected_model} from {', '.join(row.model for row in comparison.rows)}"],
        )
    baseline = reconstruct_discrete_ancestral_states(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=model,
    )
    alternative = reconstruct_discrete_ancestral_states(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=compare_model,
    )
    alternative_by_node = {estimate.node: estimate for estimate in alternative.estimates if not estimate.is_tip}
    changed = sum(
        1
        for estimate in baseline.estimates
        if not estimate.is_tip
        and estimate.node in alternative_by_node
        and estimate.most_likely_state != alternative_by_node[estimate.node].most_likely_state
    )
    return AncestralSensitivitySummary(
        label="model",
        changed_node_count=changed,
        notes=[f"compared discrete models {model} and {compare_model}"],
    )


def _tree_sensitivity_summary(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    reconstruction_kind: str,
    model: str,
    taxon_column: str | None,
    alpha: float,
    compare_tree_path: Path | None,
) -> AncestralSensitivitySummary | None:
    if compare_tree_path is None:
        return None
    if reconstruction_kind == "continuous":
        comparison = compare_continuous_ancestral_trees(
            tree_path,
            compare_tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
            model=model,
            alpha=alpha,
        )
        changed = sum(1 for row in comparison.rows if abs(row.estimate_delta) > 1e-9)
    else:
        comparison = compare_discrete_ancestral_trees(
            tree_path,
            compare_tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
            model=model,
        )
        changed = sum(1 for row in comparison.rows if row.differs)
    return AncestralSensitivitySummary(
        label="tree",
        changed_node_count=changed,
        notes=[f"compared baseline tree with {compare_tree_path.name}"],
    )


def _pruning_sensitivity_summary(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    reconstruction_kind: str,
    model: str,
    taxon_column: str | None,
    alpha: float,
    drop_taxa: list[str] | None,
) -> AncestralSensitivitySummary | None:
    if not drop_taxa:
        return None
    tree = load_tree(tree_path)
    kept_taxa = [taxon for taxon in tree.tip_names if taxon not in set(drop_taxa)]
    pruned_tree, pruning_report = prune_tree_to_requested_taxa(tree_path, kept_taxa)
    pruned_tree_path = Path(tempfile.mkstemp(prefix="bijux-ancestral-pruned-", suffix=".nwk")[1])
    pruned_traits_path = Path(tempfile.mkstemp(prefix="bijux-ancestral-pruned-", suffix=Path(traits_path).suffix)[1])
    try:
        write_newick(pruned_tree_path, pruned_tree)
        _write_filtered_trait_table(traits_path, pruned_traits_path, kept_taxa=pruning_report.kept_taxa, taxon_column=taxon_column)
        if reconstruction_kind == "continuous":
            comparison = compare_continuous_ancestral_trees(
                tree_path,
                pruned_tree_path,
                pruned_traits_path,
                trait=trait,
                taxon_column=taxon_column,
                model=model,
                alpha=alpha,
            )
            changed = sum(1 for row in comparison.rows if abs(row.estimate_delta) > 1e-9)
        else:
            comparison = compare_discrete_ancestral_trees(
                tree_path,
                pruned_tree_path,
                pruned_traits_path,
                trait=trait,
                taxon_column=taxon_column,
                model=model,
            )
            changed = sum(1 for row in comparison.rows if row.differs)
    finally:
        pruned_tree_path.unlink(missing_ok=True)
        pruned_traits_path.unlink(missing_ok=True)
    return AncestralSensitivitySummary(
        label="pruning",
        changed_node_count=changed,
        notes=[f"dropped taxa: {', '.join(sorted(drop_taxa))}"],
    )


def _trait_coding_sensitivity_summary(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    reconstruction_kind: str,
    model: str,
    taxon_column: str | None,
    coding_map: dict[str, str] | None,
) -> AncestralSensitivitySummary | None:
    if reconstruction_kind != "discrete" or not coding_map:
        return None
    recoded_path = Path(tempfile.mkstemp(prefix="bijux-ancestral-coding-", suffix=Path(traits_path).suffix)[1])
    try:
        _write_recoded_trait_table(
            traits_path,
            recoded_path,
            trait=trait,
            coding_map=coding_map,
            taxon_column=taxon_column,
        )
        baseline = reconstruct_discrete_ancestral_states(
            tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
            model=model,
        )
        recoded = reconstruct_discrete_ancestral_states(
            tree_path,
            recoded_path,
            trait=trait,
            taxon_column=taxon_column,
            model=model,
        )
    finally:
        recoded_path.unlink(missing_ok=True)
    recoded_by_node = {estimate.node: estimate for estimate in recoded.estimates if not estimate.is_tip}
    changed = sum(
        1
        for estimate in baseline.estimates
        if not estimate.is_tip
        and estimate.node in recoded_by_node
        and estimate.most_likely_state != recoded_by_node[estimate.node].most_likely_state
    )
    return AncestralSensitivitySummary(
        label="trait-coding",
        changed_node_count=changed,
        notes=[f"recoded {len(coding_map)} raw states into {len(set(coding_map.values()))} labels"],
    )


def _write_filtered_trait_table(
    source_path: Path,
    out_path: Path,
    *,
    kept_taxa: list[str],
    taxon_column: str | None,
) -> Path:
    table = load_taxon_table(source_path, taxon_column=taxon_column)
    keep = set(kept_taxa)
    return write_taxon_rows(
        out_path,
        columns=table.columns,
        rows=[row for row in table.rows if row[table.taxon_column] in keep],
    )


def _write_recoded_trait_table(
    source_path: Path,
    out_path: Path,
    *,
    trait: str,
    coding_map: dict[str, str],
    taxon_column: str | None,
) -> Path:
    table = load_taxon_table(source_path, taxon_column=taxon_column)
    rows = []
    for row in table.rows:
        updated = dict(row)
        raw_value = updated.get(trait, "")
        if raw_value in coding_map:
            updated[trait] = coding_map[raw_value]
        rows.append(updated)
    return write_taxon_rows(out_path, columns=table.columns, rows=rows)
