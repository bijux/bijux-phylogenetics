from __future__ import annotations

import math
from pathlib import Path
import tempfile

from bijux_phylogenetics.ancestral.common import stable_value
from bijux_phylogenetics.ancestral.discrete import (
    reconstruct_discrete_ancestral_states,
)
from bijux_phylogenetics.ancestral.discrete.policy import (
    resolve_discrete_model_name,
)
from bijux_phylogenetics.ancestral.discrete.state_resolution import (
    resolve_clade_consensus_state,
)
from bijux_phylogenetics.datasets.study_inputs import load_taxon_table
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.runtime.errors import AncestralReconstructionError

from .models import (
    AncestralTreeSetTreeRow,
    DiscreteAncestralTreeSetCladeSummaryRow,
    DiscreteAncestralTreeSetNodeRow,
    DiscreteAncestralTreeSetReport,
)
from .preparation import (
    load_tree_set_trees,
    prepare_analysis_tree_set,
    shared_taxa,
    validate_burnin_fraction,
)


def summarize_discrete_ancestral_tree_set(
    tree_set_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    model: str = "fitch",
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
    burnin_fraction: float = 0.0,
) -> DiscreteAncestralTreeSetReport:
    """Run discrete ancestral reconstruction across a retained posterior or bootstrap tree set."""
    if model == "meristic":
        resolve_discrete_model_name(model)
    if model not in {"fitch", "equal-rates", "symmetric", "all-rates-different"}:
        raise AncestralReconstructionError(
            f"unsupported discrete ancestral tree-set model: {model}"
        )
    if model == "fitch" and state_ordering != "unordered":
        raise AncestralReconstructionError(
            "ordered discrete ancestral tree-set reconstruction requires a likelihood model"
        )
    validate_burnin_fraction(burnin_fraction)
    _source_format, trees = load_tree_set_trees(tree_set_path)
    total_tree_count = len(trees)
    burnin_tree_count = math.floor(total_tree_count * burnin_fraction)
    kept_tree_entries = [
        (source_tree_index, tree)
        for source_tree_index, tree in enumerate(trees, start=1)
    ][burnin_tree_count:]
    if not kept_tree_entries:
        raise AncestralReconstructionError(
            "ancestral tree-set analysis retains no trees after burn-in removal"
        )
    kept_trees = [tree for _, tree in kept_tree_entries]
    shared_tree_taxa = sorted(shared_taxa(kept_trees))
    warnings: list[str] = []
    if any(set(tree.tip_names) != set(shared_tree_taxa) for tree in kept_trees):
        warnings.append(
            "retained trees do not share identical tip sets and were reduced to their shared taxa"
        )
    (
        analysis_trees,
        topology_summary,
        analysis_taxa,
        exclusions,
        dataset_warnings,
        resolved_taxon_column,
    ) = prepare_analysis_tree_set(
        traits_path=traits_path,
        taxon_column=taxon_column,
        trait=trait,
        kept_tree_entries=kept_tree_entries,
        shared_tree_taxa=shared_tree_taxa,
        dataset_kind="discrete",
    )
    warnings.extend(dataset_warnings)
    taxon_table = load_taxon_table(traits_path, taxon_column=resolved_taxon_column)
    analysis_taxon_set = set(analysis_taxa)
    observed_states_by_taxon = {
        row[resolved_taxon_column]: row[trait].strip()
        for row in taxon_table.rows
        if row[resolved_taxon_column] in analysis_taxon_set and row[trait].strip()
    }
    tree_rows = [
        AncestralTreeSetTreeRow(
            source_tree_index=source_tree_index,
            post_burnin_index=record.index,
            rooted_topology_id=record.rooted_topology_id,
            unrooted_topology_id=record.unrooted_topology_id,
            internal_clade_count=max(len(analysis_tree.tip_names) - 1, 0),
        )
        for (source_tree_index, analysis_tree), record in zip(
            analysis_trees, topology_summary.records, strict=True
        )
    ]
    node_rows: list[DiscreteAncestralTreeSetNodeRow] = []
    with tempfile.TemporaryDirectory(
        prefix="bijux-phylogenetics-ancestral-tree-set-discrete-"
    ) as tmp_dir:
        current_tree_path = Path(tmp_dir) / "ancestral-tree-set-current-tree.nwk"
        for tree_row, (_source_tree_index, analysis_tree) in zip(
            tree_rows,
            analysis_trees,
            strict=True,
        ):
            current_tree_path.write_text(
                dumps_newick(analysis_tree) + "\n",
                encoding="utf-8",
            )
            report = reconstruct_discrete_ancestral_states(
                current_tree_path,
                traits_path,
                trait=trait,
                taxon_column=taxon_column,
                model=model,
                state_ordering=state_ordering,
                ordered_states=ordered_states,
            )
            for estimate in report.estimates:
                if estimate.is_tip:
                    continue
                node_rows.append(
                    DiscreteAncestralTreeSetNodeRow(
                        source_tree_index=tree_row.source_tree_index,
                        post_burnin_index=tree_row.post_burnin_index,
                        rooted_topology_id=tree_row.rooted_topology_id,
                        unrooted_topology_id=tree_row.unrooted_topology_id,
                        clade_id=estimate.node,
                        clade_taxa=estimate.descendant_taxa,
                        most_likely_state=estimate.most_likely_state,
                        state_set=estimate.state_set,
                        confidence=estimate.confidence,
                        ambiguous=estimate.ambiguous,
                        unstable=estimate.unstable,
                    )
                )
    clade_summaries = summarize_discrete_clades(
        node_rows,
        kept_tree_count=len(tree_rows),
        observed_states_by_taxon=observed_states_by_taxon,
    )
    if any(row.tree_presence_fraction < 1.0 for row in clade_summaries):
        warnings.append(
            "one or more comparable ancestral clades are absent from some retained trees"
        )
    if any(row.stability_class != "stable" for row in clade_summaries):
        warnings.append(
            "one or more discrete ancestral clades change state or support profile across retained trees"
        )
    return DiscreteAncestralTreeSetReport(
        tree_set_path=tree_set_path,
        traits_path=traits_path,
        trait=trait,
        taxon_column=resolved_taxon_column,
        model=model,
        state_ordering=state_ordering,
        ordered_states=list(ordered_states or []),
        burnin_fraction=burnin_fraction,
        total_tree_count=total_tree_count,
        burnin_tree_count=burnin_tree_count,
        kept_tree_count=len(tree_rows),
        shared_tree_taxa=shared_tree_taxa,
        analysis_taxa=analysis_taxa,
        rooted_topology_count=topology_summary.rooted_topology_count,
        unrooted_topology_count=topology_summary.unrooted_topology_count,
        tree_rows=tree_rows,
        node_rows=node_rows,
        clade_summaries=clade_summaries,
        exclusions=exclusions,
        warnings=warnings,
    )


def summarize_discrete_clades(
    rows: list[DiscreteAncestralTreeSetNodeRow],
    *,
    kept_tree_count: int,
    observed_states_by_taxon: dict[str, str],
) -> list[DiscreteAncestralTreeSetCladeSummaryRow]:
    """Summarize comparable discrete clades across retained trees."""
    grouped: dict[str, list[DiscreteAncestralTreeSetNodeRow]] = {}
    for row in rows:
        grouped.setdefault(row.clade_id, []).append(row)
    summaries: list[DiscreteAncestralTreeSetCladeSummaryRow] = []
    for clade_id, clade_rows in sorted(grouped.items()):
        presence_fraction = stable_value(len(clade_rows) / kept_tree_count)
        state_counts: dict[str, int] = {}
        for row in clade_rows:
            resolved_state = resolve_clade_consensus_state(
                clade_taxa=row.clade_taxa,
                candidate_states=row.state_set,
                observed_states_by_taxon=observed_states_by_taxon,
                fallback_state=row.most_likely_state,
            )
            state_counts[resolved_state] = state_counts.get(resolved_state, 0) + 1
        dominant_state = max(
            sorted(state_counts),
            key=lambda state: (state_counts[state], state),
        )
        dominant_state_tree_count = state_counts[dominant_state]
        dominant_state_fraction = stable_value(
            dominant_state_tree_count / len(clade_rows)
        )
        ambiguous_tree_count = sum(row.ambiguous for row in clade_rows)
        unstable_tree_count = sum(row.unstable for row in clade_rows)
        ambiguous_tree_fraction = stable_value(ambiguous_tree_count / len(clade_rows))
        unstable_tree_fraction = stable_value(unstable_tree_count / len(clade_rows))
        instability_score = stable_value(
            (1.0 - presence_fraction)
            + (1.0 - dominant_state_fraction)
            + ambiguous_tree_fraction
            + unstable_tree_fraction
        )
        if presence_fraction < 1.0:
            stability_class = "topology_sensitive"
        elif len(state_counts) > 1:
            stability_class = "state_conflict"
        elif ambiguous_tree_count > 0 or unstable_tree_count > 0:
            stability_class = "low_confidence"
        else:
            stability_class = "stable"
        summaries.append(
            DiscreteAncestralTreeSetCladeSummaryRow(
                clade_id=clade_id,
                clade_taxa=clade_rows[0].clade_taxa,
                tree_presence_count=len(clade_rows),
                tree_presence_fraction=presence_fraction,
                unique_state_count=len(state_counts),
                dominant_state=dominant_state,
                dominant_state_tree_count=dominant_state_tree_count,
                dominant_state_fraction=dominant_state_fraction,
                ambiguous_tree_count=ambiguous_tree_count,
                ambiguous_tree_fraction=ambiguous_tree_fraction,
                unstable_tree_count=unstable_tree_count,
                unstable_tree_fraction=unstable_tree_fraction,
                state_distribution=dict(sorted(state_counts.items())),
                instability_score=instability_score,
                stability_class=stability_class,
            )
        )
    return summaries
