from __future__ import annotations

from pathlib import Path
import tempfile

from bijux_phylogenetics.ancestral.discrete import reconstruct_discrete_ancestral_states
from bijux_phylogenetics.ancestral.tree_set.preparation import (
    load_tree_set_trees,
    prepare_analysis_tree_set,
    shared_taxa,
    validate_burnin_fraction,
)
from bijux_phylogenetics.datasets.study_inputs import load_taxon_table
from bijux_phylogenetics.io.newick import dumps_newick

from .aggregation import (
    _summarize_transition_rows,
    _summarize_transition_rows_across_trees,
)
from .branch_analysis import _build_transition_branch_rows
from .contracts import (
    AncestralTransitionCountRow,
    AncestralTransitionExclusion,
    AncestralTransitionSummary,
    AncestralTransitionTreeRow,
    AncestralTransitionTreeSetBranchRow,
    AncestralTransitionTreeSetReport,
    AncestralTransitionTreeSetSummary,
)


def summarize_ancestral_transition_tree_set(
    tree_set_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    model: str = "fitch",
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
    burnin_fraction: float = 0.0,
) -> AncestralTransitionTreeSetReport:
    """Count inferred discrete ancestral transitions across a retained tree set."""
    validate_burnin_fraction(burnin_fraction)
    _source_format, trees = load_tree_set_trees(tree_set_path)
    total_tree_count = len(trees)
    burnin_tree_count = int(total_tree_count * burnin_fraction)
    kept_tree_entries = [
        (source_tree_index, tree)
        for source_tree_index, tree in enumerate(trees, start=1)
    ][burnin_tree_count:]
    if not kept_tree_entries:
        raise ValueError(
            "ancestral transition tree-set analysis retains no trees after burn-in removal"
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
        raw_exclusions,
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
    exclusions = [
        AncestralTransitionExclusion(taxon=row.taxon, reason=row.reason)
        for row in raw_exclusions
    ]
    tree_rows: list[AncestralTransitionTreeRow] = []
    branch_rows: list[AncestralTransitionTreeSetBranchRow] = []
    per_tree_transition_counts: dict[str, list[AncestralTransitionCountRow]] = {}
    with tempfile.TemporaryDirectory(
        prefix="bijux-phylogenetics-ancestral-transition-tree-set-"
    ) as tmp_dir:
        current_tree_path = Path(tmp_dir) / "ancestral-transition-current-tree.nwk"
        for (source_tree_index, analysis_tree), topology_record in zip(
            analysis_trees,
            topology_summary.records,
            strict=True,
        ):
            current_tree_path.write_text(
                dumps_newick(analysis_tree) + "\n",
                encoding="utf-8",
            )
            reconstruction = reconstruct_discrete_ancestral_states(
                current_tree_path,
                traits_path,
                trait=trait,
                taxon_column=taxon_column,
                model=model,
                state_ordering=state_ordering,
                ordered_states=ordered_states,
            )
            resolved_branch_rows = _build_transition_branch_rows(
                reconstruction,
                observed_states_by_taxon=observed_states_by_taxon,
            )
            resolved_transition_rows = _summarize_transition_rows(resolved_branch_rows)
            changed_rows = [row for row in resolved_branch_rows if row.changed]
            summary = AncestralTransitionSummary(
                trait=reconstruction.trait,
                taxon_column=reconstruction.taxon_column,
                model=reconstruction.model,
                state_ordering=reconstruction.state_ordering,
                analyzed_taxon_count=reconstruction.taxon_count,
                excluded_taxon_count=len(reconstruction.dropped_missing_taxa),
                total_branch_count=len(resolved_branch_rows),
                changed_branch_count=len(changed_rows),
                certain_change_count=sum(
                    row.certainty_class == "certain_change" for row in changed_rows
                ),
                uncertain_change_count=sum(
                    row.certainty_class == "uncertain_change" for row in changed_rows
                ),
                unique_transition_count=len(resolved_transition_rows),
                warning_count=len(reconstruction.warnings),
            )
            tree_rows.append(
                AncestralTransitionTreeRow(
                    source_tree_index=source_tree_index,
                    post_burnin_index=topology_record.index,
                    rooted_topology_id=topology_record.rooted_topology_id,
                    unrooted_topology_id=topology_record.unrooted_topology_id,
                    total_branch_count=summary.total_branch_count,
                    changed_branch_count=summary.changed_branch_count,
                    certain_change_count=summary.certain_change_count,
                    uncertain_change_count=summary.uncertain_change_count,
                )
            )
            for row in resolved_branch_rows:
                branch_rows.append(
                    AncestralTransitionTreeSetBranchRow(
                        source_tree_index=source_tree_index,
                        post_burnin_index=topology_record.index,
                        rooted_topology_id=topology_record.rooted_topology_id,
                        unrooted_topology_id=topology_record.unrooted_topology_id,
                        parent_node=row.parent_node,
                        child_node=row.child_node,
                        child_descendant_taxa=row.child_descendant_taxa,
                        branch_length=row.branch_length,
                        parent_most_likely_state=row.parent_most_likely_state,
                        child_most_likely_state=row.child_most_likely_state,
                        parent_state_set=row.parent_state_set,
                        child_state_set=row.child_state_set,
                        overlapping_states=row.overlapping_states,
                        changed=row.changed,
                        transition=row.transition,
                        certainty_class=row.certainty_class,
                    )
                )
            per_tree_transition_counts[str(source_tree_index)] = (
                resolved_transition_rows
            )
    transition_rows = _summarize_transition_rows_across_trees(
        tree_rows=tree_rows,
        per_tree_transition_counts=per_tree_transition_counts,
    )
    if any(row.stability_class == "topology_sensitive" for row in transition_rows):
        warnings.append(
            "one or more inferred transition pairs are absent from some retained trees"
        )
    if any(row.stability_class == "uncertainty_sensitive" for row in transition_rows):
        warnings.append(
            "one or more inferred transition pairs depend on uncertain branchwise ancestral changes"
        )
    return AncestralTransitionTreeSetReport(
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
        branch_rows=branch_rows,
        transition_rows=transition_rows,
        exclusions=exclusions,
        warnings=warnings,
    )


def summarize_ancestral_transition_tree_set_report(
    report: AncestralTransitionTreeSetReport,
) -> AncestralTransitionTreeSetSummary:
    """Summarize the main review facts for one ancestral transition tree-set report."""
    return AncestralTransitionTreeSetSummary(
        trait=report.trait,
        taxon_column=report.taxon_column,
        model=report.model,
        state_ordering=report.state_ordering,
        total_tree_count=report.total_tree_count,
        burnin_tree_count=report.burnin_tree_count,
        kept_tree_count=report.kept_tree_count,
        shared_tree_taxon_count=len(report.shared_tree_taxa),
        analysis_taxon_count=len(report.analysis_taxa),
        rooted_topology_count=report.rooted_topology_count,
        unrooted_topology_count=report.unrooted_topology_count,
        transition_pair_count=len(report.transition_rows),
        topology_sensitive_transition_pair_count=sum(
            row.stability_class == "topology_sensitive"
            for row in report.transition_rows
        ),
        uncertainty_sensitive_transition_pair_count=sum(
            row.stability_class == "uncertainty_sensitive"
            for row in report.transition_rows
        ),
        warning_count=len(report.warnings),
    )
