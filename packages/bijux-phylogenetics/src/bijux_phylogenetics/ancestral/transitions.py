from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import statistics
import tempfile

from bijux_phylogenetics.ancestral.common import write_ancestral_rows
from bijux_phylogenetics.ancestral.discrete import (
    DiscreteAncestralReport,
    reconstruct_discrete_ancestral_states,
)
from bijux_phylogenetics.ancestral.tree_set import (
    _load_tree_set_trees,
    _prepare_analysis_tree_set,
    _shared_taxa,
    _validate_burnin_fraction,
)
from bijux_phylogenetics.core.tree import TreeNode
from bijux_phylogenetics.io.newick import dumps_newick, loads_newick


@dataclass(frozen=True, slots=True)
class AncestralTransitionBranchRow:
    """One branchwise ancestral transition review row."""

    parent_node: str
    child_node: str
    child_descendant_taxa: list[str]
    branch_length: float | None
    parent_most_likely_state: str
    child_most_likely_state: str
    parent_state_set: list[str]
    child_state_set: list[str]
    overlapping_states: list[str]
    changed: bool
    transition: str
    certainty_class: str


@dataclass(frozen=True, slots=True)
class AncestralTransitionCountRow:
    """One transition-pair count summary for a discrete ancestral reconstruction."""

    transition: str
    source_state: str
    target_state: str
    certain_change_count: int
    uncertain_change_count: int
    total_change_count: int


@dataclass(frozen=True, slots=True)
class AncestralTransitionExclusion:
    """One excluded tip from ancestral transition counting."""

    taxon: str
    reason: str


@dataclass(frozen=True, slots=True)
class AncestralTransitionSummary:
    """Reviewer-facing summary for one ancestral transition count report."""

    trait: str
    taxon_column: str
    model: str
    state_ordering: str
    analyzed_taxon_count: int
    excluded_taxon_count: int
    total_branch_count: int
    changed_branch_count: int
    certain_change_count: int
    uncertain_change_count: int
    unique_transition_count: int
    warning_count: int


@dataclass(slots=True)
class AncestralTransitionReport:
    """Discrete ancestral transition counts for one analyzed rooted tree."""

    tree_path: Path
    traits_path: Path
    trait: str
    taxon_column: str
    model: str
    state_ordering: str
    ordered_states: list[str]
    taxon_count: int
    branch_rows: list[AncestralTransitionBranchRow]
    transition_rows: list[AncestralTransitionCountRow]
    exclusions: list[AncestralTransitionExclusion]
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class AncestralTransitionTreeRow:
    """One retained tree summary from ancestral transition counting over a tree set."""

    source_tree_index: int
    post_burnin_index: int
    rooted_topology_id: str
    unrooted_topology_id: str
    total_branch_count: int
    changed_branch_count: int
    certain_change_count: int
    uncertain_change_count: int


@dataclass(frozen=True, slots=True)
class AncestralTransitionTreeSetBranchRow:
    """One branchwise ancestral transition row from one retained tree."""

    source_tree_index: int
    post_burnin_index: int
    rooted_topology_id: str
    unrooted_topology_id: str
    parent_node: str
    child_node: str
    child_descendant_taxa: list[str]
    branch_length: float | None
    parent_most_likely_state: str
    child_most_likely_state: str
    parent_state_set: list[str]
    child_state_set: list[str]
    overlapping_states: list[str]
    changed: bool
    transition: str
    certainty_class: str


@dataclass(frozen=True, slots=True)
class AncestralTransitionTreeSetCountRow:
    """One transition-pair summary across retained trees."""

    transition: str
    source_state: str
    target_state: str
    tree_presence_count: int
    tree_presence_fraction: float
    mean_certain_change_count: float
    mean_uncertain_change_count: float
    mean_total_change_count: float
    minimum_total_change_count: int
    maximum_total_change_count: int
    lower_95_total_change_count: float
    upper_95_total_change_count: float
    stability_class: str


@dataclass(frozen=True, slots=True)
class AncestralTransitionTreeSetSummary:
    """Reviewer-facing summary for ancestral transition counts across a tree set."""

    trait: str
    taxon_column: str
    model: str
    state_ordering: str
    total_tree_count: int
    burnin_tree_count: int
    kept_tree_count: int
    shared_tree_taxon_count: int
    analysis_taxon_count: int
    rooted_topology_count: int
    unrooted_topology_count: int
    transition_pair_count: int
    topology_sensitive_transition_pair_count: int
    uncertainty_sensitive_transition_pair_count: int
    warning_count: int


@dataclass(slots=True)
class AncestralTransitionTreeSetReport:
    """Discrete ancestral transition counts aggregated across a posterior/bootstrap tree set."""

    tree_set_path: Path
    traits_path: Path
    trait: str
    taxon_column: str
    model: str
    state_ordering: str
    ordered_states: list[str]
    burnin_fraction: float
    total_tree_count: int
    burnin_tree_count: int
    kept_tree_count: int
    shared_tree_taxa: list[str]
    analysis_taxa: list[str]
    rooted_topology_count: int
    unrooted_topology_count: int
    tree_rows: list[AncestralTransitionTreeRow]
    branch_rows: list[AncestralTransitionTreeSetBranchRow]
    transition_rows: list[AncestralTransitionTreeSetCountRow]
    exclusions: list[AncestralTransitionExclusion]
    warnings: list[str]


def summarize_ancestral_transitions(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    model: str = "fitch",
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
) -> AncestralTransitionReport:
    """Count inferred discrete ancestral transitions for one rooted tree."""
    reconstruction = reconstruct_discrete_ancestral_states(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=model,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
    )
    branch_rows = _build_transition_branch_rows(reconstruction)
    transition_rows = _summarize_transition_rows(branch_rows)
    exclusions = [
        AncestralTransitionExclusion(
            taxon=taxon,
            reason="missing_discrete_trait_state",
        )
        for taxon in reconstruction.dropped_missing_taxa
    ]
    warnings = list(reconstruction.warnings)
    if any(row.uncertain_change_count > 0 for row in transition_rows):
        warnings.append(
            "one or more inferred ancestral transitions remain uncertain because parent and child state sets overlap"
        )
    return AncestralTransitionReport(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        taxon_column=reconstruction.taxon_column,
        model=reconstruction.model,
        state_ordering=reconstruction.state_ordering,
        ordered_states=reconstruction.ordered_states,
        taxon_count=reconstruction.taxon_count,
        branch_rows=branch_rows,
        transition_rows=transition_rows,
        exclusions=exclusions,
        warnings=warnings,
    )


def summarize_ancestral_transition_report(
    report: AncestralTransitionReport,
) -> AncestralTransitionSummary:
    """Summarize the main review facts for one ancestral transition report."""
    changed_rows = [row for row in report.branch_rows if row.changed]
    certain_change_count = sum(
        row.certainty_class == "certain_change" for row in changed_rows
    )
    uncertain_change_count = sum(
        row.certainty_class == "uncertain_change" for row in changed_rows
    )
    return AncestralTransitionSummary(
        trait=report.trait,
        taxon_column=report.taxon_column,
        model=report.model,
        state_ordering=report.state_ordering,
        analyzed_taxon_count=report.taxon_count,
        excluded_taxon_count=len(report.exclusions),
        total_branch_count=len(report.branch_rows),
        changed_branch_count=len(changed_rows),
        certain_change_count=certain_change_count,
        uncertain_change_count=uncertain_change_count,
        unique_transition_count=len(report.transition_rows),
        warning_count=len(report.warnings),
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
    _validate_burnin_fraction(burnin_fraction)
    _source_format, trees = _load_tree_set_trees(tree_set_path)
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
    shared_tree_taxa = sorted(_shared_taxa(kept_trees))
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
    ) = _prepare_analysis_tree_set(
        tree_set_path=tree_set_path,
        traits_path=traits_path,
        taxon_column=taxon_column,
        trait=trait,
        kept_tree_entries=kept_tree_entries,
        shared_tree_taxa=shared_tree_taxa,
        dataset_kind="discrete",
    )
    warnings.extend(dataset_warnings)
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
            report = summarize_ancestral_transitions(
                current_tree_path,
                traits_path,
                trait=trait,
                taxon_column=taxon_column,
                model=model,
                state_ordering=state_ordering,
                ordered_states=ordered_states,
            )
            summary = summarize_ancestral_transition_report(report)
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
            for row in report.branch_rows:
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
            per_tree_transition_counts[str(source_tree_index)] = report.transition_rows
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


def write_ancestral_transition_summary_table(
    path: Path,
    report: AncestralTransitionReport,
) -> Path:
    """Write one summary ledger for ancestral transition counting."""
    summary = summarize_ancestral_transition_report(report)
    return write_ancestral_rows(
        path,
        columns=[
            "trait",
            "taxon_column",
            "model",
            "state_ordering",
            "analyzed_taxon_count",
            "excluded_taxon_count",
            "total_branch_count",
            "changed_branch_count",
            "certain_change_count",
            "uncertain_change_count",
            "unique_transition_count",
            "warning_count",
        ],
        rows=[
            {
                "trait": summary.trait,
                "taxon_column": summary.taxon_column,
                "model": summary.model,
                "state_ordering": summary.state_ordering,
                "analyzed_taxon_count": str(summary.analyzed_taxon_count),
                "excluded_taxon_count": str(summary.excluded_taxon_count),
                "total_branch_count": str(summary.total_branch_count),
                "changed_branch_count": str(summary.changed_branch_count),
                "certain_change_count": str(summary.certain_change_count),
                "uncertain_change_count": str(summary.uncertain_change_count),
                "unique_transition_count": str(summary.unique_transition_count),
                "warning_count": str(summary.warning_count),
            }
        ],
    )


def write_ancestral_transition_branch_table(
    path: Path,
    report: AncestralTransitionReport,
) -> Path:
    """Write one branchwise ancestral transition ledger."""
    return write_ancestral_rows(
        path,
        columns=[
            "parent_node",
            "child_node",
            "child_descendant_taxa",
            "branch_length",
            "parent_most_likely_state",
            "child_most_likely_state",
            "parent_state_set",
            "child_state_set",
            "overlapping_states",
            "changed",
            "transition",
            "certainty_class",
        ],
        rows=[
            {
                "parent_node": row.parent_node,
                "child_node": row.child_node,
                "child_descendant_taxa": ",".join(row.child_descendant_taxa),
                "branch_length": _format_optional_float(row.branch_length),
                "parent_most_likely_state": row.parent_most_likely_state,
                "child_most_likely_state": row.child_most_likely_state,
                "parent_state_set": ",".join(row.parent_state_set),
                "child_state_set": ",".join(row.child_state_set),
                "overlapping_states": ",".join(row.overlapping_states),
                "changed": str(row.changed).lower(),
                "transition": row.transition,
                "certainty_class": row.certainty_class,
            }
            for row in report.branch_rows
        ],
    )


def write_ancestral_transition_count_table(
    path: Path,
    report: AncestralTransitionReport,
) -> Path:
    """Write one transition-pair count ledger."""
    return write_ancestral_rows(
        path,
        columns=[
            "transition",
            "source_state",
            "target_state",
            "certain_change_count",
            "uncertain_change_count",
            "total_change_count",
        ],
        rows=[
            {
                "transition": row.transition,
                "source_state": row.source_state,
                "target_state": row.target_state,
                "certain_change_count": str(row.certain_change_count),
                "uncertain_change_count": str(row.uncertain_change_count),
                "total_change_count": str(row.total_change_count),
            }
            for row in report.transition_rows
        ],
    )


def write_ancestral_transition_exclusion_table(
    path: Path,
    report: AncestralTransitionReport,
) -> Path:
    """Write one excluded-tip ledger for ancestral transition counting."""
    return write_ancestral_rows(
        path,
        columns=["taxon", "reason"],
        rows=[
            {
                "taxon": row.taxon,
                "reason": row.reason,
            }
            for row in report.exclusions
        ],
    )


def write_ancestral_transition_tree_set_summary_table(
    path: Path,
    report: AncestralTransitionTreeSetReport,
) -> Path:
    """Write one summary ledger for ancestral transition counting across a tree set."""
    summary = summarize_ancestral_transition_tree_set_report(report)
    return write_ancestral_rows(
        path,
        columns=[
            "trait",
            "taxon_column",
            "model",
            "state_ordering",
            "total_tree_count",
            "burnin_tree_count",
            "kept_tree_count",
            "shared_tree_taxon_count",
            "analysis_taxon_count",
            "rooted_topology_count",
            "unrooted_topology_count",
            "transition_pair_count",
            "topology_sensitive_transition_pair_count",
            "uncertainty_sensitive_transition_pair_count",
            "warning_count",
        ],
        rows=[
            {
                "trait": summary.trait,
                "taxon_column": summary.taxon_column,
                "model": summary.model,
                "state_ordering": summary.state_ordering,
                "total_tree_count": str(summary.total_tree_count),
                "burnin_tree_count": str(summary.burnin_tree_count),
                "kept_tree_count": str(summary.kept_tree_count),
                "shared_tree_taxon_count": str(summary.shared_tree_taxon_count),
                "analysis_taxon_count": str(summary.analysis_taxon_count),
                "rooted_topology_count": str(summary.rooted_topology_count),
                "unrooted_topology_count": str(summary.unrooted_topology_count),
                "transition_pair_count": str(summary.transition_pair_count),
                "topology_sensitive_transition_pair_count": str(
                    summary.topology_sensitive_transition_pair_count
                ),
                "uncertainty_sensitive_transition_pair_count": str(
                    summary.uncertainty_sensitive_transition_pair_count
                ),
                "warning_count": str(summary.warning_count),
            }
        ],
    )


def write_ancestral_transition_tree_set_tree_table(
    path: Path,
    report: AncestralTransitionTreeSetReport,
) -> Path:
    """Write one retained-tree transition summary ledger."""
    return write_ancestral_rows(
        path,
        columns=[
            "source_tree_index",
            "post_burnin_index",
            "rooted_topology_id",
            "unrooted_topology_id",
            "total_branch_count",
            "changed_branch_count",
            "certain_change_count",
            "uncertain_change_count",
        ],
        rows=[
            {
                "source_tree_index": str(row.source_tree_index),
                "post_burnin_index": str(row.post_burnin_index),
                "rooted_topology_id": row.rooted_topology_id,
                "unrooted_topology_id": row.unrooted_topology_id,
                "total_branch_count": str(row.total_branch_count),
                "changed_branch_count": str(row.changed_branch_count),
                "certain_change_count": str(row.certain_change_count),
                "uncertain_change_count": str(row.uncertain_change_count),
            }
            for row in report.tree_rows
        ],
    )


def write_ancestral_transition_tree_set_branch_table(
    path: Path,
    report: AncestralTransitionTreeSetReport,
) -> Path:
    """Write one per-branch ancestral transition ledger across retained trees."""
    return write_ancestral_rows(
        path,
        columns=[
            "source_tree_index",
            "post_burnin_index",
            "rooted_topology_id",
            "unrooted_topology_id",
            "parent_node",
            "child_node",
            "child_descendant_taxa",
            "branch_length",
            "parent_most_likely_state",
            "child_most_likely_state",
            "parent_state_set",
            "child_state_set",
            "overlapping_states",
            "changed",
            "transition",
            "certainty_class",
        ],
        rows=[
            {
                "source_tree_index": str(row.source_tree_index),
                "post_burnin_index": str(row.post_burnin_index),
                "rooted_topology_id": row.rooted_topology_id,
                "unrooted_topology_id": row.unrooted_topology_id,
                "parent_node": row.parent_node,
                "child_node": row.child_node,
                "child_descendant_taxa": ",".join(row.child_descendant_taxa),
                "branch_length": _format_optional_float(row.branch_length),
                "parent_most_likely_state": row.parent_most_likely_state,
                "child_most_likely_state": row.child_most_likely_state,
                "parent_state_set": ",".join(row.parent_state_set),
                "child_state_set": ",".join(row.child_state_set),
                "overlapping_states": ",".join(row.overlapping_states),
                "changed": str(row.changed).lower(),
                "transition": row.transition,
                "certainty_class": row.certainty_class,
            }
            for row in report.branch_rows
        ],
    )


def write_ancestral_transition_tree_set_count_table(
    path: Path,
    report: AncestralTransitionTreeSetReport,
) -> Path:
    """Write one transition-pair summary ledger across retained trees."""
    return write_ancestral_rows(
        path,
        columns=[
            "transition",
            "source_state",
            "target_state",
            "tree_presence_count",
            "tree_presence_fraction",
            "mean_certain_change_count",
            "mean_uncertain_change_count",
            "mean_total_change_count",
            "minimum_total_change_count",
            "maximum_total_change_count",
            "lower_95_total_change_count",
            "upper_95_total_change_count",
            "stability_class",
        ],
        rows=[
            {
                "transition": row.transition,
                "source_state": row.source_state,
                "target_state": row.target_state,
                "tree_presence_count": str(row.tree_presence_count),
                "tree_presence_fraction": str(row.tree_presence_fraction),
                "mean_certain_change_count": str(row.mean_certain_change_count),
                "mean_uncertain_change_count": str(row.mean_uncertain_change_count),
                "mean_total_change_count": str(row.mean_total_change_count),
                "minimum_total_change_count": str(row.minimum_total_change_count),
                "maximum_total_change_count": str(row.maximum_total_change_count),
                "lower_95_total_change_count": str(row.lower_95_total_change_count),
                "upper_95_total_change_count": str(row.upper_95_total_change_count),
                "stability_class": row.stability_class,
            }
            for row in report.transition_rows
        ],
    )


def _build_transition_branch_rows(
    report: DiscreteAncestralReport,
) -> list[AncestralTransitionBranchRow]:
    tree = loads_newick(report.analysis_tree_newick)
    estimate_by_node = {estimate.node: estimate for estimate in report.estimates}
    branch_rows: list[AncestralTransitionBranchRow] = []

    def visit(node: TreeNode) -> None:
        parent_estimate = estimate_by_node[_node_signature(node)]
        for child in node.children:
            child_estimate = estimate_by_node[_node_signature(child)]
            overlapping_states = sorted(
                set(parent_estimate.state_set) & set(child_estimate.state_set)
            )
            changed = (
                parent_estimate.most_likely_state != child_estimate.most_likely_state
            )
            branch_rows.append(
                AncestralTransitionBranchRow(
                    parent_node=parent_estimate.node,
                    child_node=child_estimate.node,
                    child_descendant_taxa=child_estimate.descendant_taxa,
                    branch_length=child.branch_length,
                    parent_most_likely_state=parent_estimate.most_likely_state,
                    child_most_likely_state=child_estimate.most_likely_state,
                    parent_state_set=parent_estimate.state_set,
                    child_state_set=child_estimate.state_set,
                    overlapping_states=overlapping_states,
                    changed=changed,
                    transition=(
                        f"{parent_estimate.most_likely_state}->{child_estimate.most_likely_state}"
                        if changed
                        else ""
                    ),
                    certainty_class=_transition_certainty_class(
                        changed=changed,
                        overlapping_states=overlapping_states,
                        parent_state_set=parent_estimate.state_set,
                        child_state_set=child_estimate.state_set,
                    ),
                )
            )
            visit(child)

    visit(tree.root)
    return branch_rows


def _summarize_transition_rows(
    branch_rows: list[AncestralTransitionBranchRow],
) -> list[AncestralTransitionCountRow]:
    grouped: dict[str, list[AncestralTransitionBranchRow]] = {}
    for row in branch_rows:
        if row.changed:
            grouped.setdefault(row.transition, []).append(row)
    transition_rows: list[AncestralTransitionCountRow] = []
    for transition in sorted(grouped):
        transition_branch_rows = grouped[transition]
        source_state, target_state = transition.split("->", maxsplit=1)
        certain_change_count = sum(
            row.certainty_class == "certain_change" for row in transition_branch_rows
        )
        uncertain_change_count = sum(
            row.certainty_class == "uncertain_change" for row in transition_branch_rows
        )
        transition_rows.append(
            AncestralTransitionCountRow(
                transition=transition,
                source_state=source_state,
                target_state=target_state,
                certain_change_count=certain_change_count,
                uncertain_change_count=uncertain_change_count,
                total_change_count=len(transition_branch_rows),
            )
        )
    return transition_rows


def _summarize_transition_rows_across_trees(
    *,
    tree_rows: list[AncestralTransitionTreeRow],
    per_tree_transition_counts: dict[str, list[AncestralTransitionCountRow]],
) -> list[AncestralTransitionTreeSetCountRow]:
    counts_by_transition: dict[str, dict[int, AncestralTransitionCountRow]] = {}
    tree_count = len(tree_rows)
    for tree_row in tree_rows:
        transition_rows = per_tree_transition_counts[str(tree_row.source_tree_index)]
        for transition_row in transition_rows:
            counts_by_transition.setdefault(transition_row.transition, {})[
                tree_row.source_tree_index
            ] = transition_row
    summary_rows: list[AncestralTransitionTreeSetCountRow] = []
    for transition in sorted(counts_by_transition):
        present_rows = counts_by_transition[transition]
        source_state, target_state = transition.split("->", maxsplit=1)
        certain_counts = [
            float(row.certain_change_count) for row in present_rows.values()
        ]
        uncertain_counts = [
            float(row.uncertain_change_count) for row in present_rows.values()
        ]
        total_counts = [float(row.total_change_count) for row in present_rows.values()]
        tree_presence_fraction = len(present_rows) / tree_count
        mean_uncertain_change_count = statistics.fmean(uncertain_counts)
        if tree_presence_fraction < 1.0:
            stability_class = "topology_sensitive"
        elif mean_uncertain_change_count > 0.0:
            stability_class = "uncertainty_sensitive"
        else:
            stability_class = "stable"
        summary_rows.append(
            AncestralTransitionTreeSetCountRow(
                transition=transition,
                source_state=source_state,
                target_state=target_state,
                tree_presence_count=len(present_rows),
                tree_presence_fraction=tree_presence_fraction,
                mean_certain_change_count=statistics.fmean(certain_counts),
                mean_uncertain_change_count=mean_uncertain_change_count,
                mean_total_change_count=statistics.fmean(total_counts),
                minimum_total_change_count=int(min(total_counts)),
                maximum_total_change_count=int(max(total_counts)),
                lower_95_total_change_count=_empirical_quantile(total_counts, 0.025),
                upper_95_total_change_count=_empirical_quantile(total_counts, 0.975),
                stability_class=stability_class,
            )
        )
    return summary_rows


def _transition_certainty_class(
    *,
    changed: bool,
    overlapping_states: list[str],
    parent_state_set: list[str],
    child_state_set: list[str],
) -> str:
    if changed:
        if not overlapping_states:
            return "certain_change"
        return "uncertain_change"
    if len(parent_state_set) == 1 and len(child_state_set) == 1:
        return "certain_no_change"
    return "uncertain_no_change"


def _node_signature(node: TreeNode) -> str:
    if node.is_leaf():
        return node.name or "<unnamed>"
    descendant_taxa: list[str] = []
    for child in node.children:
        descendant_taxa.extend(_node_signature_taxa(child))
    return "|".join(sorted(descendant_taxa))


def _node_signature_taxa(node: TreeNode) -> list[str]:
    if node.is_leaf():
        return [node.name] if node.name is not None else []
    descendant_taxa: list[str] = []
    for child in node.children:
        descendant_taxa.extend(_node_signature_taxa(child))
    return descendant_taxa


def _format_optional_float(value: float | None) -> str:
    if value is None:
        return ""
    return format(value, ".15g")


def _empirical_quantile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    index = (len(ordered) - 1) * probability
    lower_index = int(index)
    upper_index = min(lower_index + 1, len(ordered) - 1)
    if lower_index == upper_index:
        return ordered[lower_index]
    fraction = index - lower_index
    return ordered[lower_index] + (
        (ordered[upper_index] - ordered[lower_index]) * fraction
    )
