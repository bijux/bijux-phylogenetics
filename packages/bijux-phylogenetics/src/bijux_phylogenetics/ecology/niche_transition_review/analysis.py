from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.ancestral.discrete import DiscreteAncestralReport
from bijux_phylogenetics.comparative.discrete_evolution import (
    audit_discrete_state_coding,
)
from bijux_phylogenetics.io.newick import loads_newick

from .contracts import (
    NicheStateNodeRow,
    NicheTransitionBranchRow,
    NicheTransitionCladeRow,
    NicheTransitionCountRow,
    NicheTransitionExclusionRow,
    NicheTransitionRateRow,
    NicheTransitionSummary,
)
from .shared import (
    node_signature,
    stable_float,
    transition_certainty_class,
)


def build_node_rows(report: DiscreteAncestralReport) -> list[NicheStateNodeRow]:
    root_node = report.estimates[0].node
    return [
        NicheStateNodeRow(
            node=estimate.node,
            node_name=estimate.node_name,
            descendant_taxa=list(estimate.descendant_taxa),
            most_likely_niche=estimate.most_likely_state,
            niche_probabilities=dict(sorted(estimate.state_probabilities.items())),
            confidence=stable_float(estimate.confidence),
            ambiguous=estimate.ambiguous,
            is_root=estimate.node == root_node,
        )
        for estimate in report.estimates
        if not estimate.is_tip
    ]


def build_rate_rows(report: DiscreteAncestralReport) -> list[NicheTransitionRateRow]:
    return [
        NicheTransitionRateRow(
            source_niche=row.source_state,
            target_niche=row.target_state,
            transition_allowed=row.transition_allowed,
            step_distance=row.step_distance,
            rate=stable_float(row.rate),
        )
        for row in report.transition_rate_rows
    ]


def build_branch_rows(
    report: DiscreteAncestralReport,
) -> list[NicheTransitionBranchRow]:
    tree = loads_newick(report.analysis_tree_newick)
    estimate_by_node = {estimate.node: estimate for estimate in report.estimates}
    branch_rows: list[NicheTransitionBranchRow] = []

    def visit(node) -> None:
        parent_estimate = estimate_by_node[node_signature(node)]
        for child in node.children:
            child_estimate = estimate_by_node[node_signature(child)]
            overlapping_niches = sorted(
                set(parent_estimate.state_set) & set(child_estimate.state_set)
            )
            changed = (
                parent_estimate.most_likely_state != child_estimate.most_likely_state
            )
            support = stable_float(
                min(parent_estimate.confidence, child_estimate.confidence)
            )
            branch_rows.append(
                NicheTransitionBranchRow(
                    branch_id=child_estimate.node,
                    parent_node=parent_estimate.node,
                    child_node=child_estimate.node,
                    child_descendant_taxa=list(child_estimate.descendant_taxa),
                    branch_length=child.branch_length,
                    parent_most_likely_niche=parent_estimate.most_likely_state,
                    child_most_likely_niche=child_estimate.most_likely_state,
                    parent_niche_set=list(parent_estimate.state_set),
                    child_niche_set=list(child_estimate.state_set),
                    overlapping_niches=overlapping_niches,
                    changed=changed,
                    transition=(
                        f"{parent_estimate.most_likely_state}->{child_estimate.most_likely_state}"
                        if changed
                        else ""
                    ),
                    certainty_class=transition_certainty_class(
                        changed=changed,
                        overlapping_niches=overlapping_niches,
                        parent_niche_set=parent_estimate.state_set,
                        child_niche_set=child_estimate.state_set,
                    ),
                    support=support,
                    strongly_supported=support >= 0.9 and changed,
                    parent_confidence=stable_float(parent_estimate.confidence),
                    child_confidence=stable_float(child_estimate.confidence),
                )
            )
            visit(child)

    visit(tree.root)
    return branch_rows


def build_count_rows(
    branch_rows: list[NicheTransitionBranchRow],
) -> list[NicheTransitionCountRow]:
    grouped: dict[str, list[NicheTransitionBranchRow]] = {}
    for row in branch_rows:
        if row.changed:
            grouped.setdefault(row.transition, []).append(row)
    count_rows: list[NicheTransitionCountRow] = []
    for transition in sorted(grouped):
        rows = grouped[transition]
        source_niche, target_niche = transition.split("->", maxsplit=1)
        count_rows.append(
            NicheTransitionCountRow(
                transition=transition,
                source_niche=source_niche,
                target_niche=target_niche,
                certain_transition_count=sum(
                    row.certainty_class == "certain_transition" for row in rows
                ),
                uncertain_transition_count=sum(
                    row.certainty_class == "uncertain_transition" for row in rows
                ),
                total_transition_count=len(rows),
                strongly_supported_transition_count=sum(
                    row.strongly_supported for row in rows
                ),
            )
        )
    return count_rows


def build_clade_rows(
    report: DiscreteAncestralReport,
    branch_rows: list[NicheTransitionBranchRow],
) -> list[NicheTransitionCladeRow]:
    node_rows = build_node_rows(report)
    non_root_rows = [row for row in node_rows if not row.is_root]
    unsorted_rows: list[NicheTransitionCladeRow] = []
    for node_row in non_root_rows:
        descendant_taxa = set(node_row.descendant_taxa)
        clade_branch_rows = [
            row
            for row in branch_rows
            if set(row.child_descendant_taxa).issubset(descendant_taxa)
        ]
        changed_rows = [row for row in clade_branch_rows if row.changed]
        transition_counts: dict[str, int] = {}
        for row in changed_rows:
            transition_counts[row.transition] = (
                transition_counts.get(row.transition, 0) + 1
            )
        dominant_transition = ""
        dominant_transition_count = 0
        if transition_counts:
            dominant_transition, dominant_transition_count = sorted(
                transition_counts.items(),
                key=lambda item: (-item[1], item[0]),
            )[0]
        certain_transition_count = sum(
            row.certainty_class == "certain_transition" for row in changed_rows
        )
        uncertain_transition_count = sum(
            row.certainty_class == "uncertain_transition" for row in changed_rows
        )
        strongly_supported_transition_count = sum(
            row.strongly_supported for row in changed_rows
        )
        descendant_internal_node_count = sum(
            other.node != node_row.node
            and set(other.descendant_taxa).issubset(descendant_taxa)
            for other in non_root_rows
        )
        shift_burden_score = stable_float(
            certain_transition_count
            + 0.5 * uncertain_transition_count
            + 0.25 * strongly_supported_transition_count
        )
        unsorted_rows.append(
            NicheTransitionCladeRow(
                node=node_row.node,
                node_name=node_row.node_name,
                descendant_taxa=node_row.descendant_taxa,
                descendant_taxon_count=len(node_row.descendant_taxa),
                descendant_internal_node_count=descendant_internal_node_count,
                changed_branch_count=len(changed_rows),
                certain_transition_count=certain_transition_count,
                uncertain_transition_count=uncertain_transition_count,
                strongly_supported_transition_count=strongly_supported_transition_count,
                transition_diversity=len(transition_counts),
                dominant_transition=dominant_transition,
                dominant_transition_count=dominant_transition_count,
                shift_burden_score=shift_burden_score,
                contains_repeated_shifts=len(changed_rows) >= 2,
                rank=0,
            )
        )
    ranked_rows = sorted(
        unsorted_rows,
        key=lambda row: (
            -row.shift_burden_score,
            -row.changed_branch_count,
            ",".join(row.descendant_taxa),
        ),
    )
    return [
        NicheTransitionCladeRow(
            node=row.node,
            node_name=row.node_name,
            descendant_taxa=row.descendant_taxa,
            descendant_taxon_count=row.descendant_taxon_count,
            descendant_internal_node_count=row.descendant_internal_node_count,
            changed_branch_count=row.changed_branch_count,
            certain_transition_count=row.certain_transition_count,
            uncertain_transition_count=row.uncertain_transition_count,
            strongly_supported_transition_count=row.strongly_supported_transition_count,
            transition_diversity=row.transition_diversity,
            dominant_transition=row.dominant_transition,
            dominant_transition_count=row.dominant_transition_count,
            shift_burden_score=row.shift_burden_score,
            contains_repeated_shifts=row.contains_repeated_shifts,
            rank=index,
        )
        for index, row in enumerate(ranked_rows, start=1)
    ]


def build_exclusion_rows(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None,
) -> list[NicheTransitionExclusionRow]:
    audit = audit_discrete_state_coding(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    return [
        NicheTransitionExclusionRow(
            taxon=row.taxon,
            raw_niche=row.raw_state,
            normalized_niche=row.normalized_state,
            reason=row.issue_code or "excluded",
            note=row.note,
        )
        for row in audit.rows
        if not row.included
    ]


def build_summary(
    *,
    reconstruction: DiscreteAncestralReport,
    model: str,
    node_rows: list[NicheStateNodeRow],
    rate_rows: list[NicheTransitionRateRow],
    branch_rows: list[NicheTransitionBranchRow],
    clade_rows: list[NicheTransitionCladeRow],
    exclusion_rows: list[NicheTransitionExclusionRow],
    count_rows: list[NicheTransitionCountRow],
) -> NicheTransitionSummary:
    root_estimate = next(
        estimate
        for estimate in reconstruction.estimates
        if estimate.node == reconstruction.estimates[0].node
    )
    return NicheTransitionSummary(
        trait=reconstruction.trait,
        taxon_column=reconstruction.taxon_column,
        model=model,
        internal_model=reconstruction.model,
        analyzed_taxon_count=reconstruction.taxon_count,
        excluded_taxon_count=len(exclusion_rows),
        observed_niche_count=len(reconstruction.observed_states),
        internal_node_count=len(node_rows),
        ambiguous_internal_node_count=sum(row.ambiguous for row in node_rows),
        log_likelihood=stable_float(reconstruction.log_likelihood or 0.0),
        parameter_count=reconstruction.parameter_count or 0,
        aic=stable_float(reconstruction.aic or 0.0),
        transition_rate_row_count=len(rate_rows),
        changed_branch_count=sum(row.total_transition_count for row in count_rows),
        certain_transition_count=sum(
            row.certain_transition_count for row in count_rows
        ),
        uncertain_transition_count=sum(
            row.uncertain_transition_count for row in count_rows
        ),
        strongly_supported_transition_count=sum(
            row.strongly_supported_transition_count for row in count_rows
        ),
        clade_shift_row_count=len(clade_rows),
        repeated_shift_clade_count=sum(
            row.contains_repeated_shifts for row in clade_rows
        ),
        root_niche=root_estimate.most_likely_state,
        root_confidence=stable_float(root_estimate.confidence),
        warning_count=len(reconstruction.warnings),
    )
