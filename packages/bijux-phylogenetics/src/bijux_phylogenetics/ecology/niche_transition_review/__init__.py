from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.ancestral.discrete import (
    DiscreteAncestralReport,
    reconstruct_discrete_ancestral_states,
)
from bijux_phylogenetics.comparative.discrete_evolution import (
    audit_discrete_state_coding,
)
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.io.newick import loads_newick
from bijux_phylogenetics.runtime.errors import AncestralReconstructionError
from .contracts import (
    NicheStateNodeRow,
    NicheTransitionBranchRow,
    NicheTransitionCladeRow,
    NicheTransitionCountRow,
    NicheTransitionExclusionRow,
    NicheTransitionRateRow,
    NicheTransitionReport,
    NicheTransitionSummary,
)

_MODEL_ALIAS_TO_INTERNAL = {
    "er": "equal-rates",
    "sym": "symmetric",
    "ard": "all-rates-different",
}


def summarize_niche_transitions(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    model: str = "er",
) -> NicheTransitionReport:
    """Fit one ecological niche transition model and summarize clade shifts."""
    internal_model = _resolve_internal_model(model)
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
    node_rows = _build_node_rows(reconstruction)
    rate_rows = _build_rate_rows(reconstruction)
    branch_rows = _build_branch_rows(reconstruction)
    count_rows = _build_count_rows(branch_rows)
    clade_rows = _build_clade_rows(reconstruction, branch_rows)
    exclusion_rows = _build_exclusion_rows(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    summary = _build_summary(
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


def write_niche_transition_summary_table(
    path: Path,
    report: NicheTransitionReport,
) -> Path:
    """Write one overall ecological niche transition summary ledger."""
    summary = report.summary
    return write_taxon_rows(
        path,
        columns=[
            "trait",
            "taxon_column",
            "model",
            "internal_model",
            "analyzed_taxon_count",
            "excluded_taxon_count",
            "observed_niche_count",
            "internal_node_count",
            "ambiguous_internal_node_count",
            "log_likelihood",
            "parameter_count",
            "aic",
            "transition_rate_row_count",
            "changed_branch_count",
            "certain_transition_count",
            "uncertain_transition_count",
            "strongly_supported_transition_count",
            "clade_shift_row_count",
            "repeated_shift_clade_count",
            "root_niche",
            "root_confidence",
            "warning_count",
        ],
        rows=[
            {
                "trait": summary.trait,
                "taxon_column": summary.taxon_column,
                "model": summary.model,
                "internal_model": summary.internal_model,
                "analyzed_taxon_count": str(summary.analyzed_taxon_count),
                "excluded_taxon_count": str(summary.excluded_taxon_count),
                "observed_niche_count": str(summary.observed_niche_count),
                "internal_node_count": str(summary.internal_node_count),
                "ambiguous_internal_node_count": str(
                    summary.ambiguous_internal_node_count
                ),
                "log_likelihood": str(summary.log_likelihood),
                "parameter_count": str(summary.parameter_count),
                "aic": str(summary.aic),
                "transition_rate_row_count": str(summary.transition_rate_row_count),
                "changed_branch_count": str(summary.changed_branch_count),
                "certain_transition_count": str(summary.certain_transition_count),
                "uncertain_transition_count": str(summary.uncertain_transition_count),
                "strongly_supported_transition_count": str(
                    summary.strongly_supported_transition_count
                ),
                "clade_shift_row_count": str(summary.clade_shift_row_count),
                "repeated_shift_clade_count": str(summary.repeated_shift_clade_count),
                "root_niche": summary.root_niche,
                "root_confidence": str(summary.root_confidence),
                "warning_count": str(summary.warning_count),
            }
        ],
    )


def write_niche_state_node_table(
    path: Path,
    report: NicheTransitionReport,
) -> Path:
    """Write one internal-node ecological niche probability ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "node",
            "node_name",
            "descendant_taxa",
            "most_likely_niche",
            "niche_probabilities",
            "confidence",
            "ambiguous",
            "is_root",
        ],
        rows=[
            {
                "node": row.node,
                "node_name": row.node_name or "",
                "descendant_taxa": ",".join(row.descendant_taxa),
                "most_likely_niche": row.most_likely_niche,
                "niche_probabilities": json.dumps(
                    row.niche_probabilities,
                    sort_keys=True,
                ),
                "confidence": str(row.confidence),
                "ambiguous": str(row.ambiguous).lower(),
                "is_root": str(row.is_root).lower(),
            }
            for row in report.node_rows
        ],
    )


def write_niche_transition_rate_table(
    path: Path,
    report: NicheTransitionReport,
) -> Path:
    """Write one ecological niche transition-rate ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "source_niche",
            "target_niche",
            "transition_allowed",
            "step_distance",
            "rate",
        ],
        rows=[
            {
                "source_niche": row.source_niche,
                "target_niche": row.target_niche,
                "transition_allowed": str(row.transition_allowed).lower(),
                "step_distance": str(row.step_distance),
                "rate": str(row.rate),
            }
            for row in report.rate_rows
        ],
    )


def write_niche_transition_branch_table(
    path: Path,
    report: NicheTransitionReport,
) -> Path:
    """Write one branchwise ecological niche transition ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "branch_id",
            "parent_node",
            "child_node",
            "child_descendant_taxa",
            "branch_length",
            "parent_most_likely_niche",
            "child_most_likely_niche",
            "parent_niche_set",
            "child_niche_set",
            "overlapping_niches",
            "changed",
            "transition",
            "certainty_class",
            "support",
            "strongly_supported",
            "parent_confidence",
            "child_confidence",
        ],
        rows=[
            {
                "branch_id": row.branch_id,
                "parent_node": row.parent_node,
                "child_node": row.child_node,
                "child_descendant_taxa": ",".join(row.child_descendant_taxa),
                "branch_length": _format_optional_float(row.branch_length),
                "parent_most_likely_niche": row.parent_most_likely_niche,
                "child_most_likely_niche": row.child_most_likely_niche,
                "parent_niche_set": ",".join(row.parent_niche_set),
                "child_niche_set": ",".join(row.child_niche_set),
                "overlapping_niches": ",".join(row.overlapping_niches),
                "changed": str(row.changed).lower(),
                "transition": row.transition,
                "certainty_class": row.certainty_class,
                "support": str(row.support),
                "strongly_supported": str(row.strongly_supported).lower(),
                "parent_confidence": str(row.parent_confidence),
                "child_confidence": str(row.child_confidence),
            }
            for row in report.branch_rows
        ],
    )


def write_niche_transition_count_table(
    path: Path,
    report: NicheTransitionReport,
) -> Path:
    """Write one aggregated ecological niche transition-count ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "transition",
            "source_niche",
            "target_niche",
            "certain_transition_count",
            "uncertain_transition_count",
            "total_transition_count",
            "strongly_supported_transition_count",
        ],
        rows=[
            {
                "transition": row.transition,
                "source_niche": row.source_niche,
                "target_niche": row.target_niche,
                "certain_transition_count": str(row.certain_transition_count),
                "uncertain_transition_count": str(row.uncertain_transition_count),
                "total_transition_count": str(row.total_transition_count),
                "strongly_supported_transition_count": str(
                    row.strongly_supported_transition_count
                ),
            }
            for row in report.count_rows
        ],
    )


def write_niche_transition_clade_table(
    path: Path,
    report: NicheTransitionReport,
) -> Path:
    """Write one internal-clade ecological niche shift ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "node",
            "node_name",
            "descendant_taxa",
            "descendant_taxon_count",
            "descendant_internal_node_count",
            "changed_branch_count",
            "certain_transition_count",
            "uncertain_transition_count",
            "strongly_supported_transition_count",
            "transition_diversity",
            "dominant_transition",
            "dominant_transition_count",
            "shift_burden_score",
            "contains_repeated_shifts",
            "rank",
        ],
        rows=[
            {
                "node": row.node,
                "node_name": row.node_name or "",
                "descendant_taxa": ",".join(row.descendant_taxa),
                "descendant_taxon_count": str(row.descendant_taxon_count),
                "descendant_internal_node_count": str(
                    row.descendant_internal_node_count
                ),
                "changed_branch_count": str(row.changed_branch_count),
                "certain_transition_count": str(row.certain_transition_count),
                "uncertain_transition_count": str(row.uncertain_transition_count),
                "strongly_supported_transition_count": str(
                    row.strongly_supported_transition_count
                ),
                "transition_diversity": str(row.transition_diversity),
                "dominant_transition": row.dominant_transition,
                "dominant_transition_count": str(row.dominant_transition_count),
                "shift_burden_score": str(row.shift_burden_score),
                "contains_repeated_shifts": str(row.contains_repeated_shifts).lower(),
                "rank": str(row.rank),
            }
            for row in report.clade_rows
        ],
    )


def write_niche_transition_exclusion_table(
    path: Path,
    report: NicheTransitionReport,
) -> Path:
    """Write one excluded ecological niche row ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "taxon",
            "raw_niche",
            "normalized_niche",
            "reason",
            "note",
        ],
        rows=[
            {
                "taxon": row.taxon,
                "raw_niche": row.raw_niche,
                "normalized_niche": row.normalized_niche or "",
                "reason": row.reason,
                "note": row.note,
            }
            for row in report.exclusion_rows
        ],
    )


def _build_node_rows(report: DiscreteAncestralReport) -> list[NicheStateNodeRow]:
    root_node = report.estimates[0].node
    return [
        NicheStateNodeRow(
            node=estimate.node,
            node_name=estimate.node_name,
            descendant_taxa=list(estimate.descendant_taxa),
            most_likely_niche=estimate.most_likely_state,
            niche_probabilities=dict(sorted(estimate.state_probabilities.items())),
            confidence=_stable_float(estimate.confidence),
            ambiguous=estimate.ambiguous,
            is_root=estimate.node == root_node,
        )
        for estimate in report.estimates
        if not estimate.is_tip
    ]


def _build_rate_rows(report: DiscreteAncestralReport) -> list[NicheTransitionRateRow]:
    return [
        NicheTransitionRateRow(
            source_niche=row.source_state,
            target_niche=row.target_state,
            transition_allowed=row.transition_allowed,
            step_distance=row.step_distance,
            rate=_stable_float(row.rate),
        )
        for row in report.transition_rate_rows
    ]


def _build_branch_rows(
    report: DiscreteAncestralReport,
) -> list[NicheTransitionBranchRow]:
    tree = loads_newick(report.analysis_tree_newick)
    estimate_by_node = {estimate.node: estimate for estimate in report.estimates}
    branch_rows: list[NicheTransitionBranchRow] = []

    def visit(node) -> None:
        parent_estimate = estimate_by_node[_node_signature(node)]
        for child in node.children:
            child_estimate = estimate_by_node[_node_signature(child)]
            overlapping_niches = sorted(
                set(parent_estimate.state_set) & set(child_estimate.state_set)
            )
            changed = (
                parent_estimate.most_likely_state != child_estimate.most_likely_state
            )
            support = _stable_float(
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
                    certainty_class=_transition_certainty_class(
                        changed=changed,
                        overlapping_niches=overlapping_niches,
                        parent_niche_set=parent_estimate.state_set,
                        child_niche_set=child_estimate.state_set,
                    ),
                    support=support,
                    strongly_supported=support >= 0.9 and changed,
                    parent_confidence=_stable_float(parent_estimate.confidence),
                    child_confidence=_stable_float(child_estimate.confidence),
                )
            )
            visit(child)

    visit(tree.root)
    return branch_rows


def _build_count_rows(
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


def _build_clade_rows(
    report: DiscreteAncestralReport,
    branch_rows: list[NicheTransitionBranchRow],
) -> list[NicheTransitionCladeRow]:
    node_rows = _build_node_rows(report)
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
        shift_burden_score = _stable_float(
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


def _build_exclusion_rows(
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


def _build_summary(
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
        log_likelihood=_stable_float(reconstruction.log_likelihood or 0.0),
        parameter_count=reconstruction.parameter_count or 0,
        aic=_stable_float(reconstruction.aic or 0.0),
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
        root_confidence=_stable_float(root_estimate.confidence),
        warning_count=len(reconstruction.warnings),
    )


def _resolve_internal_model(model: str) -> str:
    try:
        return _MODEL_ALIAS_TO_INTERNAL[model]
    except KeyError as error:
        raise ValueError(
            "ecological niche transition model must be one of: er, sym, ard"
        ) from error


def _transition_certainty_class(
    *,
    changed: bool,
    overlapping_niches: list[str],
    parent_niche_set: list[str],
    child_niche_set: list[str],
) -> str:
    if not changed:
        return "no_transition"
    if overlapping_niches:
        return "uncertain_transition"
    if len(parent_niche_set) == 1 and len(child_niche_set) == 1:
        return "certain_transition"
    return "uncertain_transition"


def _format_optional_float(value: float | None) -> str:
    if value is None:
        return ""
    return str(_stable_float(value))


def _stable_float(value: float) -> float:
    normalized = round(float(value), 8)
    return 0.0 if normalized == -0.0 else normalized


def _node_signature(node) -> str:
    if node.is_leaf():
        return node.name
    return "|".join(sorted(leaf.name for leaf in node.iter_leaves()))
