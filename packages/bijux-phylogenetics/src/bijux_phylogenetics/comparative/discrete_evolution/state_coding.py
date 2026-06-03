from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import load_taxon_table, load_tsv_summary
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.runtime.errors import AncestralReconstructionError

from .models import (
    DominantStateBiasReport,
    SparseStateInstabilityReport,
    StateCodingAuditReport,
    StateCodingAuditRow,
    StateCodingIssue,
    StateCodingValidationReport,
    StateImbalanceReport,
    StateImbalanceWarning,
)

_DEFAULT_ALLOWED_STATE_PATTERN_BLOCKLIST = ("|", "/", ";")


def _detect_sparse_state_instability(
    *,
    state_counts: dict[str, int],
    count_matrix: dict[str, dict[str, float]],
) -> SparseStateInstabilityReport:
    sparse_states = sorted(state for state, count in state_counts.items() if count < 2)
    zero_support_transitions = sorted(
        f"{source}->{target}"
        for source, targets in count_matrix.items()
        for target, value in targets.items()
        if value == 0.0
    )
    warning_count = int(bool(sparse_states)) + int(bool(zero_support_transitions))
    return SparseStateInstabilityReport(
        sparse_states=sparse_states,
        zero_support_transitions=zero_support_transitions,
        warning_count=warning_count,
        unstable=bool(sparse_states or zero_support_transitions),
    )


def _summarize_dominant_state_bias(
    state_counts: dict[str, int],
) -> DominantStateBiasReport:
    if not state_counts:
        return DominantStateBiasReport(
            dominant_states=[],
            dominant_fraction=0.0,
            biased=False,
            message=None,
        )
    total = sum(state_counts.values())
    max_count = max(state_counts.values())
    dominant_states = sorted(
        state for state, count in state_counts.items() if count == max_count
    )
    dominant_fraction = float(format(max_count / max(total, 1), ".15g"))
    biased = dominant_fraction >= 0.8
    return DominantStateBiasReport(
        dominant_states=dominant_states,
        dominant_fraction=dominant_fraction,
        biased=biased,
        message=(
            "one state dominates most taxa and may compress minority-state transition evidence"
            if biased
            else None
        ),
    )


def validate_discrete_state_coding(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    allowed_states: list[str] | None = None,
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
) -> StateCodingValidationReport:
    """Detect impossible or unsupported discrete-state labels."""
    tree = load_tree(tree_path)
    table = (
        load_tsv_summary(traits_path)
        if taxon_column is None
        else load_taxon_table(traits_path, taxon_column=taxon_column)
    )
    if trait not in table.columns:
        raise AncestralReconstructionError(
            f"trait table does not contain column '{trait}'"
        )
    tree_taxa = set(tree.tip_names)
    if state_ordering not in {"unordered", "ordered"}:
        raise ValueError(f"unsupported state ordering: {state_ordering}")
    ordered = list(ordered_states or [])
    if ordered and len(set(ordered)) != len(ordered):
        raise AncestralReconstructionError(
            "ordered state vocabulary contains duplicate labels"
        )
    allowed = list(allowed_states or ordered or [])
    allowed_set = set(allowed)
    issues: list[StateCodingIssue] = []
    usable_taxa: list[str] = []
    observed_states: set[str] = set()
    for row in table.rows:
        taxon = row[table.taxon_column]
        if taxon not in tree_taxa:
            continue
        raw_state = row[trait].strip()
        if not raw_state:
            continue
        invalid_delimiter = next(
            (
                token
                for token in _DEFAULT_ALLOWED_STATE_PATTERN_BLOCKLIST
                if token in raw_state
            ),
            None,
        )
        if invalid_delimiter is not None:
            issues.append(
                StateCodingIssue(
                    taxon=taxon,
                    raw_state=raw_state,
                    code="unsupported-state-delimiter",
                    message=f"state label contains reserved delimiter '{invalid_delimiter}'",
                )
            )
            continue
        if allowed and raw_state not in allowed_set:
            issues.append(
                StateCodingIssue(
                    taxon=taxon,
                    raw_state=raw_state,
                    code="unordered-state-vocabulary"
                    if state_ordering == "ordered" and ordered
                    else "unsupported-state-label",
                    message=(
                        "state label is not present in the declared ordered vocabulary"
                        if state_ordering == "ordered" and ordered
                        else "state label is not present in the allowed state vocabulary"
                    ),
                )
            )
            continue
        observed_states.add(raw_state)
        usable_taxa.append(taxon)
    if state_ordering == "ordered" and ordered:
        missing_from_order = sorted(observed_states - set(ordered))
        issues.extend(
            StateCodingIssue(
                taxon="",
                raw_state=state,
                code="unordered-state-vocabulary",
                message="observed state is missing from the declared ordered vocabulary",
            )
            for state in missing_from_order
        )
    return StateCodingValidationReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=table.taxon_column,
        trait=trait,
        allowed_states=allowed,
        state_ordering=state_ordering,
        ordered_states=ordered,
        valid=not issues,
        issues=issues,
        observed_states=sorted(observed_states),
        usable_taxa=sorted(usable_taxa),
    )


def audit_discrete_state_coding(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    allowed_states: list[str] | None = None,
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
    coding_map: dict[str, str] | None = None,
) -> StateCodingAuditReport:
    """Show how raw metadata states become model states or get excluded."""
    tree = load_tree(tree_path)
    table = (
        load_tsv_summary(traits_path)
        if taxon_column is None
        else load_taxon_table(traits_path, taxon_column=taxon_column)
    )
    if trait not in table.columns:
        raise AncestralReconstructionError(
            f"trait table does not contain column '{trait}'"
        )
    tree_taxa = set(tree.tip_names)
    mapping = dict(coding_map or {})
    ordered = list(ordered_states or [])
    allowed = list(allowed_states or ordered or [])
    allowed_set = set(allowed)
    rows: list[StateCodingAuditRow] = []
    for row in table.rows:
        taxon = row[table.taxon_column]
        raw_state = row[trait].strip()
        in_tree = taxon in tree_taxa
        if not raw_state:
            rows.append(
                StateCodingAuditRow(
                    taxon=taxon,
                    raw_state=raw_state,
                    normalized_state=None,
                    in_tree=in_tree,
                    included=False,
                    issue_code="missing-state",
                    note="state is blank and cannot be used",
                )
            )
            continue
        normalized_state = mapping.get(raw_state, raw_state)
        invalid_delimiter = next(
            (
                token
                for token in _DEFAULT_ALLOWED_STATE_PATTERN_BLOCKLIST
                if token in normalized_state
            ),
            None,
        )
        if not in_tree:
            rows.append(
                StateCodingAuditRow(
                    taxon=taxon,
                    raw_state=raw_state,
                    normalized_state=normalized_state,
                    in_tree=False,
                    included=False,
                    issue_code="taxon-not-in-tree",
                    note="taxon does not overlap the tree and is excluded before state modeling",
                )
            )
            continue
        if invalid_delimiter is not None:
            rows.append(
                StateCodingAuditRow(
                    taxon=taxon,
                    raw_state=raw_state,
                    normalized_state=normalized_state,
                    in_tree=True,
                    included=False,
                    issue_code="unsupported-state-delimiter",
                    note=f"normalized state contains reserved delimiter '{invalid_delimiter}'",
                )
            )
            continue
        if allowed and normalized_state not in allowed_set:
            rows.append(
                StateCodingAuditRow(
                    taxon=taxon,
                    raw_state=raw_state,
                    normalized_state=normalized_state,
                    in_tree=True,
                    included=False,
                    issue_code="unsupported-state-label"
                    if state_ordering != "ordered"
                    else "unordered-state-vocabulary",
                    note=(
                        "normalized state is absent from the declared ordered vocabulary"
                        if state_ordering == "ordered" and ordered
                        else "normalized state is absent from the allowed state vocabulary"
                    ),
                )
            )
            continue
        rows.append(
            StateCodingAuditRow(
                taxon=taxon,
                raw_state=raw_state,
                normalized_state=normalized_state,
                in_tree=True,
                included=True,
                issue_code=None,
                note="state is retained for discrete-state modeling",
            )
        )
    included_row_count = sum(1 for row in rows if row.included)
    return StateCodingAuditReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=table.taxon_column,
        trait=trait,
        state_ordering=state_ordering,
        ordered_states=ordered,
        coding_map=mapping,
        row_count=len(rows),
        included_row_count=included_row_count,
        excluded_row_count=len(rows) - included_row_count,
        rows=rows,
    )


def detect_state_imbalance_problems(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
) -> StateImbalanceReport:
    """Flag rare or degenerate discrete-state inputs over the tree overlap."""
    tree = load_tree(tree_path)
    table = (
        load_tsv_summary(traits_path)
        if taxon_column is None
        else load_taxon_table(traits_path, taxon_column=taxon_column)
    )
    if trait not in table.columns:
        raise AncestralReconstructionError(
            f"trait table does not contain column '{trait}'"
        )
    tree_taxa = set(tree.tip_names)
    state_counts: dict[str, int] = {}
    for row in table.rows:
        taxon = row[table.taxon_column]
        state = row[trait].strip()
        if taxon in tree_taxa and state:
            state_counts[state] = state_counts.get(state, 0) + 1
    observed_states = sorted(state_counts)
    warnings: list[StateImbalanceWarning] = []
    rare_states = sorted(state for state, count in state_counts.items() if count < 2)
    if len(observed_states) < 2:
        warnings.append(
            StateImbalanceWarning(
                code="single-state-dataset",
                message="only one observed state remains after pruning to usable tree taxa",
                affected_states=observed_states,
            )
        )
    if rare_states:
        warnings.append(
            StateImbalanceWarning(
                code="rare-states",
                message="one or more states are represented by fewer than two taxa",
                affected_states=rare_states,
            )
        )
    dominant_fraction = (
        max(state_counts.values()) / max(sum(state_counts.values()), 1)
        if state_counts
        else 0.0
    )
    if dominant_fraction >= 0.8 and observed_states:
        dominant_states = [
            state
            for state, count in state_counts.items()
            if count == max(state_counts.values())
        ]
        warnings.append(
            StateImbalanceWarning(
                code="dominant-state-skew",
                message="one state dominates most observed taxa and may overwhelm transition inference",
                affected_states=sorted(dominant_states),
            )
        )
    return StateImbalanceReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=table.taxon_column,
        trait=trait,
        taxon_count=sum(state_counts.values()),
        observed_states=observed_states,
        state_counts=state_counts,
        warnings=warnings,
    )
