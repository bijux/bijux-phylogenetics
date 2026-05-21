from __future__ import annotations

from pathlib import Path

from ..registry import ApeParityCase
from .mismatch_policy import (
    _consensus_tree_mismatch_reason,
    _drop_tip_tree_mismatch_reason,
    _extract_clade_tree_mismatch_reason,
    _keep_tip_tree_mismatch_reason,
    _neighbor_joining_tree_mismatch_reason,
    _root_tree_outgroup_mismatch_reason,
    _transition_rate_rows_match,
    _tree_set_structure_mismatch_reason,
    _tree_structure_mismatch_reason,
    _unroot_tree_mismatch_reason,
)
from .normalization import _compare_json

_EXPECTED_ERROR_CONTRACTS: dict[str, tuple[str, str, str]] = {
    "parse-error": (
        "bijux_expected_parse_error_missing",
        "reference_expected_parse_error_missing",
        "parse_error_message_missing",
    ),
    "rooting-error": (
        "bijux_expected_rooting_error_missing",
        "reference_expected_rooting_error_missing",
        "rooting_error_message_missing",
    ),
    "clade-extraction-error": (
        "bijux_expected_clade_extraction_error_missing",
        "reference_expected_clade_extraction_error_missing",
        "clade_extraction_error_message_missing",
    ),
    "mrca-error": (
        "bijux_expected_mrca_error_missing",
        "reference_expected_mrca_error_missing",
        "mrca_error_message_missing",
    ),
    "monophyly-error": (
        "bijux_expected_monophyly_error_missing",
        "reference_expected_monophyly_error_missing",
        "monophyly_error_message_missing",
    ),
    "consensus-error": (
        "bijux_expected_consensus_error_missing",
        "reference_expected_consensus_error_missing",
        "consensus_error_message_missing",
    ),
    "prop-clades-error": (
        "bijux_expected_prop_clades_error_missing",
        "reference_expected_prop_clades_error_missing",
        "prop_clades_error_message_missing",
    ),
    "dna-distance-error": (
        "bijux_expected_dna_distance_error_missing",
        "reference_expected_dna_distance_error_missing",
        "dna_distance_error_message_missing",
    ),
}


def _determine_reference_mismatch_reason(
    *,
    case: ApeParityCase,
    execution_root: Path,
    reference_summary: dict[str, object] | None,
    bijux_summary: dict[str, object] | None,
    reference_rows: list[dict[str, object]] | None,
    bijux_rows: list[dict[str, object]] | None,
    reference_normalized_text: str | None,
    bijux_normalized_text: str | None,
) -> str | None:
    if case.operation in {"read-tree-structure", "write-tree-structure"}:
        return _tree_structure_mismatch_reason(case, execution_root)
    if case.operation == "root-tree-outgroup":
        return _root_tree_outgroup_mismatch_reason(case, execution_root)
    if case.operation == "unroot-tree":
        return _unroot_tree_mismatch_reason(case, execution_root)
    if case.operation == "drop-tree-taxa":
        mismatch_reason = _drop_tip_tree_mismatch_reason(case, execution_root)
        if mismatch_reason is None and not _compare_json(
            reference_summary,
            bijux_summary,
            tolerance=case.tolerance,
        ):
            return "summary_mismatch"
        return mismatch_reason
    if case.operation == "keep-tree-taxa":
        mismatch_reason = _keep_tip_tree_mismatch_reason(case, execution_root)
        if mismatch_reason is None and not _compare_json(
            reference_summary,
            bijux_summary,
            tolerance=case.tolerance,
        ):
            return "summary_mismatch"
        return mismatch_reason
    if case.operation == "extract-tree-clade":
        mismatch_reason = _extract_clade_tree_mismatch_reason(case, execution_root)
        if mismatch_reason is None and not _compare_json(
            reference_summary,
            bijux_summary,
            tolerance=case.tolerance,
        ):
            return "summary_mismatch"
        return mismatch_reason
    if case.operation in {"get-tree-mrca", "assess-tree-monophyly"}:
        if not _compare_json(
            reference_summary,
            bijux_summary,
            tolerance=case.tolerance,
        ):
            return "summary_mismatch"
        return None
    if case.operation in {"read-tree-set-structure", "write-tree-set-structure"}:
        return _tree_set_structure_mismatch_reason(case, execution_root)
    if case.operation == "tree-consensus":
        mismatch_reason = _consensus_tree_mismatch_reason(case, execution_root)
        if mismatch_reason is None and not _compare_json(
            reference_summary,
            bijux_summary,
            tolerance=case.tolerance,
        ):
            return "summary_mismatch"
        if mismatch_reason is None and not _compare_json(
            reference_rows,
            bijux_rows,
            tolerance=case.tolerance,
        ):
            return "rows_mismatch"
        return mismatch_reason
    if case.operation == "distance-matrix-neighbor-joining":
        mismatch_reason = _neighbor_joining_tree_mismatch_reason(case, execution_root)
        if mismatch_reason is None and not _compare_json(
            reference_summary,
            bijux_summary,
            tolerance=case.tolerance,
        ):
            return "summary_mismatch"
        return mismatch_reason
    if case.operation == "tree-discrete-ancestral-states":
        reference_transition_rows = (
            []
            if reference_summary is None
            else reference_summary.get("transition_rate_rows", [])
        )
        bijux_transition_rows = (
            []
            if bijux_summary is None
            else bijux_summary.get("transition_rate_rows", [])
        )
        reference_summary_without_transition_rows = (
            {}
            if reference_summary is None
            else {
                key: value
                for key, value in reference_summary.items()
                if key != "transition_rate_rows"
            }
        )
        bijux_summary_without_transition_rows = (
            {}
            if bijux_summary is None
            else {
                key: value
                for key, value in bijux_summary.items()
                if key != "transition_rate_rows"
            }
        )
        if not _compare_json(
            reference_summary_without_transition_rows,
            bijux_summary_without_transition_rows,
            tolerance=case.tolerance,
        ):
            return "summary_mismatch"
        if not _transition_rate_rows_match(
            reference_rows=reference_transition_rows,
            bijux_rows=bijux_transition_rows,
            reference_summary=reference_summary,
            bijux_summary=bijux_summary,
            tolerance=(
                case.transition_rate_tolerance
                if case.transition_rate_tolerance is not None
                else case.tolerance
            ),
        ):
            return "transition_rate_rows_mismatch"
        return None
    if not _compare_json(reference_summary, bijux_summary, tolerance=case.tolerance):
        return "summary_mismatch"
    if not _compare_json(reference_rows, bijux_rows, tolerance=case.tolerance):
        return "rows_mismatch"
    if reference_normalized_text != bijux_normalized_text:
        return "normalized_text_mismatch"
    return None


def _apply_expected_status_contract(
    *,
    case: ApeParityCase,
    bijux_error: dict[str, object] | None,
    reference_error: dict[str, object] | None,
    status: str,
    mismatch_reason: str | None,
) -> tuple[str, str | None]:
    contract = _EXPECTED_ERROR_CONTRACTS.get(case.expected_status)
    if contract is None:
        return status, mismatch_reason
    missing_bijux_reason, missing_reference_reason, missing_message_reason = contract
    if bijux_error is None:
        return "failed", missing_bijux_reason
    if reference_error is None:
        return "failed", missing_reference_reason
    if not bijux_error.get("message") or not reference_error.get("message"):
        return "failed", missing_message_reason
    return "passed", None
