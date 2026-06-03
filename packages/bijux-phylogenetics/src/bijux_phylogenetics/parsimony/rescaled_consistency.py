from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.runtime.errors import ParsimonyAnalysisError

from .consistency import consistency_index
from .matrix import load_parsimony_character_matrix
from .models import (
    FitchCharacterMatrix,
    ParsimonyRescaledConsistencyCharacterIndex,
    ParsimonyRescaledConsistencyIndexReport,
)
from .retention import retention_index

_RESCALED_CONSISTENCY_METHODS = frozenset({"fitch", "acctran", "deltran"})


def rescaled_consistency_index(
    tree: Path,
    matrix: FitchCharacterMatrix | Path,
    *,
    method: str,
) -> ParsimonyRescaledConsistencyIndexReport:
    """Compute rescaled consistency index from tested CI and RI surfaces."""
    resolved_method = _resolve_rescaled_consistency_method(method)
    resolved_matrix = (
        matrix
        if isinstance(matrix, FitchCharacterMatrix)
        else load_parsimony_character_matrix(matrix)
    )
    ci_report = consistency_index(tree, resolved_matrix, method=resolved_method)
    ri_report = retention_index(tree, resolved_matrix, method=resolved_method)
    ri_rows_by_character = {row.character_id: row for row in ri_report.character_rows}

    character_rows: list[ParsimonyRescaledConsistencyCharacterIndex] = []
    for ci_row in ci_report.character_rows:
        try:
            ri_row = ri_rows_by_character[ci_row.character_id]
        except KeyError as error:
            raise ParsimonyAnalysisError(
                "rescaled consistency index requires aligned CI and RI character surfaces",
                code="parsimony_rescaled_consistency_character_mismatch",
                details={"character_id": ci_row.character_id},
            ) from error
        undefined_reason = _combine_undefined_reasons(
            ci_row.undefined_reason,
            ri_row.undefined_reason,
        )
        rc = None
        if ci_row.consistency_index is not None and ri_row.retention_index is not None:
            rc = ci_row.consistency_index * ri_row.retention_index
        character_rows.append(
            ParsimonyRescaledConsistencyCharacterIndex(
                character_id=ci_row.character_id,
                ci=ci_row.consistency_index,
                ri=ri_row.retention_index,
                rc=rc,
                undefined_reason=undefined_reason,
            )
        )

    undefined_reason = _combine_undefined_reasons(
        ci_report.undefined_reason,
        ri_report.undefined_reason,
    )
    rc = None
    if (
        ci_report.consistency_index is not None
        and ri_report.retention_index is not None
    ):
        rc = ci_report.consistency_index * ri_report.retention_index
    return ParsimonyRescaledConsistencyIndexReport(
        algorithm="parsimony-rescaled-consistency-index",
        method=resolved_method,
        tree_path=ci_report.tree_path,
        matrix_path=ci_report.matrix_path,
        taxon_column=ci_report.taxon_column,
        taxon_count=ci_report.taxon_count,
        character_count=ci_report.character_count,
        ci=ci_report.consistency_index,
        ri=ri_report.retention_index,
        rc=rc,
        undefined_reason=undefined_reason,
        character_rows=character_rows,
    )


def _resolve_rescaled_consistency_method(method: str) -> str:
    resolved_method = method.strip().lower()
    if resolved_method not in _RESCALED_CONSISTENCY_METHODS:
        raise ParsimonyAnalysisError(
            "rescaled consistency index currently supports only methods with owned CI and RI surfaces",
            code="parsimony_rescaled_consistency_method_unsupported",
            details={
                "method": method,
                "supported_methods": sorted(_RESCALED_CONSISTENCY_METHODS),
            },
        )
    return resolved_method


def _combine_undefined_reasons(
    ci_reason: str | None,
    ri_reason: str | None,
) -> str | None:
    if ci_reason is None and ri_reason is None:
        return None
    reasons = [reason for reason in (ci_reason, ri_reason) if reason is not None]
    if len(reasons) == 1:
        return reasons[0]
    return "|".join(reasons)
