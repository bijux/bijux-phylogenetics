from __future__ import annotations

from ..models import (
    RabiesMethodSensitivityCladeRow,
    RabiesMethodSensitivityVariantRun,
)


def _aggregate_clades(
    variant_runs: list[RabiesMethodSensitivityVariantRun], *, stable_only: bool
) -> list[RabiesMethodSensitivityCladeRow]:
    counts: dict[tuple[str, str], tuple[str, str, int]] = {}
    for variant in variant_runs:
        for row in variant.inference_comparison.conclusion_rows:
            is_stable = row.conclusion_class == "stable_clade"
            if is_stable != stable_only:
                continue
            key = (row.split_id, row.conclusion_class)
            evidence_class, detail, count = counts.get(
                key, (row.evidence_class, row.detail, 0)
            )
            counts[key] = (evidence_class, detail, count + 1)
    return [
        RabiesMethodSensitivityCladeRow(
            split_id=split_id,
            conclusion_class=conclusion_class,
            evidence_class=evidence_class,
            occurrence_count=count,
            variant_count=len(variant_runs),
            detail=detail,
        )
        for (split_id, conclusion_class), (evidence_class, detail, count) in sorted(
            counts.items()
        )
    ]
