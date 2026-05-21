from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .contracts import RabiesMethodSensitivitySlurmTreeRetentionCheckRow
from .shared import (
    _CONFIG_FILENAME,
    _SLURM_OUTPUT_EXPLOSION_SUMMARY_FILENAME,
    _SLURM_STORAGE_CATEGORIES_FILENAME,
    _load_json,
    _read_tsv_rows,
)


@dataclass(frozen=True, slots=True)
class TreeRetentionInputs:
    bundle_root: Path
    config: dict[str, object]
    configured_variant_ids: list[str]
    output_explosion_summary: dict[str, object]
    storage_category_by_id: dict[str, dict[str, str]]
    checks: tuple[RabiesMethodSensitivitySlurmTreeRetentionCheckRow, ...]


def load_tree_retention_inputs(bundle_root: Path) -> TreeRetentionInputs:
    bundle_root = bundle_root.resolve()
    config = _load_json(bundle_root / _CONFIG_FILENAME)
    storage_category_rows = _read_tsv_rows(
        bundle_root / _SLURM_STORAGE_CATEGORIES_FILENAME
    )
    output_explosion_summary = _load_json(
        bundle_root / _SLURM_OUTPUT_EXPLOSION_SUMMARY_FILENAME
    )
    checks: list[RabiesMethodSensitivitySlurmTreeRetentionCheckRow] = []

    configured_variant_ids = sorted(
        str(row["variant_id"]) for row in list(config.get("variants", []))
    )
    storage_category_by_id = {
        str(row["category_id"]): row for row in storage_category_rows
    }
    _add_check(
        checks,
        "storage:category-coverage",
        surface="storage",
        condition={"logs", "outputs", "posterior_samples", "reports", "trees"}
        == set(storage_category_by_id),
        expected=sorted(["logs", "outputs", "posterior_samples", "reports", "trees"]),
        observed=sorted(storage_category_by_id),
        detail="tree-retention policy reads the explicit storage-category surface",
    )

    return TreeRetentionInputs(
        bundle_root=bundle_root,
        config=config,
        configured_variant_ids=configured_variant_ids,
        output_explosion_summary=output_explosion_summary,
        storage_category_by_id=storage_category_by_id,
        checks=tuple(checks),
    )


def _add_check(
    checks: list[RabiesMethodSensitivitySlurmTreeRetentionCheckRow],
    check_id: str,
    *,
    surface: str,
    condition: bool,
    expected: object,
    observed: object,
    detail: str,
) -> None:
    checks.append(
        RabiesMethodSensitivitySlurmTreeRetentionCheckRow(
            check_id=check_id,
            surface=surface,
            status="passed" if condition else "failed",
            expected="" if expected is None else str(expected),
            observed="" if observed is None else str(observed),
            detail=detail,
        )
    )
