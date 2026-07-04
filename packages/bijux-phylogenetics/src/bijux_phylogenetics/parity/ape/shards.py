from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .registry import list_ape_parity_cases


@dataclass(frozen=True, slots=True)
class ApeParityShard:
    """One governed `ape` parity shard grouped by owned function surface."""

    shard_id: str
    function_name: str
    python_function_name: str
    case_ids: tuple[str, ...]

    @property
    def case_count(self) -> int:
        return len(self.case_ids)


def list_ape_parity_shards(fixtures_root: Path | None = None) -> list[ApeParityShard]:
    """Return governed `ape` parity shards grouped by owned function surface."""
    cases = list_ape_parity_cases(fixtures_root=fixtures_root)
    shard_map: dict[tuple[str, str], list[str]] = {}
    for case in cases:
        shard_map.setdefault(
            (case.python_function_name, case.function_name),
            [],
        ).append(case.case_id)
    return [
        ApeParityShard(
            shard_id=python_function_name,
            function_name=function_name,
            python_function_name=python_function_name,
            case_ids=tuple(case_ids),
        )
        for (python_function_name, function_name), case_ids in sorted(
            shard_map.items(),
            key=lambda item: item[0][0],
        )
    ]


def resolve_ape_parity_case_ids(
    *,
    case_ids: list[str] | None,
    shard_ids: list[str] | None,
    fixtures_root: Path | None = None,
) -> list[str] | None:
    """Resolve one optional case-and-shard selection into concrete case ids."""
    if case_ids is None and shard_ids is None:
        return None
    cases = list_ape_parity_cases(fixtures_root=fixtures_root)
    case_order = [case.case_id for case in cases]
    case_index = set(case_order)
    if case_ids is not None:
        missing_case_ids = [
            case_id for case_id in case_ids if case_id not in case_index
        ]
        if missing_case_ids:
            raise ValueError(
                f"unknown ape parity case ids: {', '.join(sorted(missing_case_ids))}"
            )
    selected_case_ids = set(case_ids) if case_ids is not None else set(case_order)
    if shard_ids is not None:
        shard_map = {
            shard.shard_id: set(shard.case_ids)
            for shard in list_ape_parity_shards(fixtures_root=fixtures_root)
        }
        missing_shard_ids = [
            shard_id for shard_id in shard_ids if shard_id not in shard_map
        ]
        if missing_shard_ids:
            raise ValueError(
                f"unknown ape parity shard ids: {', '.join(sorted(missing_shard_ids))}"
            )
        shard_case_ids = set().union(*(shard_map[shard_id] for shard_id in shard_ids))
        selected_case_ids &= shard_case_ids
    if not selected_case_ids:
        raise ValueError("ape parity selection resolved to zero cases")
    return [case_id for case_id in case_order if case_id in selected_case_ids]


__all__ = [
    "ApeParityShard",
    "list_ape_parity_shards",
    "resolve_ape_parity_case_ids",
]
