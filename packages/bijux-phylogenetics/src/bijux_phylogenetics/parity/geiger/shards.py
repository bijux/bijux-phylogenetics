from __future__ import annotations

from dataclasses import dataclass
import re

from .registry import list_geiger_parity_cases

_NON_ALPHANUMERIC_PATTERN = re.compile(r"[^a-z0-9]+")


@dataclass(frozen=True, slots=True)
class GeigerParityShard:
    """One governed `geiger` parity shard grouped by one governed reference call."""

    shard_id: str
    function_name: str
    python_function_name: str
    model_name: str
    case_ids: tuple[str, ...]

    @property
    def case_count(self) -> int:
        return len(self.case_ids)


def _normalize_shard_token(value: str) -> str:
    return _NON_ALPHANUMERIC_PATTERN.sub("-", value.strip().lower()).strip("-")


def _build_geiger_shard_id(*, function_name: str) -> str:
    return _normalize_shard_token(function_name)


def list_geiger_parity_shards() -> list[GeigerParityShard]:
    """Return governed `geiger` parity shards grouped by one governed reference call."""
    shard_map: dict[tuple[str, str, str], list[str]] = {}
    for case in list_geiger_parity_cases():
        shard_map.setdefault(
            (
                case.python_function_name,
                case.function_name,
                case.model_name,
            ),
            [],
        ).append(case.case_id)
    return [
        GeigerParityShard(
            shard_id=_build_geiger_shard_id(function_name=function_name),
            function_name=function_name,
            python_function_name=python_function_name,
            model_name=model_name,
            case_ids=tuple(case_ids),
        )
        for (python_function_name, function_name, model_name), case_ids in sorted(
            shard_map.items(),
            key=lambda item: (
                item[0][0],
                _normalize_shard_token(item[0][1]),
            ),
        )
    ]


def resolve_geiger_parity_case_ids(
    *,
    case_ids: list[str] | None,
    shard_ids: list[str] | None,
) -> list[str] | None:
    """Resolve one optional case-and-shard selection into concrete case ids."""
    if case_ids is None and shard_ids is None:
        return None
    cases = list_geiger_parity_cases()
    case_order = [case.case_id for case in cases]
    case_index = set(case_order)
    if case_ids is not None:
        missing_case_ids = [case_id for case_id in case_ids if case_id not in case_index]
        if missing_case_ids:
            raise ValueError(
                f"unknown geiger parity case ids: {', '.join(sorted(missing_case_ids))}"
            )
    selected_case_ids = set(case_ids) if case_ids is not None else set(case_order)
    if shard_ids is not None:
        shard_map = {
            shard.shard_id: set(shard.case_ids)
            for shard in list_geiger_parity_shards()
        }
        missing_shard_ids = [
            shard_id for shard_id in shard_ids if shard_id not in shard_map
        ]
        if missing_shard_ids:
            raise ValueError(
                "unknown geiger parity shard ids: "
                + ", ".join(sorted(missing_shard_ids))
            )
        shard_case_ids = set().union(*(shard_map[shard_id] for shard_id in shard_ids))
        selected_case_ids &= shard_case_ids
    if not selected_case_ids:
        raise ValueError("geiger parity selection resolved to zero cases")
    return [case_id for case_id in case_order if case_id in selected_case_ids]


__all__ = [
    "GeigerParityShard",
    "list_geiger_parity_shards",
    "resolve_geiger_parity_case_ids",
]
