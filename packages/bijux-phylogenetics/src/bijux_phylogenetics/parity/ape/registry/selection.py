from __future__ import annotations

from pathlib import Path

from .models import ApeParityCase


def build_case_lookup(fixtures_root: Path | None = None) -> dict[str, ApeParityCase]:
    from . import list_ape_parity_cases

    return {case.case_id: case for case in list_ape_parity_cases(fixtures_root)}


def select_cases(
    *,
    case_ids: list[str] | None,
    fixtures_root: Path | None = None,
) -> list[ApeParityCase]:
    cases = build_case_lookup(fixtures_root)
    if case_ids is None:
        return list(cases.values())
    selected: list[ApeParityCase] = []
    for case_id in case_ids:
        try:
            selected.append(cases[case_id])
        except KeyError as error:
            supported = ", ".join(sorted(cases))
            raise ValueError(
                f"unsupported ape parity case '{case_id}'; expected one of: {supported}"
            ) from error
    return selected
