from __future__ import annotations

from .models import PhytoolsParityCase


def build_case_lookup() -> dict[str, PhytoolsParityCase]:
    from . import list_phytools_parity_cases

    return {case.case_id: case for case in list_phytools_parity_cases()}


def select_cases(*, case_ids: list[str] | None) -> list[PhytoolsParityCase]:
    cases = build_case_lookup()
    if case_ids is None:
        return list(cases.values())
    missing = [case_id for case_id in case_ids if case_id not in cases]
    if missing:
        missing_text = ", ".join(sorted(missing))
        raise ValueError(f"unknown phytools parity case id(s): {missing_text}")
    return [cases[case_id] for case_id in case_ids]
