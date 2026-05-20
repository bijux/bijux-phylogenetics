from __future__ import annotations

from ..registry import PhytoolsParityCase
from .comparative_payloads import build_comparative_case_payload
from .continuous_payloads import build_continuous_case_payload
from .discrete_payloads import build_discrete_case_payload
from .signal_payloads import build_signal_case_payload


def build_bijux_case_payload(
    case: PhytoolsParityCase,
) -> tuple[dict[str, object], list[dict[str, object]] | None]:
    tree_path = case.input_fixtures[0]
    traits_path = case.input_fixtures[1] if len(case.input_fixtures) > 1 else None
    comparative_payload = build_comparative_case_payload(
        case,
        tree_path=tree_path,
        traits_path=traits_path,
    )
    if comparative_payload is not None:
        return comparative_payload
    discrete_payload = build_discrete_case_payload(
        case,
        tree_path=tree_path,
        traits_path=traits_path,
    )
    if discrete_payload is not None:
        return discrete_payload
    continuous_payload = build_continuous_case_payload(
        case,
        tree_path=tree_path,
        traits_path=traits_path,
    )
    if continuous_payload is not None:
        return continuous_payload
    signal_payload = build_signal_case_payload(
        case,
        tree_path=tree_path,
        traits_path=traits_path,
    )
    if signal_payload is not None:
        return signal_payload
    raise ValueError(f"unsupported phytools parity operation: {case.operation}")
