"""Discrete trait and history simulation workflows."""

from __future__ import annotations

from importlib import import_module
from typing import TYPE_CHECKING, Any

_PUBLIC_NAME_TO_MODULE = {
    "simulate_discrete_histories": ".histories",
    "simulate_discrete_traits": ".traits",
    "write_discrete_history_branch_truth_table": ".histories",
    "write_discrete_history_event_table": ".histories",
    "write_discrete_history_node_truth_table": ".histories",
    "write_discrete_history_segment_table": ".histories",
    "write_discrete_history_summary_table": ".histories",
    "write_discrete_history_tip_truth_table": ".histories",
    "write_discrete_trait_table": ".traits",
}

__all__ = sorted(_PUBLIC_NAME_TO_MODULE)

if TYPE_CHECKING:
    from .histories import (
        simulate_discrete_histories,
        write_discrete_history_branch_truth_table,
        write_discrete_history_event_table,
        write_discrete_history_node_truth_table,
        write_discrete_history_segment_table,
        write_discrete_history_summary_table,
        write_discrete_history_tip_truth_table,
    )
    from .traits import (
        simulate_discrete_traits,
        write_discrete_trait_table,
    )


def __getattr__(name: str) -> Any:
    if name not in _PUBLIC_NAME_TO_MODULE:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module = import_module(_PUBLIC_NAME_TO_MODULE[name], __name__)
    value = getattr(module, name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
