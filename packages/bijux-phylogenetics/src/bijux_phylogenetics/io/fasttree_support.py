from __future__ import annotations

from dataclasses import dataclass
import math

__all__ = ["FastTreeBranchSupportLabel", "parse_fasttree_branch_support_label"]


@dataclass(frozen=True, slots=True)
class FastTreeBranchSupportLabel:
    """One parsed FastTree SH-like local support label."""

    raw_label: str
    local_support: float


def parse_fasttree_branch_support_label(
    value: str | None,
) -> FastTreeBranchSupportLabel | None:
    """Parse one FastTree internal support label on the documented 0-1 scale."""
    if value is None:
        return None
    text = value.strip()
    if not text:
        return None
    try:
        local_support = float(text)
    except ValueError:
        return None
    if not math.isfinite(local_support):
        return None
    if local_support < 0.0 or local_support > 1.0:
        return None
    return FastTreeBranchSupportLabel(raw_label=text, local_support=local_support)
