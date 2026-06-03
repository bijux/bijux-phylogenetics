from __future__ import annotations

from dataclasses import dataclass
import math

__all__ = [
    "IqtreeBranchSupportLabel",
    "parse_iqtree_branch_support_label",
    "support_fraction",
]


@dataclass(frozen=True, slots=True)
class IqtreeBranchSupportLabel:
    """One parsed IQ-TREE internal support label."""

    raw_label: str
    sh_alrt_support: float | None
    ufboot_support: float | None


def parse_iqtree_branch_support_label(
    value: str | None,
) -> IqtreeBranchSupportLabel | None:
    """Parse one IQ-TREE support label from an internal branch."""
    if value is None:
        return None
    text = value.strip()
    if not text:
        return None
    if "/" not in text:
        try:
            ufboot_support = float(text)
            if not math.isfinite(ufboot_support):
                return None
            return IqtreeBranchSupportLabel(
                raw_label=text,
                sh_alrt_support=None,
                ufboot_support=ufboot_support,
            )
        except ValueError:
            return None
    parts = [part.strip() for part in text.split("/")]
    if len(parts) != 2:
        return None
    try:
        sh_alrt_support = float(parts[0])
        ufboot_support = float(parts[1])
    except ValueError:
        return None
    if not math.isfinite(sh_alrt_support) or not math.isfinite(ufboot_support):
        return None
    return IqtreeBranchSupportLabel(
        raw_label=text,
        sh_alrt_support=sh_alrt_support,
        ufboot_support=ufboot_support,
    )


def support_fraction(value: float | None) -> float | None:
    """Normalize a support value to a 0-1 fraction when expressed as a percent."""
    if value is None:
        return None
    if not math.isfinite(value):
        return None
    return value / 100.0 if value > 1.0 else value
