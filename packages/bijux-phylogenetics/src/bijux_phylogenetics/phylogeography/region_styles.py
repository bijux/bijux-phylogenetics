from __future__ import annotations

from collections.abc import Iterable

_GEOGRAPHIC_STATE_PALETTE = (
    "#0f766e",
    "#1d4ed8",
    "#c2410c",
    "#7c3aed",
    "#b91c1c",
    "#047857",
    "#a16207",
    "#0f172a",
)

_TRANSITION_SUPPORT_COLORS = {
    "strong": "#166534",
    "moderate": "#b45309",
    "weak": "#7f1d1d",
}

_TRANSITION_SUPPORT_DETAILS = {
    "strong": "strongly supported geographic transition",
    "moderate": "moderately supported geographic transition",
    "weak": "weakly supported geographic transition",
}


def build_geographic_state_color_map(
    state_labels: Iterable[str],
) -> dict[str, str]:
    """Assign stable publication colors to geographic state labels."""
    ordered = sorted({label for label in state_labels if label})
    return {
        label: _GEOGRAPHIC_STATE_PALETTE[index % len(_GEOGRAPHIC_STATE_PALETTE)]
        for index, label in enumerate(ordered)
    }


def geographic_transition_support_colors() -> dict[str, str]:
    """Return stable publication colors for transition-support classes."""
    return dict(_TRANSITION_SUPPORT_COLORS)


def geographic_transition_support_details() -> dict[str, str]:
    """Return reviewer-facing descriptions for transition-support classes."""
    return dict(_TRANSITION_SUPPORT_DETAILS)
