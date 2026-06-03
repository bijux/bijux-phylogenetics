from __future__ import annotations

_MODEL_ALIAS_TO_INTERNAL = {
    "er": "equal-rates",
    "sym": "symmetric",
    "ard": "all-rates-different",
}


def resolve_internal_model(model: str) -> str:
    try:
        return _MODEL_ALIAS_TO_INTERNAL[model]
    except KeyError as error:
        raise ValueError(
            "ecological niche transition model must be one of: er, sym, ard"
        ) from error


def transition_certainty_class(
    *,
    changed: bool,
    overlapping_niches: list[str],
    parent_niche_set: list[str],
    child_niche_set: list[str],
) -> str:
    if not changed:
        return "no_transition"
    if overlapping_niches:
        return "uncertain_transition"
    if len(parent_niche_set) == 1 and len(child_niche_set) == 1:
        return "certain_transition"
    return "uncertain_transition"


def format_optional_float(value: float | None) -> str:
    if value is None:
        return ""
    return str(stable_float(value))


def stable_float(value: float) -> float:
    normalized = round(float(value), 8)
    return 0.0 if normalized == -0.0 else normalized


def node_signature(node) -> str:
    if node.is_leaf():
        return node.name
    return "|".join(sorted(leaf.name for leaf in node.iter_leaves()))
