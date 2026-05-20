from __future__ import annotations

from importlib import import_module

_PUBLIC_SURFACES = (
    (
        "figure_bundle",
        (
            "AncestralFigurePackageResult",
            "build_ancestral_figure_package",
        ),
    ),
    (
        "visualization",
        (
            "AncestralVisualizationResult",
            "render_ancestral_state_visualization",
        ),
    ),
)

__all__ = [
    name
    for _, names in _PUBLIC_SURFACES
    for name in names
]

_NAME_TO_MODULE = {
    name: module_name
    for module_name, names in _PUBLIC_SURFACES
    for name in names
}


def __getattr__(name: str):
    """Resolve presentation exports lazily from their owning submodules."""
    try:
        module_name = _NAME_TO_MODULE[name]
    except KeyError as error:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from error
    value = getattr(import_module(f"{__name__}.{module_name}"), name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
