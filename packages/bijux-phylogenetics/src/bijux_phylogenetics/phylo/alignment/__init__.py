from __future__ import annotations

from importlib import import_module
from typing import Any

_EXPORT_MODULES = ("concatenation", "models", "occupancy", "partitions")


def __getattr__(name: str) -> Any:
    for module_name in _EXPORT_MODULES:
        module = import_module(f".{module_name}", __name__)
        if not hasattr(module, name):
            continue
        value = getattr(module, name)
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    names = set(globals())
    for module_name in _EXPORT_MODULES:
        module = import_module(f".{module_name}", __name__)
        names.update(name for name in dir(module) if not name.startswith("_"))
    return sorted(names)
