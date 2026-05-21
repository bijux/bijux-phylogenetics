"""Comparative-analysis methods and helpers."""

from __future__ import annotations

from .public_api import EXPORT_MODULES as _EXPORT_MODULES
from .public_api import EXPORT_NAMES as _EXPORT_NAMES

for _export_module in _EXPORT_MODULES:
    globals().update(
        {
            export_name: getattr(_export_module, export_name)
            for export_name in _export_module.__all__
        }
    )

__all__ = list(_EXPORT_NAMES)

del _export_module
del _EXPORT_MODULES
del _EXPORT_NAMES
