"""Comparative-analysis methods and helpers."""

from __future__ import annotations

from .public_api import EXPORT_MODULES, EXPORT_NAMES


def _install_public_exports() -> None:
    """Expose the comparative root contract from the owned public API modules."""

    for export_module in EXPORT_MODULES:
        globals().update(
            {
                export_name: getattr(export_module, export_name)
                for export_name in export_module.__all__
            }
        )


_install_public_exports()

__all__ = list(EXPORT_NAMES)

del EXPORT_MODULES
del EXPORT_NAMES
del _install_public_exports
