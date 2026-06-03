"""Owned public export modules for comparative root imports."""

from __future__ import annotations

from . import (
    clades,
    continuous_traits,
    diversification,
    recovery_validation,
    regression_pgls,
    reporting,
    signal_modes,
)

EXPORT_MODULES = (
    continuous_traits,
    clades,
    diversification,
    recovery_validation,
    regression_pgls,
    reporting,
    signal_modes,
)

EXPORT_NAMES = tuple(
    export_name
    for export_module in EXPORT_MODULES
    for export_name in export_module.__all__
)

__all__ = [
    "EXPORT_MODULES",
    "EXPORT_NAMES",
    "clades",
    "continuous_traits",
    "diversification",
    "recovery_validation",
    "regression_pgls",
    "reporting",
    "signal_modes",
]
