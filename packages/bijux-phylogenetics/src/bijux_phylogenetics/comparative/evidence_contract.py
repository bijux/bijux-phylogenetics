from __future__ import annotations

from importlib import import_module
from typing import Any

SUPPORTED_EVIDENCE_API_MODULES = ("bijux_phylogenetics.comparative",)
SUPPORTED_EVIDENCE_API_LOCATORS = (
    "bijux_phylogenetics.comparative:inspect_pgls_inputs",
    "bijux_phylogenetics.comparative:run_pgls",
    "bijux_phylogenetics.comparative:run_pgls_multiple_testing",
    "bijux_phylogenetics.comparative:compute_phylogenetic_signal_test",
    "bijux_phylogenetics.comparative:estimate_pagels_lambda",
    "bijux_phylogenetics.comparative:audit_ou_identifiability_reference_examples",
    "bijux_phylogenetics.comparative:compare_brownian_and_ou_models",
    "bijux_phylogenetics.comparative:fit_brownian_motion_model",
    "bijux_phylogenetics.comparative:fit_ornstein_uhlenbeck_model",
    "bijux_phylogenetics.comparative:transform_tree_for_evolutionary_mode",
)


def resolve_supported_evidence_api(locator: str) -> Any:
    """Resolve one governed evidence-consumer locator against the public runtime API."""
    if locator not in SUPPORTED_EVIDENCE_API_LOCATORS:
        raise KeyError(locator)
    module_name, export_name = locator.split(":", maxsplit=1)
    module = import_module(module_name)
    return getattr(module, export_name)
