from __future__ import annotations

from .continuous_comparison_cases import build_continuous_comparison_cases
from .continuous_fit_cases import build_continuous_fit_cases
from .discrete_fit_cases import build_discrete_fit_cases
from .fixtures import build_geiger_registry_fixture_catalog
from .models import GeigerParityCase


def list_geiger_parity_cases() -> list[GeigerParityCase]:
    """Return the governed live `geiger` parity cases."""

    fixture_catalog = build_geiger_registry_fixture_catalog()
    return [
        *build_continuous_fit_cases(fixture_catalog),
        *build_continuous_comparison_cases(fixture_catalog),
        *build_discrete_fit_cases(fixture_catalog),
    ]
