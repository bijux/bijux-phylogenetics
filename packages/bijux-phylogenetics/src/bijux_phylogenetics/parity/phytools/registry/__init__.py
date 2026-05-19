from __future__ import annotations

from .comparative_cases import build_comparative_cases
from .continuous_cases import build_continuous_cases
from .fixtures import build_phytools_registry_fixture_catalog
from .models import PhytoolsParityCase
from .discrete_history_cases import build_discrete_history_cases
from .discrete_model_cases import build_discrete_model_cases
from .signal_cases import build_signal_cases
from .stochastic_map_cases import build_stochastic_map_cases


def list_phytools_parity_cases() -> list[PhytoolsParityCase]:
    """Return the governed live `phytools` parity cases."""
    fixture_catalog = build_phytools_registry_fixture_catalog()
    return [
        *build_signal_cases(fixture_catalog),
        *build_discrete_model_cases(fixture_catalog),
        *build_stochastic_map_cases(fixture_catalog),
        *build_discrete_history_cases(fixture_catalog),
        *build_continuous_cases(fixture_catalog),
        *build_comparative_cases(fixture_catalog),
    ]
