from __future__ import annotations

from pathlib import Path

from .comparative_cases import build_comparative_cases
from .continuous_cases import build_continuous_cases
from .discrete_history_cases import build_discrete_history_cases
from .discrete_model_cases import build_discrete_model_cases
from .fixtures import build_phytools_registry_fixture_catalog
from .models import PhytoolsParityCase
from .selection import select_cases
from .serialization import write_case_file
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


def _selected_cases(case_ids: list[str] | None) -> list[PhytoolsParityCase]:
    return select_cases(case_ids=case_ids)


def _write_case_file(path: Path, case: PhytoolsParityCase) -> Path:
    return write_case_file(path, case)
