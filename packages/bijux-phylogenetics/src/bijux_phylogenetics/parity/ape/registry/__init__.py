from __future__ import annotations

from pathlib import Path

from .comparative_cases import build_comparative_case_catalog
from .fixtures import ApeParityFixtureResolver
from .models import ApeParityCase
from .selection import select_cases
from .sequence_cases import build_sequence_case_catalog
from .serialization import write_case_file
from .tree_cases import build_tree_case_catalog


def list_ape_parity_cases(fixtures_root: Path | None = None) -> list[ApeParityCase]:
    """Return the governed live `ape` parity cases."""
    resolver = ApeParityFixtureResolver(fixtures_root)

    return [
        *build_tree_case_catalog(resolver),
        *build_comparative_case_catalog(resolver),
        *build_sequence_case_catalog(resolver),
    ]


def _selected_cases(
    *,
    case_ids: list[str] | None,
    fixtures_root: Path | None = None,
) -> list[ApeParityCase]:
    return select_cases(case_ids=case_ids, fixtures_root=fixtures_root)


def _write_case_file(path: Path, case: ApeParityCase) -> Path:
    return write_case_file(path, case)
