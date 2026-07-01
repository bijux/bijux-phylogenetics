from __future__ import annotations

from bijux_phylogenetics.datasets.shared_fixtures.geiger_discrete import (
    list_shared_geiger_discrete_fixtures,
)
from bijux_phylogenetics.parity.geiger.registry import list_geiger_parity_cases
from bijux_phylogenetics.parity.geiger.registry.models import GeigerParityCase

from .contracts import (
    SharedDiscreteFixtureCatalogEntry,
    SharedDiscreteFixtureCatalogReport,
)


def build_shared_discrete_fixture_catalog() -> SharedDiscreteFixtureCatalogReport:
    """Summarize the governed discrete fixtures consumed by live `geiger` parity."""

    live_cases = [
        case
        for case in list_geiger_parity_cases()
        if case.function_name.startswith("geiger::fitDiscrete")
    ]
    case_groups: dict[str, list[GeigerParityCase]] = {}
    for case in live_cases:
        case_groups.setdefault(case.fixture_id, []).append(case)

    fixtures = sorted(
        [
            fixture
            for fixture in list_shared_geiger_discrete_fixtures()
            if fixture.fixture_id in case_groups
        ],
        key=lambda fixture: fixture.fixture_id,
    )
    rows = []
    for fixture in fixtures:
        grouped_cases = case_groups.get(fixture.fixture_id, [])
        rows.append(
            SharedDiscreteFixtureCatalogEntry(
                fixture_id=fixture.fixture_id,
                tree_fixture_id=fixture.tree_fixture.fixture_id,
                trait_table_fixture_id=fixture.trait_table_fixture_id,
                trait_name=fixture.trait_name,
                taxon_column=fixture.taxon_column,
                trait_kind=fixture.trait_kind,
                taxon_count=fixture.tree_fixture.tip_count,
                trait_row_count=fixture.trait_table_fixture.row_count,
                supported_model_names=list(fixture.supported_model_names),
                covered_model_names=sorted({case.model_name for case in grouped_cases}),
                live_case_count=len(grouped_cases),
                validation_expectation=fixture.validation_expectation,
                geiger_reference_expectation=fixture.geiger_reference_expectation,
                feature_tags=list(fixture.feature_tags),
                transition_matrix_metadata=fixture.transition_matrix_metadata,
                notes=fixture.notes,
            )
        )
    return SharedDiscreteFixtureCatalogReport(
        goal_id=253,
        fixture_count=len(rows),
        live_case_count=len(live_cases),
        fixtures=rows,
        limitations=[
            "the catalog covers the governed discrete fixtures consumed by the live fitDiscrete parity harness; it does not enumerate governed fixtures that are not yet part of that live surface",
            "the catalog proves which governed fixtures are consumed by the live fitDiscrete parity harness; it does not by itself prove parity passed for every covered model",
            "fixture membership is governed by the shared fixture catalog and live geiger registry, so coverage changes only when one of those owned sources changes",
        ],
    )
