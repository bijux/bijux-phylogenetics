from __future__ import annotations

import json
from pathlib import Path

from .fixtures import ApeParityFixtureResolver
from .models import ApeParityCase
from .tree_cases import build_tree_case_catalog
from .comparative_cases import build_comparative_case_catalog
from .sequence_cases import build_sequence_case_catalog


def list_ape_parity_cases(fixtures_root: Path | None = None) -> list[ApeParityCase]:
    """Return the governed live `ape` parity cases."""
    resolver = ApeParityFixtureResolver(fixtures_root)

    return [
        *build_tree_case_catalog(resolver),
        *build_comparative_case_catalog(resolver),
        *build_sequence_case_catalog(resolver),
    ]


def _build_case_lookup(fixtures_root: Path | None = None) -> dict[str, ApeParityCase]:
    return {case.case_id: case for case in list_ape_parity_cases(fixtures_root)}


def _selected_cases(
    *,
    case_ids: list[str] | None,
    fixtures_root: Path | None = None,
) -> list[ApeParityCase]:
    cases = _build_case_lookup(fixtures_root)
    if case_ids is None:
        return list(cases.values())
    selected: list[ApeParityCase] = []
    for case_id in case_ids:
        try:
            selected.append(cases[case_id])
        except KeyError as error:
            supported = ", ".join(sorted(cases))
            raise ValueError(
                f"unsupported ape parity case '{case_id}'; expected one of: {supported}"
            ) from error
    return selected


def _write_case_file(path: Path, case: ApeParityCase) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            {
                "case_id": case.case_id,
                "fixture_kind": case.fixture_kind,
                "fixture_id": case.fixture_id,
                "function_name": case.function_name,
                "operation": case.operation,
                "input_fixture": str(case.input_fixture),
                "tolerance": case.tolerance,
                "expected_status": case.expected_status,
                "pairwise_deletion": case.pairwise_deletion,
                "distance_model": case.distance_model,
                "genetic_code_id": case.genetic_code_id,
                "outgroup_taxa": list(case.outgroup_taxa),
                "excluded_taxa": list(case.excluded_taxa),
                "requested_taxa": list(case.requested_taxa),
                "node_id": case.node_id,
                "mrca_taxa": list(case.mrca_taxa),
                "monophyly_reroot": case.monophyly_reroot,
                "ultrametric_option": case.ultrametric_option,
                "rf_mode": case.rf_mode,
                "consensus_method": case.consensus_method,
                "reference_tree_path": (
                    None
                    if case.reference_tree_path is None
                    else str(case.reference_tree_path)
                ),
                "ancestral_model": case.ancestral_model,
                "trait_fixture_id": case.trait_fixture_id,
                "trait_table_path": (
                    None
                    if case.trait_table_path is None
                    else str(case.trait_table_path)
                ),
                "trait_name": case.trait_name,
                "trait_taxon_column": case.trait_taxon_column,
                "transition_rate_tolerance": case.transition_rate_tolerance,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return path

