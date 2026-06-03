from __future__ import annotations

import json
from pathlib import Path

from .models import ApeParityCase


def write_case_file(path: Path, case: ApeParityCase) -> Path:
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
