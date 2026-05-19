from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class ApeParityCase:
    """One governed live `ape` parity case."""

    case_id: str
    fixture_kind: str
    fixture_id: str
    function_name: str
    python_function_name: str
    operation: str
    input_fixture: Path
    tolerance: float
    expected_status: str = "ok"
    pairwise_deletion: bool | None = None
    distance_model: str | None = None
    genetic_code_id: int | None = None
    outgroup_taxa: tuple[str, ...] = ()
    excluded_taxa: tuple[str, ...] = ()
    requested_taxa: tuple[str, ...] = ()
    node_id: int | None = None
    mrca_taxa: tuple[str, ...] = ()
    monophyly_reroot: bool | None = None
    ultrametric_option: int | None = None
    rf_mode: str | None = None
    consensus_method: str | None = None
    reference_tree_path: Path | None = None
    ancestral_model: str | None = None
    trait_fixture_id: str | None = None
    trait_table_path: Path | None = None
    trait_name: str | None = None
    trait_taxon_column: str | None = None
    transition_rate_tolerance: float | None = None
