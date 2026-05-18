from __future__ import annotations

import bijux_phylogenetics.parity as parity_api
from bijux_phylogenetics.parity import (
    list_geiger_parity_cases,
    run_geiger_parity_cases,
    write_geiger_parity_observation_table,
    write_geiger_parity_summary_table,
)


def test_public_runtime_exports_include_geiger_parity_surface() -> None:
    assert parity_api.list_geiger_parity_cases is list_geiger_parity_cases
    assert parity_api.run_geiger_parity_cases is run_geiger_parity_cases
    assert (
        parity_api.write_geiger_parity_summary_table
        is write_geiger_parity_summary_table
    )
    assert (
        parity_api.write_geiger_parity_observation_table
        is write_geiger_parity_observation_table
    )
