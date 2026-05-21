from __future__ import annotations

import bijux_phylogenetics.parity as parity_api
from bijux_phylogenetics.parity import (
    GeneratedGeigerParityReport,
    build_generated_geiger_parity_report,
    list_geiger_parity_cases,
    run_geiger_parity_cases,
    write_geiger_boundary_warning_table,
    write_geiger_likelihood_policy_table,
    write_geiger_model_confidence_table,
    write_geiger_optimizer_triage_table,
    write_geiger_parameterization_registry_table,
    write_geiger_parity_observation_table,
    write_geiger_parity_summary_table,
    write_generated_geiger_parity_report_json,
    write_generated_geiger_parity_report_markdown,
)


def test_public_runtime_exports_include_geiger_parity_surface() -> None:
    assert parity_api.GeneratedGeigerParityReport is GeneratedGeigerParityReport
    assert parity_api.list_geiger_parity_cases is list_geiger_parity_cases
    assert (
        parity_api.build_generated_geiger_parity_report
        is build_generated_geiger_parity_report
    )
    assert parity_api.run_geiger_parity_cases is run_geiger_parity_cases
    assert (
        parity_api.write_geiger_parity_summary_table
        is write_geiger_parity_summary_table
    )
    assert (
        parity_api.write_geiger_parity_observation_table
        is write_geiger_parity_observation_table
    )
    assert (
        parity_api.write_geiger_optimizer_triage_table
        is write_geiger_optimizer_triage_table
    )
    assert (
        parity_api.write_geiger_boundary_warning_table
        is write_geiger_boundary_warning_table
    )
    assert (
        parity_api.write_geiger_model_confidence_table
        is write_geiger_model_confidence_table
    )
    assert (
        parity_api.write_geiger_likelihood_policy_table
        is write_geiger_likelihood_policy_table
    )
    assert (
        parity_api.write_geiger_parameterization_registry_table
        is write_geiger_parameterization_registry_table
    )
    assert (
        parity_api.write_generated_geiger_parity_report_json
        is write_generated_geiger_parity_report_json
    )
    assert (
        parity_api.write_generated_geiger_parity_report_markdown
        is write_generated_geiger_parity_report_markdown
    )
