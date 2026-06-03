from __future__ import annotations

import bijux_phylogenetics.comparative as comparative_api
from bijux_phylogenetics.comparative import (
    summarize_comparative_covariance_audit,
    write_comparative_covariance_audit_candidate_table,
    write_comparative_covariance_audit_excluded_taxa_table,
    write_comparative_covariance_audit_summary_table,
)


def test_comparative_covariance_audit_surfaces_export_publicly() -> None:
    assert (
        comparative_api.summarize_comparative_covariance_audit
        is summarize_comparative_covariance_audit
    )
    assert (
        comparative_api.write_comparative_covariance_audit_summary_table
        is write_comparative_covariance_audit_summary_table
    )
    assert (
        comparative_api.write_comparative_covariance_audit_candidate_table
        is write_comparative_covariance_audit_candidate_table
    )
    assert (
        comparative_api.write_comparative_covariance_audit_excluded_taxa_table
        is write_comparative_covariance_audit_excluded_taxa_table
    )
