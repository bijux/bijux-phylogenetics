from __future__ import annotations

from dataclasses import asdict
from typing import Any

from ..artifacts import section
from .linked_evidence import TaxonLinkedEvidence

TAXON_REPORT_TITLE = "Bijux Taxon Audit Report"


def build_taxon_reviewer_summary(
    *, audit: Any, linked_evidence: TaxonLinkedEvidence
) -> list[str]:
    """Build the reviewer summary for a taxon audit report."""
    reviewer_summary = [
        f"taxon audit status: {audit.status}",
        f"tree tip count: {audit.tree_tip_count}",
        *audit.summary,
    ]
    if linked_evidence.taxon_crosswalk is not None:
        reviewer_summary.append(
            f"crosswalk rows: {len(linked_evidence.taxon_crosswalk.rows)} across linked dataset surfaces"
        )
    if linked_evidence.taxon_exclusions is not None:
        reviewer_summary.append(
            f"excluded taxa with explicit causes: {len(linked_evidence.taxon_exclusions.rows)}"
        )
    if linked_evidence.taxon_workflow_loss is not None:
        reviewer_summary.append(
            f"workflow loss stages observed: {len(linked_evidence.taxon_workflow_loss.loss_stage_counts)}"
        )
    if linked_evidence.taxon_stability is not None:
        reviewer_summary.append(
            f"unstable taxa across linked sources: {len(linked_evidence.taxon_stability.unstable_taxa)}"
        )
    return reviewer_summary


def build_taxon_sections(
    *,
    audit: Any,
    linked_evidence: TaxonLinkedEvidence,
    reviewer_summary: list[str],
    limitations: list[str],
) -> list[tuple[str, Any]]:
    """Build taxon report sections for reviewer-facing HTML and sidecars."""
    taxon_workflow_loss = linked_evidence.taxon_workflow_loss
    return [
        section("reviewer-summary", reviewer_summary),
        section("taxon-audit", asdict(audit)),
        section("taxon-identity", asdict(audit.identity_audit)),
        section("taxon-safety", asdict(audit.safety_report)),
        section("taxon-namespaces", asdict(audit.namespace_report)),
        section("taxon-rank-consistency", asdict(audit.rank_consistency)),
        *(
            [section("taxon-synonyms", asdict(audit.synonym_audit))]
            if audit.synonym_audit is not None
            else []
        ),
        section("taxon-duplicate-identities", asdict(audit.duplicate_identities)),
        section("taxon-mapping-conflicts", asdict(audit.mapping_conflicts)),
        *(
            [section("taxon-accepted-names", asdict(audit.accepted_name_export))]
            if audit.accepted_name_export is not None
            else []
        ),
        *(
            [section("taxon-crosswalk", asdict(linked_evidence.taxon_crosswalk))]
            if linked_evidence.taxon_crosswalk is not None
            else []
        ),
        *(
            [section("taxon-exclusions", asdict(linked_evidence.taxon_exclusions))]
            if linked_evidence.taxon_exclusions is not None
            else []
        ),
        *(
            [section("taxon-loss", asdict(taxon_workflow_loss))]
            if taxon_workflow_loss is not None
            else []
        ),
        *(
            [
                section(
                    "taxon-loss-events",
                    [
                        {
                            "taxon": row.taxon,
                            "first_loss_stage": row.first_loss_stage,
                            "loss_events": [asdict(event) for event in row.loss_events],
                        }
                        for row in taxon_workflow_loss.rows
                        if row.loss_events
                    ],
                )
            ]
            if taxon_workflow_loss is not None
            else []
        ),
        *(
            [section("taxon-stability", asdict(linked_evidence.taxon_stability))]
            if linked_evidence.taxon_stability is not None
            else []
        ),
        section("limitations", limitations),
    ]
