from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

from bijux_phylogenetics.io.fasta.cleaning import (
    compute_pairwise_sequence_identity_matrix,
    list_alignment_filter_profiles,
)
from bijux_phylogenetics.io.fasta.coding import inspect_coding_alignment
from bijux_phylogenetics.io.fasta.quality import (
    assess_alignment_low_information,
    build_alignment_forensic_report,
    build_alignment_quality_report,
    build_ambiguous_alignment_column_report,
    build_duplicate_sequence_policy_report,
    build_sequence_quality_ranking,
)
from bijux_phylogenetics.io.fasta.records import summarise_fasta
from bijux_phylogenetics.render.html import write_html_report

from ..artifacts import report_sidecar_path, section, write_machine_manifest
from ..ledger import sha256
from ..models import AlignmentReportBuildResult


def render_alignment_report(
    *, alignment_path: Path, out_path: Path
) -> AlignmentReportBuildResult:
    """Build a reviewer-facing alignment-only report."""
    alignment = summarise_fasta(alignment_path)
    alignment_quality = build_alignment_quality_report(alignment_path)
    alignment_low_information = assess_alignment_low_information(alignment_path)
    alignment_duplicate_policy = build_duplicate_sequence_policy_report(alignment_path)
    alignment_ambiguous_columns = build_ambiguous_alignment_column_report(
        alignment_path
    )
    alignment_sequence_ranking = build_sequence_quality_ranking(alignment_path)
    alignment_forensic = build_alignment_forensic_report(alignment_path)
    alignment_coding = (
        inspect_coding_alignment(alignment_path)
        if alignment.inferred_alphabet in {"dna", "rna"}
        else None
    )
    alignment_identity_matrix = compute_pairwise_sequence_identity_matrix(
        alignment_path
    )
    title = "Bijux Alignment Report"
    reviewer_summary = [
        f"alignment quality score: {alignment_quality.quality_score}",
        (
            "alignment suspicious diagnostics: flagged"
            if alignment_quality.suspicious_alignment
            else "alignment suspicious diagnostics: clear"
        ),
        (
            "alignment remains suitable for at least one inference family"
            if any(
                (
                    alignment_forensic.safe_for_distance_analysis,
                    alignment_forensic.safe_for_maximum_likelihood,
                    alignment_forensic.safe_for_bayesian_inference,
                    alignment_forensic.safe_for_coding_analysis,
                )
            )
            else "alignment is currently blocked for the main inference families reviewed here"
        ),
        f"reviewer-facing warnings: {len(alignment_forensic.warnings)}",
    ]
    if alignment_quality.suspicious_alignment:
        reviewer_summary.append(
            "longest concentrated missing-data run: "
            f"{alignment_quality.missing_data_concentration.longest_concentrated_run}"
        )
    limitations = sorted(
        dict.fromkeys([*alignment_forensic.limitations, *alignment_forensic.warnings])
    )
    sections = [
        section("reviewer-summary", reviewer_summary),
        section("alignment-summary", asdict(alignment)),
        section("alignment-quality", asdict(alignment_quality)),
        section("alignment-readiness", asdict(alignment_forensic.readiness)),
        section("alignment-low-information", asdict(alignment_low_information)),
        section("alignment-duplicate-policy", asdict(alignment_duplicate_policy)),
        section("alignment-ambiguous-columns", asdict(alignment_ambiguous_columns)),
        section("alignment-sequence-ranking", asdict(alignment_sequence_ranking)),
        section(
            "alignment-filter-profiles",
            [asdict(profile) for profile in list_alignment_filter_profiles()],
        ),
        section(
            "alignment-suspicious-windows",
            {
                "over_aligned_regions": [
                    asdict(row) for row in alignment_forensic.over_aligned_regions
                ],
                "under_aligned_regions": [
                    asdict(row) for row in alignment_forensic.under_aligned_regions
                ],
            },
        ),
        section("alignment-forensic", asdict(alignment_forensic)),
        *(
            [section("alignment-coding", asdict(alignment_coding))]
            if alignment_coding is not None
            else []
        ),
        section("alignment-identity-matrix", asdict(alignment_identity_matrix)),
        section("limitations", limitations),
    ]
    machine_manifest = {
        "report_kind": "alignment",
        "title": title,
        "input_paths": [str(alignment_path)],
        "input_checksums": {str(alignment_path): sha256(alignment_path)},
        "sections": [name for name, _ in sections],
        "metrics": {
            "sequence_count": alignment.sequence_count,
            "alignment_length": alignment.alignment_length,
            "quality_score": alignment_quality.quality_score,
        },
        "reviewer_summary": reviewer_summary,
        "limitations": limitations,
    }
    machine_manifest_path = write_machine_manifest(
        report_sidecar_path(out_path), machine_manifest
    )
    write_html_report(
        title=title,
        sections=sections,
        out_path=out_path,
        embedded_json=machine_manifest,
    )
    return AlignmentReportBuildResult(
        output_path=out_path,
        machine_manifest_path=machine_manifest_path,
        report_kind="alignment",
        title=title,
        alignment=alignment,
        alignment_quality=alignment_quality,
        alignment_forensic=alignment_forensic,
        alignment_low_information=alignment_low_information,
        alignment_duplicate_policy=alignment_duplicate_policy,
        alignment_ambiguous_columns=alignment_ambiguous_columns,
        alignment_sequence_ranking=alignment_sequence_ranking,
        alignment_coding=alignment_coding,
        alignment_identity_matrix=alignment_identity_matrix,
        machine_manifest=machine_manifest,
    )
