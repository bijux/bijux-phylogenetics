from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.phylo.alignment import AlignmentForensicReport

from ..coding import inspect_coding_alignment
from ..records import summarise_fasta
from .quality_report import build_alignment_quality_report
from .readiness import summarize_alignment_readiness
from .sequence_review import (
    build_duplicate_sequence_policy_report,
    build_sequence_quality_ranking,
)
from .site_diagnostics import (
    assess_alignment_low_information,
    build_ambiguous_alignment_column_report,
)
from .window_diagnostics import (
    detect_over_aligned_regions,
    detect_under_aligned_regions,
)


def build_alignment_forensic_report(path: Path) -> AlignmentForensicReport:
    """Integrate alignment quality, readiness, coding, and suspicious-region diagnostics."""
    quality = build_alignment_quality_report(path)
    readiness = summarize_alignment_readiness(path)
    summary = summarise_fasta(path)
    low_information = assess_alignment_low_information(path)
    coding = (
        inspect_coding_alignment(path)
        if summary.inferred_alphabet in {"dna", "rna"}
        else None
    )
    duplicate_policy = build_duplicate_sequence_policy_report(path)
    ambiguous_columns = build_ambiguous_alignment_column_report(path)
    sequence_ranking = build_sequence_quality_ranking(path)
    over_aligned = detect_over_aligned_regions(path)
    under_aligned = detect_under_aligned_regions(path)
    method_by_name = {method.analysis: method for method in readiness.methods}
    warnings = list(
        dict.fromkeys(
            [
                *quality.warnings,
                *readiness.warnings,
                *(
                    ["alignment mixes coding-like and noncoding-like sequences"]
                    if coding is not None and coding.mixed_coding_signals
                    else []
                ),
            ]
        )
    )
    limitations = [
        "alignment forensics summarize data quality and readiness but do not replace model-based tree inference validation",
        "publication readiness still depends on explicit method reporting, biological context, and reviewer inspection",
    ]
    return AlignmentForensicReport(
        path=path,
        quality=quality,
        readiness=readiness,
        low_information=low_information,
        coding=coding,
        duplicate_policy=duplicate_policy,
        ambiguous_columns=ambiguous_columns,
        sequence_ranking=sequence_ranking,
        over_aligned_regions=over_aligned,
        under_aligned_regions=under_aligned,
        safe_for_distance_analysis=method_by_name["distance"].ready,
        safe_for_maximum_likelihood=method_by_name["maximum_likelihood"].ready,
        safe_for_bayesian_inference=method_by_name["bayesian"].ready,
        safe_for_coding_analysis=method_by_name["coding"].ready,
        safe_for_publication=(
            quality.quality_score >= 75.0
            and not quality.suspicious_alignment
            and not warnings
        ),
        warnings=warnings,
        limitations=limitations,
    )
