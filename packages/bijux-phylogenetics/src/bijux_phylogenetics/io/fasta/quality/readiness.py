from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.phylo.alignment import (
    AlignmentMethodReadiness,
    AlignmentReadinessReport,
)

from ..coding import inspect_coding_alignment
from ..core import infer_alignment_alphabet, load_fasta_records
from ..records import (
    classify_alignment_sequences,
    detect_sequence_length_outliers,
    summarise_fasta,
)
from .quality_report import build_alignment_quality_report
from .site_diagnostics import assess_alignment_low_information
from .window_diagnostics import (
    detect_over_aligned_regions,
    detect_under_aligned_regions,
)


def summarize_alignment_readiness(path: Path) -> AlignmentReadinessReport:
    """Classify whether an input alignment is ready for key downstream analysis families."""
    sequence_kind = classify_alignment_sequences(path)
    records = load_fasta_records(path)
    inferred_alphabet = infer_alignment_alphabet(records)
    sequence_count = len(records)
    alignment_length = (
        sequence_kind.max_sequence_length
        if sequence_kind.state != "raw_sequence_fasta"
        else None
    )
    methods: list[AlignmentMethodReadiness] = []
    warnings: list[str] = []
    length_outliers = detect_sequence_length_outliers(path)

    if sequence_kind.state == "raw_sequence_fasta":
        common_blocker = ["input sequences are not yet aligned"]
        methods.extend(
            [
                AlignmentMethodReadiness(
                    analysis="distance",
                    ready=False,
                    blockers=common_blocker,
                    warnings=[],
                ),
                AlignmentMethodReadiness(
                    analysis="maximum_likelihood",
                    ready=False,
                    blockers=common_blocker,
                    warnings=[],
                ),
                AlignmentMethodReadiness(
                    analysis="bayesian",
                    ready=False,
                    blockers=common_blocker,
                    warnings=[],
                ),
                AlignmentMethodReadiness(
                    analysis="coding", ready=False, blockers=common_blocker, warnings=[]
                ),
                AlignmentMethodReadiness(
                    analysis="protein",
                    ready=inferred_alphabet == "protein",
                    blockers=[] if inferred_alphabet == "protein" else common_blocker,
                    warnings=[],
                ),
            ]
        )
        if length_outliers:
            warnings.append(
                "raw sequences include substantial length outliers before alignment"
            )
        return AlignmentReadinessReport(
            path=path,
            sequence_kind=sequence_kind,
            inferred_alphabet=inferred_alphabet,
            sequence_count=sequence_count,
            alignment_length=alignment_length,
            methods=methods,
            warnings=warnings,
        )

    summary = summarise_fasta(path)
    quality = build_alignment_quality_report(path)
    low_information = assess_alignment_low_information(path)
    over_aligned = detect_over_aligned_regions(path)
    under_aligned = detect_under_aligned_regions(path)
    coding = (
        inspect_coding_alignment(path) if inferred_alphabet in {"dna", "rna"} else None
    )

    def _method(
        analysis: str, ready: bool, blockers: list[str], extra_warnings: list[str]
    ) -> AlignmentMethodReadiness:
        return AlignmentMethodReadiness(
            analysis=analysis, ready=ready, blockers=blockers, warnings=extra_warnings
        )

    generic_warnings: list[str] = []
    if sequence_kind.state == "ambiguous_equal_length_fasta":
        generic_warnings.append(
            "equal-length ungapped FASTA may be aligned, but prior alignment cannot be proved from sequence shape alone"
        )
    if over_aligned:
        generic_warnings.append("one or more windows look suspiciously over-aligned")
    if under_aligned:
        generic_warnings.append("one or more windows look suspiciously under-aligned")
    if quality.composition_outliers:
        generic_warnings.append("composition outliers may bias downstream inference")
    if quality.sequence_length_outliers:
        generic_warnings.append(
            "sequence length outliers suggest problematic trimming, concatenation, or mixed loci"
        )
    if coding is not None and coding.mixed_coding_signals:
        generic_warnings.append(
            "alignment mixes coding-like and noncoding-like sequence behavior"
        )
    if low_information.low_information:
        generic_warnings.extend(low_information.reasons)
    warnings.extend(generic_warnings)

    base_alignment_blockers = (
        ["alignment contains invalid characters for the inferred alphabet"]
        if quality.invalid_characters
        else []
    )
    distance_blockers = list(base_alignment_blockers)
    if summary.variable_site_count == 0:
        distance_blockers.append("alignment has no variable sites")
    if low_information.low_information:
        distance_blockers.append(
            "alignment has too few parsimony-informative sites for defensible inference"
        )
    methods.append(
        _method(
            "distance",
            ready=not distance_blockers,
            blockers=distance_blockers,
            extra_warnings=generic_warnings,
        )
    )

    model_blockers = list(base_alignment_blockers)
    if summary.variable_site_count == 0:
        model_blockers.append("alignment has no variable sites")
    if low_information.low_information:
        model_blockers.append(
            "alignment has too few parsimony-informative sites for defensible inference"
        )
    methods.append(
        _method(
            "maximum_likelihood",
            ready=not model_blockers,
            blockers=model_blockers,
            extra_warnings=generic_warnings,
        )
    )
    methods.append(
        _method(
            "bayesian",
            ready=not model_blockers,
            blockers=model_blockers,
            extra_warnings=generic_warnings,
        )
    )

    coding_blockers: list[str] = []
    coding_warnings = list(generic_warnings)
    if inferred_alphabet not in {"dna", "rna"}:
        coding_blockers.append("coding analysis requires a nucleotide alignment")
    elif coding is not None:
        if not coding.alignment_length_multiple_of_three:
            coding_blockers.append("alignment length is not divisible by three")
        if coding.frameshift_like_sequences:
            coding_blockers.append(
                "one or more sequences contain partial codons after gaps and missing data are removed"
            )
        if any(not stop.terminal for stop in coding.stop_codons):
            coding_blockers.append(
                "one or more sequences contain premature stop codons"
            )
        if any(stop.terminal for stop in coding.stop_codons):
            coding_warnings.append(
                "terminal stop codons were detected and should be verified against coding conventions"
            )
    methods.append(
        _method(
            "coding",
            ready=not coding_blockers,
            blockers=coding_blockers,
            extra_warnings=coding_warnings,
        )
    )

    protein_blockers = (
        []
        if inferred_alphabet == "protein"
        else ["protein analysis requires an amino-acid alignment"]
    )
    methods.append(
        _method(
            "protein",
            ready=not protein_blockers,
            blockers=protein_blockers,
            extra_warnings=generic_warnings,
        )
    )

    return AlignmentReadinessReport(
        path=path,
        sequence_kind=sequence_kind,
        inferred_alphabet=inferred_alphabet,
        sequence_count=summary.sequence_count,
        alignment_length=summary.alignment_length,
        methods=methods,
        warnings=warnings,
    )
