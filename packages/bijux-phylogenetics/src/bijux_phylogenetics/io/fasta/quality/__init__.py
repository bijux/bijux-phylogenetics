from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.phylo.alignment import (
    AlignmentAmbiguousColumnReport,
    AlignmentForensicReport,
    AlignmentLowInformationReport,
    AlignmentMethodReadiness,
    AlignmentMissingDataConcentration,
    AlignmentQualityReport,
    AlignmentReadinessReport,
    AlignmentRecord,
    AlignmentSummary,
    AlignmentSuspiciousRegion,
    AlignmentWindowSummary,
    AmbiguousAlignmentColumn,
    DuplicateSequencePolicyAction,
    DuplicateSequencePolicyReport,
    SequenceQualityRankingReport,
    SequenceQualityRankingRow,
)

from ..cleaning import (
    detect_identical_duplicate_sequences,
    detect_near_duplicate_sequences,
)
from ..coding import (
    inspect_coding_alignment,
)
from ..core import (
    _GAP_CHARACTERS,
    _detect_sequence_length_outlier_rows,
    _is_ambiguity_character,
    _is_explicit_missing,
    _validate_fraction_threshold,
    infer_alignment_alphabet,
    load_fasta_alignment,
    load_fasta_records,
)
from ..records import (
    classify_alignment_sequences,
    detect_sequence_length_outliers,
    summarise_fasta,
    summarise_records_as_alignment_summary,
)

_LOW_INFORMATION_SITE_THRESHOLD = 1
_LOW_INFORMATION_FRACTION_THRESHOLD = 0.0


def _summarize_alignment_windows_from_records(
    summary: AlignmentSummary,
    records: list[AlignmentRecord],
    *,
    window_size: int,
    step_size: int,
) -> list[AlignmentWindowSummary]:
    windows: list[AlignmentWindowSummary] = []
    for start_index in range(0, summary.alignment_length, step_size):
        end_index = min(start_index + window_size, summary.alignment_length)
        if end_index <= start_index:
            continue
        columns = list(
            zip(
                *(record.sequence[start_index:end_index] for record in records),
                strict=True,
            )
        )
        if not columns:
            continue

        gap_count = 0
        missing_count = 0
        ambiguity_count = 0
        variable_sites = 0
        disagreement_fractions: list[float] = []
        comparable_sites = 0
        total_residues = len(records) * len(columns)
        alphabet = summary.inferred_alphabet
        for column in columns:
            gap_count += sum(1 for residue in column if residue in _GAP_CHARACTERS)
            missing_count += sum(
                1 for residue in column if _is_explicit_missing(residue)
            )
            ambiguity_count += sum(
                1
                for residue in column
                if _is_ambiguity_character(residue, alphabet=alphabet)
            )
            comparable = [
                residue.upper()
                for residue in column
                if residue not in _GAP_CHARACTERS
                and not _is_explicit_missing(residue)
                and not _is_ambiguity_character(residue, alphabet=alphabet)
            ]
            if comparable:
                comparable_sites += len(comparable)
                state_counts = {
                    state: comparable.count(state) for state in set(comparable)
                }
                if len(state_counts) > 1:
                    variable_sites += 1
                disagreement_fractions.append(
                    1.0 - (max(state_counts.values()) / len(comparable))
                )
            else:
                disagreement_fractions.append(0.0)

        windows.append(
            AlignmentWindowSummary(
                start=start_index + 1,
                end=end_index,
                site_count=len(columns),
                gap_fraction=round(gap_count / total_residues, 15),
                missing_fraction=round(missing_count / total_residues, 15),
                ambiguity_fraction=round(ambiguity_count / total_residues, 15),
                variable_fraction=round(variable_sites / len(columns), 15),
                disagreement_fraction=round(
                    sum(disagreement_fractions) / len(disagreement_fractions),
                    15,
                ),
                comparable_fraction=round(comparable_sites / total_residues, 15),
            )
        )
        if end_index == summary.alignment_length:
            break
    return windows


def _detect_over_aligned_regions_from_windows(
    windows: list[AlignmentWindowSummary],
) -> list[AlignmentSuspiciousRegion]:
    regions: list[AlignmentSuspiciousRegion] = []
    for window in windows:
        uncertainty_fraction = (
            window.gap_fraction + window.missing_fraction + window.ambiguity_fraction
        )
        if uncertainty_fraction >= 0.4 and window.variable_fraction <= 0.2:
            regions.append(
                AlignmentSuspiciousRegion(
                    start=window.start,
                    end=window.end,
                    kind="over_aligned",
                    score=round(uncertainty_fraction - window.variable_fraction, 15),
                    note="gap- or ambiguity-heavy window with little residual variability; review for aggressive or artifactual alignment",
                )
            )
    return regions


def _detect_under_aligned_regions_from_windows(
    windows: list[AlignmentWindowSummary],
) -> list[AlignmentSuspiciousRegion]:
    regions: list[AlignmentSuspiciousRegion] = []
    for window in windows:
        if window.variable_fraction >= 0.7 and window.disagreement_fraction >= 0.35:
            regions.append(
                AlignmentSuspiciousRegion(
                    start=window.start,
                    end=window.end,
                    kind="under_aligned",
                    score=round(
                        window.variable_fraction + window.disagreement_fraction,
                        15,
                    ),
                    note="high local mismatch and disagreement suggest the region may require realignment or masking",
                )
            )
    return regions


def _assess_alignment_low_information_from_summary(
    summary: AlignmentSummary,
    *,
    minimum_informative_sites: int,
    minimum_informative_fraction: float,
) -> AlignmentLowInformationReport:
    informative_fraction = (
        0.0
        if summary.alignment_length == 0
        else round(
            summary.parsimony_informative_site_count / summary.alignment_length,
            15,
        )
    )
    reasons: list[str] = []
    if (
        summary.parsimony_informative_site_count < minimum_informative_sites
        and summary.variable_site_count == 0
    ):
        reasons.append(
            "alignment has fewer parsimony-informative sites than the minimum threshold for defensible inference"
        )
    if (
        summary.variable_site_count == 0
        and informative_fraction < minimum_informative_fraction
    ):
        reasons.append(
            "alignment has a very low parsimony-informative-site fraction and may not support stable inference"
        )
    return AlignmentLowInformationReport(
        sequence_count=summary.sequence_count,
        alignment_length=summary.alignment_length,
        parsimony_informative_site_count=summary.parsimony_informative_site_count,
        parsimony_informative_fraction=informative_fraction,
        threshold_site_count=minimum_informative_sites,
        threshold_fraction=minimum_informative_fraction,
        low_information=bool(reasons),
        reasons=reasons,
    )


def _build_ambiguous_alignment_column_report_from_summary(
    path: Path,
    summary: AlignmentSummary,
    *,
    threshold: float,
) -> AlignmentAmbiguousColumnReport:
    uncertainty_by_position = {
        row.position: row for row in summary.per_site_uncertainty
    }
    rows: list[AmbiguousAlignmentColumn] = []
    for position in range(1, summary.alignment_length + 1):
        uncertainty = uncertainty_by_position[position]
        ambiguity_burden = (
            uncertainty.ambiguity_fraction
            + uncertainty.missing_fraction
            + uncertainty.gap_fraction
        )
        if ambiguity_burden < threshold:
            continue
        rows.append(
            AmbiguousAlignmentColumn(
                position=position,
                ambiguity_fraction=uncertainty.ambiguity_fraction,
                missing_fraction=uncertainty.missing_fraction,
                gap_fraction=uncertainty.gap_fraction,
                comparable_fraction=round(max(0.0, 1.0 - ambiguity_burden), 15),
                note="column is dominated by ambiguity, missing data, or gaps and should be reviewed for masking",
            )
        )
    warnings = (
        [
            "alignment contains ambiguity-heavy columns that may be unsuitable for inference without masking"
        ]
        if rows
        else []
    )
    return AlignmentAmbiguousColumnReport(
        path=path,
        threshold=threshold,
        rows=rows,
        warnings=warnings,
    )


def _alignment_quality_components(summary: AlignmentSummary) -> dict[str, float]:
    informative_density = (
        0.0
        if summary.alignment_length == 0
        else summary.parsimony_informative_site_count / summary.alignment_length
    )
    return {
        "missingness": round(
            max(0.0, 1.0 - min(summary.missing_data_fraction, 1.0)),
            15,
        ),
        "gap_burden": round(max(0.0, 1.0 - min(summary.gap_fraction, 1.0)), 15),
        "composition_outliers": round(
            max(
                0.0,
                1.0
                - min(
                    len(summary.composition_outliers) / max(summary.sequence_count, 1),
                    1.0,
                ),
            ),
            15,
        ),
        "duplicates": round(
            max(
                0.0,
                1.0
                - min(
                    len(summary.duplicate_sequence_groups)
                    / max(summary.sequence_count, 1),
                    1.0,
                ),
            ),
            15,
        ),
        "informative_density": round(min(informative_density / 0.2, 1.0), 15),
    }


def _alignment_quality_score(components: dict[str, float]) -> float:
    return round(sum(components.values()) / max(len(components), 1) * 100.0, 3)


def _summarize_missing_data_concentration(
    summary: AlignmentSummary,
    *,
    threshold: float = 0.5,
) -> AlignmentMissingDataConcentration:
    concentrated_columns = [
        row.position
        for row in summary.per_site_missingness
        if row.missing_fraction >= threshold
    ]
    longest_run = 0
    longest_start: int | None = None
    longest_end: int | None = None
    current_run = 0
    current_start: int | None = None
    previous_position: int | None = None
    for position in concentrated_columns:
        if previous_position is None or position != previous_position + 1:
            current_run = 1
            current_start = position
        else:
            current_run += 1
        if current_run > longest_run:
            longest_run = current_run
            longest_start = current_start
            longest_end = position
        previous_position = position
    maximum_missing_fraction = max(
        (row.missing_fraction for row in summary.per_site_missingness),
        default=0.0,
    )
    maximum_missing_positions = [
        row.position
        for row in summary.per_site_missingness
        if row.missing_fraction == maximum_missing_fraction
    ]
    return AlignmentMissingDataConcentration(
        threshold=threshold,
        concentrated_column_count=len(concentrated_columns),
        concentrated_column_fraction=(
            0.0
            if summary.alignment_length == 0
            else len(concentrated_columns) / summary.alignment_length
        ),
        longest_concentrated_run=longest_run,
        longest_concentrated_run_start=longest_start,
        longest_concentrated_run_end=longest_end,
        maximum_missing_fraction=maximum_missing_fraction,
        maximum_missing_positions=maximum_missing_positions,
    )


def _alignment_suspicion_reasons(
    *,
    low_information: AlignmentLowInformationReport,
    missing_data_concentration: AlignmentMissingDataConcentration,
    ambiguous_column_count: int,
    over_aligned_count: int,
    under_aligned_count: int,
    invalid_character_count: int,
) -> list[str]:
    reasons: list[str] = []
    if low_information.low_information:
        reasons.append("alignment has low information content for defensible inference")
    if missing_data_concentration.longest_concentrated_run >= 2:
        reasons.append("alignment concentrates missing data into adjacent columns")
    elif missing_data_concentration.concentrated_column_count > 0:
        reasons.append("alignment contains one or more highly missing columns")
    if ambiguous_column_count > 0:
        reasons.append("alignment contains ambiguity-heavy columns")
    if over_aligned_count > 0:
        reasons.append("alignment contains suspiciously over-aligned windows")
    if under_aligned_count > 0:
        reasons.append("alignment contains suspiciously under-aligned windows")
    if invalid_character_count > 0:
        reasons.append(
            "alignment contains invalid characters for the inferred alphabet"
        )
    return reasons


def summarize_alignment_windows(
    path: Path,
    *,
    window_size: int = 30,
    step_size: int = 10,
) -> list[AlignmentWindowSummary]:
    """Summarize an aligned FASTA file in sliding windows."""
    if window_size <= 0:
        raise ValueError(f"window_size must be positive, got {window_size}")
    if step_size <= 0:
        raise ValueError(f"step_size must be positive, got {step_size}")

    summary = summarise_fasta(path)
    records = load_fasta_alignment(path)
    return _summarize_alignment_windows_from_records(
        summary,
        records,
        window_size=window_size,
        step_size=step_size,
    )


def detect_over_aligned_regions(
    path: Path,
    *,
    window_size: int = 30,
    step_size: int = 10,
) -> list[AlignmentSuspiciousRegion]:
    """Flag suspicious windows that look excessively gap-heavy or over-regularized."""
    return _detect_over_aligned_regions_from_windows(
        summarize_alignment_windows(path, window_size=window_size, step_size=step_size)
    )


def detect_under_aligned_regions(
    path: Path,
    *,
    window_size: int = 30,
    step_size: int = 10,
) -> list[AlignmentSuspiciousRegion]:
    """Flag suspicious windows with strong local mismatch or gap disorder."""
    return _detect_under_aligned_regions_from_windows(
        summarize_alignment_windows(path, window_size=window_size, step_size=step_size)
    )


def assess_alignment_low_information(
    path: Path,
    *,
    minimum_informative_sites: int = _LOW_INFORMATION_SITE_THRESHOLD,
    minimum_informative_fraction: float = _LOW_INFORMATION_FRACTION_THRESHOLD,
) -> AlignmentLowInformationReport:
    """Assess whether an alignment carries enough signal for defensible inference."""
    _validate_fraction_threshold(minimum_informative_fraction)
    if minimum_informative_sites < 0:
        raise ValueError(
            f"minimum_informative_sites must be non-negative, got {minimum_informative_sites}"
        )
    summary = summarise_fasta(path)
    return _assess_alignment_low_information_from_summary(
        summary,
        minimum_informative_sites=minimum_informative_sites,
        minimum_informative_fraction=minimum_informative_fraction,
    )


def build_duplicate_sequence_policy_report(
    path: Path,
    *,
    near_duplicate_threshold: float = 0.99,
) -> DuplicateSequencePolicyReport:
    """Build reviewer-facing policy recommendations for duplicate sequences."""
    exact_duplicates = detect_identical_duplicate_sequences(path)
    near_duplicates = detect_near_duplicate_sequences(
        path, identity_threshold=near_duplicate_threshold
    )
    actions: list[DuplicateSequencePolicyAction] = []
    warnings: list[str] = []
    if exact_duplicates:
        warnings.append(
            "exact duplicate sequences should be deduplicated or explicitly justified before inference"
        )
        for group in exact_duplicates:
            actions.append(
                DuplicateSequencePolicyAction(
                    action="collapse_exact_duplicates",
                    rationale=(
                        f"retain one representative such as {group.identifiers[0]} unless metadata show that "
                        "the duplicated labels represent distinct biological samples"
                    ),
                    affected_identifiers=group.identifiers,
                )
            )
    if near_duplicates:
        warnings.append(
            "near-duplicate sequences should be checked for replicate samples, contamination, or oversampling bias"
        )
        seen_pairs: set[tuple[str, ...]] = set()
        for pair in near_duplicates:
            identifiers = tuple(sorted((pair.left_identifier, pair.right_identifier)))
            if identifiers in seen_pairs:
                continue
            seen_pairs.add(identifiers)
            actions.append(
                DuplicateSequencePolicyAction(
                    action="review_near_duplicates",
                    rationale="inspect metadata and sampling provenance before keeping highly similar sequences together in inference",
                    affected_identifiers=list(identifiers),
                )
            )
    if not warnings:
        warnings.append("no duplicate-sequence policy actions are currently required")
    return DuplicateSequencePolicyReport(
        path=path,
        exact_duplicate_groups=exact_duplicates,
        near_duplicate_pairs=near_duplicates,
        policy_actions=actions,
        warnings=warnings,
    )


def build_ambiguous_alignment_column_report(
    path: Path,
    *,
    threshold: float = 0.5,
) -> AlignmentAmbiguousColumnReport:
    """Report columns dominated by ambiguity, explicit missingness, or gaps."""
    _validate_fraction_threshold(threshold)
    summary = summarise_fasta(path)
    return _build_ambiguous_alignment_column_report_from_summary(
        path,
        summary,
        threshold=threshold,
    )


def build_sequence_quality_ranking(path: Path) -> SequenceQualityRankingReport:
    """Rank aligned sequences by transparent quality burdens."""
    summary = summarise_fasta(path)
    composition_outlier_ids = {row.identifier for row in summary.composition_outliers}
    exact_duplicate_ids = {
        identifier
        for group in summary.duplicate_sequence_groups
        for identifier in group.identifiers
    }
    near_duplicate_ids = {
        identifier
        for pair in summary.near_duplicate_pairs
        for identifier in (pair.left_identifier, pair.right_identifier)
    }
    uncertainty_by_id = {
        row.identifier: row for row in summary.per_sequence_uncertainty
    }
    ranked: list[tuple[float, str, SequenceQualityRankingRow]] = []
    for identifier in summary.ids:
        uncertainty = uncertainty_by_id[identifier]
        composition_outlier = identifier in composition_outlier_ids
        if identifier in exact_duplicate_ids:
            duplicate_status = "exact_duplicate"
        elif identifier in near_duplicate_ids:
            duplicate_status = "near_duplicate"
        else:
            duplicate_status = "unique"
        penalty = (
            uncertainty.missing_fraction * 40.0
            + uncertainty.gap_fraction * 25.0
            + uncertainty.ambiguity_fraction * 20.0
            + (10.0 if composition_outlier else 0.0)
            + (
                10.0
                if duplicate_status == "exact_duplicate"
                else 5.0
                if duplicate_status == "near_duplicate"
                else 0.0
            )
        )
        score = round(max(0.0, 100.0 - penalty), 3)
        note_parts: list[str] = []
        if uncertainty.missing_fraction > 0.0:
            note_parts.append("missing data")
        if uncertainty.gap_fraction > 0.0:
            note_parts.append("gaps")
        if uncertainty.ambiguity_fraction > 0.0:
            note_parts.append("ambiguity codes")
        if composition_outlier:
            note_parts.append("composition outlier")
        if duplicate_status != "unique":
            note_parts.append(duplicate_status.replace("_", " "))
        note = (
            "quality burdens: " + ", ".join(note_parts)
            if note_parts
            else "no major quality burdens detected"
        )
        ranked.append(
            (
                score,
                identifier,
                SequenceQualityRankingRow(
                    identifier=identifier,
                    rank=0,
                    score=score,
                    missing_fraction=uncertainty.missing_fraction,
                    gap_fraction=uncertainty.gap_fraction,
                    ambiguity_fraction=uncertainty.ambiguity_fraction,
                    composition_outlier=composition_outlier,
                    duplicate_status=duplicate_status,
                    note=note,
                ),
            )
        )
    ranked.sort(key=lambda item: (item[0], item[1]))
    rows = [
        SequenceQualityRankingRow(
            identifier=row.identifier,
            rank=index,
            score=row.score,
            missing_fraction=row.missing_fraction,
            gap_fraction=row.gap_fraction,
            ambiguity_fraction=row.ambiguity_fraction,
            composition_outlier=row.composition_outlier,
            duplicate_status=row.duplicate_status,
            note=row.note,
        )
        for index, (_, _, row) in enumerate(ranked, start=1)
    ]
    warnings = (
        ["lower-ranked sequences should be reviewed before publication or inference"]
        if rows and any(row.score < 85.0 for row in rows)
        else []
    )
    if not summary.near_duplicate_scan_performed:
        warnings.append(
            "near-duplicate sequence ranking was skipped because the alignment exceeds the governed pairwise review threshold"
        )
    return SequenceQualityRankingReport(path=path, rows=rows, warnings=warnings)


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


def build_alignment_quality_report(path: Path) -> AlignmentQualityReport:
    """Generate a higher-level alignment quality report from composition and identity diagnostics."""
    records = load_fasta_alignment(path)
    summary = summarise_records_as_alignment_summary(path=path, records=records)
    inferred_alphabet = summary.inferred_alphabet
    invalid_characters = summary.invalid_characters
    composition_outliers = summary.composition_outliers
    sequence_length_outliers = _detect_sequence_length_outlier_rows(
        [(record.identifier, len(record.sequence)) for record in records]
    )
    duplicate_sequence_groups = summary.duplicate_sequence_groups
    near_duplicate_pairs = summary.near_duplicate_pairs
    low_information = _assess_alignment_low_information_from_summary(
        summary,
        minimum_informative_sites=_LOW_INFORMATION_SITE_THRESHOLD,
        minimum_informative_fraction=_LOW_INFORMATION_FRACTION_THRESHOLD,
    )
    ambiguous_columns = _build_ambiguous_alignment_column_report_from_summary(
        path,
        summary,
        threshold=0.5,
    )
    windows = _summarize_alignment_windows_from_records(
        summary,
        records,
        window_size=30,
        step_size=10,
    )
    over_aligned = _detect_over_aligned_regions_from_windows(windows)
    under_aligned = _detect_under_aligned_regions_from_windows(windows)
    missing_data_concentration = _summarize_missing_data_concentration(summary)
    suspicious_reasons = _alignment_suspicion_reasons(
        low_information=low_information,
        missing_data_concentration=missing_data_concentration,
        ambiguous_column_count=len(ambiguous_columns.rows),
        over_aligned_count=len(over_aligned),
        under_aligned_count=len(under_aligned),
        invalid_character_count=len(invalid_characters),
    )
    quality_components = _alignment_quality_components(summary)
    warnings: list[str] = []
    if invalid_characters:
        warnings.append(
            "alignment contains characters invalid for the inferred alphabet"
        )
    if composition_outliers:
        warnings.append("alignment contains composition outlier sequences")
    if sequence_length_outliers:
        warnings.append("alignment contains raw-sequence length outliers")
    if duplicate_sequence_groups:
        warnings.append("alignment contains identical duplicate sequences")
    if near_duplicate_pairs:
        warnings.append("alignment contains near-duplicate sequences")
    if not summary.near_duplicate_scan_performed:
        warnings.append(
            "near-duplicate sequence scan was skipped because the alignment exceeds the governed pairwise review threshold"
        )
    warnings.extend(reason for reason in suspicious_reasons if reason not in warnings)
    return AlignmentQualityReport(
        path=path,
        sequence_count=summary.sequence_count,
        alignment_length=summary.alignment_length,
        invariant_site_count=summary.constant_site_count,
        missing_data_fraction=summary.missing_data_fraction,
        gap_fraction=summary.gap_fraction,
        ambiguity_fraction=summary.ambiguity_fraction,
        variable_site_count=summary.variable_site_count,
        parsimony_informative_site_count=summary.parsimony_informative_site_count,
        per_sequence_uncertainty=summary.per_sequence_uncertainty,
        per_site_uncertainty=summary.per_site_uncertainty,
        inferred_alphabet=inferred_alphabet,
        invalid_characters=invalid_characters,
        composition_outliers=composition_outliers,
        sequence_length_outliers=sequence_length_outliers,
        duplicate_sequence_groups=duplicate_sequence_groups,
        near_duplicate_pairs=near_duplicate_pairs,
        missing_data_concentration=missing_data_concentration,
        suspicious_alignment=bool(suspicious_reasons),
        suspicious_reasons=suspicious_reasons,
        quality_score=_alignment_quality_score(quality_components),
        quality_components=quality_components,
        warnings=warnings,
        near_duplicate_scan_performed=summary.near_duplicate_scan_performed,
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
