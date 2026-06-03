from __future__ import annotations

from ._shared import (
    _TABULAR_WARNING_PREFIX_PATTERN,
    DEFAULT_BURNIN_FRACTIONS,
    Path,
    TraceConvergenceReport,
    _beast_artifact_error,
    _beast_state_field,
    _classify_beast_parameter,
    _mean_beast_parameter,
    _normalize_tabular_field,
    _summary_parameters_by_category,
    compute_clade_frequency_table,
    compute_consensus_tree,
    csv,
    normalize_burnin_fractions,
    summarize_burnin_clade_shifts,
    summarize_burnin_parameter_shifts,
    summarize_maximum_clade_credibility_tree,
    summarize_trace_convergence,
    write_taxon_rows,
)
from .models import (
    BeastBurninSensitivityReport,
    BeastBurninSensitivitySlice,
    BeastChainMixingIssue,
    BeastChainMixingReport,
    BeastConvergenceReport,
    BeastLogParameterSummary,
    BeastLogReport,
    BeastLogRow,
    BeastLogSummaryReport,
    BeastLogValidationIssue,
    BeastPosteriorClade,
    BeastPosteriorDecompositionReport,
    BeastPosteriorDecompositionRow,
    BeastPosteriorLogValidationReport,
)

_POSTERIOR_DECOMPOSITION_TOLERANCE = 1e-9


def parse_beast_log(path: Path) -> BeastLogReport:
    """Parse a BEAST-style log table into deterministic numeric rows."""
    if not path.exists():
        raise _beast_artifact_error(
            f"BEAST log file was not found: {path}",
            code="beast_log_missing_file",
            path=path,
            artifact_kind="beast-log",
            details={"expected_section": "posterior log file"},
        )
    logical_lines = _read_tabular_artifact_lines(path)
    header_fields, data_lines = _split_beast_log_table(logical_lines, path=path)
    state_field = _beast_state_field(header_fields)
    if state_field is None:
        raise _beast_artifact_error(
            f"BEAST log lacks a state column: {path}",
            code="beast_log_missing_state_column",
            path=path,
            artifact_kind="beast-log",
            details={"columns": list(header_fields)},
        )
    columns = [field for field in header_fields if field and field != state_field]
    rows: list[BeastLogRow] = []
    for row_number, fields in enumerate(data_lines, start=2):
        row = _build_tabular_row(header_fields, fields)
        raw_state = _normalize_tabular_field(row.get(state_field))
        if raw_state in {None, ""}:
            raise _beast_artifact_error(
                f"BEAST log contains an empty state value on row {row_number}: {path}",
                code="beast_log_missing_state_value",
                path=path,
                artifact_kind="beast-log",
                details={"row_number": row_number, "expected_section": "state column"},
            )
        try:
            state = int(float(raw_state))
        except ValueError as error:
            raise _beast_artifact_error(
                f"BEAST log contains a non-numeric state value on row {row_number}: {path}",
                code="beast_log_invalid_state_value",
                path=path,
                artifact_kind="beast-log",
                details={
                    "row_number": row_number,
                    "value": raw_state,
                    "expected_section": "state column",
                },
            ) from error
        values: dict[str, float] = {}
        for column in columns:
            raw_value = _normalize_tabular_field(row.get(column))
            if raw_value in {None, ""}:
                raise _beast_artifact_error(
                    f"BEAST log is missing a sampled value for '{column}' on row {row_number}: {path}",
                    code="beast_log_missing_parameter_value",
                    path=path,
                    artifact_kind="beast-log",
                    details={
                        "row_number": row_number,
                        "column": column,
                        "expected_section": "sampled parameter row",
                    },
                )
            try:
                values[column] = float(raw_value)
            except ValueError as error:
                raise _beast_artifact_error(
                    f"BEAST log contains a non-numeric value for '{column}' on row {row_number}: {path}",
                    code="beast_log_invalid_parameter_value",
                    path=path,
                    artifact_kind="beast-log",
                    details={
                        "row_number": row_number,
                        "column": column,
                        "value": raw_value,
                        "expected_section": "sampled parameter row",
                    },
                ) from error
        rows.append(BeastLogRow(state=state, values=values))
    if not rows:
        raise _beast_artifact_error(
            f"BEAST log contains no sampled rows: {path}",
            code="beast_log_missing_rows",
            path=path,
            artifact_kind="beast-log",
            details={"expected_section": "sampled parameter rows"},
        )
    return BeastLogReport(path=path, row_count=len(rows), columns=columns, rows=rows)


def _read_tabular_artifact_lines(path: Path) -> list[str]:
    lines: list[str] = []
    with path.open(encoding="utf-8", newline="") as handle:
        for raw_line in handle:
            if not raw_line.strip() or raw_line.lstrip().startswith("#"):
                continue
            lines.append(raw_line.rstrip("\r\n"))
    return lines


def _split_tabular_fields(line: str) -> list[str]:
    return [_normalize_tabular_field(field) or "" for field in line.split("\t")]


def _is_tabular_warning_line(line: str) -> bool:
    stripped = line.strip()
    if not stripped or "\t" in stripped:
        return False
    return _TABULAR_WARNING_PREFIX_PATTERN.match(stripped) is not None


def _trim_trailing_empty_fields(
    fields: list[str],
    *,
    expected_count: int,
) -> list[str]:
    trimmed = list(fields)
    while len(trimmed) > expected_count and trimmed[-1] == "":
        trimmed.pop()
    return trimmed


def _build_tabular_row(
    header_fields: list[str],
    fields: list[str],
) -> dict[str, str | None]:
    return {
        header_fields[index]: (fields[index] if index < len(fields) else None)
        for index in range(len(header_fields))
    }


def _split_beast_log_table(
    logical_lines: list[str],
    *,
    path: Path,
) -> tuple[list[str], list[list[str]]]:
    if not logical_lines:
        raise _beast_artifact_error(
            f"BEAST log contains no header row: {path}",
            code="beast_log_missing_header",
            path=path,
            artifact_kind="beast-log",
            details={"expected_section": "tabular header"},
        )
    header_fields: list[str] | None = None
    header_index = -1
    for index, line in enumerate(logical_lines):
        if _is_tabular_warning_line(line):
            continue
        fields = _split_tabular_fields(line)
        if _beast_state_field(fields) is not None:
            header_fields = fields
            header_index = index
            break
        if len(fields) > 1:
            raise _beast_artifact_error(
                f"BEAST log lacks a state column: {path}",
                code="beast_log_missing_state_column",
                path=path,
                artifact_kind="beast-log",
                details={
                    "columns": fields,
                    "expected_section": "state column",
                },
            )
    if header_fields is None:
        raise _beast_artifact_error(
            f"BEAST log contains no header row: {path}",
            code="beast_log_missing_header",
            path=path,
            artifact_kind="beast-log",
            details={"expected_section": "tabular header"},
        )
    data_lines: list[list[str]] = []
    for line in logical_lines[header_index + 1 :]:
        if _is_tabular_warning_line(line):
            continue
        fields = _trim_trailing_empty_fields(
            _split_tabular_fields(line),
            expected_count=len(header_fields),
        )
        if len(fields) > len(header_fields):
            raise _beast_artifact_error(
                f"BEAST log contains more fields than its header: {path}",
                code="beast_log_unexpected_field_count",
                path=path,
                artifact_kind="beast-log",
                details={
                    "expected_field_count": len(header_fields),
                    "observed_field_count": len(fields),
                    "expected_section": "sampled parameter row",
                },
            )
        data_lines.append(fields)
    return header_fields, data_lines


def summarize_beast_log(
    path: Path,
    *,
    burnin_fraction: float = 0.0,
) -> BeastLogSummaryReport:
    """Summarize a BEAST posterior log after discarding an optional burn-in fraction."""
    report = parse_beast_log(path)
    burnin_row_count, kept_rows = _split_beast_log_rows(
        report, burnin_fraction=burnin_fraction
    )
    convergence = summarize_trace_convergence(
        path=path,
        rows=[row.values for row in kept_rows],
        columns=report.columns,
        ess_threshold=0.0,
        mean_shift_threshold=float("inf"),
    )
    parameter_summaries = [
        BeastLogParameterSummary(
            parameter=summary.parameter,
            parameter_category=_classify_beast_parameter(summary.parameter),
            sample_count=summary.sample_count,
            effective_sample_size=summary.effective_sample_size,
            mean=summary.mean,
            median=summary.median,
            standard_deviation=summary.standard_deviation,
            minimum=summary.minimum,
            maximum=summary.maximum,
            hpd_95_lower=summary.hpd_95_lower,
            hpd_95_upper=summary.hpd_95_upper,
            first_half_mean=summary.first_half_mean,
            second_half_mean=summary.second_half_mean,
            standardized_mean_shift=summary.standardized_mean_shift,
        )
        for summary in convergence.series
    ]
    return BeastLogSummaryReport(
        path=path,
        burnin_fraction=burnin_fraction,
        burnin_row_count=burnin_row_count,
        kept_row_count=len(kept_rows),
        first_kept_state=kept_rows[0].state,
        last_kept_state=kept_rows[-1].state,
        posterior_parameters=_summary_parameters_by_category(
            parameter_summaries, category="posterior"
        ),
        likelihood_parameters=_summary_parameters_by_category(
            parameter_summaries, category="likelihood"
        ),
        prior_parameters=_summary_parameters_by_category(
            parameter_summaries, category="prior"
        ),
        clock_parameters=_summary_parameters_by_category(
            parameter_summaries, category="clock"
        ),
        tree_parameters=_summary_parameters_by_category(
            parameter_summaries, category="tree"
        ),
        other_parameters=_summary_parameters_by_category(
            parameter_summaries, category="other"
        ),
        parameter_summaries=parameter_summaries,
    )


def write_beast_log_summary_table(path: Path, report: BeastLogSummaryReport) -> Path:
    """Write a reviewer-facing TSV summary of one BEAST posterior log."""
    rows = [
        {
            "parameter_category": summary.parameter_category,
            "parameter": summary.parameter,
            "sample_count": str(summary.sample_count),
            "effective_sample_size": format(summary.effective_sample_size, ".15g"),
            "mean": format(summary.mean, ".15g"),
            "median": format(summary.median, ".15g"),
            "standard_deviation": format(summary.standard_deviation, ".15g"),
            "minimum": format(summary.minimum, ".15g"),
            "maximum": format(summary.maximum, ".15g"),
            "hpd_95_lower": format(summary.hpd_95_lower, ".15g"),
            "hpd_95_upper": format(summary.hpd_95_upper, ".15g"),
            "first_half_mean": format(summary.first_half_mean, ".15g"),
            "second_half_mean": format(summary.second_half_mean, ".15g"),
            "standardized_mean_shift": format(summary.standardized_mean_shift, ".15g"),
            "burnin_fraction": format(report.burnin_fraction, ".15g"),
            "burnin_row_count": str(report.burnin_row_count),
            "kept_row_count": str(report.kept_row_count),
            "first_kept_state": str(report.first_kept_state),
            "last_kept_state": str(report.last_kept_state),
        }
        for summary in report.parameter_summaries
    ]
    return write_taxon_rows(
        path,
        columns=[
            "parameter_category",
            "parameter",
            "sample_count",
            "effective_sample_size",
            "mean",
            "median",
            "standard_deviation",
            "minimum",
            "maximum",
            "hpd_95_lower",
            "hpd_95_upper",
            "first_half_mean",
            "second_half_mean",
            "standardized_mean_shift",
            "burnin_fraction",
            "burnin_row_count",
            "kept_row_count",
            "first_kept_state",
            "last_kept_state",
        ],
        rows=rows,
    )


def summarize_beast_posterior_decomposition(
    path: Path,
    *,
    burnin_fraction: float = 0.0,
    identity_tolerance: float = _POSTERIOR_DECOMPOSITION_TOLERANCE,
) -> BeastPosteriorDecompositionReport:
    """Separate posterior, likelihood, and prior terms from a BEAST log."""
    report = parse_beast_log(path)
    burnin_row_count, kept_rows = _split_beast_log_rows(
        report, burnin_fraction=burnin_fraction
    )
    if "posterior" not in report.columns or "likelihood" not in report.columns:
        raise _beast_artifact_error(
            f"BEAST posterior decomposition requires posterior and likelihood columns: {path}",
            code="beast_log_missing_posterior_terms",
            path=path,
            artifact_kind="beast-log",
            details={
                "required_columns": ["posterior", "likelihood"],
                "columns": list(report.columns),
            },
        )
    prior_logged = "prior" in report.columns
    rows: list[BeastPosteriorDecompositionRow] = []
    maximum_absolute_delta = 0.0
    verified = True
    for row in kept_rows:
        log_posterior = row.values["posterior"]
        log_likelihood = row.values["likelihood"]
        log_prior = (
            row.values["prior"] if prior_logged else log_posterior - log_likelihood
        )
        decomposition_delta = log_posterior - (log_likelihood + log_prior)
        decomposition_valid = abs(decomposition_delta) <= identity_tolerance
        maximum_absolute_delta = max(maximum_absolute_delta, abs(decomposition_delta))
        verified = verified and decomposition_valid
        rows.append(
            BeastPosteriorDecompositionRow(
                state=row.state,
                log_posterior=log_posterior,
                log_likelihood=log_likelihood,
                log_prior=log_prior,
                decomposition_delta=decomposition_delta,
                decomposition_valid=decomposition_valid,
            )
        )
    return BeastPosteriorDecompositionReport(
        path=path,
        burnin_fraction=burnin_fraction,
        burnin_row_count=burnin_row_count,
        kept_row_count=len(kept_rows),
        first_kept_state=kept_rows[0].state,
        last_kept_state=kept_rows[-1].state,
        posterior_term_source="logged",
        likelihood_term_source="logged",
        prior_term_source=(
            "logged" if prior_logged else "derived_from_posterior_and_likelihood"
        ),
        identity_tolerance=identity_tolerance,
        verified=verified,
        maximum_absolute_delta=maximum_absolute_delta,
        rows=rows,
    )


def write_beast_posterior_decomposition_table(
    path: Path,
    report: BeastPosteriorDecompositionReport,
) -> Path:
    """Write one row per retained BEAST sample with posterior decomposition terms."""
    return write_taxon_rows(
        path,
        columns=[
            "state",
            "log_posterior",
            "log_likelihood",
            "log_prior",
            "decomposition_delta",
            "decomposition_valid",
            "posterior_term_source",
            "likelihood_term_source",
            "prior_term_source",
            "identity_tolerance",
        ],
        rows=[
            {
                "state": str(row.state),
                "log_posterior": format(row.log_posterior, ".15g"),
                "log_likelihood": format(row.log_likelihood, ".15g"),
                "log_prior": format(row.log_prior, ".15g"),
                "decomposition_delta": format(row.decomposition_delta, ".15g"),
                "decomposition_valid": str(row.decomposition_valid).lower(),
                "posterior_term_source": report.posterior_term_source,
                "likelihood_term_source": report.likelihood_term_source,
                "prior_term_source": report.prior_term_source,
                "identity_tolerance": format(report.identity_tolerance, ".15g"),
            }
            for row in report.rows
        ],
    )


def validate_beast_posterior_log(
    path: Path,
    *,
    required_columns: tuple[str, ...] = ("posterior", "likelihood"),
) -> BeastPosteriorLogValidationReport:
    """Validate that a BEAST posterior log contains the required fields and monotonic sampled states."""
    issues: list[BeastLogValidationIssue] = []
    with path.open(encoding="utf-8", newline="") as handle:
        filtered_lines = [
            line
            for line in handle
            if line.strip() and not line.lstrip().startswith("#")
        ]
    reader = csv.DictReader(filtered_lines, delimiter="\t")
    if reader.fieldnames is None:
        issues.append(
            BeastLogValidationIssue(
                code="missing-header", message="BEAST log contains no header row"
            )
        )
        return BeastPosteriorLogValidationReport(
            path=path,
            row_count=0,
            state_count=0,
            required_columns=list(required_columns),
            observed_columns=[],
            missing_columns=list(required_columns),
            issues=issues,
            valid=False,
        )
    state_field = _beast_state_field(reader.fieldnames)
    observed_columns = [field for field in reader.fieldnames if field]
    if state_field is None:
        issues.append(
            BeastLogValidationIssue(
                code="missing-state-column", message="BEAST log lacks a state column"
            )
        )
    missing_columns = [
        column for column in required_columns if column not in observed_columns
    ]
    for column in missing_columns:
        issues.append(
            BeastLogValidationIssue(
                code="missing-required-column",
                message=f"missing required BEAST log column '{column}'",
                column=column,
            )
        )
    previous_state: int | None = None
    row_count = 0
    state_count = 0
    for row_index, row in enumerate(reader, start=2):
        row_count += 1
        if state_field is not None:
            raw_state = row.get(state_field, "")
            if not raw_state:
                issues.append(
                    BeastLogValidationIssue(
                        code="missing-state-value",
                        message="row is missing a sampled state",
                        row=row_index,
                        column=state_field,
                    )
                )
            else:
                try:
                    state_value = int(float(raw_state))
                    state_count += 1
                    if previous_state is not None and state_value <= previous_state:
                        issues.append(
                            BeastLogValidationIssue(
                                code="nonmonotonic-state",
                                message="sampled states must increase strictly through the log",
                                row=row_index,
                                column=state_field,
                            )
                        )
                    previous_state = state_value
                except ValueError:
                    issues.append(
                        BeastLogValidationIssue(
                            code="invalid-state-value",
                            message="sampled state must be numeric",
                            row=row_index,
                            column=state_field,
                        )
                    )
        for column in observed_columns:
            if column == state_field:
                continue
            raw_value = row.get(column, "")
            if raw_value in {None, ""}:
                issues.append(
                    BeastLogValidationIssue(
                        code="missing-parameter-value",
                        message=f"missing sampled value for '{column}'",
                        row=row_index,
                        column=column,
                    )
                )
                continue
            try:
                float(raw_value)
            except ValueError:
                issues.append(
                    BeastLogValidationIssue(
                        code="invalid-parameter-value",
                        message=f"sampled value for '{column}' must be numeric",
                        row=row_index,
                        column=column,
                    )
                )
    if row_count == 0:
        issues.append(
            BeastLogValidationIssue(
                code="missing-rows", message="BEAST log contains no sampled rows"
            )
        )
    return BeastPosteriorLogValidationReport(
        path=path,
        row_count=row_count,
        state_count=state_count,
        required_columns=list(required_columns),
        observed_columns=observed_columns,
        missing_columns=missing_columns,
        issues=issues,
        valid=not issues,
    )


def write_beast_burnin_sensitivity_slice_table(
    path: Path,
    report: BeastBurninSensitivityReport,
) -> Path:
    """Write one row per tested BEAST burn-in fraction."""
    return write_taxon_rows(
        path,
        columns=[
            "burnin_fraction",
            "burnin_tree_count",
            "kept_tree_count",
            "rooted_topology_count",
            "selected_tree_index",
            "clade_credibility_score",
            "clade_frequency_count",
            "kept_row_count",
            "first_kept_state",
            "last_kept_state",
            "posterior_mean",
            "likelihood_mean",
            "tree_height_mean",
            "consensus_newick",
        ],
        rows=[
            {
                "burnin_fraction": format(row.burnin_fraction, ".15g"),
                "burnin_tree_count": str(row.burnin_tree_count),
                "kept_tree_count": str(row.kept_tree_count),
                "rooted_topology_count": str(row.rooted_topology_count),
                "selected_tree_index": str(row.selected_tree_index),
                "clade_credibility_score": format(row.clade_credibility_score, ".15g"),
                "clade_frequency_count": str(row.clade_frequency_count),
                "kept_row_count": ""
                if row.kept_row_count is None
                else str(row.kept_row_count),
                "first_kept_state": ""
                if row.first_kept_state is None
                else str(row.first_kept_state),
                "last_kept_state": ""
                if row.last_kept_state is None
                else str(row.last_kept_state),
                "posterior_mean": ""
                if row.posterior_mean is None
                else format(row.posterior_mean, ".15g"),
                "likelihood_mean": ""
                if row.likelihood_mean is None
                else format(row.likelihood_mean, ".15g"),
                "tree_height_mean": ""
                if row.tree_height_mean is None
                else format(row.tree_height_mean, ".15g"),
                "consensus_newick": row.consensus_newick,
            }
            for row in report.slices
        ],
    )


def assess_beast_burnin_sensitivity(
    posterior_tree_path: Path,
    *,
    log_path: Path | None = None,
    burnin_fractions: tuple[float, ...] = DEFAULT_BURNIN_FRACTIONS,
) -> BeastBurninSensitivityReport:
    """Compare posterior summaries across multiple BEAST burn-in fractions."""
    ordered_fractions = normalize_burnin_fractions(burnin_fractions)
    slices: list[BeastBurninSensitivitySlice] = []
    previous_newick: str | None = None
    previous_consensus: str | None = None
    changed_mcc_count = 0
    changed_consensus_count = 0
    parameter_summaries_by_fraction: dict[float, list[BeastLogParameterSummary]] = {}
    clade_frequencies_by_fraction: dict[float, list[BeastPosteriorClade]] = {}
    for fraction in ordered_fractions:
        _, mcc_report = summarize_maximum_clade_credibility_tree(
            posterior_tree_path,
            burnin_fraction=fraction,
        )
        _consensus_tree, consensus_report = compute_consensus_tree(
            mcc_report.filtered_tree_set_path
        )
        clade_report = compute_clade_frequency_table(mcc_report.filtered_tree_set_path)
        posterior_mean = None
        likelihood_mean = None
        tree_height_mean = None
        kept_row_count = None
        first_kept_state = None
        last_kept_state = None
        if log_path is not None:
            log_summary = summarize_beast_log(log_path, burnin_fraction=fraction)
            parameter_summaries_by_fraction[fraction] = log_summary.parameter_summaries
            kept_row_count = log_summary.kept_row_count
            first_kept_state = log_summary.first_kept_state
            last_kept_state = log_summary.last_kept_state
            posterior_mean = _mean_beast_parameter(log_summary, "posterior")
            likelihood_mean = _mean_beast_parameter(log_summary, "likelihood")
            tree_height_mean = _mean_beast_parameter(log_summary, "treeHeight")
            if tree_height_mean is None:
                tree_height_mean = _mean_beast_parameter(log_summary, "tree.height")
        clade_frequencies_by_fraction[fraction] = [
            BeastPosteriorClade(
                clade=row.clade,
                tree_count=row.tree_count,
                frequency=row.frequency,
            )
            for row in clade_report.clade_frequencies
        ]
        slices.append(
            BeastBurninSensitivitySlice(
                burnin_fraction=fraction,
                burnin_tree_count=mcc_report.burnin_tree_count,
                kept_tree_count=mcc_report.kept_tree_count,
                rooted_topology_count=mcc_report.rooted_topology_count,
                selected_tree_index=mcc_report.selected_tree_index,
                clade_credibility_score=mcc_report.clade_credibility_score,
                consensus_newick=consensus_report.consensus_newick,
                clade_frequency_count=len(clade_report.clade_frequencies),
                kept_row_count=kept_row_count,
                first_kept_state=first_kept_state,
                last_kept_state=last_kept_state,
                posterior_mean=posterior_mean,
                likelihood_mean=likelihood_mean,
                tree_height_mean=tree_height_mean,
            )
        )
        if previous_newick is not None and previous_newick != mcc_report.mcc_newick:
            changed_mcc_count += 1
        previous_newick = mcc_report.mcc_newick
        if (
            previous_consensus is not None
            and previous_consensus != consensus_report.consensus_newick
        ):
            changed_consensus_count += 1
        previous_consensus = consensus_report.consensus_newick
    parameter_shifts = summarize_burnin_parameter_shifts(
        parameter_summaries_by_fraction
    )
    clade_shifts = summarize_burnin_clade_shifts(clade_frequencies_by_fraction)
    warnings: list[str] = []
    if changed_mcc_count:
        warnings.append(
            "maximum clade credibility topology changes across tested burn-in fractions"
        )
    if changed_consensus_count:
        warnings.append(
            "majority-rule consensus topology changes across tested burn-in fractions"
        )
    if len({row.rooted_topology_count for row in slices}) > 1:
        warnings.append(
            "rooted topology diversity changes across tested burn-in fractions"
        )
    if any(shift.unstable for shift in parameter_shifts):
        warnings.append(
            "one or more posterior parameter 95% HPD intervals do not overlap across tested burn-in fractions"
        )
    if any(shift.unstable for shift in clade_shifts):
        warnings.append(
            "one or more posterior clade probabilities cross the majority-rule threshold across tested burn-in fractions"
        )
    return BeastBurninSensitivityReport(
        posterior_tree_path=posterior_tree_path,
        log_path=log_path,
        slices=slices,
        changed_mcc_count=changed_mcc_count,
        changed_consensus_count=changed_consensus_count,
        parameter_shifts=parameter_shifts,
        clade_shifts=clade_shifts,
        unstable_parameter_count=sum(1 for shift in parameter_shifts if shift.unstable),
        unstable_clade_count=sum(1 for shift in clade_shifts if shift.unstable),
        warnings=warnings,
    )


def assess_beast_chain_mixing(
    log_paths: list[Path],
    *,
    ess_threshold: float = 200.0,
    mean_shift_threshold: float = 0.5,
    cross_chain_mean_shift_threshold: float = 0.75,
    stuck_parameter_span_threshold: float = 1e-9,
) -> BeastChainMixingReport:
    """Flag low ESS, mean drift, stuck parameters, and inconsistent chain means across BEAST logs."""
    if not log_paths:
        raise ValueError("assess_beast_chain_mixing requires at least one log path")
    chain_summaries = [
        assess_beast_convergence(
            path,
            ess_threshold=ess_threshold,
            mean_shift_threshold=mean_shift_threshold,
        )
        for path in log_paths
    ]
    issues: list[BeastChainMixingIssue] = []
    for summary in chain_summaries:
        for warning in summary.warnings:
            issues.append(
                BeastChainMixingIssue(
                    path=summary.path,
                    parameter=str(warning["parameter"]),
                    code=str(warning["code"]),
                    message=str(warning["message"]),
                    observed_value=float(warning["observed_value"]),
                    threshold=float(warning["threshold"]),
                )
            )
        for parameter_summary in summary.parameter_summaries:
            span = float(parameter_summary["maximum"]) - float(
                parameter_summary["minimum"]
            )
            if span <= stuck_parameter_span_threshold:
                issues.append(
                    BeastChainMixingIssue(
                        path=summary.path,
                        parameter=str(parameter_summary["parameter"]),
                        code="stuck-parameter",
                        message="parameter shows effectively no movement across the sampled chain",
                        observed_value=span,
                        threshold=stuck_parameter_span_threshold,
                    )
                )
    parameter_to_means: dict[str, list[tuple[Path, float]]] = {}
    for summary in chain_summaries:
        for parameter_summary in summary.parameter_summaries:
            parameter_to_means.setdefault(
                str(parameter_summary["parameter"]), []
            ).append((summary.path, float(parameter_summary["mean"])))
    for parameter, chain_means in parameter_to_means.items():
        if len(chain_means) < 2:
            continue
        mean_values = [value for _, value in chain_means]
        span = max(mean_values) - min(mean_values)
        if span > cross_chain_mean_shift_threshold:
            issues.append(
                BeastChainMixingIssue(
                    path=None,
                    parameter=parameter,
                    code="inconsistent-chains",
                    message="independent chains disagree more than the allowed mean-shift threshold",
                    observed_value=span,
                    threshold=cross_chain_mean_shift_threshold,
                )
            )
    return BeastChainMixingReport(
        log_paths=log_paths,
        chain_count=len(log_paths),
        converged=not issues,
        issues=issues,
        chain_summaries=chain_summaries,
    )


def assess_beast_convergence(
    path: Path,
    *,
    burnin_fraction: float = 0.0,
    ess_threshold: float = 200.0,
    mean_shift_threshold: float = 0.5,
) -> BeastConvergenceReport:
    """Flag low-ESS or unstable BEAST trace parameters."""
    report = parse_beast_log(path)
    burnin_row_count, kept_rows = _split_beast_log_rows(
        report, burnin_fraction=burnin_fraction
    )
    convergence = summarize_trace_convergence(
        path=path,
        rows=[row.values for row in kept_rows],
        columns=report.columns,
        ess_threshold=ess_threshold,
        mean_shift_threshold=mean_shift_threshold,
    )
    return _build_beast_convergence_report(
        convergence,
        burnin_fraction=burnin_fraction,
        burnin_row_count=burnin_row_count,
    )


def _build_beast_convergence_report(
    convergence: TraceConvergenceReport,
    *,
    burnin_fraction: float,
    burnin_row_count: int,
) -> BeastConvergenceReport:
    return BeastConvergenceReport(
        path=convergence.path,
        burnin_fraction=burnin_fraction,
        burnin_row_count=burnin_row_count,
        sample_count=convergence.sample_count,
        converged=convergence.converged,
        ess_threshold=convergence.ess_threshold,
        mean_shift_threshold=convergence.mean_shift_threshold,
        warnings=[
            {
                "parameter": warning.parameter,
                "code": warning.code,
                "message": warning.message,
                "observed_value": warning.observed_value,
                "threshold": warning.threshold,
            }
            for warning in convergence.warnings
        ],
        parameter_summaries=[
            {
                "parameter": summary.parameter,
                "sample_count": summary.sample_count,
                "effective_sample_size": summary.effective_sample_size,
                "mean": summary.mean,
                "median": summary.median,
                "standard_deviation": summary.standard_deviation,
                "minimum": summary.minimum,
                "maximum": summary.maximum,
                "hpd_95_lower": summary.hpd_95_lower,
                "hpd_95_upper": summary.hpd_95_upper,
                "first_half_mean": summary.first_half_mean,
                "second_half_mean": summary.second_half_mean,
                "standardized_mean_shift": summary.standardized_mean_shift,
            }
            for summary in convergence.series
        ],
    )


def _split_beast_log_rows(
    report: BeastLogReport,
    *,
    burnin_fraction: float,
) -> tuple[int, list[BeastLogRow]]:
    if not 0.0 <= burnin_fraction < 1.0:
        raise ValueError(
            f"burnin_fraction must be between 0 and 1, got {burnin_fraction}"
        )
    burnin_row_count = int(report.row_count * burnin_fraction)
    kept_rows = report.rows[burnin_row_count:]
    if not kept_rows:
        raise ValueError(
            f"BEAST log contains no sampled rows after burn-in filtering: {report.path}"
        )
    return burnin_row_count, kept_rows
