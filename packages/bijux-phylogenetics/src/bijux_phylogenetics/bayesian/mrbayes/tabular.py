from __future__ import annotations

from pathlib import Path
import re

from .artifacts import _mrbayes_artifact_error
from .models import (
    MrBayesMcmcReport,
    MrBayesMcmcRow,
    MrBayesTraceReport,
    MrBayesTraceRow,
)

_TABULAR_WARNING_PREFIX_PATTERN = re.compile(
    r"^(warning|warn|caution|note|info)\b",
    flags=re.IGNORECASE,
)


def parse_mrbayes_parameter_traces(path: Path) -> MrBayesTraceReport:
    """Parse a MrBayes parameter trace table into deterministic numeric rows."""
    if not path.exists():
        raise _mrbayes_artifact_error(
            f"MrBayes trace file was not found: {path}",
            code="mrbayes_trace_missing_file",
            path=path,
            artifact_kind="mrbayes-trace",
            details={"expected_section": "parameter trace file"},
        )
    rows: list[MrBayesTraceRow] = []
    logical_lines = _read_mrbayes_tabular_lines(path)
    header_fields, data_lines = _split_mrbayes_tabular_table(
        logical_lines,
        path=path,
        artifact_kind="mrbayes-trace",
        missing_header_code="mrbayes_trace_missing_header",
        missing_generation_code="mrbayes_trace_missing_generation_column",
        unexpected_field_code="mrbayes_trace_unexpected_field_count",
    )
    generation_field = _mrbayes_generation_field(header_fields)
    columns = [field for field in header_fields if field and field != generation_field]
    for row_number, fields in enumerate(data_lines, start=2):
        raw_row = _build_mrbayes_tabular_row(header_fields, fields)
        generation_text = _normalize_tabular_field(raw_row.get(generation_field))
        if generation_text is None:
            raise _mrbayes_artifact_error(
                f"MrBayes trace file lacks a Gen column: {path}",
                code="mrbayes_trace_missing_generation_column",
                path=path,
                artifact_kind="mrbayes-trace",
                details={"expected_section": "Gen column"},
            )
        try:
            generation = int(float(generation_text))
        except ValueError as error:
            raise _mrbayes_artifact_error(
                f"MrBayes trace file contains a non-numeric generation value on row {row_number}: {path}",
                code="mrbayes_trace_invalid_generation_value",
                path=path,
                artifact_kind="mrbayes-trace",
                details={"row_number": row_number, "expected_section": "Gen column"},
            ) from error
        values: dict[str, float] = {}
        for column in columns:
            raw_value = _normalize_tabular_field(raw_row.get(column))
            if raw_value in {None, ""}:
                raise _mrbayes_artifact_error(
                    f"MrBayes trace file is missing a sampled value for '{column}' on row {row_number}: {path}",
                    code="mrbayes_trace_missing_parameter_value",
                    path=path,
                    artifact_kind="mrbayes-trace",
                    details={
                        "row_number": row_number,
                        "column": column,
                        "expected_section": "sampled parameter row",
                    },
                )
            try:
                values[column] = float(raw_value)
            except ValueError as error:
                raise _mrbayes_artifact_error(
                    f"MrBayes trace file contains a non-numeric value for '{column}' on row {row_number}: {path}",
                    code="mrbayes_trace_invalid_parameter_value",
                    path=path,
                    artifact_kind="mrbayes-trace",
                    details={
                        "row_number": row_number,
                        "column": column,
                        "expected_section": "sampled parameter row",
                    },
                ) from error
        rows.append(MrBayesTraceRow(generation=generation, values=values))
    if not rows:
        raise _mrbayes_artifact_error(
            f"MrBayes trace file contains no sampled rows: {path}",
            code="mrbayes_trace_missing_rows",
            path=path,
            artifact_kind="mrbayes-trace",
            details={"expected_section": "sampled parameter rows"},
        )
    return MrBayesTraceReport(
        path=path, row_count=len(rows), columns=columns, rows=rows
    )


def parse_mrbayes_mcmc_diagnostics(path: Path) -> MrBayesMcmcReport:
    """Parse a MrBayes .mcmc diagnostics table into deterministic rows."""
    if not path.exists():
        raise _mrbayes_artifact_error(
            f"MrBayes MCMC diagnostics file was not found: {path}",
            code="mrbayes_mcmc_missing_file",
            path=path,
            artifact_kind="mrbayes-mcmc",
            details={"expected_section": "MCMC diagnostics file"},
        )
    comment_lines: list[str] = []
    logical_lines = _read_mrbayes_tabular_lines(path, comment_lines=comment_lines)
    header_fields, data_lines = _split_mrbayes_tabular_table(
        logical_lines,
        path=path,
        artifact_kind="mrbayes-mcmc",
        missing_header_code="mrbayes_mcmc_missing_header",
        missing_generation_code="mrbayes_mcmc_missing_generation_column",
        unexpected_field_code="mrbayes_mcmc_unexpected_field_count",
    )
    generation_field = _mrbayes_generation_field(header_fields)
    columns = [field for field in header_fields if field and field != generation_field]
    rows: list[MrBayesMcmcRow] = []
    for row_number, fields in enumerate(data_lines, start=2):
        raw_row = _build_mrbayes_tabular_row(header_fields, fields)
        generation_text = _normalize_tabular_field(raw_row.get(generation_field))
        if generation_text is None:
            raise _mrbayes_artifact_error(
                f"MrBayes MCMC diagnostics file lacks a Gen column: {path}",
                code="mrbayes_mcmc_missing_generation_column",
                path=path,
                artifact_kind="mrbayes-mcmc",
                details={"expected_section": "Gen column"},
            )
        try:
            generation = int(float(generation_text))
        except ValueError as error:
            raise _mrbayes_artifact_error(
                f"MrBayes MCMC diagnostics file contains a non-numeric generation value on row {row_number}: {path}",
                code="mrbayes_mcmc_invalid_generation_value",
                path=path,
                artifact_kind="mrbayes-mcmc",
                details={"row_number": row_number, "expected_section": "Gen column"},
            ) from error
        values: dict[str, float | None] = {}
        for column in columns:
            raw_value = _normalize_tabular_field(raw_row.get(column))
            if raw_value in {None, ""}:
                raise _mrbayes_artifact_error(
                    f"MrBayes MCMC diagnostics file is missing a sampled value for '{column}' on row {row_number}: {path}",
                    code="mrbayes_mcmc_missing_parameter_value",
                    path=path,
                    artifact_kind="mrbayes-mcmc",
                    details={
                        "row_number": row_number,
                        "column": column,
                        "expected_section": "sampled diagnostics row",
                    },
                )
            normalized = raw_value.strip()
            if normalized.lower() in {"na", "nan"}:
                values[column] = None
            else:
                try:
                    values[column] = float(normalized)
                except ValueError as error:
                    raise _mrbayes_artifact_error(
                        f"MrBayes MCMC diagnostics file contains a non-numeric value for '{column}' on row {row_number}: {path}",
                        code="mrbayes_mcmc_invalid_parameter_value",
                        path=path,
                        artifact_kind="mrbayes-mcmc",
                        details={
                            "row_number": row_number,
                            "column": column,
                            "expected_section": "sampled diagnostics row",
                        },
                    ) from error
        rows.append(MrBayesMcmcRow(generation=generation, values=values))
    if not rows:
        raise _mrbayes_artifact_error(
            f"MrBayes MCMC diagnostics file contains no sampled rows: {path}",
            code="mrbayes_mcmc_missing_rows",
            path=path,
            artifact_kind="mrbayes-mcmc",
            details={"expected_section": "sampled diagnostics rows"},
        )
    return MrBayesMcmcReport(
        path=path,
        row_count=len(rows),
        columns=columns,
        comment_lines=comment_lines,
        rows=rows,
    )


def _read_mrbayes_tabular_lines(
    path: Path,
    *,
    comment_lines: list[str] | None = None,
) -> list[str]:
    lines: list[str] = []
    with path.open(encoding="utf-8", newline="") as handle:
        for raw_line in handle:
            stripped = raw_line.strip()
            if not stripped:
                continue
            if stripped.startswith("["):
                if comment_lines is not None:
                    comment_lines.append(stripped)
                continue
            lines.append(raw_line.rstrip("\r\n"))
    return lines


def _normalize_tabular_field(value: str | None) -> str | None:
    if value is None:
        return None
    return value.lstrip("\ufeff").strip()


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


def _mrbayes_generation_field(fieldnames: list[str]) -> str:
    for fieldname in fieldnames:
        normalized = (_normalize_tabular_field(fieldname) or "").lower()
        if normalized == "gen":
            return fieldname
    return "Gen"


def _build_mrbayes_tabular_row(
    header_fields: list[str],
    fields: list[str],
) -> dict[str, str | None]:
    return {
        header_fields[index]: (fields[index] if index < len(fields) else None)
        for index in range(len(header_fields))
    }


def _split_mrbayes_tabular_table(
    logical_lines: list[str],
    *,
    path: Path,
    artifact_kind: str,
    missing_header_code: str,
    missing_generation_code: str,
    unexpected_field_code: str,
) -> tuple[list[str], list[list[str]]]:
    if not logical_lines:
        raise _mrbayes_artifact_error(
            f"MrBayes {artifact_kind} file contains no header: {path}",
            code=missing_header_code,
            path=path,
            artifact_kind=artifact_kind,
            details={"expected_section": "tabular header"},
        )
    header_fields: list[str] | None = None
    header_index = -1
    for index, line in enumerate(logical_lines):
        if _is_tabular_warning_line(line):
            continue
        fields = _split_tabular_fields(line)
        if any(
            (_normalize_tabular_field(field) or "").lower() == "gen" for field in fields
        ):
            header_fields = fields
            header_index = index
            break
        if len(fields) > 1:
            raise _mrbayes_artifact_error(
                f"MrBayes {artifact_kind} file lacks a Gen column: {path}",
                code=missing_generation_code,
                path=path,
                artifact_kind=artifact_kind,
                details={"columns": fields, "expected_section": "Gen column"},
            )
    if header_fields is None:
        raise _mrbayes_artifact_error(
            f"MrBayes {artifact_kind} file contains no header: {path}",
            code=missing_header_code,
            path=path,
            artifact_kind=artifact_kind,
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
            raise _mrbayes_artifact_error(
                f"MrBayes {artifact_kind} file contains more fields than its header: {path}",
                code=unexpected_field_code,
                path=path,
                artifact_kind=artifact_kind,
                details={
                    "expected_field_count": len(header_fields),
                    "observed_field_count": len(fields),
                    "expected_section": "sampled row",
                },
            )
        data_lines.append(fields)
    return header_fields, data_lines
