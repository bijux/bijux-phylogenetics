from __future__ import annotations

from dataclasses import MISSING, dataclass, fields, is_dataclass
import json
from pathlib import Path
import types
from typing import Any, Union, get_args, get_origin, get_type_hints

from bijux_phylogenetics.phylo.likelihood.models import (
    LikelihoodOptimizationBoundaryWarning,
    LikelihoodWrapperCorrespondenceReport,
    LikelihoodWrapperCorrespondenceSummaryRow,
    NucleotideLikelihoodBootstrapTreeInferenceReport,
    NucleotideLikelihoodMultiStartRunSummary,
    NucleotideLikelihoodNniSearchReport,
    NucleotideLikelihoodNniTraceRow,
    NucleotideLikelihoodSprSearchReport,
    NucleotideLikelihoodSprTraceRow,
    NucleotideLikelihoodTbrSearchReport,
    NucleotideLikelihoodTbrTraceRow,
    NucleotideLikelihoodTreeInferenceReport,
    NucleotideShLikeBranchSupportReport,
    SubstitutionModelSelectionRow,
)


@dataclass(slots=True)
class NucleotideMaximumLikelihoodSupportReports:
    """Optional branch-support reports attached to one native ML tree result."""

    bootstrap: NucleotideLikelihoodBootstrapTreeInferenceReport | None = None
    sh_like: NucleotideShLikeBranchSupportReport | None = None


@dataclass(slots=True)
class NucleotideMaximumLikelihoodBenchmarkMetadata:
    """Stable benchmark metadata attached to one native ML tree result."""

    wrapper_correspondence_case_count: int
    wrapper_correspondence_supported_case_count: int
    wrapper_correspondence_blocking_case_count: int
    wrapper_correspondence_all_supported_cases_clear: bool
    wrapper_correspondence_summary_rows: list[LikelihoodWrapperCorrespondenceSummaryRow]


@dataclass(slots=True)
class NucleotideMaximumLikelihoodResult:
    """Unified native ML tree result with inference, traces, support, and warnings."""

    inference_report: NucleotideLikelihoodTreeInferenceReport
    final_tree_newick: str
    final_topology_fingerprint: str
    final_log_likelihood: float
    model_name: str
    model_selection_strategy: str
    model_selection_criterion: str | None
    parameter_values: dict[str, float]
    search_method: str
    run_summaries: list[NucleotideLikelihoodMultiStartRunSummary]
    search_trace_rows: list[
        NucleotideLikelihoodNniTraceRow
        | NucleotideLikelihoodSprTraceRow
        | NucleotideLikelihoodTbrTraceRow
    ]
    support_reports: NucleotideMaximumLikelihoodSupportReports
    warning_messages: list[str]
    benchmark_metadata: NucleotideMaximumLikelihoodBenchmarkMetadata | None = None

    def to_dict(self) -> dict[str, object]:
        """Serialize the result into one JSON-ready payload."""
        payload = _serialize_value(self)
        if not isinstance(payload, dict):
            raise TypeError("maximum-likelihood result payload must serialize as a dict")
        payload["result_type"] = self.__class__.__name__
        return payload

    @classmethod
    def from_dict(cls, payload: dict[str, object]) -> NucleotideMaximumLikelihoodResult:
        """Rebuild one ML result from a `to_dict()` payload."""
        result_type = payload.get("result_type")
        if result_type != cls.__name__:
            raise ValueError(
                f"expected result_type {cls.__name__!r}, received {result_type!r}"
            )
        body = dict(payload)
        body.pop("result_type", None)
        return _restore_dataclass(cls, body)

    def write_json(self, path: Path) -> Path:
        """Write the ML result payload to one governed JSON path."""
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            json.dumps(self.to_dict(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return path

    @classmethod
    def load_json(cls, path: Path) -> NucleotideMaximumLikelihoodResult:
        """Load one ML result from a governed JSON path."""
        payload = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(payload, dict):
            raise ValueError("maximum-likelihood result JSON must contain one object")
        return cls.from_dict(payload)


def build_nucleotide_maximum_likelihood_result(
    inference_report: NucleotideLikelihoodTreeInferenceReport,
    *,
    support_reports: NucleotideMaximumLikelihoodSupportReports | None = None,
    benchmark_metadata: NucleotideMaximumLikelihoodBenchmarkMetadata | None = None,
) -> NucleotideMaximumLikelihoodResult:
    """Build one unified native ML result from a completed tree-inference report."""
    selected_model_row = _selected_model_row(inference_report)
    resolved_support_reports = (
        NucleotideMaximumLikelihoodSupportReports()
        if support_reports is None
        else _validate_support_reports(inference_report, support_reports)
    )
    return NucleotideMaximumLikelihoodResult(
        inference_report=inference_report,
        final_tree_newick=inference_report.best_final_tree_newick,
        final_topology_fingerprint=inference_report.best_final_topology_fingerprint,
        final_log_likelihood=inference_report.best_final_log_likelihood,
        model_name=inference_report.selected_model_name,
        model_selection_strategy=inference_report.model_selection_strategy,
        model_selection_criterion=inference_report.model_selection_criterion,
        parameter_values=dict(selected_model_row.parameter_values),
        search_method=inference_report.search_method,
        run_summaries=list(inference_report.run_summaries),
        search_trace_rows=list(inference_report.best_search_report.trace_rows),
        support_reports=resolved_support_reports,
        warning_messages=_collect_inference_warning_messages(
            inference_report,
            selected_model_row=selected_model_row,
            support_reports=resolved_support_reports,
        ),
        benchmark_metadata=benchmark_metadata,
    )


def build_nucleotide_maximum_likelihood_support_reports(
    inference_report: NucleotideLikelihoodTreeInferenceReport,
    *,
    bootstrap_report: NucleotideLikelihoodBootstrapTreeInferenceReport | None = None,
    sh_like_report: NucleotideShLikeBranchSupportReport | None = None,
) -> NucleotideMaximumLikelihoodSupportReports:
    """Build one validated support bundle for a completed native ML tree result."""
    return _validate_support_reports(
        inference_report,
        NucleotideMaximumLikelihoodSupportReports(
            bootstrap=bootstrap_report,
            sh_like=sh_like_report,
        ),
    )


def build_nucleotide_maximum_likelihood_benchmark_metadata(
    correspondence_report: LikelihoodWrapperCorrespondenceReport,
) -> NucleotideMaximumLikelihoodBenchmarkMetadata:
    """Condense one wrapper-correspondence report into stable ML benchmark metadata."""
    return NucleotideMaximumLikelihoodBenchmarkMetadata(
        wrapper_correspondence_case_count=correspondence_report.case_count,
        wrapper_correspondence_supported_case_count=(
            correspondence_report.supported_case_count
        ),
        wrapper_correspondence_blocking_case_count=(
            correspondence_report.blocking_case_count
        ),
        wrapper_correspondence_all_supported_cases_clear=(
            correspondence_report.all_supported_cases_clear
        ),
        wrapper_correspondence_summary_rows=list(
            correspondence_report.summary_rows
        ),
    )


def write_nucleotide_maximum_likelihood_result_json(
    path: Path,
    result: NucleotideMaximumLikelihoodResult,
) -> Path:
    """Write one machine-readable native maximum-likelihood result payload."""
    return result.write_json(path)


def load_nucleotide_maximum_likelihood_result_json(
    path: Path,
) -> NucleotideMaximumLikelihoodResult:
    """Load one machine-readable native maximum-likelihood result payload."""
    return NucleotideMaximumLikelihoodResult.load_json(path)


def _selected_model_row(
    inference_report: NucleotideLikelihoodTreeInferenceReport,
) -> SubstitutionModelSelectionRow:
    for row in inference_report.model_selection_report.rows:
        if row.model_name == inference_report.selected_model_name:
            return row
    raise ValueError(
        "native maximum-likelihood result could not resolve the selected model row"
    )


def _collect_inference_warning_messages(
    inference_report: NucleotideLikelihoodTreeInferenceReport,
    *,
    selected_model_row: SubstitutionModelSelectionRow,
    support_reports: NucleotideMaximumLikelihoodSupportReports,
) -> list[str]:
    warnings: list[str] = []
    warnings.extend(inference_report.model_selection_report.warnings)
    warnings.extend(selected_model_row.warnings)
    warnings.extend(
        warning.message
        for warning in selected_model_row.boundary_warning_rows
        if isinstance(warning, LikelihoodOptimizationBoundaryWarning)
    )
    warnings.extend(inference_report.best_search_report.substitution_parameter_warnings)
    for row in inference_report.best_search_report.trace_rows:
        warnings.extend(row.boundary_warning_messages)
    if support_reports.sh_like is not None:
        warnings.append(support_reports.sh_like.caution_label)
    return _unique_in_order(warnings)


def _validate_support_reports(
    inference_report: NucleotideLikelihoodTreeInferenceReport,
    support_reports: NucleotideMaximumLikelihoodSupportReports,
) -> NucleotideMaximumLikelihoodSupportReports:
    if support_reports.bootstrap is not None:
        _validate_bootstrap_report(inference_report, support_reports.bootstrap)
    if support_reports.sh_like is not None:
        _validate_sh_like_report(inference_report, support_reports.sh_like)
    return support_reports


def _validate_bootstrap_report(
    inference_report: NucleotideLikelihoodTreeInferenceReport,
    bootstrap_report: NucleotideLikelihoodBootstrapTreeInferenceReport,
) -> None:
    if bootstrap_report.reference_tree_newick != inference_report.best_final_tree_newick:
        raise ValueError(
            "bootstrap support report tree does not match the inference result tree"
        )
    if (
        bootstrap_report.selected_reference_model_name
        != inference_report.selected_model_name
    ):
        raise ValueError(
            "bootstrap support report model does not match the inference result model"
        )


def _validate_sh_like_report(
    inference_report: NucleotideLikelihoodTreeInferenceReport,
    sh_like_report: NucleotideShLikeBranchSupportReport,
) -> None:
    if sh_like_report.reference_tree_newick != inference_report.best_final_tree_newick:
        raise ValueError(
            "SH-like support report tree does not match the inference result tree"
        )
    if sh_like_report.model_name.upper() != inference_report.selected_model_name.upper():
        raise ValueError(
            "SH-like support report model does not match the inference result model"
        )


def _unique_in_order(values: list[str]) -> list[str]:
    seen: set[str] = set()
    ordered: list[str] = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        ordered.append(value)
    return ordered


def _serialize_value(value: object) -> object:
    if isinstance(value, Path):
        return str(value)
    if is_dataclass(value):
        return {
            field.name: _serialize_value(getattr(value, field.name))
            for field in fields(value)
        }
    if isinstance(value, dict):
        return {str(key): _serialize_value(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_serialize_value(item) for item in value]
    if isinstance(value, tuple):
        return [_serialize_value(item) for item in value]
    return value


def _restore_dataclass(dataclass_type: type[Any], payload: object) -> Any:
    if not isinstance(payload, dict):
        raise ValueError(f"{dataclass_type.__name__} payload must be a dict")
    dataclass_fields = fields(dataclass_type)
    field_names = {field.name for field in dataclass_fields}
    unexpected_names = sorted(set(payload) - field_names)
    if unexpected_names:
        raise ValueError(
            f"{dataclass_type.__name__} payload contained unexpected fields: "
            + ", ".join(unexpected_names)
        )
    missing_names = [
        field.name
        for field in dataclass_fields
        if field.name not in payload
        and field.default is MISSING
        and field.default_factory is MISSING
    ]
    if missing_names:
        raise ValueError(
            f"{dataclass_type.__name__} payload is missing required fields: "
            + ", ".join(sorted(missing_names))
        )
    type_hints = _type_hints(dataclass_type)
    restored_fields = {
        field.name: _restore_value(type_hints[field.name], payload[field.name])
        for field in dataclass_fields
        if field.name in payload
    }
    return dataclass_type(**restored_fields)


def _restore_value(type_hint: Any, payload: object) -> object:
    origin = get_origin(type_hint)
    if origin is list:
        (inner_type,) = get_args(type_hint)
        if not isinstance(payload, list):
            raise ValueError("list payload must be a JSON array")
        return [_restore_value(inner_type, item) for item in payload]
    if origin is dict:
        _, value_type = get_args(type_hint)
        if not isinstance(payload, dict):
            raise ValueError("dict payload must be a JSON object")
        return {
            str(key): _restore_value(value_type, value)
            for key, value in payload.items()
        }
    if origin in (Union, types.UnionType):
        return _restore_union(get_args(type_hint), payload)
    if isinstance(type_hint, type) and is_dataclass(type_hint):
        return _restore_dataclass(type_hint, payload)
    return payload


def _restore_union(candidate_types: tuple[Any, ...], payload: object) -> object:
    if payload is None:
        if type(None) in candidate_types:
            return None
        raise ValueError("received null for non-optional union value")
    last_error: ValueError | None = None
    for candidate_type in candidate_types:
        if candidate_type is type(None):
            continue
        try:
            return _restore_value(candidate_type, payload)
        except ValueError as error:
            last_error = error
            continue
    if last_error is not None:
        raise last_error
    raise ValueError("could not restore union payload")


def _type_hints(dataclass_type: type[Any]) -> dict[str, Any]:
    return get_type_hints(dataclass_type)
