# ruff: noqa: F401
from __future__ import annotations

from copy import deepcopy
import csv
from dataclasses import dataclass
from pathlib import Path
import re
from typing import Any, TypeAlias
from xml.etree import ElementTree  # nosec B405

from defusedxml import ElementTree as SafeXmlET
from defusedxml.common import DefusedXmlException

from bijux_phylogenetics.bayesian.posterior_sets.burnin import (
    DEFAULT_BURNIN_FRACTIONS,
    BurninSensitivityCladeShift,
    BurninSensitivityParameterShift,
    normalize_burnin_fractions,
    summarize_burnin_clade_shifts,
    summarize_burnin_parameter_shifts,
)
from bijux_phylogenetics.bayesian.posterior_sets.diagnostics import (
    TraceConvergenceReport,
    summarize_trace_convergence,
)
from bijux_phylogenetics.bayesian.posterior_sets.tree_sets import (
    summarize_maximum_clade_credibility_tree,
)
from bijux_phylogenetics.datasets.study_inputs import load_taxon_table, write_taxon_rows
from bijux_phylogenetics.diagnostics.validation import validate_tree_path
from bijux_phylogenetics.engines.common import (
    build_file_checksums,
    execute_engine_command,
    read_engine_version,
    resolve_engine_executable,
    validate_timeout_seconds,
)
from bijux_phylogenetics.engines.workflows.models import EngineWorkflowReport
from bijux_phylogenetics.engines.workflows.state import (
    _ensure_inference_ready_alignment,
    _persist_workflow_report,
    _record_output_validation_failure,
    _resolve_incomplete_workflow_state,
    _resume_existing_workflow,
)
from bijux_phylogenetics.io.biopython import loads_biophylo
from bijux_phylogenetics.io.fasta import (
    infer_alignment_alphabet,
    load_fasta_alignment,
)
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, descendant_taxa
from bijux_phylogenetics.runtime.errors import (
    EngineWorkflowError,
    InvalidAlignmentError,
    PhylogeneticsError,
    TreeParseError,
)
from bijux_phylogenetics.trees import (
    compute_clade_frequency_table,
    compute_consensus_tree,
    load_tree_set,
    summarize_posterior_topology_diversity,
)

from .models import (
    BeastLogParameterSummary,
    BeastLogSummaryReport,
    CalibrationValidationIssue,
    TipDatingValidationReport,
)

_BEAST_TREE_PATTERN = re.compile(
    r"tree\s+([^\s=]+)\s*=\s*(.+?);", flags=re.IGNORECASE | re.DOTALL
)
_BEAST_TREE_STATE_PATTERN = re.compile(r"STATE_(\d+)$", flags=re.IGNORECASE)
_TABULAR_WARNING_PREFIX_PATTERN = re.compile(
    r"^(warning|warn|caution|note|info)\b",
    flags=re.IGNORECASE,
)
XmlElement: TypeAlias = Any


def _beast_artifact_error(
    message: str,
    *,
    code: str,
    path: Path,
    artifact_kind: str,
    details: dict[str, object] | None = None,
) -> EngineWorkflowError:
    payload: dict[str, object] = {
        "path": str(path),
        "artifact_kind": artifact_kind,
    }
    if details is not None:
        payload.update(details)
    return EngineWorkflowError(message, code=code, details=payload)


def _beast_output_path(xml_path: Path, *, seed: int, suffix: str) -> Path:
    return xml_path.with_name(f"{xml_path.stem}.{seed}.{suffix}")


def _read_delimited_rows(path: Path) -> list[dict[str, str]]:
    delimiter = "," if path.suffix.lower() == ".csv" else "\t"
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ValueError(f"calibration table contains no header: {path}")
        return [
            {key: value or "" for key, value in row.items() if key is not None}
            for row in reader
        ]


def _named_clades(tree: PhyloTree) -> dict[str, list[str]]:
    named: dict[str, list[str]] = {}
    for node in tree.iter_nodes():
        if node.name:
            named[node.name] = descendant_taxa(node)
    return named


def _clade_taxon_sets(tree: PhyloTree) -> set[frozenset[str]]:
    clades: set[frozenset[str]] = set()
    for node in tree.iter_nodes():
        taxa = descendant_taxa(node)
        if taxa:
            clades.add(frozenset(taxa))
    return clades


def _parse_target_taxa(raw: str) -> list[str]:
    if not raw.strip():
        return []
    normalized = raw.replace(",", "|").replace(";", "|")
    return sorted({token.strip() for token in normalized.split("|") if token.strip()})


def _parse_age(
    raw: str,
    *,
    calibration_id: str,
    field_name: str,
    issues: list[CalibrationValidationIssue],
) -> float | None:
    if not raw.strip():
        return None
    try:
        value = float(raw)
    except ValueError:
        issues.append(
            CalibrationValidationIssue(
                calibration_id=calibration_id,
                code="invalid-age",
                message=f"{field_name} must be numeric when provided",
            )
        )
        return None
    return value


_XML_IDENTIFIER_PATTERN = re.compile(r"[^0-9A-Za-z._-]+")


def _format_decimal(value: float) -> str:
    return format(value, ".15g")


def _normalize_tabular_field(value: str | None) -> str | None:
    if value is None:
        return None
    return value.lstrip("\ufeff").strip()


def _beast_state_field(fieldnames: list[str]) -> str | None:
    for fieldname in fieldnames:
        normalized = (_normalize_tabular_field(fieldname) or "").lower()
        if normalized in {"state", "sample"}:
            return fieldname
    return None


def _tip_date_trait_value(report: TipDatingValidationReport) -> str:
    parts = [
        f"{tip.taxon}={_format_decimal(tip.date)}"
        for tip in report.tip_dates
        if tip.valid and tip.date is not None
    ]
    return ",".join(parts)


def _tree_root_age(tree_path: Path) -> float:
    tree = load_tree(tree_path)
    lengths = [length for length in tree.root_to_tip_lengths() if length is not None]
    if not lengths:
        return 0.0
    return round(max(float(length) for length in lengths), 15)


def _mean_or_none(values: list[float]) -> float | None:
    if not values:
        return None
    return round(sum(values) / len(values), 6)


def _mean_beast_parameter(
    report: BeastLogSummaryReport,
    parameter: str,
) -> float | None:
    for summary in report.parameter_summaries:
        if summary.parameter == parameter:
            return summary.mean
    return None


def _classify_beast_parameter(parameter: str) -> str:
    normalized = parameter.lower()
    if normalized == "posterior" or normalized.endswith(".posterior"):
        return "posterior"
    if normalized == "likelihood" or "likelihood" in normalized:
        return "likelihood"
    if normalized == "prior" or "prior" in normalized:
        return "prior"
    if any(
        token in normalized
        for token in (
            "clock",
            "clockrate",
            "ucld",
            "branchrate",
            "mutationrate",
            "meanrate",
        )
    ):
        return "clock"
    if any(
        token in normalized
        for token in (
            "tree",
            "birthrate",
            "deathrate",
            "popsize",
            "coalescent",
            "tmrca",
            "origin",
        )
    ):
        return "tree"
    return "other"


def _summary_parameters_by_category(
    parameter_summaries: list[BeastLogParameterSummary],
    *,
    category: str,
) -> list[str]:
    return [
        summary.parameter
        for summary in parameter_summaries
        if summary.parameter_category == category
    ]
