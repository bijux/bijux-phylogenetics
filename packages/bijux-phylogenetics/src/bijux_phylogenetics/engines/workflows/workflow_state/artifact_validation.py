from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.diagnostics.validation import validate_tree_path
from bijux_phylogenetics.io.fasta import load_fasta_alignment
from bijux_phylogenetics.io.fasta.records import summarise_fasta
from bijux_phylogenetics.io.newick import loads_newick
from bijux_phylogenetics.phylo.alignment import AlignmentSummary
from bijux_phylogenetics.runtime.errors import EngineWorkflowError
from bijux_phylogenetics.trees import load_tree_set

from ...common import build_engine_output_error


def _require_nonempty_text_output(
    path: Path,
    *,
    engine_name: str,
    workflow: str,
    output_name: str,
    artifact_kind: str,
) -> str:
    if not path.exists():
        raise build_engine_output_error(
            f"{engine_name} {workflow} did not produce required output '{output_name}': {path}",
            code="engine_required_output_missing",
            engine_name=engine_name,
            workflow=workflow,
            path=path,
            output_name=output_name,
            artifact_kind=artifact_kind,
        )
    text = path.read_text(encoding="utf-8", errors="replace")
    if not text.strip():
        raise build_engine_output_error(
            f"{engine_name} {workflow} produced an empty required output '{output_name}': {path}",
            code="engine_output_empty",
            engine_name=engine_name,
            workflow=workflow,
            path=path,
            output_name=output_name,
            artifact_kind=artifact_kind,
        )
    return text


def _validate_alignment_output(
    path: Path,
    *,
    engine_name: str,
    workflow: str,
    output_name: str,
    artifact_kind: str,
) -> AlignmentSummary:
    _require_nonempty_text_output(
        path,
        engine_name=engine_name,
        workflow=workflow,
        output_name=output_name,
        artifact_kind=artifact_kind,
    )
    records = load_fasta_alignment(path)
    if not records or len(records[0].sequence) < 1:
        raise build_engine_output_error(
            f"{engine_name} {workflow} produced an empty alignment after validation: {path}",
            code="engine_output_empty",
            engine_name=engine_name,
            workflow=workflow,
            path=path,
            output_name=output_name,
            artifact_kind=artifact_kind,
        )
    return summarise_fasta(path)


def _validate_tree_output(
    path: Path,
    *,
    engine_name: str,
    workflow: str,
    output_name: str,
    artifact_kind: str,
) -> None:
    tree_text = _require_nonempty_text_output(
        path,
        engine_name=engine_name,
        workflow=workflow,
        output_name=output_name,
        artifact_kind=artifact_kind,
    )
    loads_newick(tree_text)
    validate_tree_path(path)


def _validate_tree_set_output(
    path: Path,
    *,
    engine_name: str,
    workflow: str,
    output_name: str,
    artifact_kind: str,
) -> None:
    _require_nonempty_text_output(
        path,
        engine_name=engine_name,
        workflow=workflow,
        output_name=output_name,
        artifact_kind=artifact_kind,
    )
    load_tree_set(path)


def _validate_iqtree_required_artifacts(
    prefix_path: Path,
    *,
    workflow: str,
) -> None:
    _require_nonempty_text_output(
        prefix_path.with_suffix(".iqtree"),
        engine_name="iqtree",
        workflow=workflow,
        output_name="iqtree_report",
        artifact_kind="iqtree-report",
    )
    _require_nonempty_text_output(
        prefix_path.with_suffix(".log"),
        engine_name="iqtree",
        workflow=workflow,
        output_name="iqtree_log",
        artifact_kind="iqtree-log",
    )


def _validate_support_value_count(
    *,
    engine_name: str,
    workflow: str,
    path: Path,
    output_name: str,
    artifact_kind: str,
    support_value_count: int,
    support_kind: str,
) -> None:
    if support_value_count > 0:
        return
    raise build_engine_output_error(
        f"{engine_name} {workflow} did not expose any parsable {support_kind} values in '{output_name}': {path}",
        code="engine_support_values_missing",
        engine_name=engine_name,
        workflow=workflow,
        path=path,
        output_name=output_name,
        artifact_kind=artifact_kind,
        details={"support_kind": support_kind},
    )


def _ensure_inference_ready_alignment(path: Path) -> None:
    records = load_fasta_alignment(path)
    if not records or len(records[0].sequence) < 1:
        raise EngineWorkflowError(
            f"inference alignment is empty after filtering: {path}"
        )
    summarise_fasta(path)
