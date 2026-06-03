from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.diagnostics.validation import validate_tree_path
from bijux_phylogenetics.io.fasta import load_fasta_alignment
from bijux_phylogenetics.io.fasta.records import summarise_fasta
from bijux_phylogenetics.io.newick import loads_newick
from bijux_phylogenetics.io.trees import load_tree
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


def _validate_complete_support_coverage(
    *,
    engine_name: str,
    workflow: str,
    path: Path,
    output_name: str,
    artifact_kind: str,
    annotated_branch_count: int,
    support_kind: str,
) -> None:
    expected_branch_count = _count_supported_internal_branches(path)
    if annotated_branch_count == expected_branch_count:
        return
    raise build_engine_output_error(
        f"{engine_name} {workflow} exposed incomplete {support_kind} coverage in '{output_name}': {path}",
        code="engine_support_values_incomplete",
        engine_name=engine_name,
        workflow=workflow,
        path=path,
        output_name=output_name,
        artifact_kind=artifact_kind,
        details={
            "support_kind": support_kind,
            "expected_supported_branch_count": expected_branch_count,
            "observed_supported_branch_count": annotated_branch_count,
        },
    )


def _validate_matching_tree_taxa(
    *,
    engine_name: str,
    workflow: str,
    reference_tree_path: Path,
    comparison_tree_set_path: Path,
    reference_output_name: str,
    comparison_output_name: str,
    artifact_kind: str,
) -> None:
    reference_tree = load_tree(reference_tree_path)
    tree_set_report = load_tree_set(comparison_tree_set_path)
    expected_taxa = sorted(reference_tree.tip_names)
    if (
        tree_set_report.shared_taxa == expected_taxa
        and tree_set_report.taxa_union == expected_taxa
    ):
        return
    raise build_engine_output_error(
        f"{engine_name} {workflow} produced inconsistent taxa across '{reference_output_name}' and '{comparison_output_name}'",
        code="engine_output_taxa_mismatch",
        engine_name=engine_name,
        workflow=workflow,
        path=comparison_tree_set_path,
        output_name=comparison_output_name,
        artifact_kind=artifact_kind,
        details={
            "reference_tree_path": str(reference_tree_path),
            "reference_output_name": reference_output_name,
            "expected_taxa": expected_taxa,
            "shared_tree_set_taxa": tree_set_report.shared_taxa,
            "tree_set_taxa_union": tree_set_report.taxa_union,
        },
    )


def _count_supported_internal_branches(path: Path) -> int:
    tree = load_tree(path)
    return sum(
        1 for node in tree.iter_nodes() if node is not tree.root and not node.is_leaf()
    )


def _ensure_inference_ready_alignment(path: Path) -> None:
    records = load_fasta_alignment(path)
    if not records or len(records[0].sequence) < 1:
        raise EngineWorkflowError(
            f"inference alignment is empty after filtering: {path}"
        )
    summarise_fasta(path)
