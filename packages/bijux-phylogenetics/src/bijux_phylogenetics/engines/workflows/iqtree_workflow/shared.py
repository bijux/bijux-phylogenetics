from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.io.fasta import infer_alignment_alphabet, load_fasta_alignment
from bijux_phylogenetics.phylo.alignment import AlignmentAlphabet
from bijux_phylogenetics.runtime.errors import EngineWorkflowError

from ...artifacts.iqtree import (
    IqtreeModelSelectionSummary,
    parse_best_model_file,
    parse_iqtree_model_selection_summary,
    parse_log_likelihood_file,
    resolve_iqtree_model_sidecar,
)
from ...common import build_engine_output_error
from ...validation import summarize_bootstrap_support_distribution
from ..models import IqtreeSupportValue, IqtreeWorkflowSummary

_MINIMUM_UFBOOT_REPLICATES = 1000


def _iqtree_partition_supports_fixed_model(
    *,
    model: str,
    mixed_data_types: bool,
) -> bool:
    if not mixed_data_types:
        return True
    normalized = model.strip().upper()
    return normalized in {
        "TEST",
        "TESTONLY",
        "TESTNEWONLY",
        "MF",
        "MFP",
        "TESTMERGE",
        "TESTMERGEONLY",
        "MF+MERGE",
        "MFP+MERGE",
    }


def _iqtree_sequence_type_flag(
    path: Path, sequence_type: AlignmentAlphabet | None
) -> list[str]:
    detected = sequence_type
    if detected is None:
        detected = infer_alignment_alphabet(load_fasta_alignment(path))
    if detected in {"dna", "rna"}:
        return ["-st", "DNA"]
    if detected == "protein":
        return ["-st", "AA"]
    return []


def _iqtree_execution_controls(*, seed: int, threads: int) -> list[str]:
    if seed < 1:
        raise ValueError(f"iqtree seed must be positive, got {seed}")
    if threads < 1:
        raise ValueError(f"iqtree threads must be positive, got {threads}")
    return ["-seed", str(seed), "-nt", str(threads)]


def _validate_ufboot_replicates(replicates: int) -> None:
    if replicates < _MINIMUM_UFBOOT_REPLICATES:
        raise EngineWorkflowError(
            "iqtree ultrafast bootstrap requires at least "
            f"{_MINIMUM_UFBOOT_REPLICATES} replicates, got {replicates}"
        )


def _validate_sh_alrt_replicates(replicates: int) -> None:
    if replicates < 1:
        raise ValueError(f"sh-alrt replicates must be positive, got {replicates}")


def _parse_best_model_artifact(prefix_path: Path) -> str | None:
    for candidate in (
        prefix_path.with_suffix(".iqtree"),
        resolve_iqtree_model_sidecar(prefix_path),
    ):
        if candidate is None:
            continue
        model = parse_best_model_file(candidate)
        if model is not None:
            return model
    return None


def _parse_log_likelihood_artifact(prefix_path: Path) -> float | None:
    for candidate in (
        prefix_path.with_suffix(".iqtree"),
        prefix_path.with_suffix(".log"),
    ):
        if not candidate.exists():
            continue
        log_likelihood = parse_log_likelihood_file(candidate)
        if log_likelihood is not None:
            return log_likelihood
    return None


def _validate_iqtree_model_result(
    prefix_path: Path,
    *,
    workflow: str,
    default_selected_model: str | None = None,
) -> str:
    selected_model = _parse_best_model_artifact(prefix_path)
    if selected_model is None and default_selected_model is not None:
        return default_selected_model
    if selected_model is None:
        raise build_engine_output_error(
            f"iqtree {workflow} did not expose a parsable best-fit model result",
            code="engine_model_result_missing",
            engine_name="iqtree",
            workflow=workflow,
            path=prefix_path.with_suffix(".iqtree"),
            output_name="iqtree_report",
            artifact_kind="iqtree-model-result",
            details={
                "model_sidecar_path": (
                    None
                    if resolve_iqtree_model_sidecar(prefix_path) is None
                    else str(resolve_iqtree_model_sidecar(prefix_path))
                )
            },
        )
    return selected_model


def _existing_iqtree_outputs(
    prefix_path: Path,
    *,
    include_tree: bool = False,
    include_bootstrap: bool = False,
    include_consensus: bool = False,
) -> dict[str, Path]:
    outputs: dict[str, Path] = {}
    for key, candidate in (
        ("iqtree_report", prefix_path.with_suffix(".iqtree")),
        ("iqtree_log", prefix_path.with_suffix(".log")),
    ):
        if candidate.exists():
            outputs[key] = candidate
    model_sidecar = resolve_iqtree_model_sidecar(prefix_path)
    if model_sidecar is not None:
        outputs["model_selection_sidecar"] = model_sidecar
    tree_candidate = prefix_path.with_suffix(".treefile")
    if include_tree and tree_candidate.exists():
        outputs["tree"] = tree_candidate
    bootstrap_candidate = prefix_path.with_suffix(".ufboot")
    if include_bootstrap and bootstrap_candidate.exists():
        outputs["bootstrap_trees"] = bootstrap_candidate
    consensus_candidate = prefix_path.with_suffix(".contree")
    if include_consensus and consensus_candidate.exists():
        outputs["consensus_tree"] = consensus_candidate
    return outputs


def _build_iqtree_summary(
    prefix_path: Path,
    *,
    default_selected_model: str | None,
    support_tree_path: Path | None = None,
) -> IqtreeWorkflowSummary:
    selected_model = _parse_best_model_artifact(prefix_path) or default_selected_model
    log_likelihood = _parse_log_likelihood_artifact(prefix_path)
    support_values: list[IqtreeSupportValue] = []
    minimum_support: float | None = None
    maximum_support: float | None = None
    if support_tree_path is not None and support_tree_path.exists():
        support_summary = summarize_bootstrap_support_distribution(support_tree_path)
        support_values = [
            IqtreeSupportValue(
                node=node.node,
                descendant_taxa=list(node.descendant_taxa),
                support=node.support,
                support_fraction=node.support_fraction,
                is_backbone=node.is_backbone,
            )
            for node in support_summary.nodes
        ]
        minimum_support = support_summary.minimum_support
        maximum_support = support_summary.maximum_support
    return IqtreeWorkflowSummary(
        selected_model=selected_model,
        log_likelihood=log_likelihood,
        support_value_count=len(support_values),
        minimum_support=minimum_support,
        maximum_support=maximum_support,
        support_values=support_values,
    )


def _build_iqtree_model_selection_summary(
    prefix_path: Path,
) -> IqtreeModelSelectionSummary | None:
    return parse_iqtree_model_selection_summary(
        iqtree_report_path=prefix_path.with_suffix(".iqtree"),
        model_sidecar_path=resolve_iqtree_model_sidecar(prefix_path),
    )
