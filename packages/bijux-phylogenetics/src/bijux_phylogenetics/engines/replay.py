from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import shutil
from typing import Any

from bijux_phylogenetics.compare.topology import (
    compare_branch_lengths,
    compare_support_values,
    compare_tree_paths,
)
from bijux_phylogenetics.runtime.errors import EngineWorkflowError
from bijux_phylogenetics.io.fasta import load_fasta_alignment

from .common import build_file_checksums, load_engine_manifest, read_engine_version
from .fasta_to_tree import run_fasta_to_tree_workflow
from .inference_comparison import run_tree_inference_comparison
from .inference_reproducibility import run_inference_reproducibility_check
from .large_alignment_inference import run_large_alignment_inference
from .workflows.alignment import (
    run_alignment_trimming,
    run_codon_aware_multiple_sequence_alignment,
    run_multiple_sequence_alignment,
)
from .workflows.fasttree import run_fast_tree_inference
from .workflows.iqtree import (
    run_bootstrap_consensus_tree,
    run_bootstrap_support_estimation,
    run_maximum_likelihood_tree_inference,
    run_model_selection,
    run_sh_alrt_support_estimation,
)

__all__ = [
    "ManifestReplayComparison",
    "ManifestReplayDrift",
    "ManifestReplayReport",
    "replay_workflow_manifest",
]

_ENGINE_VERSION_ARGS: dict[str, tuple[str, ...]] = {
    "mafft": ("--version",),
    "trimal": ("--version",),
    "iqtree": ("--version",),
    "fasttree": ("-help",),
    "mrbayes": ("-v",),
    "beast": ("-version",),
}


@dataclass(slots=True)
class ManifestReplayDrift:
    kind: str
    label: str
    expected: str
    observed: str | None
    matched: bool


@dataclass(slots=True)
class ManifestReplayComparison:
    label: str
    status: str
    detail: str


@dataclass(slots=True)
class ManifestReplayReport:
    manifest_path: Path
    workflow: str
    replay_out_dir: Path
    replay_manifest_path: Path
    input_drift: list[ManifestReplayDrift]
    engine_version_drift: list[ManifestReplayDrift]
    comparisons: list[ManifestReplayComparison]
    input_drift_detected: bool
    engine_version_drift_detected: bool
    outputs_equivalent: bool
    notes: list[str]


def _payload_workflow(payload: dict[str, Any]) -> str:
    workflow = payload.get("workflow")
    if workflow is None:
        raise EngineWorkflowError(
            "manifest replay requires a workflow identifier",
            code="manifest_replay_missing_workflow",
        )
    return str(workflow)


def _path_map(values: dict[str, Any]) -> dict[str, Path]:
    return {str(key): Path(value) for key, value in values.items()}


def _recorded_input_paths(payload: dict[str, Any]) -> list[Path]:
    if "input_paths" in payload:
        return [Path(path) for path in payload["input_paths"]]
    return [Path(payload["input_path"])]


def _default_replay_out_dir(manifest_path: Path) -> Path:
    stem = manifest_path.name.removesuffix(".manifest.json")
    return manifest_path.parent / "replay" / stem


def _engine_key_from_name(engine_name: str) -> str:
    normalized = engine_name.strip().lower()
    if normalized in {"mafft", "trimal", "iqtree", "fasttree", "mrbayes", "beast"}:
        return normalized
    raise EngineWorkflowError(
        f"unsupported replay engine name: {engine_name}",
        code="manifest_replay_unknown_engine",
        details={"engine_name": engine_name},
    )


def _version_drift_for_engine_manifest(
    payload: dict[str, Any],
    *,
    executables: dict[str, str | Path | None],
) -> list[ManifestReplayDrift]:
    run_payload = dict(payload["run"])
    version_payload = dict(run_payload["version"])
    engine_key = _engine_key_from_name(str(payload["engine_name"]))
    executable = executables.get(engine_key) or str(run_payload["executable"])
    version = read_engine_version(
        str(payload["engine_name"]),
        executable,
        version_args=_ENGINE_VERSION_ARGS[engine_key],
        timeout_seconds=(
            None
            if run_payload.get("timeout_seconds") is None
            else float(run_payload["timeout_seconds"])
        ),
    )
    expected = str(version_payload["text"])
    observed = version.text
    return [
        ManifestReplayDrift(
            kind="engine-version",
            label=str(payload["engine_name"]),
            expected=expected,
            observed=observed,
            matched=expected == observed,
        )
    ]


def _version_drift_for_large_alignment(
    payload: dict[str, Any],
    *,
    executables: dict[str, str | Path | None],
) -> list[ManifestReplayDrift]:
    engine_key = _engine_key_from_name(str(payload["engine_name"]))
    command = [str(item) for item in payload["command"]]
    executable = executables.get(engine_key) or command[0]
    version = read_engine_version(
        str(payload["engine_name"]),
        executable,
        version_args=_ENGINE_VERSION_ARGS[engine_key],
    )
    expected = str(payload["engine_version_text"])
    observed = version.text
    return [
        ManifestReplayDrift(
            kind="engine-version",
            label=str(payload["engine_name"]),
            expected=expected,
            observed=observed,
            matched=expected == observed,
        )
    ]


def _collect_engine_version_drift(
    payload: dict[str, Any],
    *,
    executables: dict[str, str | Path | None],
) -> list[ManifestReplayDrift]:
    if "run" in payload:
        return _version_drift_for_engine_manifest(payload, executables=executables)
    if "step_manifests" in payload:
        drifts: list[ManifestReplayDrift] = []
        for step_manifest_path in _path_map(dict(payload["step_manifests"])).values():
            drifts.extend(
                _version_drift_for_engine_manifest(
                    load_engine_manifest(step_manifest_path),
                    executables=executables,
                )
            )
        return drifts
    if "engine_version_text" in payload and "command" in payload:
        return _version_drift_for_large_alignment(payload, executables=executables)
    raise EngineWorkflowError(
        "manifest replay could not determine engine versions from the manifest",
        code="manifest_replay_missing_engine_versions",
    )


def _collect_input_drift(payload: dict[str, Any]) -> list[ManifestReplayDrift]:
    recorded = {
        str(key): str(value)
        for key, value in dict(payload.get("input_checksums", {})).items()
    }
    observed = build_file_checksums([Path(path) for path in recorded])
    drifts: list[ManifestReplayDrift] = []
    for path_text, expected in recorded.items():
        current = observed.get(path_text)
        drifts.append(
            ManifestReplayDrift(
                kind="input-checksum",
                label=path_text,
                expected=expected,
                observed=current,
                matched=current == expected,
            )
        )
    return drifts


def _copy_input_for_inplace_engine(input_path: Path, replay_out_dir: Path) -> Path:
    replay_out_dir.mkdir(parents=True, exist_ok=True)
    replay_input = replay_out_dir / input_path.name
    shutil.copy2(input_path, replay_input)
    return replay_input


def _replay_prefix(manifest_path: Path, fallback: str) -> str:
    stem = manifest_path.name.removesuffix(".manifest.json")
    return fallback if fallback else stem


def _replay_engine_workflow(
    payload: dict[str, Any],
    *,
    replay_out_dir: Path,
    executables: dict[str, str | Path | None],
) -> Any:
    workflow = _payload_workflow(payload)
    engine_name = str(payload["engine_name"])
    engine_key = _engine_key_from_name(engine_name)
    executable = executables.get(engine_key) or str(dict(payload["run"])["executable"])
    config = dict(payload.get("config", {}))
    input_paths = _recorded_input_paths(payload)
    manifest_path = Path(payload["manifest_path"])
    output_paths = _path_map(dict(payload["output_paths"]))
    if workflow == "multiple-sequence-alignment":
        return run_multiple_sequence_alignment(
            input_paths[0],
            replay_out_dir / output_paths["alignment"].name,
            executable=executable,
            mode=str(config.get("mode", "auto")),
            extra_args=tuple(str(item) for item in config.get("extra_args", [])),
            timeout_seconds=(
                None
                if config.get("timeout_seconds") is None
                else float(config["timeout_seconds"])
            ),
        )
    if workflow == "codon-aware-multiple-sequence-alignment":
        return run_codon_aware_multiple_sequence_alignment(
            input_paths[0],
            replay_out_dir / output_paths["alignment"].name,
            executable=executable,
            mode=str(config.get("mode", "auto")),
            sequence_type=(
                None
                if config.get("sequence_type") is None
                else str(config["sequence_type"])
            ),
            genetic_code=config.get("genetic_code_id"),
            timeout_seconds=(
                None
                if config.get("timeout_seconds") is None
                else float(config["timeout_seconds"])
            ),
        )
    if workflow == "alignment-trimming":
        return run_alignment_trimming(
            input_paths[0],
            replay_out_dir / output_paths["trimmed_alignment"].name,
            executable=executable,
            mode=str(config.get("mode", "gap-threshold")),
            gap_threshold=float(config.get("gap_threshold", 0.1)),
            timeout_seconds=(
                None
                if config.get("timeout_seconds") is None
                else float(config["timeout_seconds"])
            ),
        )
    if workflow == "model-selection":
        return run_model_selection(
            input_paths[0],
            out_dir=replay_out_dir,
            prefix=_replay_prefix(manifest_path, "model-selection"),
            executable=executable,
            sequence_type=(
                None
                if config.get("sequence_type") is None
                else str(config["sequence_type"])
            ),
            partition_path=(
                None
                if config.get("partition_path") is None
                else Path(str(config["partition_path"]))
            ),
            seed=int(config.get("seed", 1)),
            threads=int(config.get("threads", 1)),
            timeout_seconds=(
                None
                if config.get("timeout_seconds") is None
                else float(config["timeout_seconds"])
            ),
        )
    if workflow == "maximum-likelihood-tree":
        return run_maximum_likelihood_tree_inference(
            input_paths[0],
            out_dir=replay_out_dir,
            prefix=_replay_prefix(manifest_path, "maximum-likelihood"),
            executable=executable,
            model=str(config["model"]),
            sequence_type=(
                None
                if config.get("sequence_type") is None
                else str(config["sequence_type"])
            ),
            partition_path=(
                None
                if config.get("partition_path") is None
                else Path(str(config["partition_path"]))
            ),
            seed=int(config.get("seed", 1)),
            threads=int(config.get("threads", 1)),
            timeout_seconds=(
                None
                if config.get("timeout_seconds") is None
                else float(config["timeout_seconds"])
            ),
        )
    if workflow == "bootstrap-support":
        return run_bootstrap_support_estimation(
            input_paths[0],
            out_dir=replay_out_dir,
            model=str(config["model"]),
            replicates=int(config.get("replicates", 1000)),
            prefix=_replay_prefix(manifest_path, "bootstrap-support"),
            executable=executable,
            sequence_type=(
                None
                if config.get("sequence_type") is None
                else str(config["sequence_type"])
            ),
            partition_path=(
                None
                if config.get("partition_path") is None
                else Path(str(config["partition_path"]))
            ),
            seed=int(config.get("seed", 1)),
            threads=int(config.get("threads", 1)),
            timeout_seconds=(
                None
                if config.get("timeout_seconds") is None
                else float(config["timeout_seconds"])
            ),
        )
    if workflow == "sh-alrt-support":
        return run_sh_alrt_support_estimation(
            input_paths[0],
            out_dir=replay_out_dir,
            model=str(config["model"]),
            sh_alrt_replicates=int(config.get("sh_alrt_replicates", 1000)),
            bootstrap_replicates=int(config.get("bootstrap_replicates", 1000)),
            prefix=_replay_prefix(manifest_path, "sh-alrt-support"),
            executable=executable,
            sequence_type=(
                None
                if config.get("sequence_type") is None
                else str(config["sequence_type"])
            ),
            partition_path=(
                None
                if config.get("partition_path") is None
                else Path(str(config["partition_path"]))
            ),
            seed=int(config.get("seed", 1)),
            threads=int(config.get("threads", 1)),
            timeout_seconds=(
                None
                if config.get("timeout_seconds") is None
                else float(config["timeout_seconds"])
            ),
        )
    if workflow == "bootstrap-consensus":
        return run_bootstrap_consensus_tree(
            input_paths[0],
            out_dir=replay_out_dir,
            prefix=_replay_prefix(manifest_path, "bootstrap-consensus"),
            executable=executable,
            minimum_support=float(config.get("minimum_support", 0.5)),
            timeout_seconds=(
                None
                if config.get("timeout_seconds") is None
                else float(config["timeout_seconds"])
            ),
        )
    if workflow == "fast-approximate-tree":
        return run_fast_tree_inference(
            input_paths[0],
            replay_out_dir / output_paths["tree"].name,
            executable=executable,
            sequence_type=(
                None
                if config.get("sequence_type") is None
                else str(config["sequence_type"])
            ),
            timeout_seconds=(
                None
                if config.get("timeout_seconds") is None
                else float(config["timeout_seconds"])
            ),
        )
    if workflow == "posterior-tree-inference" and engine_name == "MrBayes":
        from bijux_phylogenetics.bayesian.mrbayes import (
            run_mrbayes_posterior_inference,
        )

        replay_input = _copy_input_for_inplace_engine(input_paths[0], replay_out_dir)
        return run_mrbayes_posterior_inference(
            replay_input,
            executable=executable,
            timeout_seconds=(
                None
                if config.get("timeout_seconds") is None
                else float(config["timeout_seconds"])
            ),
        )
    if workflow == "posterior-tree-inference" and engine_name == "BEAST":
        from bijux_phylogenetics.bayesian.beast.execution import (
            run_beast_posterior_inference,
        )

        replay_input = _copy_input_for_inplace_engine(input_paths[0], replay_out_dir)
        return run_beast_posterior_inference(
            replay_input,
            executable=executable,
            overwrite=bool(config.get("overwrite", True)),
            threads=int(config.get("threads", 1)),
            seed=int(config.get("seed", 1)),
            timeout_seconds=(
                None
                if config.get("timeout_seconds") is None
                else float(config["timeout_seconds"])
            ),
        )
    raise EngineWorkflowError(
        f"manifest replay does not support workflow '{workflow}' for engine '{engine_name}'",
        code="manifest_replay_unsupported_workflow",
        details={"workflow": workflow, "engine_name": engine_name},
    )


def _replay_composite_workflow(
    payload: dict[str, Any],
    *,
    replay_out_dir: Path,
    executables: dict[str, str | Path | None],
) -> Any:
    workflow = _payload_workflow(payload)
    config = dict(payload.get("config", {}))
    input_path = Path(payload["input_path"])
    prefix = str(payload["prefix"])
    if workflow == "fasta-to-tree":
        return run_fasta_to_tree_workflow(
            input_path,
            out_dir=replay_out_dir,
            prefix=prefix,
            sequence_type=(
                None
                if payload.get("sequence_type") is None
                else str(payload["sequence_type"])
            ),
            mafft_executable=executables.get("mafft") or "mafft",
            alignment_mode=str(payload["alignment_mode"]),
            trimal_executable=executables.get("trimal") or "trimal",
            trimming_mode=str(payload["trimming_mode"]),
            iqtree_executable=executables.get("iqtree") or "iqtree2",
            iqtree_seed=int(payload["iqtree_seed"]),
            iqtree_threads=int(payload["iqtree_threads"]),
            trim_gap_threshold=float(payload["trim_gap_threshold"]),
            bootstrap_replicates=int(payload["bootstrap_replicates"]),
            resume=False,
            timeout_seconds=(
                None
                if config.get("timeout_seconds") is None
                else float(config["timeout_seconds"])
            ),
            incomplete_run_policy=str(config.get("incomplete_run_policy", "reject")),
        )
    if workflow == "large-alignment-inference":
        return run_large_alignment_inference(
            input_path,
            out_dir=replay_out_dir,
            prefix=prefix,
            sequence_type=(
                None
                if payload.get("sequence_type") is None
                else str(payload["sequence_type"])
            ),
            executable=executables.get("fasttree") or "FastTree",
            timeout_seconds=(
                None
                if payload.get("timeout_seconds") is None
                else float(payload["timeout_seconds"])
            ),
        )
    if workflow == "tree-inference-comparison":
        return run_tree_inference_comparison(
            input_path,
            out_dir=replay_out_dir,
            prefix=prefix,
            sequence_type=(
                None
                if payload.get("sequence_type") is None
                else str(payload["sequence_type"])
            ),
            iqtree_executable=executables.get("iqtree") or "iqtree2",
            fasttree_executable=executables.get("fasttree") or "FastTree",
            iqtree_seed=int(payload["iqtree_seed"]),
            iqtree_threads=int(payload["iqtree_threads"]),
            bootstrap_replicates=int(payload["bootstrap_replicates"]),
            timeout_seconds=(
                None
                if payload.get("timeout_seconds") is None
                else float(payload["timeout_seconds"])
            ),
            incomplete_run_policy=str(config.get("incomplete_run_policy", "reject")),
        )
    if workflow == "inference-reproducibility":
        return run_inference_reproducibility_check(
            input_path,
            out_dir=replay_out_dir,
            prefix=prefix,
            sequence_type=(
                None
                if payload.get("sequence_type") is None
                else str(payload["sequence_type"])
            ),
            executable=executables.get("iqtree") or "iqtree2",
            repeats=int(payload["repeat_count"]),
            bootstrap_replicates=int(payload["bootstrap_replicates"]),
            seed=int(payload["iqtree_seed"]),
            threads=int(payload["iqtree_threads"]),
        )
    raise EngineWorkflowError(
        f"manifest replay does not support composite workflow '{workflow}'",
        code="manifest_replay_unsupported_workflow",
        details={"workflow": workflow},
    )


def _compare_alignment_outputs(
    left_path: Path, right_path: Path
) -> ManifestReplayComparison:
    left_records = [
        (row.identifier, row.sequence) for row in load_fasta_alignment(left_path)
    ]
    right_records = [
        (row.identifier, row.sequence) for row in load_fasta_alignment(right_path)
    ]
    return ManifestReplayComparison(
        label=left_path.name,
        status="equivalent" if left_records == right_records else "different",
        detail=(
            "aligned FASTA records matched exactly"
            if left_records == right_records
            else "aligned FASTA records differed"
        ),
    )


def _compare_tree_outputs(
    left_path: Path, right_path: Path, *, label: str
) -> list[ManifestReplayComparison]:
    topology = compare_tree_paths(left_path, right_path)
    support = compare_support_values(left_path, right_path)
    branch_lengths = compare_branch_lengths(left_path, right_path)
    branch_score = branch_lengths.branch_score.branch_score_distance
    same_support = len(support.conflicting_clades) == 0
    same_branches = branch_score is None or branch_score <= 1e-6
    equivalent = topology.same_unrooted_topology and same_support
    return [
        ManifestReplayComparison(
            label=label,
            status="equivalent" if equivalent else "different",
            detail=(
                "topology and support were scientifically equivalent"
                if equivalent
                else "tree topology or support differed"
            ),
        ),
        ManifestReplayComparison(
            label=f"{label}-branch-lengths",
            status="exact" if same_branches else "different",
            detail=(
                "branch lengths matched within tolerance"
                if same_branches
                else "branch lengths differed beyond tolerance"
            ),
        ),
    ]


def _compare_posterior_outputs(
    payload: dict[str, Any],
    replay_report: Any,
) -> list[ManifestReplayComparison]:
    output_paths = _path_map(dict(payload["output_paths"]))
    if str(payload["engine_name"]) == "MrBayes":
        from bijux_phylogenetics.bayesian.mrbayes import (
            parse_mrbayes_consensus_tree,
            parse_mrbayes_parameter_traces,
            parse_mrbayes_posterior_tree_samples,
        )

        original_trace = parse_mrbayes_parameter_traces(
            output_paths["parameter_traces"]
        )
        replay_trace = parse_mrbayes_parameter_traces(
            replay_report.output_paths["parameter_traces"]
        )
        original_consensus = parse_mrbayes_consensus_tree(
            output_paths["consensus_tree"]
        )
        replay_consensus = parse_mrbayes_consensus_tree(
            replay_report.output_paths["consensus_tree"]
        )
        tree_comparisons = _compare_tree_outputs(
            output_paths["consensus_tree"],
            replay_report.output_paths["consensus_tree"],
            label="consensus-tree",
        )
        return [
            ManifestReplayComparison(
                label="parameter-traces",
                status=(
                    "equivalent"
                    if original_trace.columns == replay_trace.columns
                    and original_trace.row_count == replay_trace.row_count
                    else "different"
                ),
                detail="posterior parameter trace schema matched"
                if original_trace.columns == replay_trace.columns
                and original_trace.row_count == replay_trace.row_count
                else "posterior parameter trace schema differed",
            ),
            ManifestReplayComparison(
                label="posterior-trees",
                status=(
                    "equivalent"
                    if parse_mrbayes_posterior_tree_samples(
                        output_paths["posterior_trees"]
                    ).tree_count
                    == parse_mrbayes_posterior_tree_samples(
                        replay_report.output_paths["posterior_trees"]
                    ).tree_count
                    else "different"
                ),
                detail="posterior tree sample counts matched",
            ),
            ManifestReplayComparison(
                label="consensus-tip-set",
                status=(
                    "equivalent"
                    if original_consensus.tip_names == replay_consensus.tip_names
                    else "different"
                ),
                detail="consensus tree taxon set matched"
                if original_consensus.tip_names == replay_consensus.tip_names
                else "consensus tree taxon set differed",
            ),
            *tree_comparisons,
        ]
    from bijux_phylogenetics.bayesian.beast.logs import (
        parse_beast_log,
    )
    from bijux_phylogenetics.bayesian.beast.posterior_trees import (
        parse_beast_posterior_tree_samples,
    )

    original_log = parse_beast_log(output_paths["posterior_log"])
    replay_log = parse_beast_log(replay_report.output_paths["posterior_log"])
    original_trees = parse_beast_posterior_tree_samples(
        output_paths["posterior_trees"],
        burnin_fraction=0.0,
    )
    replay_trees = parse_beast_posterior_tree_samples(
        replay_report.output_paths["posterior_trees"],
        burnin_fraction=0.0,
    )
    return [
        ManifestReplayComparison(
            label="posterior-log",
            status=(
                "equivalent"
                if original_log.columns == replay_log.columns
                and original_log.row_count == replay_log.row_count
                else "different"
            ),
            detail="posterior log schema matched"
            if original_log.columns == replay_log.columns
            and original_log.row_count == replay_log.row_count
            else "posterior log schema differed",
        ),
        ManifestReplayComparison(
            label="posterior-trees",
            status=(
                "equivalent"
                if original_trees.kept_tree_count == replay_trees.kept_tree_count
                and original_trees.tip_names == replay_trees.tip_names
                else "different"
            ),
            detail="posterior tree sample counts and taxa matched"
            if original_trees.kept_tree_count == replay_trees.kept_tree_count
            and original_trees.tip_names == replay_trees.tip_names
            else "posterior tree sample counts or taxa differed",
        ),
    ]


def _compare_outputs(
    payload: dict[str, Any], replay_report: Any
) -> list[ManifestReplayComparison]:
    workflow = _payload_workflow(payload)
    output_paths = _path_map(dict(payload["output_paths"]))
    if workflow in {
        "multiple-sequence-alignment",
        "codon-aware-multiple-sequence-alignment",
        "alignment-trimming",
    }:
        key = "alignment" if "alignment" in output_paths else "trimmed_alignment"
        return [
            _compare_alignment_outputs(
                output_paths[key], replay_report.output_paths[key]
            )
        ]
    if workflow == "model-selection":
        original_model = str(payload.get("selected_model"))
        replay_model = getattr(replay_report, "selected_model", None)
        return [
            ManifestReplayComparison(
                label="selected-model",
                status="equivalent" if original_model == replay_model else "different",
                detail="selected substitution model matched"
                if original_model == replay_model
                else "selected substitution model differed",
            )
        ]
    if workflow == "maximum-likelihood-tree":
        return _compare_tree_outputs(
            output_paths["tree"],
            replay_report.output_paths["tree"],
            label="maximum-likelihood-tree",
        )
    if workflow == "bootstrap-support":
        return _compare_tree_outputs(
            output_paths["support_tree"],
            replay_report.output_paths["support_tree"],
            label="bootstrap-support-tree",
        )
    if workflow == "sh-alrt-support":
        return _compare_tree_outputs(
            output_paths["support_tree"],
            replay_report.output_paths["support_tree"],
            label="sh-alrt-support-tree",
        )
    if workflow == "bootstrap-consensus":
        return _compare_tree_outputs(
            output_paths["consensus_tree"],
            replay_report.output_paths["consensus_tree"],
            label="bootstrap-consensus-tree",
        )
    if workflow == "fast-approximate-tree":
        return _compare_tree_outputs(
            output_paths["tree"],
            replay_report.output_paths["tree"],
            label="fast-approximate-tree",
        )
    if workflow == "posterior-tree-inference":
        return _compare_posterior_outputs(payload, replay_report)
    if workflow == "fasta-to-tree":
        comparisons = [
            ManifestReplayComparison(
                label="selected-model",
                status=(
                    "equivalent"
                    if str(payload["selected_model"])
                    == getattr(replay_report, "selected_model", None)
                    else "different"
                ),
                detail="selected substitution model matched"
                if str(payload["selected_model"])
                == getattr(replay_report, "selected_model", None)
                else "selected substitution model differed",
            )
        ]
        comparisons.extend(
            _compare_tree_outputs(
                output_paths["tree"],
                replay_report.output_paths["tree"],
                label="fasta-to-tree",
            )
        )
        return comparisons
    if workflow == "large-alignment-inference":
        return _compare_tree_outputs(
            output_paths["tree"],
            replay_report.output_paths["tree"],
            label="large-alignment-inference",
        )
    if workflow == "tree-inference-comparison":
        return [
            *_compare_tree_outputs(
                output_paths["fasttree_tree"],
                replay_report.output_paths["fasttree_tree"],
                label="comparison-fasttree",
            ),
            *_compare_tree_outputs(
                output_paths["iqtree_support_tree"],
                replay_report.output_paths["iqtree_support_tree"],
                label="comparison-iqtree",
            ),
        ]
    if workflow == "inference-reproducibility":
        original_status = str(payload["overall_status"])
        replay_status = str(replay_report.overall_status)
        return [
            ManifestReplayComparison(
                label="overall-status",
                status="equivalent"
                if original_status == replay_status
                else "different",
                detail="reproducibility classification matched"
                if original_status == replay_status
                else "reproducibility classification differed",
            )
        ]
    raise EngineWorkflowError(
        f"manifest replay does not define output comparison for workflow '{workflow}'",
        code="manifest_replay_missing_comparison",
        details={"workflow": workflow},
    )


def replay_workflow_manifest(
    manifest_path: Path,
    *,
    out_dir: Path | None = None,
    executables: dict[str, str | Path | None] | None = None,
) -> ManifestReplayReport:
    payload = load_engine_manifest(manifest_path)
    workflow = _payload_workflow(payload)
    replay_out_dir = (
        _default_replay_out_dir(manifest_path) if out_dir is None else out_dir
    )
    replay_out_dir.mkdir(parents=True, exist_ok=True)
    executable_overrides = {} if executables is None else dict(executables)

    input_drift = _collect_input_drift(payload)
    input_drift_detected = any(not drift.matched for drift in input_drift)
    if input_drift_detected:
        raise EngineWorkflowError(
            "manifest replay detected changed inputs and refused to rerun",
            code="manifest_replay_input_changed",
            details={
                "manifest_path": str(manifest_path),
                "workflow": workflow,
                "changed_inputs": [
                    drift.label for drift in input_drift if not drift.matched
                ],
            },
        )

    engine_version_drift = _collect_engine_version_drift(
        payload,
        executables=executable_overrides,
    )
    if "run" in payload:
        replay_report = _replay_engine_workflow(
            payload,
            replay_out_dir=replay_out_dir,
            executables=executable_overrides,
        )
    else:
        replay_report = _replay_composite_workflow(
            payload,
            replay_out_dir=replay_out_dir,
            executables=executable_overrides,
        )
    comparisons = _compare_outputs(payload, replay_report)
    return ManifestReplayReport(
        manifest_path=manifest_path,
        workflow=workflow,
        replay_out_dir=replay_out_dir,
        replay_manifest_path=Path(replay_report.manifest_path),
        input_drift=input_drift,
        engine_version_drift=engine_version_drift,
        comparisons=comparisons,
        input_drift_detected=False,
        engine_version_drift_detected=any(
            not drift.matched for drift in engine_version_drift
        ),
        outputs_equivalent=all(
            comparison.status in {"exact", "equivalent"} for comparison in comparisons
        ),
        notes=[
            "replay reran the governed workflow from its recorded manifest inputs and configuration"
        ],
    )
