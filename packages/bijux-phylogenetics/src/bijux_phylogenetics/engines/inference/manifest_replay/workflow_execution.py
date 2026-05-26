from __future__ import annotations

from pathlib import Path
import shutil
from typing import Any

from bijux_phylogenetics.runtime.errors import EngineWorkflowError

from ...workflows.alignment import (
    run_alignment_trimming,
    run_codon_aware_multiple_sequence_alignment,
    run_multiple_sequence_alignment,
)
from ...workflows.fasttree import run_fast_tree_inference
from ...workflows.iqtree import (
    run_bootstrap_consensus_tree,
    run_bootstrap_support_estimation,
    run_maximum_likelihood_tree_inference,
    run_model_selection,
    run_sh_alrt_support_estimation,
)
from ..comparison import run_tree_inference_comparison
from ..fasta_to_tree import run_fasta_to_tree_workflow
from ..large_alignment import run_large_alignment_inference
from ..reproducibility import run_inference_reproducibility_check
from .manifest_policy import (
    engine_key_from_name,
    path_map,
    payload_workflow,
    recorded_command_executable,
    recorded_composite_executable,
    recorded_input_paths,
)


def copy_input_for_inplace_engine(input_path: Path, replay_out_dir: Path) -> Path:
    replay_out_dir.mkdir(parents=True, exist_ok=True)
    replay_input = replay_out_dir / input_path.name
    shutil.copy2(input_path, replay_input)
    return replay_input


def replay_prefix(manifest_path: Path, fallback: str) -> str:
    stem = manifest_path.name.removesuffix(".manifest.json")
    return fallback if fallback else stem


def replay_engine_workflow(
    payload: dict[str, Any],
    *,
    replay_out_dir: Path,
    executables: dict[str, str | Path | None],
) -> Any:
    workflow = payload_workflow(payload)
    engine_name = str(payload["engine_name"])
    engine_key = engine_key_from_name(engine_name)
    executable = executables.get(engine_key) or str(dict(payload["run"])["executable"])
    config = dict(payload.get("config", {}))
    input_paths = recorded_input_paths(payload)
    manifest_path = Path(payload["manifest_path"])
    output_paths = path_map(dict(payload["output_paths"]))
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
            prefix=replay_prefix(manifest_path, "model-selection"),
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
            prefix=replay_prefix(manifest_path, "maximum-likelihood"),
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
            prefix=replay_prefix(manifest_path, "bootstrap-support"),
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
            prefix=replay_prefix(manifest_path, "sh-alrt-support"),
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
            prefix=replay_prefix(manifest_path, "bootstrap-consensus"),
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

        replay_input = copy_input_for_inplace_engine(input_paths[0], replay_out_dir)
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

        replay_input = copy_input_for_inplace_engine(input_paths[0], replay_out_dir)
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


def replay_composite_workflow(
    payload: dict[str, Any],
    *,
    replay_out_dir: Path,
    executables: dict[str, str | Path | None],
) -> Any:
    workflow = payload_workflow(payload)
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
            mafft_executable=executables.get("mafft")
            or recorded_composite_executable(payload, engine_key="mafft"),
            alignment_mode=str(payload["alignment_mode"]),
            trimal_executable=executables.get("trimal")
            or recorded_composite_executable(payload, engine_key="trimal"),
            trimming_mode=str(payload["trimming_mode"]),
            iqtree_executable=executables.get("iqtree")
            or recorded_composite_executable(payload, engine_key="iqtree"),
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
            executable=executables.get("fasttree")
            or recorded_command_executable(payload),
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
            iqtree_executable=executables.get("iqtree")
            or recorded_composite_executable(payload, engine_key="iqtree"),
            fasttree_executable=executables.get("fasttree")
            or recorded_composite_executable(payload, engine_key="fasttree"),
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
            executable=executables.get("iqtree")
            or recorded_composite_executable(payload, engine_key="iqtree"),
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
