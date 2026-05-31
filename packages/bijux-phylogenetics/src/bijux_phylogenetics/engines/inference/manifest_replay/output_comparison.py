from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.compare.topology import (
    compare_branch_lengths,
    compare_support_values,
    compare_tree_paths,
)
from bijux_phylogenetics.io.fasta import load_fasta_alignment
from bijux_phylogenetics.runtime.errors import EngineWorkflowError

from .contracts import ManifestReplayComparison
from .manifest_policy import path_map, payload_workflow


def compare_alignment_outputs(
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


def compare_tree_outputs(
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


def compare_posterior_outputs(
    payload: dict[str, Any],
    replay_report: Any,
) -> list[ManifestReplayComparison]:
    output_paths = path_map(dict(payload["output_paths"]))
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
        tree_comparisons = compare_tree_outputs(
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


def compare_outputs(
    payload: dict[str, Any], replay_report: Any
) -> list[ManifestReplayComparison]:
    workflow = payload_workflow(payload)
    output_paths = path_map(dict(payload["output_paths"]))
    if workflow in {
        "multiple-sequence-alignment",
        "codon-aware-multiple-sequence-alignment",
        "alignment-trimming",
    }:
        key = "alignment" if "alignment" in output_paths else "trimmed_alignment"
        return [
            compare_alignment_outputs(
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
        return compare_tree_outputs(
            output_paths["tree"],
            replay_report.output_paths["tree"],
            label="maximum-likelihood-tree",
        )
    if workflow == "bootstrap-support":
        return compare_tree_outputs(
            output_paths["support_tree"],
            replay_report.output_paths["support_tree"],
            label="bootstrap-support-tree",
        )
    if workflow == "sh-alrt-support":
        return compare_tree_outputs(
            output_paths["support_tree"],
            replay_report.output_paths["support_tree"],
            label="sh-alrt-support-tree",
        )
    if workflow == "bootstrap-consensus":
        return compare_tree_outputs(
            output_paths["consensus_tree"],
            replay_report.output_paths["consensus_tree"],
            label="bootstrap-consensus-tree",
        )
    if workflow == "fast-approximate-tree":
        return compare_tree_outputs(
            output_paths["tree"],
            replay_report.output_paths["tree"],
            label="fast-approximate-tree",
        )
    if workflow == "posterior-tree-inference":
        return compare_posterior_outputs(payload, replay_report)
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
            compare_tree_outputs(
                output_paths["tree"],
                replay_report.output_paths["tree"],
                label="fasta-to-tree",
            )
        )
        return comparisons
    if workflow == "large-alignment-inference":
        return compare_tree_outputs(
            output_paths["tree"],
            replay_report.output_paths["tree"],
            label="large-alignment-inference",
        )
    if workflow == "tree-inference-comparison":
        return [
            *compare_tree_outputs(
                output_paths["fasttree_tree"],
                replay_report.output_paths["fasttree_tree"],
                label="comparison-fasttree",
            ),
            *compare_tree_outputs(
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
