from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import (
    _add_manifest_argument,
    _add_preflight_executable_arguments,
)
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.engines import (
    export_workflow_result_bundle,
    inspect_external_engine_preflight,
    list_external_engine_workflows,
    replay_workflow_manifest,
    require_preflight_workflow,
    run_phylo_workflow_config,
    validate_workflow_result_bundle,
)
from bijux_phylogenetics.parsimony import (
    bootstrap_parsimony,
    build_parsimony_stepwise_addition_tree,
    compute_parsimony_bremer_support,
    consistency_index,
    jackknife_parsimony,
    load_fitch_character_matrix,
    load_parsimony_character_matrix,
    load_parsimony_character_weights,
    load_sankoff_cost_matrix,
    place_parsimony_queries,
    reconstruct_acctran,
    reconstruct_deltran,
    rescaled_consistency_index,
    retention_index,
    run_parsimony_ratchet,
    score_camin_sokal,
    score_dollo,
    score_fitch,
    score_sankoff,
    score_wagner,
    search_parsimony_nni,
    search_parsimony_spr,
    summarize_equal_best_parsimony_trees,
    tree_length,
    write_camin_sokal_artifacts,
    write_dollo_artifacts,
    write_fitch_artifacts,
    write_parsimony_bootstrap_artifacts,
    write_parsimony_bremer_support_artifacts,
    write_parsimony_consistency_artifacts,
    write_parsimony_equal_best_consensus_artifacts,
    write_parsimony_jackknife_artifacts,
    write_parsimony_nni_artifacts,
    write_parsimony_placement_artifacts,
    write_parsimony_ratchet_artifacts,
    write_parsimony_reconstruction_artifacts,
    write_parsimony_rescaled_consistency_artifacts,
    write_parsimony_retention_artifacts,
    write_parsimony_spr_artifacts,
    write_parsimony_tree_length_artifacts,
    write_sankoff_artifacts,
    write_wagner_artifacts,
)
from bijux_phylogenetics.phylo.dating import (
    cross_validate_penalized_likelihood_smoothing_from_metadata,
    fit_least_squares_dating_from_metadata,
    fit_penalized_likelihood_dating_from_metadata,
    solve_dating_calibration_constraints,
    summarize_relaxed_rate_branches_from_paths,
    write_dating_calibration_constraint_artifacts,
    write_least_squares_dating_artifacts,
    write_penalized_likelihood_cross_validation_artifacts,
    write_penalized_likelihood_dating_artifacts,
    write_relaxed_rate_branch_summary_artifacts,
)
from bijux_phylogenetics.phylo.likelihood import (
    bootstrap_nucleotide_likelihood_tree_inference_from_alignment,
    fit_local_clock_likelihood_from_alignment,
    fit_strict_clock_likelihood_from_alignment,
    infer_nucleotide_likelihood_tree_from_alignment,
    place_queries_by_likelihood_from_alignment,
    write_likelihood_placement_artifacts,
    write_local_clock_likelihood_artifacts,
    write_nucleotide_likelihood_bootstrap_artifacts,
    write_nucleotide_likelihood_tree_inference_artifacts,
    write_strict_clock_likelihood_artifacts,
)
from bijux_phylogenetics.phylo.topology import write_stepwise_addition_artifacts
from bijux_phylogenetics.runtime.errors import EngineWorkflowError
from bijux_phylogenetics.runtime.results import build_command_result


def add_phylo_commands(subparsers: Any) -> None:
    phylo = subparsers.add_parser(
        get_command_spec("phylo").name,
        help=get_command_spec("phylo").summary,
    )
    phylo_subparsers = phylo.add_subparsers(dest="phylo_command", required=True)

    phylo_preflight = phylo_subparsers.add_parser(
        "preflight",
        help="Inspect external engine availability, version support, and workflow readiness.",
    )
    phylo_preflight.add_argument(
        "--workflow",
        choices=list_external_engine_workflows(),
        help="Require one selected external-engine workflow to be runnable in the current environment.",
    )
    _add_preflight_executable_arguments(phylo_preflight)
    phylo_preflight.add_argument(
        "--json", action="store_true", help="Emit the preflight report as JSON."
    )
    _add_manifest_argument(phylo_preflight)

    phylo_run = phylo_subparsers.add_parser(
        "run",
        help="Run one governed workflow from one YAML or JSON config file and export a validated result bundle.",
    )
    phylo_run.add_argument("config_path", type=Path)
    phylo_run.add_argument(
        "--json", action="store_true", help="Emit the config-run report as JSON."
    )
    _add_manifest_argument(phylo_run)

    phylo_replay = phylo_subparsers.add_parser(
        "replay",
        help="Rerun one governed phylogenetics workflow from its manifest and compare the replayed outputs.",
    )
    phylo_replay.add_argument("manifest_path", type=Path)
    phylo_replay.add_argument("--out-dir", type=Path)
    _add_preflight_executable_arguments(phylo_replay)
    phylo_replay.add_argument(
        "--json", action="store_true", help="Emit the replay report as JSON."
    )
    _add_manifest_argument(phylo_replay)

    phylo_bundle = phylo_subparsers.add_parser(
        "bundle",
        help="Export one portable workflow-result bundle from a governed workflow manifest.",
    )
    phylo_bundle.add_argument("manifest_path", type=Path)
    phylo_bundle.add_argument("--out-dir", required=True, type=Path)
    phylo_bundle.add_argument(
        "--json", action="store_true", help="Emit the bundle report as JSON."
    )
    _add_manifest_argument(phylo_bundle)

    phylo_validate_bundle = phylo_subparsers.add_parser(
        "validate-bundle",
        help="Validate one workflow-result bundle for checksum integrity and required workflow contents.",
    )
    phylo_validate_bundle.add_argument("bundle_root", type=Path)
    phylo_validate_bundle.add_argument(
        "--json", action="store_true", help="Emit the validation report as JSON."
    )
    _add_manifest_argument(phylo_validate_bundle)

    phylo_parsimony = phylo_subparsers.add_parser(
        "parsimony",
        help="Score governed character matrices on a tree with explicit parsimony methods.",
    )
    phylo_parsimony_subparsers = phylo_parsimony.add_subparsers(
        dest="phylo_parsimony_command",
        required=True,
    )
    phylo_parsimony_fitch = phylo_parsimony_subparsers.add_parser(
        "fitch",
        help="Score one unordered discrete character matrix on one tree with Fitch parsimony.",
    )
    phylo_parsimony_fitch.add_argument("tree_path", type=Path)
    phylo_parsimony_fitch.add_argument("matrix_path", type=Path)
    phylo_parsimony_fitch.add_argument("--taxon-column")
    _add_parsimony_character_weights_argument(phylo_parsimony_fitch)
    phylo_parsimony_fitch.add_argument("--out-dir", required=True, type=Path)
    phylo_parsimony_fitch.add_argument(
        "--json", action="store_true", help="Emit the parsimony report as JSON."
    )
    _add_manifest_argument(phylo_parsimony_fitch)
    phylo_parsimony_placement = phylo_parsimony_subparsers.add_parser(
        "placement",
        help="Place one or more query character vectors on every reference-tree edge by additional Fitch parsimony steps.",
    )
    phylo_parsimony_placement.add_argument("tree_path", type=Path)
    phylo_parsimony_placement.add_argument("matrix_path", type=Path)
    phylo_parsimony_placement.add_argument("query_matrix_path", type=Path)
    phylo_parsimony_placement.add_argument(
        "--method",
        required=True,
        choices=["fitch"],
        help="Parsimony scoring method for query placement.",
    )
    phylo_parsimony_placement.add_argument("--taxon-column")
    _add_parsimony_character_weights_argument(phylo_parsimony_placement)
    phylo_parsimony_placement.add_argument("--out-dir", required=True, type=Path)
    phylo_parsimony_placement.add_argument(
        "--json", action="store_true", help="Emit the placement report as JSON."
    )
    _add_manifest_argument(phylo_parsimony_placement)
    phylo_parsimony_wagner = phylo_parsimony_subparsers.add_parser(
        "wagner",
        help="Score one ordered discrete character matrix on one tree with Wagner parsimony.",
    )
    phylo_parsimony_wagner.add_argument("tree_path", type=Path)
    phylo_parsimony_wagner.add_argument("matrix_path", type=Path)
    phylo_parsimony_wagner.add_argument("--taxon-column")
    phylo_parsimony_wagner.add_argument(
        "--state-order",
        help="Comma-separated explicit ordered state labels such as low,medium,high.",
    )
    _add_parsimony_character_weights_argument(phylo_parsimony_wagner)
    phylo_parsimony_wagner.add_argument("--out-dir", required=True, type=Path)
    phylo_parsimony_wagner.add_argument(
        "--json", action="store_true", help="Emit the parsimony report as JSON."
    )
    _add_manifest_argument(phylo_parsimony_wagner)
    phylo_parsimony_sankoff = phylo_parsimony_subparsers.add_parser(
        "sankoff",
        help="Score one discrete character matrix on one tree with a user-supplied Sankoff cost matrix.",
    )
    phylo_parsimony_sankoff.add_argument("tree_path", type=Path)
    phylo_parsimony_sankoff.add_argument("matrix_path", type=Path)
    phylo_parsimony_sankoff.add_argument("cost_matrix_path", type=Path)
    phylo_parsimony_sankoff.add_argument("--taxon-column")
    phylo_parsimony_sankoff.add_argument(
        "--allow-asymmetric-costs",
        action="store_true",
        help="Allow asymmetric Sankoff transition costs instead of requiring symmetry.",
    )
    _add_parsimony_character_weights_argument(phylo_parsimony_sankoff)
    phylo_parsimony_sankoff.add_argument("--out-dir", required=True, type=Path)
    phylo_parsimony_sankoff.add_argument(
        "--json", action="store_true", help="Emit the parsimony report as JSON."
    )
    _add_manifest_argument(phylo_parsimony_sankoff)
    phylo_parsimony_dollo = phylo_parsimony_subparsers.add_parser(
        "dollo",
        help="Score one binary character matrix on one tree with Dollo parsimony.",
    )
    phylo_parsimony_dollo.add_argument("tree_path", type=Path)
    phylo_parsimony_dollo.add_argument("matrix_path", type=Path)
    phylo_parsimony_dollo.add_argument("--taxon-column")
    _add_parsimony_character_weights_argument(phylo_parsimony_dollo)
    phylo_parsimony_dollo.add_argument("--out-dir", required=True, type=Path)
    phylo_parsimony_dollo.add_argument(
        "--json", action="store_true", help="Emit the parsimony report as JSON."
    )
    _add_manifest_argument(phylo_parsimony_dollo)
    phylo_parsimony_camin_sokal = phylo_parsimony_subparsers.add_parser(
        "camin-sokal",
        help="Score one binary character matrix on one tree with irreversible Camin-Sokal parsimony.",
    )
    phylo_parsimony_camin_sokal.add_argument("tree_path", type=Path)
    phylo_parsimony_camin_sokal.add_argument("matrix_path", type=Path)
    phylo_parsimony_camin_sokal.add_argument("--taxon-column")
    _add_parsimony_character_weights_argument(phylo_parsimony_camin_sokal)
    phylo_parsimony_camin_sokal.add_argument("--out-dir", required=True, type=Path)
    phylo_parsimony_camin_sokal.add_argument(
        "--json", action="store_true", help="Emit the parsimony report as JSON."
    )
    _add_manifest_argument(phylo_parsimony_camin_sokal)
    phylo_parsimony_acctran = phylo_parsimony_subparsers.add_parser(
        "acctran",
        help="Resolve one unordered parsimony reconstruction toward earlier branches with ACCTRAN.",
    )
    phylo_parsimony_acctran.add_argument("tree_path", type=Path)
    phylo_parsimony_acctran.add_argument("matrix_path", type=Path)
    phylo_parsimony_acctran.add_argument("--taxon-column")
    _add_parsimony_character_weights_argument(phylo_parsimony_acctran)
    phylo_parsimony_acctran.add_argument("--out-dir", required=True, type=Path)
    phylo_parsimony_acctran.add_argument(
        "--json", action="store_true", help="Emit the parsimony reconstruction as JSON."
    )
    _add_manifest_argument(phylo_parsimony_acctran)
    phylo_parsimony_deltran = phylo_parsimony_subparsers.add_parser(
        "deltran",
        help="Resolve one unordered parsimony reconstruction toward later branches with DELTRAN.",
    )
    phylo_parsimony_deltran.add_argument("tree_path", type=Path)
    phylo_parsimony_deltran.add_argument("matrix_path", type=Path)
    phylo_parsimony_deltran.add_argument("--taxon-column")
    _add_parsimony_character_weights_argument(phylo_parsimony_deltran)
    phylo_parsimony_deltran.add_argument("--out-dir", required=True, type=Path)
    phylo_parsimony_deltran.add_argument(
        "--json", action="store_true", help="Emit the parsimony reconstruction as JSON."
    )
    _add_manifest_argument(phylo_parsimony_deltran)
    phylo_parsimony_bootstrap = phylo_parsimony_subparsers.add_parser(
        "bootstrap",
        help="Resample one character matrix with replacement, infer exact small-taxon replicate trees, and map clade support by clade identity.",
    )
    phylo_parsimony_bootstrap.add_argument("matrix_path", type=Path)
    phylo_parsimony_bootstrap.add_argument(
        "--method",
        required=True,
        choices=[
            "fitch",
            "wagner",
            "sankoff",
            "dollo",
            "camin-sokal",
            "acctran",
            "deltran",
        ],
    )
    phylo_parsimony_bootstrap.add_argument("--taxon-column")
    phylo_parsimony_bootstrap.add_argument(
        "--replicate-count",
        required=True,
        type=int,
        help="Number of bootstrap replicate trees to infer.",
    )
    phylo_parsimony_bootstrap.add_argument(
        "--seed",
        required=True,
        type=int,
        help="Deterministic random seed for character resampling.",
    )
    phylo_parsimony_bootstrap.add_argument(
        "--cost-matrix",
        dest="cost_matrix_path",
        type=Path,
        help="Required when --method sankoff is selected.",
    )
    phylo_parsimony_bootstrap.add_argument(
        "--allow-asymmetric-costs",
        action="store_true",
        help="Allow asymmetric Sankoff transition costs when --method sankoff is selected.",
    )
    phylo_parsimony_bootstrap.add_argument(
        "--state-order",
        help="Comma-separated explicit ordered state labels for Wagner bootstrap scoring.",
    )
    _add_parsimony_character_weights_argument(phylo_parsimony_bootstrap)
    phylo_parsimony_bootstrap.add_argument("--out-dir", required=True, type=Path)
    phylo_parsimony_bootstrap.add_argument(
        "--json", action="store_true", help="Emit the bootstrap report as JSON."
    )
    _add_manifest_argument(phylo_parsimony_bootstrap)
    phylo_parsimony_jackknife = phylo_parsimony_subparsers.add_parser(
        "jackknife",
        help="Subsample one character matrix without replacement, infer exact small-taxon replicate trees, and map clade support by clade identity.",
    )
    phylo_parsimony_jackknife.add_argument("matrix_path", type=Path)
    phylo_parsimony_jackknife.add_argument(
        "--method",
        required=True,
        choices=[
            "fitch",
            "wagner",
            "sankoff",
            "dollo",
            "camin-sokal",
            "acctran",
            "deltran",
        ],
    )
    phylo_parsimony_jackknife.add_argument("--taxon-column")
    phylo_parsimony_jackknife.add_argument(
        "--replicate-count",
        required=True,
        type=int,
        help="Number of jackknife replicate trees to infer.",
    )
    phylo_parsimony_jackknife.add_argument(
        "--seed",
        required=True,
        type=int,
        help="Deterministic random seed for without-replacement character subsampling.",
    )
    phylo_parsimony_jackknife.add_argument(
        "--retain-probability",
        type=float,
        default=0.75,
        help="Per-character retention probability for jackknife subsampling without replacement.",
    )
    phylo_parsimony_jackknife.add_argument(
        "--cost-matrix",
        dest="cost_matrix_path",
        type=Path,
        help="Required when --method sankoff is selected.",
    )
    phylo_parsimony_jackknife.add_argument(
        "--allow-asymmetric-costs",
        action="store_true",
        help="Allow asymmetric Sankoff transition costs when --method sankoff is selected.",
    )
    phylo_parsimony_jackknife.add_argument(
        "--state-order",
        help="Comma-separated explicit ordered state labels for Wagner jackknife scoring.",
    )
    _add_parsimony_character_weights_argument(phylo_parsimony_jackknife)
    phylo_parsimony_jackknife.add_argument("--out-dir", required=True, type=Path)
    phylo_parsimony_jackknife.add_argument(
        "--json", action="store_true", help="Emit the jackknife report as JSON."
    )
    _add_manifest_argument(phylo_parsimony_jackknife)
    phylo_parsimony_bremer_support = phylo_parsimony_subparsers.add_parser(
        "bremer-support",
        help="Compute exact small-taxon Bremer support for every informative clade on one rooted reference tree.",
    )
    phylo_parsimony_bremer_support.add_argument("tree_path", type=Path)
    phylo_parsimony_bremer_support.add_argument("matrix_path", type=Path)
    phylo_parsimony_bremer_support.add_argument(
        "--method",
        required=True,
        choices=[
            "fitch",
            "wagner",
            "sankoff",
            "dollo",
            "camin-sokal",
            "acctran",
            "deltran",
        ],
    )
    phylo_parsimony_bremer_support.add_argument("--taxon-column")
    phylo_parsimony_bremer_support.add_argument(
        "--cost-matrix",
        dest="cost_matrix_path",
        type=Path,
        help="Required when --method sankoff is selected.",
    )
    phylo_parsimony_bremer_support.add_argument(
        "--allow-asymmetric-costs",
        action="store_true",
        help="Allow asymmetric Sankoff transition costs when --method sankoff is selected.",
    )
    phylo_parsimony_bremer_support.add_argument(
        "--state-order",
        help="Comma-separated explicit ordered state labels for Wagner Bremer scoring.",
    )
    _add_parsimony_character_weights_argument(phylo_parsimony_bremer_support)
    phylo_parsimony_bremer_support.add_argument(
        "--out-dir",
        required=True,
        type=Path,
    )
    phylo_parsimony_bremer_support.add_argument(
        "--json",
        action="store_true",
        help="Emit the Bremer support report as JSON.",
    )
    _add_manifest_argument(phylo_parsimony_bremer_support)
    phylo_parsimony_equal_best_consensus = phylo_parsimony_subparsers.add_parser(
        "equal-best-consensus",
        help="Enumerate the exact equal-best parsimony tree set for a small matrix and summarize strict and majority consensus when the full set is retained.",
    )
    phylo_parsimony_equal_best_consensus.add_argument("matrix_path", type=Path)
    phylo_parsimony_equal_best_consensus.add_argument(
        "--method",
        required=True,
        choices=[
            "fitch",
            "wagner",
            "sankoff",
            "dollo",
            "camin-sokal",
            "acctran",
            "deltran",
        ],
    )
    phylo_parsimony_equal_best_consensus.add_argument("--taxon-column")
    phylo_parsimony_equal_best_consensus.add_argument(
        "--max-retained-equal-best-trees",
        type=int,
        default=128,
        help="Maximum number of equally optimal trees to retain; consensus summaries are omitted if the exact-best set exceeds this cap.",
    )
    phylo_parsimony_equal_best_consensus.add_argument(
        "--cost-matrix",
        dest="cost_matrix_path",
        type=Path,
        help="Required when --method sankoff is selected.",
    )
    phylo_parsimony_equal_best_consensus.add_argument(
        "--allow-asymmetric-costs",
        action="store_true",
        help="Allow asymmetric Sankoff transition costs when --method sankoff is selected.",
    )
    phylo_parsimony_equal_best_consensus.add_argument(
        "--state-order",
        help="Comma-separated explicit ordered state labels for Wagner exact consensus scoring.",
    )
    _add_parsimony_character_weights_argument(phylo_parsimony_equal_best_consensus)
    phylo_parsimony_equal_best_consensus.add_argument(
        "--out-dir",
        required=True,
        type=Path,
    )
    phylo_parsimony_equal_best_consensus.add_argument(
        "--json",
        action="store_true",
        help="Emit the equal-best consensus report as JSON.",
    )
    _add_manifest_argument(phylo_parsimony_equal_best_consensus)
    phylo_parsimony_stepwise_addition = phylo_parsimony_subparsers.add_parser(
        "stepwise-addition",
        help="Build one rooted tree by greedy taxon insertion under governed parsimony tree-length scoring.",
    )
    phylo_parsimony_stepwise_addition.add_argument("matrix_path", type=Path)
    phylo_parsimony_stepwise_addition.add_argument(
        "--method",
        required=True,
        choices=[
            "fitch",
            "wagner",
            "sankoff",
            "dollo",
            "camin-sokal",
            "acctran",
            "deltran",
        ],
    )
    phylo_parsimony_stepwise_addition.add_argument("--taxon-column")
    phylo_parsimony_stepwise_addition.add_argument(
        "--insertion-order",
        help="Optional comma-separated taxon order for greedy insertion; defaults to matrix row order.",
    )
    phylo_parsimony_stepwise_addition.add_argument(
        "--cost-matrix",
        dest="cost_matrix_path",
        type=Path,
        help="Required when --method sankoff is selected.",
    )
    phylo_parsimony_stepwise_addition.add_argument(
        "--allow-asymmetric-costs",
        action="store_true",
        help="Allow asymmetric Sankoff transition costs when --method sankoff is selected.",
    )
    phylo_parsimony_stepwise_addition.add_argument(
        "--state-order",
        help="Comma-separated explicit ordered state labels for Wagner stepwise scoring.",
    )
    _add_parsimony_character_weights_argument(phylo_parsimony_stepwise_addition)
    phylo_parsimony_stepwise_addition.add_argument("--out-dir", required=True, type=Path)
    phylo_parsimony_stepwise_addition.add_argument(
        "--json", action="store_true", help="Emit the stepwise-addition report as JSON."
    )
    _add_manifest_argument(phylo_parsimony_stepwise_addition)
    phylo_parsimony_nni_search = phylo_parsimony_subparsers.add_parser(
        "nni-search",
        help="Hill-climb one rooted binary starting tree by accepting score-improving rooted NNI moves.",
    )
    phylo_parsimony_nni_search.add_argument("tree_path", type=Path)
    phylo_parsimony_nni_search.add_argument("matrix_path", type=Path)
    phylo_parsimony_nni_search.add_argument(
        "--method",
        required=True,
        choices=[
            "fitch",
            "wagner",
            "sankoff",
            "dollo",
            "camin-sokal",
            "acctran",
            "deltran",
        ],
    )
    phylo_parsimony_nni_search.add_argument("--taxon-column")
    phylo_parsimony_nni_search.add_argument(
        "--cost-matrix",
        dest="cost_matrix_path",
        type=Path,
        help="Required when --method sankoff is selected.",
    )
    phylo_parsimony_nni_search.add_argument(
        "--allow-asymmetric-costs",
        action="store_true",
        help="Allow asymmetric Sankoff transition costs when --method sankoff is selected.",
    )
    phylo_parsimony_nni_search.add_argument(
        "--state-order",
        help="Comma-separated explicit ordered state labels for Wagner NNI scoring.",
    )
    _add_parsimony_character_weights_argument(phylo_parsimony_nni_search)
    phylo_parsimony_nni_search.add_argument("--out-dir", required=True, type=Path)
    phylo_parsimony_nni_search.add_argument(
        "--json", action="store_true", help="Emit the NNI search report as JSON."
    )
    _add_manifest_argument(phylo_parsimony_nni_search)
    phylo_parsimony_spr_search = phylo_parsimony_subparsers.add_parser(
        "spr-search",
        help="Hill-climb one rooted binary starting tree by accepting score-improving rooted subtree-prune-regraft moves.",
    )
    phylo_parsimony_spr_search.add_argument("tree_path", type=Path)
    phylo_parsimony_spr_search.add_argument("matrix_path", type=Path)
    phylo_parsimony_spr_search.add_argument(
        "--method",
        required=True,
        choices=[
            "fitch",
            "wagner",
            "sankoff",
            "dollo",
            "camin-sokal",
            "acctran",
            "deltran",
        ],
    )
    phylo_parsimony_spr_search.add_argument("--taxon-column")
    phylo_parsimony_spr_search.add_argument(
        "--cost-matrix",
        dest="cost_matrix_path",
        type=Path,
        help="Required when --method sankoff is selected.",
    )
    phylo_parsimony_spr_search.add_argument(
        "--allow-asymmetric-costs",
        action="store_true",
        help="Allow asymmetric Sankoff transition costs when --method sankoff is selected.",
    )
    phylo_parsimony_spr_search.add_argument(
        "--state-order",
        help="Comma-separated explicit ordered state labels for Wagner SPR scoring.",
    )
    _add_parsimony_character_weights_argument(phylo_parsimony_spr_search)
    phylo_parsimony_spr_search.add_argument("--out-dir", required=True, type=Path)
    phylo_parsimony_spr_search.add_argument(
        "--json", action="store_true", help="Emit the SPR search report as JSON."
    )
    _add_manifest_argument(phylo_parsimony_spr_search)
    phylo_parsimony_ratchet = phylo_parsimony_subparsers.add_parser(
        "ratchet",
        help="Run one deterministic parsimony ratchet with temporary character reweighting and restored rooted SPR search.",
    )
    phylo_parsimony_ratchet.add_argument("tree_path", type=Path)
    phylo_parsimony_ratchet.add_argument("matrix_path", type=Path)
    phylo_parsimony_ratchet.add_argument(
        "--method",
        required=True,
        choices=[
            "fitch",
            "wagner",
            "sankoff",
            "dollo",
            "camin-sokal",
            "acctran",
            "deltran",
        ],
    )
    phylo_parsimony_ratchet.add_argument("--taxon-column")
    phylo_parsimony_ratchet.add_argument(
        "--cycle-count",
        required=True,
        type=int,
        help="Number of ratchet perturbation cycles to run.",
    )
    phylo_parsimony_ratchet.add_argument(
        "--seed",
        required=True,
        type=int,
        help="Deterministic random seed for ratchet character perturbations.",
    )
    phylo_parsimony_ratchet.add_argument(
        "--perturbed-character-count",
        required=True,
        type=int,
        help="Number of characters to upweight temporarily per cycle.",
    )
    phylo_parsimony_ratchet.add_argument(
        "--perturbation-factor",
        type=float,
        default=2.0,
        help="Temporary multiplicative weight factor applied to the perturbed characters.",
    )
    phylo_parsimony_ratchet.add_argument(
        "--cost-matrix",
        dest="cost_matrix_path",
        type=Path,
        help="Required when --method sankoff is selected.",
    )
    phylo_parsimony_ratchet.add_argument(
        "--allow-asymmetric-costs",
        action="store_true",
        help="Allow asymmetric Sankoff transition costs when --method sankoff is selected.",
    )
    phylo_parsimony_ratchet.add_argument(
        "--state-order",
        help="Comma-separated explicit ordered state labels for Wagner ratchet scoring.",
    )
    _add_parsimony_character_weights_argument(phylo_parsimony_ratchet)
    phylo_parsimony_ratchet.add_argument("--out-dir", required=True, type=Path)
    phylo_parsimony_ratchet.add_argument(
        "--json", action="store_true", help="Emit the ratchet report as JSON."
    )
    _add_manifest_argument(phylo_parsimony_ratchet)
    phylo_parsimony_tree_length = phylo_parsimony_subparsers.add_parser(
        "tree-length",
        help="Summarize per-character and total tree length for one parsimony scoring method.",
    )
    phylo_parsimony_tree_length.add_argument("tree_path", type=Path)
    phylo_parsimony_tree_length.add_argument("matrix_path", type=Path)
    phylo_parsimony_tree_length.add_argument(
        "--method",
        required=True,
        choices=[
            "fitch",
            "wagner",
            "sankoff",
            "dollo",
            "camin-sokal",
            "acctran",
            "deltran",
        ],
    )
    phylo_parsimony_tree_length.add_argument("--taxon-column")
    phylo_parsimony_tree_length.add_argument(
        "--cost-matrix",
        dest="cost_matrix_path",
        type=Path,
        help="Required when --method sankoff is selected.",
    )
    phylo_parsimony_tree_length.add_argument(
        "--allow-asymmetric-costs",
        action="store_true",
        help="Allow asymmetric Sankoff transition costs when --method sankoff is selected.",
    )
    phylo_parsimony_tree_length.add_argument(
        "--state-order",
        help="Comma-separated explicit ordered state labels for Wagner tree length.",
    )
    phylo_parsimony_tree_length.add_argument(
        "--character-weights",
        dest="character_weights_path",
        type=Path,
        help="Optional TSV with character_id and weight columns.",
    )
    phylo_parsimony_tree_length.add_argument("--out-dir", required=True, type=Path)
    phylo_parsimony_tree_length.add_argument(
        "--json", action="store_true", help="Emit the tree-length report as JSON."
    )
    _add_manifest_argument(phylo_parsimony_tree_length)
    phylo_parsimony_consistency = phylo_parsimony_subparsers.add_parser(
        "consistency-index",
        help="Compute per-character and aggregate consistency index for one supported parsimony method.",
    )
    phylo_parsimony_consistency.add_argument("tree_path", type=Path)
    phylo_parsimony_consistency.add_argument("matrix_path", type=Path)
    phylo_parsimony_consistency.add_argument(
        "--method",
        required=True,
        choices=[
            "fitch",
            "wagner",
            "dollo",
            "camin-sokal",
            "acctran",
            "deltran",
        ],
    )
    phylo_parsimony_consistency.add_argument("--taxon-column")
    phylo_parsimony_consistency.add_argument(
        "--state-order",
        help="Comma-separated explicit ordered state labels for Wagner consistency index.",
    )
    phylo_parsimony_consistency.add_argument("--out-dir", required=True, type=Path)
    phylo_parsimony_consistency.add_argument(
        "--json", action="store_true", help="Emit the consistency-index report as JSON."
    )
    _add_manifest_argument(phylo_parsimony_consistency)
    phylo_parsimony_retention = phylo_parsimony_subparsers.add_parser(
        "retention-index",
        help="Compute per-character and aggregate retention index for unordered Fitch-style methods.",
    )
    phylo_parsimony_retention.add_argument("tree_path", type=Path)
    phylo_parsimony_retention.add_argument("matrix_path", type=Path)
    phylo_parsimony_retention.add_argument(
        "--method",
        required=True,
        choices=["fitch", "acctran", "deltran"],
    )
    phylo_parsimony_retention.add_argument("--taxon-column")
    phylo_parsimony_retention.add_argument("--out-dir", required=True, type=Path)
    phylo_parsimony_retention.add_argument(
        "--json", action="store_true", help="Emit the retention-index report as JSON."
    )
    _add_manifest_argument(phylo_parsimony_retention)
    phylo_parsimony_rescaled_consistency = phylo_parsimony_subparsers.add_parser(
        "rescaled-consistency-index",
        help="Compute per-character and aggregate rescaled consistency index from tested CI and RI surfaces.",
    )
    phylo_parsimony_rescaled_consistency.add_argument("tree_path", type=Path)
    phylo_parsimony_rescaled_consistency.add_argument("matrix_path", type=Path)
    phylo_parsimony_rescaled_consistency.add_argument(
        "--method",
        required=True,
        choices=["fitch", "acctran", "deltran"],
    )
    phylo_parsimony_rescaled_consistency.add_argument("--taxon-column")
    phylo_parsimony_rescaled_consistency.add_argument("--out-dir", required=True, type=Path)
    phylo_parsimony_rescaled_consistency.add_argument(
        "--json", action="store_true", help="Emit the rescaled-consistency report as JSON."
    )
    _add_manifest_argument(phylo_parsimony_rescaled_consistency)

    phylo_likelihood = phylo_subparsers.add_parser(
        "likelihood",
        help="Fit governed fixed-topology likelihood workflows on explicit trees and alignments.",
    )
    phylo_likelihood_subparsers = phylo_likelihood.add_subparsers(
        dest="phylo_likelihood_command",
        required=True,
    )
    phylo_likelihood_strict_clock = phylo_likelihood_subparsers.add_parser(
        "strict-clock",
        help="Fit one global JC69 strict-clock rate on one time-scaled tree and one DNA alignment.",
    )
    phylo_likelihood_strict_clock.add_argument("tree_path", type=Path)
    phylo_likelihood_strict_clock.add_argument("alignment_path", type=Path)
    phylo_likelihood_strict_clock.add_argument(
        "--model",
        default="jc69",
        choices=["jc69"],
        help="Likelihood model for the strict-clock fit.",
    )
    phylo_likelihood_strict_clock.add_argument(
        "--lower-clock-rate-bound",
        type=float,
        default=1e-6,
        help="Positive lower bound for the shared strict-clock rate search.",
    )
    phylo_likelihood_strict_clock.add_argument(
        "--upper-clock-rate-bound",
        type=float,
        default=5.0,
        help="Upper bound for the shared strict-clock rate search.",
    )
    phylo_likelihood_strict_clock.add_argument("--out-dir", required=True, type=Path)
    phylo_likelihood_strict_clock.add_argument(
        "--json", action="store_true", help="Emit the strict-clock report as JSON."
    )
    _add_manifest_argument(phylo_likelihood_strict_clock)
    phylo_likelihood_local_clock = phylo_likelihood_subparsers.add_parser(
        "local-clock",
        help="Fit user-defined branch and clade regimes as separate JC69 clock rates on one time-scaled tree and one DNA alignment.",
    )
    phylo_likelihood_local_clock.add_argument("tree_path", type=Path)
    phylo_likelihood_local_clock.add_argument("alignment_path", type=Path)
    phylo_likelihood_local_clock.add_argument("regime_path", type=Path)
    phylo_likelihood_local_clock.add_argument(
        "--model",
        default="jc69",
        choices=["jc69"],
        help="Likelihood model for the local-clock fit.",
    )
    phylo_likelihood_local_clock.add_argument(
        "--lower-clock-rate-bound",
        type=float,
        default=1e-6,
        help="Positive lower bound for every fitted local-clock rate.",
    )
    phylo_likelihood_local_clock.add_argument(
        "--upper-clock-rate-bound",
        type=float,
        default=5.0,
        help="Upper bound for every fitted local-clock rate.",
    )
    phylo_likelihood_local_clock.add_argument(
        "--max-coordinate-passes",
        type=int,
        default=12,
        help="Maximum number of coordinate-search passes over local-clock rates.",
    )
    phylo_likelihood_local_clock.add_argument("--out-dir", required=True, type=Path)
    phylo_likelihood_local_clock.add_argument(
        "--json", action="store_true", help="Emit the local-clock report as JSON."
    )
    _add_manifest_argument(phylo_likelihood_local_clock)
    phylo_likelihood_placement = phylo_likelihood_subparsers.add_parser(
        "placement",
        help="Place query DNA sequences on every reference-tree edge with JC69 likelihood.",
    )
    phylo_likelihood_placement.add_argument("reference_tree_path", type=Path)
    phylo_likelihood_placement.add_argument("reference_alignment_path", type=Path)
    phylo_likelihood_placement.add_argument("query_alignment_path", type=Path)
    phylo_likelihood_placement.add_argument(
        "--model",
        default="jc69",
        choices=["jc69"],
        help="Likelihood model for the placement fit.",
    )
    phylo_likelihood_placement.add_argument(
        "--lower-pendant-length-bound",
        type=float,
        default=1e-6,
        help="Positive lower bound for the query pendant branch-length search.",
    )
    phylo_likelihood_placement.add_argument(
        "--upper-pendant-length-bound",
        type=float,
        default=5.0,
        help="Upper bound for the query pendant branch-length search.",
    )
    phylo_likelihood_placement.add_argument(
        "--max-coordinate-passes",
        type=int,
        default=12,
        help="Maximum number of coordinate-search passes over distal and pendant placement lengths.",
    )
    phylo_likelihood_placement.add_argument("--out-dir", required=True, type=Path)
    phylo_likelihood_placement.add_argument(
        "--json", action="store_true", help="Emit the placement report as JSON."
    )
    _add_manifest_argument(phylo_likelihood_placement)
    phylo_likelihood_infer_tree = phylo_likelihood_subparsers.add_parser(
        "infer-tree",
        help="Infer one rooted nucleotide maximum-likelihood tree natively from one DNA FASTA alignment.",
    )
    phylo_likelihood_infer_tree.add_argument("alignment_path", type=Path)
    phylo_likelihood_infer_tree.add_argument(
        "--model",
        default="auto",
        choices=["auto", "jc69", "k80", "f81", "hky85", "gtr"],
        help="Fixed nucleotide model to use, or auto to select one native fixed-rate base model.",
    )
    phylo_likelihood_infer_tree.add_argument(
        "--model-selection-criterion",
        default="aic",
        choices=["aic", "aicc", "bic"],
        help="Information criterion used when --model auto selects one native base model.",
    )
    phylo_likelihood_infer_tree.add_argument(
        "--search-method",
        default="nni",
        choices=["nni", "spr", "tbr"],
        help="Native local topology search method applied after start-tree generation.",
    )
    phylo_likelihood_infer_tree.add_argument(
        "--start-tree-count",
        type=int,
        default=4,
        help="Total start trees to score, including the native stepwise-addition start tree.",
    )
    phylo_likelihood_infer_tree.add_argument(
        "--start-tree-seed",
        type=int,
        default=1,
        help="Deterministic seed for random additional start-tree generation.",
    )
    phylo_likelihood_infer_tree.add_argument(
        "--branch-reoptimization-policy",
        default="coordinate-branch-lengths",
        help="Branch-length reoptimization policy used during native local topology search.",
    )
    phylo_likelihood_infer_tree.add_argument(
        "--lower-branch-length-bound",
        type=float,
        default=0.0,
        help="Nonnegative lower branch-length bound applied during native optimization.",
    )
    phylo_likelihood_infer_tree.add_argument(
        "--upper-branch-length-bound",
        type=float,
        default=5.0,
        help="Upper branch-length bound applied during native optimization.",
    )
    phylo_likelihood_infer_tree.add_argument(
        "--max-coordinate-passes",
        type=int,
        default=12,
        help="Maximum branch-length coordinate-search passes for native optimization stages.",
    )
    phylo_likelihood_infer_tree.add_argument("--out-dir", required=True, type=Path)
    phylo_likelihood_infer_tree.add_argument(
        "--json", action="store_true", help="Emit the native ML tree-inference report as JSON."
    )
    _add_manifest_argument(phylo_likelihood_infer_tree)
    phylo_likelihood_bootstrap_tree = phylo_likelihood_subparsers.add_parser(
        "bootstrap-tree",
        help="Infer one native rooted nucleotide ML tree and bootstrap clade support by site resampling.",
    )
    phylo_likelihood_bootstrap_tree.add_argument("alignment_path", type=Path)
    phylo_likelihood_bootstrap_tree.add_argument(
        "--model",
        default="auto",
        choices=["auto", "jc69", "k80", "f81", "hky85", "gtr"],
        help="Fixed nucleotide model to use, or auto to select one native fixed-rate base model for the reference tree.",
    )
    phylo_likelihood_bootstrap_tree.add_argument(
        "--model-selection-criterion",
        default="aic",
        choices=["aic", "aicc", "bic"],
        help="Information criterion used when --model auto selects one native base model for the reference tree.",
    )
    phylo_likelihood_bootstrap_tree.add_argument(
        "--search-method",
        default="nni",
        choices=["nni", "spr", "tbr"],
        help="Native local topology search method applied to the reference tree and every bootstrap replicate.",
    )
    phylo_likelihood_bootstrap_tree.add_argument(
        "--start-tree-count",
        type=int,
        default=4,
        help="Total start trees to score for the reference tree and every bootstrap replicate.",
    )
    phylo_likelihood_bootstrap_tree.add_argument(
        "--start-tree-seed",
        type=int,
        default=1,
        help="Deterministic seed for random additional start-tree generation.",
    )
    phylo_likelihood_bootstrap_tree.add_argument(
        "--branch-reoptimization-policy",
        default="coordinate-branch-lengths",
        help="Branch-length reoptimization policy used during native local topology search.",
    )
    phylo_likelihood_bootstrap_tree.add_argument(
        "--replicate-count",
        type=int,
        default=100,
        help="Number of site-resampled native ML bootstrap replicates to infer.",
    )
    phylo_likelihood_bootstrap_tree.add_argument(
        "--bootstrap-seed",
        type=int,
        default=1,
        help="Deterministic seed for site-resampling draws.",
    )
    phylo_likelihood_bootstrap_tree.add_argument(
        "--lower-branch-length-bound",
        type=float,
        default=0.0,
        help="Nonnegative lower branch-length bound applied during native optimization.",
    )
    phylo_likelihood_bootstrap_tree.add_argument(
        "--upper-branch-length-bound",
        type=float,
        default=5.0,
        help="Upper branch-length bound applied during native optimization.",
    )
    phylo_likelihood_bootstrap_tree.add_argument(
        "--max-coordinate-passes",
        type=int,
        default=12,
        help="Maximum branch-length coordinate-search passes for native optimization stages.",
    )
    phylo_likelihood_bootstrap_tree.add_argument("--out-dir", required=True, type=Path)
    phylo_likelihood_bootstrap_tree.add_argument(
        "--json",
        action="store_true",
        help="Emit the native ML bootstrap-tree report as JSON.",
    )
    _add_manifest_argument(phylo_likelihood_bootstrap_tree)
    phylo_dating = phylo_subparsers.add_parser(
        "dating",
        help="Fit governed dated-tree workflows on rooted substitution trees and tip-date tables.",
    )
    phylo_dating_subparsers = phylo_dating.add_subparsers(
        dest="phylo_dating_command",
        required=True,
    )
    phylo_dating_least_squares = phylo_dating_subparsers.add_parser(
        "least-squares",
        help="Fit one rooted tree to fixed tip dates by closed-form least squares and emit dated-tree artifacts.",
    )
    phylo_dating_least_squares.add_argument("tree_path", type=Path)
    phylo_dating_least_squares.add_argument("metadata_path", type=Path)
    phylo_dating_least_squares.add_argument("--taxon-column")
    phylo_dating_least_squares.add_argument(
        "--date-column",
        default="date",
        help="Column containing numeric sampling dates or tip dates.",
    )
    phylo_dating_least_squares.add_argument("--out-dir", required=True, type=Path)
    phylo_dating_least_squares.add_argument(
        "--json", action="store_true", help="Emit the least-squares dating report as JSON."
    )
    _add_manifest_argument(phylo_dating_least_squares)
    phylo_dating_calibration_constraints = phylo_dating_subparsers.add_parser(
        "calibration-constraints",
        help="Resolve min, max, and fixed calibration bounds onto one rooted tree and report contradictory node windows before dating optimization.",
    )
    phylo_dating_calibration_constraints.add_argument("tree_path", type=Path)
    phylo_dating_calibration_constraints.add_argument("calibration_path", type=Path)
    phylo_dating_calibration_constraints.add_argument(
        "--out-dir",
        required=True,
        type=Path,
    )
    phylo_dating_calibration_constraints.add_argument(
        "--json",
        action="store_true",
        help="Emit the calibration-constraint report as JSON.",
    )
    _add_manifest_argument(phylo_dating_calibration_constraints)
    phylo_dating_penalized_likelihood = phylo_dating_subparsers.add_parser(
        "penalized-likelihood",
        help="Fit one rooted tree to fixed tip dates with a separate data score and rate-smoothing penalty.",
    )
    phylo_dating_penalized_likelihood.add_argument("tree_path", type=Path)
    phylo_dating_penalized_likelihood.add_argument("metadata_path", type=Path)
    phylo_dating_penalized_likelihood.add_argument("--taxon-column")
    phylo_dating_penalized_likelihood.add_argument(
        "--date-column",
        default="date",
        help="Column containing numeric sampling dates or tip dates.",
    )
    phylo_dating_penalized_likelihood.add_argument(
        "--smoothing-parameter",
        type=float,
        default=1.0,
        help="Strictly positive weight applied to the adjacent-rate smoothing penalty.",
    )
    phylo_dating_penalized_likelihood.add_argument(
        "--max-coordinate-passes",
        type=int,
        default=8,
        help="Maximum number of bounded coordinate-search passes over dated-node parameters.",
    )
    phylo_dating_penalized_likelihood.add_argument(
        "--out-dir",
        required=True,
        type=Path,
    )
    phylo_dating_penalized_likelihood.add_argument(
        "--json",
        action="store_true",
        help="Emit the penalized likelihood dating report as JSON.",
    )
    _add_manifest_argument(phylo_dating_penalized_likelihood)
    phylo_dating_penalized_likelihood_cross_validation = (
        phylo_dating_subparsers.add_parser(
            "penalized-likelihood-cross-validation",
            help="Select one penalized-dating smoothing value by held-out fixed-calibration prediction error.",
        )
    )
    phylo_dating_penalized_likelihood_cross_validation.add_argument(
        "tree_path",
        type=Path,
    )
    phylo_dating_penalized_likelihood_cross_validation.add_argument(
        "metadata_path",
        type=Path,
    )
    phylo_dating_penalized_likelihood_cross_validation.add_argument(
        "calibration_path",
        type=Path,
    )
    phylo_dating_penalized_likelihood_cross_validation.add_argument("--taxon-column")
    phylo_dating_penalized_likelihood_cross_validation.add_argument(
        "--date-column",
        default="date",
        help="Column containing numeric sampling dates or tip dates.",
    )
    phylo_dating_penalized_likelihood_cross_validation.add_argument(
        "--smoothing-parameters",
        nargs="+",
        type=float,
        default=[0.01, 0.1, 1.0, 10.0, 100.0],
        help="Positive smoothing-parameter candidates to score by held-out calibration prediction error.",
    )
    phylo_dating_penalized_likelihood_cross_validation.add_argument(
        "--max-coordinate-passes",
        type=int,
        default=8,
        help="Maximum number of bounded coordinate-search passes over dated-node parameters.",
    )
    phylo_dating_penalized_likelihood_cross_validation.add_argument(
        "--out-dir",
        required=True,
        type=Path,
    )
    phylo_dating_penalized_likelihood_cross_validation.add_argument(
        "--json",
        action="store_true",
        help="Emit the smoothing cross-validation report as JSON.",
    )
    _add_manifest_argument(phylo_dating_penalized_likelihood_cross_validation)
    phylo_dating_relaxed_rate_summary = phylo_dating_subparsers.add_parser(
        "relaxed-rate-summary",
        help="Summarize branch-specific rates by dividing substitution branch lengths by dated branch durations on one matched rooted tree pair.",
    )
    phylo_dating_relaxed_rate_summary.add_argument(
        "substitution_tree_path",
        type=Path,
    )
    phylo_dating_relaxed_rate_summary.add_argument("dated_tree_path", type=Path)
    phylo_dating_relaxed_rate_summary.add_argument(
        "--outlier-threshold",
        type=float,
        default=2.0,
        help="Absolute branch-rate z-score threshold used to flag outlier branches.",
    )
    phylo_dating_relaxed_rate_summary.add_argument(
        "--out-dir",
        required=True,
        type=Path,
    )
    phylo_dating_relaxed_rate_summary.add_argument(
        "--json",
        action="store_true",
        help="Emit the relaxed-rate branch summary as JSON.",
    )
    _add_manifest_argument(phylo_dating_relaxed_rate_summary)


def run_phylo_command(args: Any) -> int:
    executables = {
        "mafft": getattr(args, "mafft_executable", None),
        "trimal": getattr(args, "trimal_executable", None),
        "iqtree": getattr(args, "iqtree_executable", None),
        "fasttree": getattr(args, "fasttree_executable", None),
        "mrbayes": getattr(args, "mrbayes_executable", None),
        "beast": getattr(args, "beast_executable", None),
    }
    if args.phylo_command == "dating":
        if args.phylo_dating_command == "calibration-constraints":
            report = solve_dating_calibration_constraints(
                args.tree_path,
                args.calibration_path,
            )
            artifact_paths = write_dating_calibration_constraint_artifacts(
                args.out_dir,
                report,
            )
            outputs = _finalize_outputs(
                args,
                command="phylo",
                inputs=[args.tree_path, args.calibration_path],
                outputs=list(artifact_paths.values()),
            )
            _print_result(
                build_command_result(
                    command="phylo",
                    inputs=[args.tree_path, args.calibration_path],
                    outputs=outputs,
                    metrics={
                        "taxon_count": len(report.taxa),
                        "tip_count": report.tip_count,
                        "internal_node_count": report.internal_node_count,
                        "calibration_count": report.calibration_count,
                        "valid_calibration_count": report.valid_calibration_count,
                        "invalid_calibration_count": report.invalid_calibration_count,
                        "resolved_calibration_count": report.resolved_calibration_count,
                        "contradictory_calibration_count": report.contradictory_calibration_count,
                        "contradictory_node_count": report.contradictory_node_count,
                        "feasible": report.feasible,
                    },
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        if args.phylo_dating_command == "least-squares":
            report = fit_least_squares_dating_from_metadata(
                args.tree_path,
                args.metadata_path,
                taxon_column=args.taxon_column,
                date_column=args.date_column,
            )
            artifact_paths = write_least_squares_dating_artifacts(
                args.out_dir,
                report,
            )
            outputs = _finalize_outputs(
                args,
                command="phylo",
                inputs=[args.tree_path, args.metadata_path],
                outputs=list(artifact_paths.values()),
            )
            _print_result(
                build_command_result(
                    command="phylo",
                    inputs=[args.tree_path, args.metadata_path],
                    outputs=outputs,
                    metrics={
                        "taxon_count": len(report.taxa),
                        "tip_count": report.tip_count,
                        "internal_node_count": report.internal_node_count,
                        "branch_count": report.branch_count,
                        "estimated_clock_rate": report.estimated_clock_rate,
                        "root_date": report.root_date,
                        "residual_sum_squares": report.residual_sum_squares,
                        "exact_fit": report.exact_fit,
                        "converged": report.converged,
                    },
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        if args.phylo_dating_command == "penalized-likelihood":
            report = fit_penalized_likelihood_dating_from_metadata(
                args.tree_path,
                args.metadata_path,
                smoothing_parameter=args.smoothing_parameter,
                max_coordinate_passes=args.max_coordinate_passes,
                taxon_column=args.taxon_column,
                date_column=args.date_column,
            )
            artifact_paths = write_penalized_likelihood_dating_artifacts(
                args.out_dir,
                report,
            )
            outputs = _finalize_outputs(
                args,
                command="phylo",
                inputs=[args.tree_path, args.metadata_path],
                outputs=list(artifact_paths.values()),
            )
            _print_result(
                build_command_result(
                    command="phylo",
                    inputs=[args.tree_path, args.metadata_path],
                    outputs=outputs,
                    metrics={
                        "taxon_count": len(report.taxa),
                        "tip_count": report.tip_count,
                        "internal_node_count": report.internal_node_count,
                        "branch_count": report.branch_count,
                        "smoothing_parameter": report.smoothing_parameter,
                        "data_score": report.data_score,
                        "penalty_score": report.penalty_score,
                        "total_score": report.total_score,
                        "root_date": report.root_date,
                        "optimization_pass_count": report.optimization_pass_count,
                        "function_evaluation_count": report.function_evaluation_count,
                        "converged": report.converged,
                    },
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        if args.phylo_dating_command == "penalized-likelihood-cross-validation":
            report = cross_validate_penalized_likelihood_smoothing_from_metadata(
                args.tree_path,
                args.metadata_path,
                args.calibration_path,
                smoothing_parameters=args.smoothing_parameters,
                max_coordinate_passes=args.max_coordinate_passes,
                taxon_column=args.taxon_column,
                date_column=args.date_column,
            )
            artifact_paths = write_penalized_likelihood_cross_validation_artifacts(
                args.out_dir,
                report,
            )
            outputs = _finalize_outputs(
                args,
                command="phylo",
                inputs=[args.tree_path, args.metadata_path, args.calibration_path],
                outputs=list(artifact_paths.values()),
            )
            _print_result(
                build_command_result(
                    command="phylo",
                    inputs=[args.tree_path, args.metadata_path, args.calibration_path],
                    outputs=outputs,
                    metrics={
                        "taxon_count": len(report.taxa),
                        "tip_count": report.tip_count,
                        "internal_node_count": report.internal_node_count,
                        "branch_count": report.branch_count,
                        "usable_calibration_count": report.usable_calibration_count,
                        "candidate_count": report.candidate_count,
                        "selected_smoothing_parameter": report.selected_smoothing_parameter,
                        "selected_root_mean_squared_error": report.selected_root_mean_squared_error,
                        "final_total_score": report.selected_fit.total_score,
                        "final_root_date": report.selected_fit.root_date,
                        "final_converged": report.selected_fit.converged,
                    },
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        if args.phylo_dating_command == "relaxed-rate-summary":
            report = summarize_relaxed_rate_branches_from_paths(
                args.substitution_tree_path,
                args.dated_tree_path,
                outlier_threshold=args.outlier_threshold,
            )
            artifact_paths = write_relaxed_rate_branch_summary_artifacts(
                args.out_dir,
                report,
            )
            outputs = _finalize_outputs(
                args,
                command="phylo",
                inputs=[args.substitution_tree_path, args.dated_tree_path],
                outputs=list(artifact_paths.values()),
            )
            _print_result(
                build_command_result(
                    command="phylo",
                    inputs=[args.substitution_tree_path, args.dated_tree_path],
                    outputs=outputs,
                    metrics={
                        "taxon_count": len(report.taxa),
                        "tip_count": report.tip_count,
                        "internal_node_count": report.internal_node_count,
                        "branch_count": report.branch_count,
                        "outlier_threshold": report.outlier_threshold,
                        "mean_branch_rate": report.mean_branch_rate,
                        "minimum_branch_rate": report.minimum_branch_rate,
                        "maximum_branch_rate": report.maximum_branch_rate,
                        "outlier_count": report.outlier_count,
                        "top_outlier_branch_id": (
                            None
                            if not report.outlier_rows
                            else report.outlier_rows[0].branch_id
                        ),
                    },
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        raise EngineWorkflowError(
            "unknown phylo dating command",
            code="phylo_dating_command_unknown",
        )
    if args.phylo_command == "likelihood":
        if args.phylo_likelihood_command == "bootstrap-tree":
            report = bootstrap_nucleotide_likelihood_tree_inference_from_alignment(
                args.alignment_path,
                model_name=args.model,
                model_selection_criterion=args.model_selection_criterion,
                search_method=args.search_method,
                start_tree_count=args.start_tree_count,
                start_tree_seed=args.start_tree_seed,
                branch_reoptimization_policy=args.branch_reoptimization_policy,
                replicate_count=args.replicate_count,
                bootstrap_seed=args.bootstrap_seed,
                lower_branch_length_bound=args.lower_branch_length_bound,
                upper_branch_length_bound=args.upper_branch_length_bound,
                max_coordinate_passes=args.max_coordinate_passes,
            )
            artifact_paths = write_nucleotide_likelihood_bootstrap_artifacts(
                args.out_dir,
                report,
            )
            outputs = _finalize_outputs(
                args,
                command="phylo",
                inputs=[args.alignment_path],
                outputs=list(artifact_paths.values()),
            )
            _print_result(
                build_command_result(
                    command="phylo",
                    inputs=[args.alignment_path],
                    outputs=outputs,
                    metrics={
                        "selected_reference_model_name": report.selected_reference_model_name,
                        "taxon_count": report.taxon_count,
                        "site_count": report.site_count,
                        "pattern_count": report.pattern_count,
                        "replicate_count": report.replicate_count,
                        "support_row_count": len(report.clade_support_rows),
                        "reference_log_likelihood": report.reference_log_likelihood,
                        "search_method": report.search_method,
                    },
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        if args.phylo_likelihood_command == "infer-tree":
            report = infer_nucleotide_likelihood_tree_from_alignment(
                args.alignment_path,
                model_name=args.model,
                model_selection_criterion=args.model_selection_criterion,
                search_method=args.search_method,
                start_tree_count=args.start_tree_count,
                start_tree_seed=args.start_tree_seed,
                branch_reoptimization_policy=args.branch_reoptimization_policy,
                lower_branch_length_bound=args.lower_branch_length_bound,
                upper_branch_length_bound=args.upper_branch_length_bound,
                max_coordinate_passes=args.max_coordinate_passes,
            )
            artifact_paths = write_nucleotide_likelihood_tree_inference_artifacts(
                args.out_dir,
                report,
            )
            outputs = _finalize_outputs(
                args,
                command="phylo",
                inputs=[args.alignment_path],
                outputs=list(artifact_paths.values()),
            )
            _print_result(
                build_command_result(
                    command="phylo",
                    inputs=[args.alignment_path],
                    outputs=outputs,
                    metrics={
                        "selected_model_name": report.selected_model_name,
                        "taxon_count": report.taxon_count,
                        "site_count": report.site_count,
                        "pattern_count": report.pattern_count,
                        "start_tree_count": len(report.run_summaries),
                        "best_run_source_label": report.best_run_source_label,
                        "search_method": report.search_method,
                        "best_final_log_likelihood": report.best_final_log_likelihood,
                    },
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        if args.phylo_likelihood_command == "placement":
            report = place_queries_by_likelihood_from_alignment(
                args.reference_tree_path,
                args.reference_alignment_path,
                args.query_alignment_path,
                model=args.model,
                lower_pendant_length_bound=args.lower_pendant_length_bound,
                upper_pendant_length_bound=args.upper_pendant_length_bound,
                max_coordinate_passes=args.max_coordinate_passes,
            )
            artifact_paths = write_likelihood_placement_artifacts(
                args.out_dir,
                report,
            )
            outputs = _finalize_outputs(
                args,
                command="phylo",
                inputs=[
                    args.reference_tree_path,
                    args.reference_alignment_path,
                    args.query_alignment_path,
                ],
                outputs=list(artifact_paths.values()),
            )
            _print_result(
                build_command_result(
                    command="phylo",
                    inputs=[
                        args.reference_tree_path,
                        args.reference_alignment_path,
                        args.query_alignment_path,
                    ],
                    outputs=outputs,
                    metrics={
                        "model_name": report.model_name,
                        "reference_taxon_count": len(report.reference_taxa),
                        "edge_count": report.edge_count,
                        "query_count": report.query_count,
                        "site_count": report.site_count,
                        "placement_count": len(report.alternative_placements),
                        "total_function_evaluation_count": report.total_function_evaluation_count,
                    },
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        if args.phylo_likelihood_command == "local-clock":
            report = fit_local_clock_likelihood_from_alignment(
                args.tree_path,
                args.alignment_path,
                args.regime_path,
                model=args.model,
                lower_clock_rate_bound=args.lower_clock_rate_bound,
                upper_clock_rate_bound=args.upper_clock_rate_bound,
                max_coordinate_passes=args.max_coordinate_passes,
            )
            artifact_paths = write_local_clock_likelihood_artifacts(
                args.out_dir,
                report,
            )
            outputs = _finalize_outputs(
                args,
                command="phylo",
                inputs=[args.tree_path, args.alignment_path, args.regime_path],
                outputs=list(artifact_paths.values()),
            )
            _print_result(
                build_command_result(
                    command="phylo",
                    inputs=[args.tree_path, args.alignment_path, args.regime_path],
                    outputs=outputs,
                    metrics={
                        "model_name": report.model_name,
                        "taxon_count": len(report.taxa),
                        "site_count": report.site_count,
                        "pattern_count": report.pattern_count,
                        "branch_count": report.branch_count,
                        "regime_count": report.regime_count,
                        "parameter_count": report.parameter_count,
                        "optimized_log_likelihood": report.optimized_log_likelihood,
                        "strict_clock_aic": report.strict_clock_aic,
                        "local_clock_aic": report.aic,
                        "aic_delta_vs_strict_clock": report.aic_delta_vs_strict_clock,
                        "preferred_model_by_aic": report.preferred_model_by_aic,
                        "function_evaluation_count": report.function_evaluation_count,
                        "optimization_pass_count": report.optimization_pass_count,
                        "converged": report.converged,
                    },
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        if args.phylo_likelihood_command == "strict-clock":
            report = fit_strict_clock_likelihood_from_alignment(
                args.tree_path,
                args.alignment_path,
                model=args.model,
                lower_clock_rate_bound=args.lower_clock_rate_bound,
                upper_clock_rate_bound=args.upper_clock_rate_bound,
            )
            artifact_paths = write_strict_clock_likelihood_artifacts(
                args.out_dir,
                report,
            )
            outputs = _finalize_outputs(
                args,
                command="phylo",
                inputs=[args.tree_path, args.alignment_path],
                outputs=list(artifact_paths.values()),
            )
            _print_result(
                build_command_result(
                    command="phylo",
                    inputs=[args.tree_path, args.alignment_path],
                    outputs=outputs,
                    metrics={
                        "model_name": report.model_name,
                        "taxon_count": len(report.taxa),
                        "site_count": report.site_count,
                        "pattern_count": report.pattern_count,
                        "branch_count": report.branch_count,
                        "optimized_clock_rate": report.optimized_clock_rate,
                        "optimized_log_likelihood": report.optimized_log_likelihood,
                        "aic": report.aic,
                        "function_evaluation_count": report.function_evaluation_count,
                        "converged": report.converged,
                    },
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        raise EngineWorkflowError(
            "unknown phylo likelihood command",
            code="phylo_likelihood_command_unknown",
        )
    if args.phylo_command == "parsimony":
        if args.phylo_parsimony_command == "fitch":
            matrix = load_fitch_character_matrix(
                args.matrix_path,
                taxon_column=args.taxon_column,
            )
            character_weights = _load_parsimony_character_weights_argument(args)
            report = score_fitch(
                args.tree_path,
                matrix,
                character_weights=character_weights,
            )
            artifact_paths = write_fitch_artifacts(args.out_dir, report)
            metrics = {
                "algorithm": report.algorithm,
                "taxon_count": report.taxon_count,
                "character_count": report.character_count,
                "total_steps": report.total_steps,
                "total_weighted_score": report.total_weighted_score,
            }
        elif args.phylo_parsimony_command == "placement":
            matrix = load_fitch_character_matrix(
                args.matrix_path,
                taxon_column=args.taxon_column,
            )
            query_matrix = load_fitch_character_matrix(
                args.query_matrix_path,
                taxon_column=args.taxon_column,
            )
            character_weights = _load_parsimony_character_weights_argument(args)
            report = place_parsimony_queries(
                args.tree_path,
                matrix,
                query_matrix,
                method=args.method,
                character_weights=character_weights,
            )
            artifact_paths = write_parsimony_placement_artifacts(args.out_dir, report)
            metrics = {
                "algorithm": report.algorithm,
                "method": report.method,
                "reference_taxon_count": report.reference_taxon_count,
                "character_count": report.character_count,
                "edge_count": report.edge_count,
                "query_count": report.query_count,
                "reference_total_steps": report.reference_total_steps,
                "placement_count": len(report.alternative_rows),
                "equally_best_placement_count": sum(
                    1 for row in report.alternative_rows if row.is_equally_best
                ),
            }
        elif args.phylo_parsimony_command == "wagner":
            matrix = load_parsimony_character_matrix(
                args.matrix_path,
                taxon_column=args.taxon_column,
            )
            character_weights = _load_parsimony_character_weights_argument(args)
            state_order = _split_state_order(getattr(args, "state_order", None))
            report = score_wagner(
                args.tree_path,
                matrix,
                state_order=state_order,
                character_weights=character_weights,
            )
            artifact_paths = write_wagner_artifacts(args.out_dir, report)
            metrics = {
                "algorithm": report.algorithm,
                "taxon_count": report.taxon_count,
                "character_count": report.character_count,
                "total_cost": report.total_cost,
                "total_weighted_score": report.total_weighted_score,
            }
        elif args.phylo_parsimony_command == "sankoff":
            matrix = load_parsimony_character_matrix(
                args.matrix_path,
                taxon_column=args.taxon_column,
            )
            character_weights = _load_parsimony_character_weights_argument(args)
            cost_matrix = load_sankoff_cost_matrix(
                args.cost_matrix_path,
                allow_asymmetric_costs=args.allow_asymmetric_costs,
            )
            report = score_sankoff(
                args.tree_path,
                matrix,
                cost_matrix,
                allow_asymmetric_costs=args.allow_asymmetric_costs,
                character_weights=character_weights,
            )
            artifact_paths = write_sankoff_artifacts(args.out_dir, report)
            metrics = {
                "algorithm": report.algorithm,
                "taxon_count": report.taxon_count,
                "character_count": report.character_count,
                "total_cost": report.total_cost,
                "total_weighted_score": report.total_weighted_score,
                "validation_warning_count": len(report.validation_warnings),
            }
        elif args.phylo_parsimony_command == "dollo":
            matrix = load_parsimony_character_matrix(
                args.matrix_path,
                taxon_column=args.taxon_column,
            )
            character_weights = _load_parsimony_character_weights_argument(args)
            report = score_dollo(
                args.tree_path,
                matrix,
                character_weights=character_weights,
            )
            artifact_paths = write_dollo_artifacts(args.out_dir, report)
            metrics = {
                "algorithm": report.algorithm,
                "taxon_count": report.taxon_count,
                "character_count": report.character_count,
                "total_gains": report.total_gains,
                "total_losses": report.total_losses,
                "total_weighted_score": report.total_weighted_score,
            }
        elif args.phylo_parsimony_command == "camin-sokal":
            matrix = load_parsimony_character_matrix(
                args.matrix_path,
                taxon_column=args.taxon_column,
            )
            character_weights = _load_parsimony_character_weights_argument(args)
            report = score_camin_sokal(
                args.tree_path,
                matrix,
                character_weights=character_weights,
            )
            artifact_paths = write_camin_sokal_artifacts(args.out_dir, report)
            metrics = {
                "algorithm": report.algorithm,
                "taxon_count": report.taxon_count,
                "character_count": report.character_count,
                "root_state": report.root_state,
                "total_gains": report.total_gains,
                "total_weighted_score": report.total_weighted_score,
            }
        elif args.phylo_parsimony_command == "acctran":
            matrix = load_parsimony_character_matrix(
                args.matrix_path,
                taxon_column=args.taxon_column,
            )
            character_weights = _load_parsimony_character_weights_argument(args)
            report = reconstruct_acctran(
                args.tree_path,
                matrix,
                character_weights=character_weights,
            )
            artifact_paths = write_parsimony_reconstruction_artifacts(
                args.out_dir,
                report,
            )
            metrics = {
                "algorithm": report.algorithm,
                "taxon_count": report.taxon_count,
                "character_count": report.character_count,
                "total_steps": report.total_steps,
                "total_weighted_score": report.total_weighted_score,
            }
        elif args.phylo_parsimony_command == "deltran":
            matrix = load_parsimony_character_matrix(
                args.matrix_path,
                taxon_column=args.taxon_column,
            )
            character_weights = _load_parsimony_character_weights_argument(args)
            report = reconstruct_deltran(
                args.tree_path,
                matrix,
                character_weights=character_weights,
            )
            artifact_paths = write_parsimony_reconstruction_artifacts(
                args.out_dir,
                report,
            )
            metrics = {
                "algorithm": report.algorithm,
                "taxon_count": report.taxon_count,
                "character_count": report.character_count,
                "total_steps": report.total_steps,
                "total_weighted_score": report.total_weighted_score,
            }
        elif args.phylo_parsimony_command == "bootstrap":
            matrix = load_parsimony_character_matrix(
                args.matrix_path,
                taxon_column=args.taxon_column,
            )
            character_weights = _load_parsimony_character_weights_argument(args)
            cost_matrix = (
                load_sankoff_cost_matrix(
                    args.cost_matrix_path,
                    allow_asymmetric_costs=args.allow_asymmetric_costs,
                )
                if getattr(args, "cost_matrix_path", None) is not None
                else None
            )
            state_order = _split_state_order(getattr(args, "state_order", None))
            report = bootstrap_parsimony(
                matrix,
                method=args.method,
                replicate_count=args.replicate_count,
                random_seed=args.seed,
                state_order=state_order,
                cost_matrix=cost_matrix,
                allow_asymmetric_costs=args.allow_asymmetric_costs,
                character_weights=character_weights,
            )
            artifact_paths = write_parsimony_bootstrap_artifacts(args.out_dir, report)
            metrics = {
                "algorithm": report.algorithm,
                "method": report.method,
                "taxon_count": report.taxon_count,
                "character_count": report.character_count,
                "replicate_count": report.replicate_count,
                "candidate_tree_count": report.candidate_tree_count,
                "reference_score": report.reference_score,
                "support_row_count": len(report.clade_support_rows),
            }
        elif args.phylo_parsimony_command == "jackknife":
            matrix = load_parsimony_character_matrix(
                args.matrix_path,
                taxon_column=args.taxon_column,
            )
            character_weights = _load_parsimony_character_weights_argument(args)
            cost_matrix = (
                load_sankoff_cost_matrix(
                    args.cost_matrix_path,
                    allow_asymmetric_costs=args.allow_asymmetric_costs,
                )
                if getattr(args, "cost_matrix_path", None) is not None
                else None
            )
            state_order = _split_state_order(getattr(args, "state_order", None))
            report = jackknife_parsimony(
                matrix,
                method=args.method,
                replicate_count=args.replicate_count,
                random_seed=args.seed,
                retain_probability=args.retain_probability,
                state_order=state_order,
                cost_matrix=cost_matrix,
                allow_asymmetric_costs=args.allow_asymmetric_costs,
                character_weights=character_weights,
            )
            artifact_paths = write_parsimony_jackknife_artifacts(args.out_dir, report)
            metrics = {
                "algorithm": report.algorithm,
                "method": report.method,
                "taxon_count": report.taxon_count,
                "character_count": report.character_count,
                "replicate_count": report.replicate_count,
                "retain_probability": report.retain_probability,
                "candidate_tree_count": report.candidate_tree_count,
                "reference_score": report.reference_score,
                "support_row_count": len(report.clade_support_rows),
            }
        elif args.phylo_parsimony_command == "equal-best-consensus":
            matrix = load_parsimony_character_matrix(
                args.matrix_path,
                taxon_column=args.taxon_column,
            )
            character_weights = _load_parsimony_character_weights_argument(args)
            cost_matrix = (
                load_sankoff_cost_matrix(
                    args.cost_matrix_path,
                    allow_asymmetric_costs=args.allow_asymmetric_costs,
                )
                if getattr(args, "cost_matrix_path", None) is not None
                else None
            )
            state_order = _split_state_order(getattr(args, "state_order", None))
            report = summarize_equal_best_parsimony_trees(
                matrix,
                method=args.method,
                state_order=state_order,
                cost_matrix=cost_matrix,
                allow_asymmetric_costs=args.allow_asymmetric_costs,
                character_weights=character_weights,
                max_retained_equal_best_trees=args.max_retained_equal_best_trees,
            )
            artifact_paths = write_parsimony_equal_best_consensus_artifacts(
                args.out_dir,
                report,
            )
            metrics = {
                "algorithm": report.algorithm,
                "method": report.method,
                "taxon_count": report.taxon_count,
                "character_count": report.character_count,
                "candidate_tree_count": report.candidate_tree_count,
                "best_score": report.best_score,
                "equal_best_tree_count": report.equal_best_tree_count,
                "retained_equal_best_tree_count": report.retained_equal_best_tree_count,
                "retained_all_equal_best_trees": report.retained_all_equal_best_trees,
            }
        elif args.phylo_parsimony_command == "bremer-support":
            matrix = load_parsimony_character_matrix(
                args.matrix_path,
                taxon_column=args.taxon_column,
            )
            character_weights = _load_parsimony_character_weights_argument(args)
            cost_matrix = (
                load_sankoff_cost_matrix(
                    args.cost_matrix_path,
                    allow_asymmetric_costs=args.allow_asymmetric_costs,
                )
                if getattr(args, "cost_matrix_path", None) is not None
                else None
            )
            state_order = _split_state_order(getattr(args, "state_order", None))
            report = compute_parsimony_bremer_support(
                args.tree_path,
                matrix,
                method=args.method,
                state_order=state_order,
                cost_matrix=cost_matrix,
                allow_asymmetric_costs=args.allow_asymmetric_costs,
                character_weights=character_weights,
            )
            artifact_paths = write_parsimony_bremer_support_artifacts(
                args.out_dir,
                report,
            )
            metrics = {
                "algorithm": report.algorithm,
                "method": report.method,
                "taxon_count": report.taxon_count,
                "character_count": report.character_count,
                "candidate_tree_count": report.candidate_tree_count,
                "reference_tree_score": report.reference_tree_score,
                "optimal_score": report.optimal_score,
                "reference_tree_is_optimal": report.reference_tree_is_optimal,
                "bremer_row_count": len(report.bremer_rows),
            }
        elif args.phylo_parsimony_command == "stepwise-addition":
            matrix = load_parsimony_character_matrix(
                args.matrix_path,
                taxon_column=args.taxon_column,
            )
            character_weights = _load_parsimony_character_weights_argument(args)
            cost_matrix = (
                load_sankoff_cost_matrix(
                    args.cost_matrix_path,
                    allow_asymmetric_costs=args.allow_asymmetric_costs,
                )
                if getattr(args, "cost_matrix_path", None) is not None
                else None
            )
            state_order = _split_state_order(getattr(args, "state_order", None))
            insertion_order = _split_taxon_order(getattr(args, "insertion_order", None))
            _tree, report = build_parsimony_stepwise_addition_tree(
                matrix,
                method=args.method,
                insertion_order=insertion_order,
                state_order=state_order,
                cost_matrix=cost_matrix,
                allow_asymmetric_costs=args.allow_asymmetric_costs,
                character_weights=character_weights,
            )
            artifact_paths = write_stepwise_addition_artifacts(args.out_dir, report)
            metrics = {
                "algorithm": report.algorithm,
                "method": args.method,
                "objective_name": report.objective_name,
                "taxon_count": report.tip_count,
                "character_count": matrix.character_count,
                "final_score": report.final_score,
                "insertion_step_count": len(report.trace_rows),
            }
        elif args.phylo_parsimony_command == "nni-search":
            matrix = load_parsimony_character_matrix(
                args.matrix_path,
                taxon_column=args.taxon_column,
            )
            character_weights = _load_parsimony_character_weights_argument(args)
            cost_matrix = (
                load_sankoff_cost_matrix(
                    args.cost_matrix_path,
                    allow_asymmetric_costs=args.allow_asymmetric_costs,
                )
                if getattr(args, "cost_matrix_path", None) is not None
                else None
            )
            state_order = _split_state_order(getattr(args, "state_order", None))
            report = search_parsimony_nni(
                args.tree_path,
                matrix,
                method=args.method,
                state_order=state_order,
                cost_matrix=cost_matrix,
                allow_asymmetric_costs=args.allow_asymmetric_costs,
                character_weights=character_weights,
            )
            artifact_paths = write_parsimony_nni_artifacts(args.out_dir, report)
            metrics = {
                "algorithm": report.algorithm,
                "method": report.method,
                "taxon_count": report.taxon_count,
                "character_count": report.character_count,
                "start_score": report.start_score,
                "final_score": report.final_score,
                "accepted_move_count": report.accepted_move_count,
                "evaluated_neighbor_count": report.evaluated_neighbor_count,
                "stopping_reason": report.stopping_reason,
            }
        elif args.phylo_parsimony_command == "spr-search":
            matrix = load_parsimony_character_matrix(
                args.matrix_path,
                taxon_column=args.taxon_column,
            )
            character_weights = _load_parsimony_character_weights_argument(args)
            cost_matrix = (
                load_sankoff_cost_matrix(
                    args.cost_matrix_path,
                    allow_asymmetric_costs=args.allow_asymmetric_costs,
                )
                if getattr(args, "cost_matrix_path", None) is not None
                else None
            )
            state_order = _split_state_order(getattr(args, "state_order", None))
            report = search_parsimony_spr(
                args.tree_path,
                matrix,
                method=args.method,
                state_order=state_order,
                cost_matrix=cost_matrix,
                allow_asymmetric_costs=args.allow_asymmetric_costs,
                character_weights=character_weights,
            )
            artifact_paths = write_parsimony_spr_artifacts(args.out_dir, report)
            metrics = {
                "algorithm": report.algorithm,
                "method": report.method,
                "taxon_count": report.taxon_count,
                "character_count": report.character_count,
                "start_score": report.start_score,
                "final_score": report.final_score,
                "accepted_move_count": report.accepted_move_count,
                "evaluated_neighbor_count": report.evaluated_neighbor_count,
                "stopping_reason": report.stopping_reason,
            }
        elif args.phylo_parsimony_command == "ratchet":
            matrix = load_parsimony_character_matrix(
                args.matrix_path,
                taxon_column=args.taxon_column,
            )
            character_weights = _load_parsimony_character_weights_argument(args)
            cost_matrix = (
                load_sankoff_cost_matrix(
                    args.cost_matrix_path,
                    allow_asymmetric_costs=args.allow_asymmetric_costs,
                )
                if getattr(args, "cost_matrix_path", None) is not None
                else None
            )
            state_order = _split_state_order(getattr(args, "state_order", None))
            report = run_parsimony_ratchet(
                args.tree_path,
                matrix,
                method=args.method,
                cycle_count=args.cycle_count,
                random_seed=args.seed,
                perturbed_character_count=args.perturbed_character_count,
                perturbation_factor=args.perturbation_factor,
                state_order=state_order,
                cost_matrix=cost_matrix,
                allow_asymmetric_costs=args.allow_asymmetric_costs,
                character_weights=character_weights,
            )
            artifact_paths = write_parsimony_ratchet_artifacts(args.out_dir, report)
            metrics = {
                "algorithm": report.algorithm,
                "method": report.method,
                "taxon_count": report.taxon_count,
                "character_count": report.character_count,
                "cycle_count": report.cycle_count,
                "random_seed": report.random_seed,
                "perturbed_character_count": report.perturbed_character_count,
                "perturbation_factor": report.perturbation_factor,
                "start_score": report.start_score,
                "final_score": report.final_score,
                "best_score": report.best_score,
                "best_tree_history_count": len(report.best_tree_history_rows),
            }
        elif args.phylo_parsimony_command == "tree-length":
            matrix = load_parsimony_character_matrix(
                args.matrix_path,
                taxon_column=args.taxon_column,
            )
            character_weights = (
                load_parsimony_character_weights(args.character_weights_path)
                if getattr(args, "character_weights_path", None) is not None
                else None
            )
            cost_matrix = (
                load_sankoff_cost_matrix(
                    args.cost_matrix_path,
                    allow_asymmetric_costs=args.allow_asymmetric_costs,
                )
                if getattr(args, "cost_matrix_path", None) is not None
                else None
            )
            state_order = _split_state_order(getattr(args, "state_order", None))
            report = tree_length(
                args.tree_path,
                matrix,
                method=args.method,
                state_order=state_order,
                cost_matrix=cost_matrix,
                allow_asymmetric_costs=args.allow_asymmetric_costs,
                character_weights=character_weights,
            )
            artifact_paths = write_parsimony_tree_length_artifacts(args.out_dir, report)
            metrics = {
                "algorithm": report.algorithm,
                "method": report.method,
                "taxon_count": report.taxon_count,
                "character_count": report.character_count,
                "raw_total_score": report.raw_total_score,
                "total_score": report.total_score,
            }
        elif args.phylo_parsimony_command == "consistency-index":
            matrix = load_parsimony_character_matrix(
                args.matrix_path,
                taxon_column=args.taxon_column,
            )
            state_order = _split_state_order(getattr(args, "state_order", None))
            report = consistency_index(
                args.tree_path,
                matrix,
                method=args.method,
                state_order=state_order,
            )
            artifact_paths = write_parsimony_consistency_artifacts(args.out_dir, report)
            metrics = {
                "algorithm": report.algorithm,
                "method": report.method,
                "taxon_count": report.taxon_count,
                "character_count": report.character_count,
                "included_character_count": report.included_character_count,
                "excluded_character_count": report.excluded_character_count,
                "minimum_possible_steps_total": report.minimum_possible_steps_total,
                "observed_steps_total": report.observed_steps_total,
                "consistency_index": report.consistency_index,
                "undefined_reason": report.undefined_reason,
            }
        elif args.phylo_parsimony_command == "retention-index":
            matrix = load_parsimony_character_matrix(
                args.matrix_path,
                taxon_column=args.taxon_column,
            )
            report = retention_index(
                args.tree_path,
                matrix,
                method=args.method,
            )
            artifact_paths = write_parsimony_retention_artifacts(args.out_dir, report)
            metrics = {
                "algorithm": report.algorithm,
                "method": report.method,
                "taxon_count": report.taxon_count,
                "character_count": report.character_count,
                "included_character_count": report.included_character_count,
                "excluded_character_count": report.excluded_character_count,
                "minimum_possible_steps_total": report.minimum_possible_steps_total,
                "maximum_possible_steps_total": report.maximum_possible_steps_total,
                "observed_steps_total": report.observed_steps_total,
                "retention_index": report.retention_index,
                "undefined_reason": report.undefined_reason,
            }
        elif args.phylo_parsimony_command == "rescaled-consistency-index":
            matrix = load_parsimony_character_matrix(
                args.matrix_path,
                taxon_column=args.taxon_column,
            )
            report = rescaled_consistency_index(
                args.tree_path,
                matrix,
                method=args.method,
            )
            artifact_paths = write_parsimony_rescaled_consistency_artifacts(
                args.out_dir,
                report,
            )
            metrics = {
                "algorithm": report.algorithm,
                "method": report.method,
                "taxon_count": report.taxon_count,
                "character_count": report.character_count,
                "ci": report.ci,
                "ri": report.ri,
                "rc": report.rc,
                "undefined_reason": report.undefined_reason,
            }
        else:
            raise EngineWorkflowError(
                "unknown phylo parsimony command",
                code="phylo_parsimony_command_unknown",
            )
        parsimony_inputs = [
            *([args.tree_path] if hasattr(args, "tree_path") else []),
            args.matrix_path,
            *(
                [args.query_matrix_path]
                if hasattr(args, "query_matrix_path")
                else []
            ),
            *(
                [args.cost_matrix_path]
                if hasattr(args, "cost_matrix_path")
                and (
                    args.phylo_parsimony_command == "sankoff"
                    or (
                        args.phylo_parsimony_command == "jackknife"
                        and getattr(args, "cost_matrix_path", None) is not None
                    )
                    or (
                        args.phylo_parsimony_command == "bootstrap"
                        and getattr(args, "cost_matrix_path", None) is not None
                    )
                    or (
                        args.phylo_parsimony_command == "equal-best-consensus"
                        and getattr(args, "cost_matrix_path", None) is not None
                    )
                    or (
                        args.phylo_parsimony_command == "stepwise-addition"
                        and getattr(args, "cost_matrix_path", None) is not None
                    )
                    or (
                        args.phylo_parsimony_command == "tree-length"
                        and getattr(args, "cost_matrix_path", None) is not None
                    )
                    or (
                        args.phylo_parsimony_command == "nni-search"
                        and getattr(args, "cost_matrix_path", None) is not None
                    )
                    or (
                        args.phylo_parsimony_command == "spr-search"
                        and getattr(args, "cost_matrix_path", None) is not None
                    )
                    or (
                        args.phylo_parsimony_command == "ratchet"
                        and getattr(args, "cost_matrix_path", None) is not None
                    )
                )
                else []
            ),
            *(
                [args.character_weights_path]
                if hasattr(args, "character_weights_path")
                and getattr(args, "character_weights_path", None) is not None
                else []
            ),
        ]
        outputs = _finalize_outputs(
            args,
            command="phylo",
            inputs=parsimony_inputs,
            outputs=list(artifact_paths.values()),
        )
        _print_result(
            build_command_result(
                command="phylo",
                inputs=parsimony_inputs,
                outputs=outputs,
                metrics=metrics,
                data=report,
                warnings=_parsimony_warning_messages(report),
            ),
            json_output=args.json,
        )
        return 0
    if args.phylo_command == "run":
        report = run_phylo_workflow_config(args.config_path)
        outputs = _finalize_outputs(
            args,
            command="phylo",
            inputs=[args.config_path],
            outputs=[
                report.fasta_to_tree_report.manifest_path,
                report.bundle_report.bundle_root,
                report.bundle_report.bundle_manifest_path,
                report.bundle_report.report_path,
            ],
        )
        _print_result(
            build_command_result(
                command="phylo",
                inputs=[args.config_path],
                outputs=outputs,
                warnings=report.warnings + report.notes,
                metrics={
                    "workflow": report.workflow,
                    "selected_workflow_status": (
                        report.selected_workflow_status.readiness_status
                    ),
                    "metadata_present": report.workflow_config.metadata_path
                    is not None,
                    "traits_present": report.workflow_config.traits_path is not None,
                    "alignment_mode": report.workflow_config.alignment_mode,
                    "trimming_mode": report.workflow_config.trimming_mode,
                    "bootstrap_replicates": report.workflow_config.bootstrap_replicates,
                    "iqtree_seed": report.workflow_config.iqtree_seed,
                    "iqtree_threads": report.workflow_config.iqtree_threads,
                    "timeout_seconds": report.workflow_config.timeout_seconds,
                    "bundle_file_count": report.bundle_report.file_count,
                    "bundle_validation_passed": report.bundle_validation.valid,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.phylo_command == "preflight":
        report = inspect_external_engine_preflight(
            executables=executables,
            selected_workflow=args.workflow,
        )
        selected_workflow_status = None
        if args.workflow is not None:
            selected_workflow_status = require_preflight_workflow(
                report, workflow_id=args.workflow
            ).readiness_status
        inputs = [] if args.workflow is None else [args.workflow]
        outputs = _finalize_outputs(args, command="phylo", inputs=inputs)
        _print_result(
            build_command_result(
                command="phylo",
                inputs=inputs,
                outputs=outputs,
                metrics={
                    "engine_count": len(report.engines),
                    "available_engine_count": sum(
                        1 for engine in report.engines if engine.available
                    ),
                    "workflow_count": len(report.workflows),
                    "runnable_workflow_count": sum(
                        1 for workflow in report.workflows if workflow.runnable
                    ),
                    "selected_workflow": args.workflow,
                    "selected_workflow_status": selected_workflow_status,
                    "overall_status": report.overall_status,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.phylo_command == "bundle":
        report = export_workflow_result_bundle(
            args.manifest_path,
            bundle_root=args.out_dir,
        )
        validation = validate_workflow_result_bundle(report.bundle_root)
        outputs = _finalize_outputs(
            args,
            command="phylo",
            inputs=[args.manifest_path],
            outputs=[
                report.bundle_root,
                report.bundle_manifest_path,
                report.report_path,
            ],
        )
        _print_result(
            build_command_result(
                command="phylo",
                inputs=[args.manifest_path],
                outputs=outputs,
                warnings=report.notes,
                metrics={
                    "workflow": report.workflow,
                    "file_count": report.file_count,
                    "copied_input_count": report.copied_input_count,
                    "copied_output_count": report.copied_output_count,
                    "copied_step_manifest_count": report.copied_step_manifest_count,
                    "copied_step_output_count": report.copied_step_output_count,
                    "copied_report_count": report.copied_report_count,
                    "missing_input_count": len(report.missing_input_paths),
                    "validation_passed": validation.valid,
                },
                data={"bundle": report, "validation": validation},
            ),
            json_output=args.json,
        )
        return 0
    if args.phylo_command == "validate-bundle":
        report = validate_workflow_result_bundle(args.bundle_root)
        if not report.valid:
            raise EngineWorkflowError(
                "workflow result bundle validation failed",
                code="workflow_bundle_validation_failed",
                details={
                    "issue_count": len(report.issues),
                    "bundle_root": str(args.bundle_root),
                },
            )
        outputs = _finalize_outputs(
            args,
            command="phylo",
            inputs=[args.bundle_root],
        )
        _print_result(
            build_command_result(
                command="phylo",
                inputs=[args.bundle_root],
                outputs=outputs,
                metrics={
                    "workflow": report.workflow,
                    "file_count": report.file_count,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    report = replay_workflow_manifest(
        args.manifest_path,
        out_dir=args.out_dir,
        executables=executables,
    )
    outputs = _finalize_outputs(
        args,
        command="phylo",
        inputs=[args.manifest_path],
        outputs=[report.replay_manifest_path],
    )
    _print_result(
        build_command_result(
            command="phylo",
            inputs=[args.manifest_path],
            outputs=outputs,
            metrics={
                "workflow": report.workflow,
                "input_drift_count": len(report.input_drift),
                "changed_input_count": sum(
                    1 for drift in report.input_drift if not drift.matched
                ),
                "engine_version_drift_count": sum(
                    1 for drift in report.engine_version_drift if not drift.matched
                ),
                "comparison_count": len(report.comparisons),
                "outputs_equivalent": report.outputs_equivalent,
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0


def _split_state_order(raw: str | None) -> list[str] | None:
    if raw is None:
        return None
    values = [value.strip() for value in raw.split(",")]
    return [value for value in values if value]


def _split_taxon_order(raw: str | None) -> list[str] | None:
    if raw is None:
        return None
    values = [value.strip() for value in raw.split(",")]
    return [value for value in values if value]


def _add_parsimony_character_weights_argument(parser: Any) -> None:
    parser.add_argument(
        "--character-weights",
        dest="character_weights_path",
        type=Path,
        help="Optional TSV with character_id and weight columns.",
    )


def _load_parsimony_character_weights_argument(args: Any):
    weights_path = getattr(args, "character_weights_path", None)
    if weights_path is None:
        return None
    return load_parsimony_character_weights(weights_path)


def _parsimony_warning_messages(report: Any) -> list[str]:
    validation_warnings = getattr(report, "validation_warnings", [])
    return [warning.message for warning in validation_warnings]
