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
    consistency_index,
    load_fitch_character_matrix,
    load_parsimony_character_weights,
    load_parsimony_character_matrix,
    load_sankoff_cost_matrix,
    reconstruct_acctran,
    reconstruct_deltran,
    rescaled_consistency_index,
    retention_index,
    score_camin_sokal,
    score_dollo,
    score_fitch,
    score_sankoff,
    score_wagner,
    tree_length,
    write_parsimony_consistency_artifacts,
    write_parsimony_bootstrap_artifacts,
    write_parsimony_reconstruction_artifacts,
    write_parsimony_rescaled_consistency_artifacts,
    write_parsimony_retention_artifacts,
    write_parsimony_tree_length_artifacts,
    write_camin_sokal_artifacts,
    write_dollo_artifacts,
    write_fitch_artifacts,
    write_sankoff_artifacts,
    write_wagner_artifacts,
)
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


def run_phylo_command(args: Any) -> int:
    executables = {
        "mafft": getattr(args, "mafft_executable", None),
        "trimal": getattr(args, "trimal_executable", None),
        "iqtree": getattr(args, "iqtree_executable", None),
        "fasttree": getattr(args, "fasttree_executable", None),
        "mrbayes": getattr(args, "mrbayes_executable", None),
        "beast": getattr(args, "beast_executable", None),
    }
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
                [args.cost_matrix_path]
                if hasattr(args, "cost_matrix_path")
                and (
                    args.phylo_parsimony_command == "sankoff"
                    or (
                        args.phylo_parsimony_command == "bootstrap"
                        and getattr(args, "cost_matrix_path", None) is not None
                    )
                    or (
                        args.phylo_parsimony_command == "tree-length"
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
