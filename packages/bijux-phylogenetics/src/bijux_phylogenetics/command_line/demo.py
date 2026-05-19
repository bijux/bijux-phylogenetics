from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.core.demo import run_capability_demo
from bijux_phylogenetics.datasets import (
    run_avian_reproductive_trait_demo,
    run_catarrhine_mitogenome_five_locus_panel_demo,
    run_central_european_seashore_flora_demo,
    run_gnathostome_ortholog_protein_benchmark_demo,
    run_influenza_a_ha_reference_demo,
    run_pleistocene_bear_cytb_fragment_demo,
    run_primate_comparative_demo,
    run_rabies_cross_host_geography_panel_demo,
    run_rabies_cross_host_panel_demo,
    run_rabies_geographic_transition_panel_demo,
    run_rabies_method_sensitivity_panel_demo,
)
from bijux_phylogenetics.datasets.continuous_mode_recovery import (
    run_continuous_mode_recovery_panel_demo,
)
from bijux_phylogenetics.datasets.discrete_mode_recovery import (
    run_discrete_mode_recovery_panel_demo,
)
from bijux_phylogenetics.datasets.data_quality_stress import (
    run_catarrhine_data_quality_stress_panel_demo,
)
from bijux_phylogenetics.datasets.known_answer_reference import (
    run_known_answer_reference_demo,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_demo_command(subparsers: Any) -> None:
    demo = subparsers.add_parser(
        get_command_spec("demo").name, help=get_command_spec("demo").summary
    )
    demo_subparsers = demo.add_subparsers(dest="demo_command", required=True)
    demo_run = demo_subparsers.add_parser(
        "run", help="Run the repository capability demo workflow."
    )
    demo_run.add_argument("--out", required=True, type=Path)
    demo_run.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_run)
    demo_primate = demo_subparsers.add_parser(
        "primate-comparative",
        help="Materialize the packaged primate dataset and comparative workflow outputs.",
    )
    demo_primate.add_argument("--out", required=True, type=Path)
    demo_primate.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_primate)
    demo_birds = demo_subparsers.add_parser(
        "avian-reproductive-traits",
        help="Materialize the packaged avian reproductive dataset and workflow outputs.",
    )
    demo_birds.add_argument("--out", required=True, type=Path)
    demo_birds.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_birds)
    demo_plants = demo_subparsers.add_parser(
        "central-european-seashore-flora",
        help="Materialize the packaged Central European plant dataset and workflow outputs.",
    )
    demo_plants.add_argument("--out", required=True, type=Path)
    demo_plants.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_plants)
    demo_viruses = demo_subparsers.add_parser(
        "influenza-a-ha-reference-panel",
        help="Materialize the packaged influenza A HA dataset and rerun the sequence-to-tree workflow outputs.",
    )
    demo_viruses.add_argument("--out", required=True, type=Path)
    demo_viruses.add_argument("--mafft-executable", type=str)
    demo_viruses.add_argument("--trimal-executable", type=str)
    demo_viruses.add_argument("--iqtree-executable", type=str)
    demo_viruses.add_argument("--iqtree-seed", type=int, default=1)
    demo_viruses.add_argument("--iqtree-threads", type=int, default=1)
    demo_viruses.add_argument("--bootstrap-replicates", type=int, default=1000)
    demo_viruses.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_viruses)
    demo_protein_benchmark = demo_subparsers.add_parser(
        "gnathostome-ortholog-protein-benchmark",
        help="Materialize the packaged gnathostome protein benchmark and rerun the governed amino-acid sequence-to-tree outputs.",
    )
    demo_protein_benchmark.add_argument("--out", required=True, type=Path)
    demo_protein_benchmark.add_argument("--mafft-executable", type=str)
    demo_protein_benchmark.add_argument("--trimal-executable", type=str)
    demo_protein_benchmark.add_argument("--iqtree-executable", type=str)
    demo_protein_benchmark.add_argument("--iqtree-seed", type=int, default=1)
    demo_protein_benchmark.add_argument("--iqtree-threads", type=int, default=1)
    demo_protein_benchmark.add_argument(
        "--bootstrap-replicates", type=int, default=1000
    )
    demo_protein_benchmark.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_protein_benchmark)
    demo_ancient_dna = demo_subparsers.add_parser(
        "pleistocene-bear-cytb-fragments",
        help="Materialize the packaged ancient-DNA-style bear dataset and rerun the sequence-to-tree workflow outputs.",
    )
    demo_ancient_dna.add_argument("--out", required=True, type=Path)
    demo_ancient_dna.add_argument("--mafft-executable", type=str)
    demo_ancient_dna.add_argument("--trimal-executable", type=str)
    demo_ancient_dna.add_argument("--iqtree-executable", type=str)
    demo_ancient_dna.add_argument("--iqtree-seed", type=int, default=1)
    demo_ancient_dna.add_argument("--iqtree-threads", type=int, default=1)
    demo_ancient_dna.add_argument("--bootstrap-replicates", type=int, default=1000)
    demo_ancient_dna.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_ancient_dna)
    demo_rabies = demo_subparsers.add_parser(
        "rabies-cross-host-panel",
        help="Materialize the packaged rabies host-switching dataset and rerun the governed host-transition review outputs.",
    )
    demo_rabies.add_argument("--out", required=True, type=Path)
    demo_rabies.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_rabies)
    demo_rabies_geography = demo_subparsers.add_parser(
        "rabies-geographic-transition-panel",
        help="Materialize the packaged rabies geography dataset and rerun the governed geographic transition review outputs.",
    )
    demo_rabies_geography.add_argument("--out", required=True, type=Path)
    demo_rabies_geography.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_rabies_geography)
    demo_rabies_host_geography = demo_subparsers.add_parser(
        "rabies-cross-host-geography-panel",
        help="Materialize the packaged rabies integrated dataset and rerun the full sequence-to-tree, host, and geography workflow outputs.",
    )
    demo_rabies_host_geography.add_argument("--out", required=True, type=Path)
    demo_rabies_host_geography.add_argument(
        "--config",
        type=Path,
        help="Optional workflow config JSON. Defaults to the packaged dataset config.",
    )
    demo_rabies_host_geography.add_argument("--mafft-executable", type=str)
    demo_rabies_host_geography.add_argument("--trimal-executable", type=str)
    demo_rabies_host_geography.add_argument("--iqtree-executable", type=str)
    demo_rabies_host_geography.add_argument("--fasttree-executable", type=str)
    demo_rabies_host_geography.add_argument("--iqtree-seed", type=int, default=1)
    demo_rabies_host_geography.add_argument("--iqtree-threads", type=int, default=1)
    demo_rabies_host_geography.add_argument(
        "--bootstrap-replicates", type=int, default=1000
    )
    demo_rabies_host_geography.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_rabies_host_geography)
    demo_rabies_method_sensitivity = demo_subparsers.add_parser(
        "rabies-method-sensitivity-panel",
        help="Materialize the packaged rabies method-sensitivity dataset and rerun the governed preprocessing and engine-comparison workflow outputs.",
    )
    demo_rabies_method_sensitivity.add_argument("--out", required=True, type=Path)
    demo_rabies_method_sensitivity.add_argument("--mafft-executable", type=str)
    demo_rabies_method_sensitivity.add_argument("--trimal-executable", type=str)
    demo_rabies_method_sensitivity.add_argument("--iqtree-executable", type=str)
    demo_rabies_method_sensitivity.add_argument("--fasttree-executable", type=str)
    demo_rabies_method_sensitivity.add_argument("--iqtree-seed", type=int, default=1)
    demo_rabies_method_sensitivity.add_argument("--iqtree-threads", type=int, default=1)
    demo_rabies_method_sensitivity.add_argument(
        "--bootstrap-replicates", type=int, default=1000
    )
    demo_rabies_method_sensitivity.add_argument(
        "--parallel-workers",
        type=int,
        default=None,
        help="Number of isolated variant workers to run in parallel. Defaults to the packaged workflow config.",
    )
    demo_rabies_method_sensitivity.add_argument(
        "--variant-id",
        action="append",
        dest="variant_ids",
        help="Restrict the workflow to one declared variant id. Repeat to preserve a specific subset order.",
    )
    demo_rabies_method_sensitivity.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_rabies_method_sensitivity)
    demo_catarrhine_mitogenome = demo_subparsers.add_parser(
        "catarrhine-mitogenome-five-locus-panel",
        help="Materialize the packaged catarrhine multi-locus dataset and rerun the governed concatenation and partitioned inference outputs.",
    )
    demo_catarrhine_mitogenome.add_argument("--out", required=True, type=Path)
    demo_catarrhine_mitogenome.add_argument("--iqtree-executable", type=str)
    demo_catarrhine_mitogenome.add_argument("--iqtree-seed", type=int, default=1)
    demo_catarrhine_mitogenome.add_argument("--iqtree-threads", type=int, default=1)
    demo_catarrhine_mitogenome.add_argument(
        "--bootstrap-replicates", type=int, default=1000
    )
    demo_catarrhine_mitogenome.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_catarrhine_mitogenome)
    demo_catarrhine_stress = demo_subparsers.add_parser(
        "catarrhine-data-quality-stress-panel",
        help="Materialize the packaged catarrhine dirty-data stress dataset and rerun the governed audit and cleanup outputs.",
    )
    demo_catarrhine_stress.add_argument("--out", required=True, type=Path)
    demo_catarrhine_stress.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_catarrhine_stress)
    demo_continuous_mode_recovery = demo_subparsers.add_parser(
        "continuous-mode-recovery-panel",
        help="Materialize the packaged continuous-trait recovery dataset and rerun the governed simulation-recovery outputs.",
    )
    demo_continuous_mode_recovery.add_argument("--out", required=True, type=Path)
    demo_continuous_mode_recovery.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_continuous_mode_recovery)
    demo_discrete_mode_recovery = demo_subparsers.add_parser(
        "discrete-mode-recovery-panel",
        help="Materialize the packaged discrete-trait recovery dataset and rerun the governed simulation-recovery outputs.",
    )
    demo_discrete_mode_recovery.add_argument("--out", required=True, type=Path)
    demo_discrete_mode_recovery.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_discrete_mode_recovery)
    demo_known_answer = demo_subparsers.add_parser(
        "known-answer-reference-panel",
        help="Materialize the packaged known-answer simulation dataset and rerun the governed recovery outputs.",
    )
    demo_known_answer.add_argument("--out", required=True, type=Path)
    demo_known_answer.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_known_answer)


def run_demo_command(args: Any) -> int:
    if args.demo_command == "run":
        result = run_capability_demo(args.out)
        outputs = _finalize_outputs(
            args,
            command="demo",
            inputs=[],
            outputs=[
                result.tree_report,
                result.dataset_report,
                result.phylo_inputs_report,
                result.comparison_report,
                result.capability_summary,
            ],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="demo",
                    inputs=[],
                    outputs=outputs,
                    metrics={"artifact_count": 5},
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_root)
        return 0

    if args.demo_command == "primate-comparative":
        result = run_primate_comparative_demo(args.out)
        outputs = _finalize_outputs(
            args,
            command="demo",
            inputs=[],
            outputs=[
                result.dataset_export.readme_path,
                result.dataset_export.tree_path,
                result.dataset_export.traits_path,
                result.workflow_bundle.summary_path,
                result.workflow_bundle.pgls_lambda_profile_path,
                result.workflow_bundle.brownian_summary_path,
                result.workflow_bundle.ou_summary_path,
                result.workflow_bundle.signal_summary_path,
                result.workflow_bundle.signal_permutations_path,
                result.workflow_bundle.continuous_ancestral_summary_path,
                result.workflow_bundle.continuous_ancestral_uncertainty_path,
                result.workflow_bundle.discrete_ancestral_summary_path,
                result.workflow_bundle.discrete_ancestral_probability_path,
                result.overview_path,
            ],
        )
        if args.json:
            expected_output_count = len(
                list(result.dataset_export.expected_output_root.glob("*"))
            )
            _print_result(
                build_command_result(
                    command="demo",
                    inputs=[],
                    outputs=outputs,
                    metrics={
                        "artifact_count": len(outputs),
                        "dataset_taxon_count": result.dataset.taxon_count,
                        "reference_output_count": expected_output_count,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_root)
        return 0

    if args.demo_command == "avian-reproductive-traits":
        result = run_avian_reproductive_trait_demo(args.out)
        outputs = _finalize_outputs(
            args,
            command="demo",
            inputs=[],
            outputs=[
                result.dataset_export.readme_path,
                result.dataset_export.tree_path,
                result.dataset_export.traits_path,
                result.workflow_bundle.summary_path,
                result.workflow_bundle.pgls_lambda_profile_path,
                result.workflow_bundle.brownian_summary_path,
                result.workflow_bundle.ou_summary_path,
                result.workflow_bundle.signal_summary_path,
                result.workflow_bundle.signal_permutations_path,
                result.workflow_bundle.continuous_ancestral_summary_path,
                result.workflow_bundle.continuous_ancestral_uncertainty_path,
                result.workflow_bundle.discrete_ancestral_summary_path,
                result.workflow_bundle.discrete_ancestral_probability_path,
                result.workflow_bundle.clade_summary_path,
                result.workflow_bundle.clade_rows_path,
                result.overview_path,
            ],
        )
        if args.json:
            expected_output_count = len(
                list(result.dataset_export.expected_output_root.glob("*"))
            )
            _print_result(
                build_command_result(
                    command="demo",
                    inputs=[],
                    outputs=outputs,
                    metrics={
                        "artifact_count": len(outputs),
                        "dataset_taxon_count": result.dataset.taxon_count,
                        "reference_output_count": expected_output_count,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_root)
        return 0

    if args.demo_command == "central-european-seashore-flora":
        result = run_central_european_seashore_flora_demo(args.out)
        outputs = _finalize_outputs(
            args,
            command="demo",
            inputs=[],
            outputs=[
                result.dataset_export.readme_path,
                result.dataset_export.tree_path,
                result.dataset_export.traits_path,
                result.workflow_bundle.summary_path,
                result.workflow_bundle.pgls_lambda_profile_path,
                result.workflow_bundle.brownian_summary_path,
                result.workflow_bundle.ou_summary_path,
                result.workflow_bundle.signal_summary_path,
                result.workflow_bundle.signal_permutations_path,
                result.workflow_bundle.continuous_ancestral_summary_path,
                result.workflow_bundle.continuous_ancestral_uncertainty_path,
                result.workflow_bundle.discrete_ancestral_summary_path,
                result.workflow_bundle.discrete_ancestral_probability_path,
                result.workflow_bundle.clade_summary_path,
                result.workflow_bundle.clade_rows_path,
                result.overview_path,
            ],
        )
        if args.json:
            expected_output_count = len(
                list(result.dataset_export.expected_output_root.glob("*"))
            )
            _print_result(
                build_command_result(
                    command="demo",
                    inputs=[],
                    outputs=outputs,
                    metrics={
                        "artifact_count": len(outputs),
                        "dataset_taxon_count": result.dataset.taxon_count,
                        "reference_output_count": expected_output_count,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_root)
        return 0

    if args.demo_command == "influenza-a-ha-reference-panel":
        result = run_influenza_a_ha_reference_demo(
            args.out,
            mafft_executable=args.mafft_executable or "mafft",
            trimal_executable=args.trimal_executable or "trimal",
            iqtree_executable=args.iqtree_executable or "iqtree2",
            iqtree_seed=args.iqtree_seed,
            iqtree_threads=args.iqtree_threads,
            bootstrap_replicates=args.bootstrap_replicates,
        )
        outputs = _finalize_outputs(
            args,
            command="demo",
            inputs=[],
            outputs=[
                result.dataset_export.readme_path,
                result.dataset_export.sequences_path,
                result.workflow_bundle.summary_path,
                result.workflow_bundle.alignment_path,
                result.workflow_bundle.trimmed_alignment_path,
                result.workflow_bundle.tree_path,
                result.workflow_bundle.model_table_path,
                result.workflow_bundle.support_table_path,
                result.workflow_bundle.log_path,
                result.workflow_bundle.manifest_path,
                result.overview_path,
            ],
        )
        if args.json:
            expected_output_count = len(
                list(result.dataset_export.expected_output_root.glob("*"))
            )
            _print_result(
                build_command_result(
                    command="demo",
                    inputs=[],
                    outputs=outputs,
                    metrics={
                        "artifact_count": len(outputs),
                        "sequence_count": result.dataset.sequence_count,
                        "sequence_type": result.dataset.sequence_type,
                        "selected_model": result.workflow_bundle.selected_model,
                        "minimum_support": result.workflow_bundle.minimum_support,
                        "maximum_support": result.workflow_bundle.maximum_support,
                        "weakly_supported_clade_count": (
                            result.workflow_bundle.weakly_supported_clade_count
                        ),
                        "reference_output_count": expected_output_count,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_root)
        return 0

    if args.demo_command == "gnathostome-ortholog-protein-benchmark":
        result = run_gnathostome_ortholog_protein_benchmark_demo(
            args.out,
            mafft_executable=args.mafft_executable or "mafft",
            trimal_executable=args.trimal_executable or "trimal",
            iqtree_executable=args.iqtree_executable or "iqtree2",
            iqtree_seed=args.iqtree_seed,
            iqtree_threads=args.iqtree_threads,
            bootstrap_replicates=args.bootstrap_replicates,
        )
        outputs = _finalize_outputs(
            args,
            command="demo",
            inputs=[],
            outputs=[
                result.dataset_export.readme_path,
                result.dataset_export.sequences_path,
                result.workflow_bundle.summary_path,
                result.workflow_bundle.assumptions_path,
                result.workflow_bundle.alignment_path,
                result.workflow_bundle.trimmed_alignment_path,
                result.workflow_bundle.tree_path,
                result.workflow_bundle.model_table_path,
                result.workflow_bundle.support_table_path,
                result.workflow_bundle.log_path,
                result.workflow_bundle.manifest_path,
                result.overview_path,
            ],
        )
        if args.json:
            expected_output_count = len(
                list(result.dataset_export.expected_output_root.glob("*"))
            )
            _print_result(
                build_command_result(
                    command="demo",
                    inputs=[],
                    outputs=outputs,
                    metrics={
                        "artifact_count": len(outputs),
                        "sequence_count": result.dataset.sequence_count,
                        "sequence_type": result.dataset.sequence_type,
                        "selected_model": result.workflow_bundle.selected_model,
                        "alignment_length": result.workflow_bundle.alignment_length,
                        "trimmed_alignment_length": (
                            result.workflow_bundle.trimmed_alignment_length
                        ),
                        "minimum_support": result.workflow_bundle.minimum_support,
                        "maximum_support": result.workflow_bundle.maximum_support,
                        "weakly_supported_clade_count": (
                            result.workflow_bundle.weakly_supported_clade_count
                        ),
                        "state_space": "amino-acid",
                        "model_selection_scope": "protein-models-only",
                        "reference_output_count": expected_output_count,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_root)
        return 0

    if args.demo_command == "pleistocene-bear-cytb-fragments":
        result = run_pleistocene_bear_cytb_fragment_demo(
            args.out,
            mafft_executable=args.mafft_executable or "mafft",
            trimal_executable=args.trimal_executable or "trimal",
            iqtree_executable=args.iqtree_executable or "iqtree2",
            iqtree_seed=args.iqtree_seed,
            iqtree_threads=args.iqtree_threads,
            bootstrap_replicates=args.bootstrap_replicates,
        )
        outputs = _finalize_outputs(
            args,
            command="demo",
            inputs=[],
            outputs=[
                result.dataset_export.readme_path,
                result.dataset_export.sequences_path,
                result.workflow_bundle.summary_path,
                result.workflow_bundle.missingness_effects_path,
                result.workflow_bundle.alignment_path,
                result.workflow_bundle.trimmed_alignment_path,
                result.workflow_bundle.cleaned_alignment_path,
                result.workflow_bundle.tree_path,
                result.workflow_bundle.model_table_path,
                result.workflow_bundle.support_table_path,
                result.overview_path,
            ],
        )
        if args.json:
            expected_output_count = len(
                list(result.dataset_export.expected_output_root.glob("*"))
            )
            _print_result(
                build_command_result(
                    command="demo",
                    inputs=[],
                    outputs=outputs,
                    metrics={
                        "artifact_count": len(outputs),
                        "sequence_count": result.dataset.sequence_count,
                        "degraded_sequence_count": len(
                            result.dataset.degraded_sequence_ids
                        ),
                        "selected_model": result.workflow_bundle.selected_model,
                        "minimum_support": result.workflow_bundle.minimum_support,
                        "maximum_support": result.workflow_bundle.maximum_support,
                        "removed_column_count": (
                            result.workflow_bundle.removed_column_count
                        ),
                        "cleaned_missing_data_fraction": (
                            result.workflow_bundle.cleaned_missing_data_fraction
                        ),
                        "reference_output_count": expected_output_count,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_root)
        return 0

    if args.demo_command == "rabies-cross-host-panel":
        result = run_rabies_cross_host_panel_demo(args.out)
        outputs = _finalize_outputs(
            args,
            command="demo",
            inputs=[],
            outputs=[
                result.dataset_export.readme_path,
                result.dataset_export.sequences_path,
                result.dataset_export.tree_path,
                result.dataset_export.hosts_path,
                result.workflow_bundle.workflow_summary_path,
                result.workflow_bundle.host_switch_summary_path,
                result.workflow_bundle.host_state_nodes_path,
                result.workflow_bundle.host_switch_branches_path,
                result.workflow_bundle.host_switch_counts_path,
                result.workflow_bundle.host_switch_fits_path,
                result.workflow_bundle.host_switch_unsupported_path,
                result.workflow_bundle.host_switch_exclusions_path,
                result.overview_path,
            ],
        )
        if args.json:
            expected_output_count = len(
                list(result.dataset_export.expected_output_root.glob("*"))
            )
            _print_result(
                build_command_result(
                    command="demo",
                    inputs=[],
                    outputs=outputs,
                    metrics={
                        "artifact_count": len(outputs),
                        "taxon_count": result.dataset.taxon_count,
                        "workflow_trait": result.dataset.workflow_trait,
                        "observed_host_group_count": (
                            result.dataset.observed_host_group_count
                        ),
                        "analysis_constraint_mode": (
                            result.workflow_bundle.analysis_constraint_mode
                        ),
                        "root_host": result.workflow_bundle.root_host,
                        "root_confidence": result.workflow_bundle.root_confidence,
                        "host_switch_count": (
                            result.workflow_bundle.host_switch_count
                        ),
                        "certain_host_switch_count": (
                            result.workflow_bundle.certain_host_switch_count
                        ),
                        "uncertain_host_switch_count": (
                            result.workflow_bundle.uncertain_host_switch_count
                        ),
                        "reference_output_count": expected_output_count,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_root)
        return 0

    if args.demo_command == "rabies-geographic-transition-panel":
        result = run_rabies_geographic_transition_panel_demo(args.out)
        outputs = _finalize_outputs(
            args,
            command="demo",
            inputs=[],
            outputs=[
                result.dataset_export.readme_path,
                result.dataset_export.sequences_path,
                result.dataset_export.tree_path,
                result.dataset_export.regions_path,
                result.workflow_bundle.workflow_summary_path,
                result.workflow_bundle.geographic_state_summary_path,
                result.workflow_bundle.geographic_region_probability_path,
                result.workflow_bundle.geographic_transition_rate_path,
                result.workflow_bundle.geographic_transition_event_path,
                result.workflow_bundle.geographic_state_exclusion_path,
                result.workflow_bundle.geographic_migration_summary_path,
                result.workflow_bundle.geographic_migration_event_path,
                result.workflow_bundle.geographic_migration_exclusion_path,
                result.overview_path,
            ],
        )
        if args.json:
            expected_output_count = len(
                list(result.dataset_export.expected_output_root.glob("*"))
            )
            _print_result(
                build_command_result(
                    command="demo",
                    inputs=[],
                    outputs=outputs,
                    metrics={
                        "artifact_count": len(outputs),
                        "taxon_count": result.dataset.taxon_count,
                        "workflow_trait": result.dataset.workflow_trait,
                        "observed_region_group_count": (
                            result.dataset.observed_region_group_count
                        ),
                        "root_region": result.workflow_bundle.root_region,
                        "root_region_probability": (
                            result.workflow_bundle.root_region_probability
                        ),
                        "changed_branch_count": (
                            result.workflow_bundle.changed_branch_count
                        ),
                        "strongly_supported_transition_count": (
                            result.workflow_bundle.strongly_supported_transition_count
                        ),
                        "migration_event_count": (
                            result.workflow_bundle.migration_event_count
                        ),
                        "strongly_supported_migration_event_count": (
                            result.workflow_bundle.strongly_supported_migration_event_count
                        ),
                        "reference_output_count": expected_output_count,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_root)
        return 0

    if args.demo_command == "rabies-cross-host-geography-panel":
        result = run_rabies_cross_host_geography_panel_demo(
            args.out,
            config_path=args.config,
            mafft_executable=args.mafft_executable or "mafft",
            trimal_executable=args.trimal_executable or "trimal",
            iqtree_executable=args.iqtree_executable or "iqtree2",
            fasttree_executable=args.fasttree_executable or "FastTree",
            iqtree_seed=args.iqtree_seed,
            iqtree_threads=args.iqtree_threads,
            bootstrap_replicates=args.bootstrap_replicates,
        )
        outputs = _finalize_outputs(
            args,
            command="demo",
            inputs=[],
            outputs=[
                result.dataset_export.readme_path,
                result.dataset_export.workflow_config_path,
                result.dataset_export.sequences_path,
                result.dataset_export.metadata_path,
                result.dataset_export.centroids_path,
                result.dataset_export.accession_table_path,
                result.workflow_bundle.workflow_summary_path,
                result.workflow_bundle.resource_observations_path,
                result.workflow_bundle.config_audit_path,
                result.workflow_bundle.resolved_config_path,
                result.workflow_bundle.input_validation_path,
                result.workflow_bundle.alignment_quality_path,
                result.workflow_bundle.alignment_sequence_ranking_path,
                result.workflow_bundle.alignment_path,
                result.workflow_bundle.trimmed_alignment_path,
                result.workflow_bundle.tree_path,
                result.workflow_bundle.rooting_report_path,
                result.workflow_bundle.model_table_path,
                result.workflow_bundle.support_table_path,
                result.workflow_bundle.clade_table_path,
                result.workflow_bundle.bootstrap_summary_path,
                result.workflow_bundle.bootstrap_tree_comparison_summary_path,
                result.workflow_bundle.host_switch_summary_path,
                result.workflow_bundle.host_switch_counts_path,
                result.workflow_bundle.biogeography_report_path,
                result.workflow_bundle.biogeography_tree_figure_path,
                result.workflow_bundle.biogeography_map_path,
                result.workflow_bundle.comparative_report_path,
                result.workflow_bundle.comparative_summary_path,
                result.workflow_bundle.conclusion_stability_summary_path,
                result.workflow_bundle.key_clade_stability_path,
                result.workflow_bundle.support_value_stability_path,
                result.workflow_bundle.ancestral_state_stability_path,
                result.workflow_bundle.comparative_coefficient_stability_path,
                result.workflow_bundle.conclusion_stability_report_path,
                result.workflow_bundle.scientific_findings_path,
                result.workflow_bundle.final_report_path,
                result.workflow_bundle.final_manifest_path,
                result.overview_path,
                result.overview_html_path,
                result.artifact_inventory_path,
                result.reproducibility_checklist_path,
                result.package_manifest_path,
            ],
        )
        if args.json:
            expected_output_count = len(
                [
                    path
                    for path in result.dataset_export.expected_output_root.rglob("*")
                    if path.is_file()
                ]
            )
            _print_result(
                build_command_result(
                    command="demo",
                    inputs=[],
                    outputs=outputs,
                    metrics={
                        "artifact_count": len(outputs),
                        "sequence_count": result.dataset.sequence_count,
                        "config_path": str(result.dataset_export.workflow_config_path),
                        "biological_question": (
                            "Do the host-associated rabies lineages in this compact panel occupy one distinct geographic regime while retaining one coherent phylogenetic signal?"
                        ),
                        "short_answer": (
                            "The rooted panel remains anchored in bat and north_asia, and `host_group[canid]` shows a nominally supported positive longitude association under the selected comparative model, but the inference remains cautionary because the panel is intentionally compact."
                        ),
                        "host_trait": result.dataset.host_trait,
                        "geography_trait": result.dataset.geography_trait,
                        "selected_model": result.workflow_bundle.selected_model,
                        "aligned_quality_score": (
                            result.workflow_bundle.aligned_quality_score
                        ),
                        "trimmed_quality_score": (
                            result.workflow_bundle.trimmed_quality_score
                        ),
                        "minimum_support": result.workflow_bundle.minimum_support,
                        "maximum_support": result.workflow_bundle.maximum_support,
                        "root_host": result.workflow_bundle.root_host,
                        "root_region": result.workflow_bundle.root_region,
                        "host_switch_count": result.workflow_bundle.host_switch_count,
                        "migration_event_count": (
                            result.workflow_bundle.migration_event_count
                        ),
                        "clade_row_count": result.workflow_bundle.clade_row_count,
                        "bootstrap_tree_count": (
                            result.workflow_bundle.bootstrap_tree_count
                        ),
                        "timeout_seconds": result.workflow_bundle.timeout_seconds,
                        "max_bootstrap_tree_count": (
                            result.workflow_bundle.max_bootstrap_tree_count
                        ),
                        "max_report_table_rows": (
                            result.workflow_bundle.max_report_table_rows
                        ),
                        "budget_warning_count": (
                            result.workflow_bundle.budget_warning_count
                        ),
                        "bootstrap_review_runtime_seconds": (
                            result.workflow_bundle.bootstrap_review_runtime_seconds
                        ),
                        "bootstrap_review_peak_memory_bytes": (
                            result.workflow_bundle.bootstrap_review_peak_memory_bytes
                        ),
                        "bootstrap_consensus_rooted_rf_distance": (
                            result.workflow_bundle.bootstrap_consensus_rooted_rf_distance
                        ),
                        "comparative_formula": (
                            result.workflow_bundle.comparative_formula
                        ),
                        "comparative_selected_model": (
                            result.workflow_bundle.comparative_selected_model
                        ),
                        "conclusion_stable_count": (
                            result.workflow_bundle.conclusion_stable_count
                        ),
                        "conclusion_weak_count": (
                            result.workflow_bundle.conclusion_weak_count
                        ),
                        "conclusion_unstable_count": (
                            result.workflow_bundle.conclusion_unstable_count
                        ),
                        "config_check_count": result.workflow_bundle.config_check_count,
                        "scientific_finding_count": (
                            result.workflow_bundle.scientific_finding_count
                        ),
                        "package_artifact_count": (
                            sum(
                                1
                                for _ in result.artifact_inventory_path.open(
                                    "r", encoding="utf-8"
                                )
                            )
                            - 1
                        ),
                        "package_checklist_item_count": (
                            sum(
                                1
                                for _ in result.reproducibility_checklist_path.open(
                                    "r", encoding="utf-8"
                                )
                            )
                            - 1
                        ),
                        "reference_output_count": expected_output_count,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_root)
        return 0

    if args.demo_command == "rabies-method-sensitivity-panel":
        result = run_rabies_method_sensitivity_panel_demo(
            args.out,
            mafft_executable=args.mafft_executable or "mafft",
            trimal_executable=args.trimal_executable or "trimal",
            iqtree_executable=args.iqtree_executable or "iqtree2",
            fasttree_executable=args.fasttree_executable or "FastTree",
            iqtree_seed=args.iqtree_seed,
            iqtree_threads=args.iqtree_threads,
            bootstrap_replicates=args.bootstrap_replicates,
            parallel_workers=args.parallel_workers,
            variant_ids=(
                tuple(args.variant_ids) if getattr(args, "variant_ids", None) else None
            ),
        )
        outputs = _finalize_outputs(
            args,
            command="demo",
            inputs=[],
            outputs=[
                result.dataset_export.readme_path,
                result.dataset_export.config_path,
                result.dataset_export.sequences_path,
                result.dataset_export.metadata_path,
                result.workflow_bundle.workflow_summary_path,
                result.workflow_bundle.variant_summary_path,
                result.workflow_bundle.parallel_summary_path,
                result.workflow_bundle.execution_record_path,
                result.workflow_bundle.preprocessing_comparison_path,
                result.workflow_bundle.stable_clades_path,
                result.workflow_bundle.changed_clades_path,
                result.workflow_bundle.conclusion_summary_path,
                result.workflow_bundle.config_path,
                result.workflow_bundle.manifest_path,
                result.workflow_bundle.report_manifest_path,
                result.workflow_bundle.slurm_job_plan_path,
                result.workflow_bundle.slurm_assumptions_path,
                result.workflow_bundle.slurm_summary_path,
                result.workflow_bundle.slurm_array_partitions_path,
                result.workflow_bundle.slurm_array_members_path,
                result.workflow_bundle.slurm_array_strategy_path,
                result.workflow_bundle.slurm_job_evidence_index_path,
                result.workflow_bundle.slurm_job_evidence_summary_path,
                result.workflow_bundle.slurm_storage_categories_path,
                result.workflow_bundle.slurm_storage_variants_path,
                result.workflow_bundle.slurm_storage_summary_path,
                result.workflow_bundle.slurm_storage_report_path,
                result.workflow_bundle.slurm_output_explosion_checks_path,
                result.workflow_bundle.slurm_output_explosion_variants_path,
                result.workflow_bundle.slurm_output_explosion_summary_path,
                result.workflow_bundle.slurm_output_explosion_report_path,
                result.workflow_bundle.slurm_tree_retention_checks_path,
                result.workflow_bundle.slurm_tree_retention_files_path,
                result.workflow_bundle.slurm_tree_retention_summary_path,
                result.workflow_bundle.slurm_tree_retention_report_path,
                result.workflow_bundle.slurm_merge_checks_path,
                result.workflow_bundle.slurm_merge_variants_path,
                result.workflow_bundle.slurm_merge_summary_path,
                result.workflow_bundle.slurm_merge_report_path,
                result.workflow_bundle.slurm_output_freshness_path,
                result.workflow_bundle.slurm_output_freshness_checks_path,
                result.workflow_bundle.slurm_output_freshness_summary_path,
                result.workflow_bundle.slurm_job_status_path,
                result.workflow_bundle.slurm_partition_status_path,
                result.workflow_bundle.slurm_workflow_status_path,
                result.workflow_bundle.slurm_failure_recovery_jobs_path,
                result.workflow_bundle.slurm_failure_recovery_partitions_path,
                result.workflow_bundle.slurm_failure_recovery_summary_path,
                result.workflow_bundle.slurm_failure_recovery_report_path,
                result.workflow_bundle.reproducibility_checks_path,
                result.workflow_bundle.reproducibility_variant_audit_path,
                result.workflow_bundle.reproducibility_audit_path,
                result.workflow_bundle.report_path,
                result.overview_path,
            ],
        )
        if args.json:
            expected_output_count = len(
                [
                    path
                    for path in result.dataset_export.expected_output_root.rglob("*")
                    if path.is_file()
                ]
            )
            _print_result(
                build_command_result(
                    command="demo",
                    inputs=[],
                    outputs=outputs,
                    metrics={
                        "artifact_count": len(outputs),
                        "taxon_count": result.dataset.taxon_count,
                        "variant_count": result.workflow_bundle.variant_count,
                        "parallel_workers": result.workflow_bundle.parallel_workers,
                        "execution_mode": result.workflow_bundle.execution_mode,
                        "stable_clade_count": (
                            result.workflow_bundle.stable_clade_count
                        ),
                        "changed_clade_count": (
                            result.workflow_bundle.changed_clade_count
                        ),
                        "preprocessing_change_pair_count": (
                            result.workflow_bundle.preprocessing_change_pair_count
                        ),
                        "rooted_engine_change_variant_count": (
                            result.workflow_bundle.rooted_engine_change_variant_count
                        ),
                        "serious_conflict_variant_count": (
                            result.workflow_bundle.serious_conflict_variant_count
                        ),
                        "report_linked_artifact_count": (
                            result.workflow_bundle.report_linked_artifact_count
                        ),
                        "report_html_size_bytes": (
                            result.workflow_bundle.report_html_size_bytes
                        ),
                        "report_linked_artifact_bytes": (
                            result.workflow_bundle.report_linked_artifact_bytes
                        ),
                        "report_total_output_bytes": (
                            result.workflow_bundle.report_total_output_bytes
                        ),
                        "slurm_job_count": result.workflow_bundle.slurm_job_count,
                        "slurm_total_estimated_core_hours": (
                            result.workflow_bundle.slurm_total_estimated_core_hours
                        ),
                        "slurm_maximum_estimated_memory_mib": (
                            result.workflow_bundle.slurm_maximum_estimated_memory_mib
                        ),
                        "slurm_maximum_estimated_wallclock_minutes": (
                            result.workflow_bundle.slurm_maximum_estimated_wallclock_minutes
                        ),
                        "slurm_total_estimated_scratch_mib": (
                            result.workflow_bundle.slurm_total_estimated_scratch_mib
                        ),
                        "slurm_total_estimated_output_mib": (
                            result.workflow_bundle.slurm_total_estimated_output_mib
                        ),
                        "slurm_array_partition_count": (
                            result.workflow_bundle.slurm_array_partition_count
                        ),
                        "slurm_array_script_count": (
                            result.workflow_bundle.slurm_array_script_count
                        ),
                        "slurm_array_largest_partition_size": (
                            result.workflow_bundle.slurm_array_largest_partition_size
                        ),
                        "slurm_job_evidence_file_count": (
                            result.workflow_bundle.slurm_job_evidence_file_count
                        ),
                        "slurm_job_evidence_total_runtime_seconds": (
                            result.workflow_bundle.slurm_job_evidence_total_runtime_seconds
                        ),
                        "slurm_job_evidence_total_output_byte_count": (
                            result.workflow_bundle.slurm_job_evidence_total_output_byte_count
                        ),
                        "slurm_storage_total_estimated_mib": (
                            result.workflow_bundle.slurm_storage_total_estimated_mib
                        ),
                        "slurm_storage_output_byte_count": (
                            result.workflow_bundle.slurm_storage_output_byte_count
                        ),
                        "slurm_storage_log_byte_count": (
                            result.workflow_bundle.slurm_storage_log_byte_count
                        ),
                        "slurm_storage_tree_byte_count": (
                            result.workflow_bundle.slurm_storage_tree_byte_count
                        ),
                        "slurm_storage_posterior_sample_byte_count": (
                            result.workflow_bundle.slurm_storage_posterior_sample_byte_count
                        ),
                        "slurm_storage_report_byte_count": (
                            result.workflow_bundle.slurm_storage_report_byte_count
                        ),
                        "slurm_storage_largest_variant_id": (
                            result.workflow_bundle.slurm_storage_largest_variant_id
                        ),
                        "slurm_output_explosion_status": (
                            result.workflow_bundle.slurm_output_explosion_status
                        ),
                        "slurm_output_explosion_global_issue_count": (
                            result.workflow_bundle.slurm_output_explosion_global_issue_count
                        ),
                        "slurm_output_explosion_warning_variant_count": (
                            result.workflow_bundle.slurm_output_explosion_warning_variant_count
                        ),
                        "slurm_output_explosion_high_risk_variant_count": (
                            result.workflow_bundle.slurm_output_explosion_high_risk_variant_count
                        ),
                        "slurm_tree_retention_status": (
                            result.workflow_bundle.slurm_tree_retention_status
                        ),
                        "slurm_tree_set_file_count": (
                            result.workflow_bundle.slurm_tree_set_file_count
                        ),
                        "slurm_tree_posterior_sample_file_count": (
                            result.workflow_bundle.slurm_tree_posterior_sample_file_count
                        ),
                        "slurm_tree_thinning_recommended_file_count": (
                            result.workflow_bundle.slurm_tree_thinning_recommended_file_count
                        ),
                        "slurm_tree_thinning_required_file_count": (
                            result.workflow_bundle.slurm_tree_thinning_required_file_count
                        ),
                        "slurm_tree_compression_recommended_file_count": (
                            result.workflow_bundle.slurm_tree_compression_recommended_file_count
                        ),
                        "slurm_tree_compression_required_file_count": (
                            result.workflow_bundle.slurm_tree_compression_required_file_count
                        ),
                        "slurm_merge_status": (
                            result.workflow_bundle.slurm_merge_status
                        ),
                        "slurm_merge_ready": (
                            result.workflow_bundle.slurm_merge_ready
                        ),
                        "slurm_mergeable_variant_count": (
                            result.workflow_bundle.slurm_mergeable_variant_count
                        ),
                        "slurm_merge_failed_check_count": (
                            result.workflow_bundle.slurm_merge_failed_check_count
                        ),
                        "slurm_output_freshness_check_count": (
                            result.workflow_bundle.slurm_output_freshness_check_count
                        ),
                        "slurm_output_freshness_failed_check_count": (
                            result.workflow_bundle.slurm_output_freshness_failed_check_count
                        ),
                        "slurm_fresh_output_job_count": (
                            result.workflow_bundle.slurm_fresh_output_job_count
                        ),
                        "slurm_stale_output_job_count": (
                            result.workflow_bundle.slurm_stale_output_job_count
                        ),
                        "slurm_completed_job_count": (
                            result.workflow_bundle.slurm_completed_job_count
                        ),
                        "slurm_failed_job_count": (
                            result.workflow_bundle.slurm_failed_job_count
                        ),
                        "slurm_pending_job_count": (
                            result.workflow_bundle.slurm_pending_job_count
                        ),
                        "slurm_stale_job_count": (
                            result.workflow_bundle.slurm_stale_job_count
                        ),
                        "slurm_failure_recovery_status": (
                            result.workflow_bundle.slurm_failure_recovery_status
                        ),
                        "slurm_failure_recovery_rerunnable_job_count": (
                            result.workflow_bundle.slurm_failure_recovery_rerunnable_job_count
                        ),
                        "slurm_failure_recovery_blocked_job_count": (
                            result.workflow_bundle.slurm_failure_recovery_blocked_job_count
                        ),
                        "slurm_failure_recovery_partition_count": (
                            result.workflow_bundle.slurm_failure_recovery_partition_count
                        ),
                        "reproducibility_passed": (
                            result.workflow_bundle.reproducibility_passed
                        ),
                        "reproducibility_check_count": (
                            result.workflow_bundle.reproducibility_check_count
                        ),
                        "reproducibility_failed_check_count": (
                            result.workflow_bundle.reproducibility_failed_check_count
                        ),
                        "reproducibility_failed_variant_count": (
                            result.workflow_bundle.reproducibility_failed_variant_count
                        ),
                        "reference_output_count": expected_output_count,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_root)
        return 0

    if args.demo_command == "catarrhine-mitogenome-five-locus-panel":
        result = run_catarrhine_mitogenome_five_locus_panel_demo(
            args.out,
            iqtree_executable=args.iqtree_executable or "iqtree2",
            iqtree_seed=args.iqtree_seed,
            iqtree_threads=args.iqtree_threads,
            bootstrap_replicates=args.bootstrap_replicates,
        )
        outputs = _finalize_outputs(
            args,
            command="demo",
            inputs=[],
            outputs=[
                result.dataset_export.readme_path,
                result.dataset_export.taxa_path,
                *sorted(result.dataset_export.locus_alignment_root.glob("*.fasta")),
                result.workflow_bundle.workflow_summary_path,
                result.workflow_bundle.supermatrix_path,
                result.workflow_bundle.partitions_path,
                result.workflow_bundle.occupancy_taxa_path,
                result.workflow_bundle.occupancy_loci_path,
                result.workflow_bundle.occupancy_matrix_path,
                result.workflow_bundle.partition_summary_path,
                result.workflow_bundle.model_candidates_path,
                result.workflow_bundle.support_tree_path,
                result.workflow_bundle.support_table_path,
                result.overview_path,
            ],
        )
        if args.json:
            expected_output_count = len(
                list(result.dataset_export.expected_output_root.glob("*"))
            )
            _print_result(
                build_command_result(
                    command="demo",
                    inputs=[],
                    outputs=outputs,
                    metrics={
                        "artifact_count": len(outputs),
                        "taxon_count": result.dataset.taxon_count,
                        "locus_count": result.dataset.locus_count,
                        "alignment_length": result.workflow_bundle.alignment_length,
                        "partition_count": result.workflow_bundle.partition_count,
                        "selected_model": result.workflow_bundle.selected_model,
                        "minimum_support": result.workflow_bundle.minimum_support,
                        "maximum_support": result.workflow_bundle.maximum_support,
                        "weakly_supported_clade_count": (
                            result.workflow_bundle.weakly_supported_clade_count
                        ),
                        "reference_output_count": expected_output_count,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_root)
        return 0

    if args.demo_command == "catarrhine-data-quality-stress-panel":
        result = run_catarrhine_data_quality_stress_panel_demo(args.out)
        outputs = _finalize_outputs(
            args,
            command="demo",
            inputs=[],
            outputs=[
                result.dataset_export.readme_path,
                result.dataset_export.raw_alignment_path,
                result.dataset_export.raw_sequence_input_path,
                result.dataset_export.raw_coding_sequences_path,
                result.dataset_export.raw_tree_path,
                result.dataset_export.raw_traits_path,
                result.dataset_export.raw_trait_mismatch_path,
                result.workflow_bundle.workflow_summary_path,
                result.workflow_bundle.raw_sequence_findings_path,
                result.workflow_bundle.raw_sequence_repair_path,
                result.workflow_bundle.repaired_sequence_input_path,
                result.workflow_bundle.repaired_sequence_validation_path,
                result.workflow_bundle.coding_sequence_exclusions_path,
                result.workflow_bundle.prepared_coding_sequences_path,
                result.workflow_bundle.raw_trait_linkage_path,
                result.workflow_bundle.trait_duplicates_path,
                result.workflow_bundle.trait_missing_values_path,
                result.workflow_bundle.sequence_outliers_path,
                result.workflow_bundle.tree_issues_path,
                result.workflow_bundle.repair_actions_path,
                result.workflow_bundle.cleaned_traits_path,
                result.workflow_bundle.cleaned_alignment_path,
                result.workflow_bundle.cleaned_tree_path,
                result.workflow_bundle.cleaned_linkage_path,
                result.workflow_bundle.cleaned_validation_path,
                result.overview_path,
            ],
        )
        if args.json:
            expected_output_count = len(
                list(result.dataset_export.expected_output_root.glob("*"))
            )
            _print_result(
                build_command_result(
                    command="demo",
                    inputs=[],
                    outputs=outputs,
                    metrics={
                        "artifact_count": len(outputs),
                        "raw_taxon_count": result.workflow_bundle.raw_taxon_count,
                        "cleaned_taxon_count": (
                            result.workflow_bundle.cleaned_taxon_count
                        ),
                        "duplicate_sequence_identifier_count": (
                            result.workflow_bundle.duplicate_sequence_identifier_count
                        ),
                        "illegal_character_count": (
                            result.workflow_bundle.illegal_character_count
                        ),
                        "empty_sequence_count": (
                            result.workflow_bundle.empty_sequence_count
                        ),
                        "raw_sequence_length_outlier_count": (
                            result.workflow_bundle.raw_sequence_length_outlier_count
                        ),
                        "duplicate_trait_taxon_count": (
                            result.workflow_bundle.duplicate_trait_taxon_count
                        ),
                        "missing_trait_value_count": (
                            result.workflow_bundle.missing_trait_value_count
                        ),
                        "sequence_outlier_count": (
                            result.workflow_bundle.sequence_outlier_count
                        ),
                        "tree_zero_length_branch_count": (
                            result.workflow_bundle.tree_zero_length_branch_count
                        ),
                        "tree_negative_branch_count": (
                            result.workflow_bundle.tree_negative_branch_count
                        ),
                        "tree_long_branch_outlier_count": (
                            result.workflow_bundle.tree_long_branch_outlier_count
                        ),
                        "coding_frame_error_count": (
                            result.workflow_bundle.coding_frame_error_count
                        ),
                        "coding_internal_stop_count": (
                            result.workflow_bundle.coding_internal_stop_count
                        ),
                        "raw_trait_missing_from_traits_count": (
                            result.workflow_bundle.raw_trait_missing_from_traits_count
                        ),
                        "raw_trait_extra_taxon_count": (
                            result.workflow_bundle.raw_trait_extra_taxon_count
                        ),
                        "dropped_taxon_count": (
                            result.workflow_bundle.dropped_taxon_count
                        ),
                        "repaired_branch_count": (
                            result.workflow_bundle.repaired_branch_count
                        ),
                        "reference_output_count": expected_output_count,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_root)
        return 0

    if args.demo_command == "continuous-mode-recovery-panel":
        result = run_continuous_mode_recovery_panel_demo(args.out)
        outputs = _finalize_outputs(
            args,
            command="demo",
            inputs=[],
            outputs=[
                result.dataset_export.readme_path,
                result.dataset_export.default_tree_path,
                result.dataset_export.simulation_cases_path,
                result.workflow_bundle.workflow_summary_path,
                result.workflow_bundle.recovery_summary_path,
                result.workflow_bundle.parameter_recovery_path,
                result.workflow_bundle.parameter_comparison_path,
                result.workflow_bundle.model_choice_path,
                result.workflow_bundle.execution_review_path,
                result.workflow_bundle.warning_review_path,
                result.workflow_bundle.geiger_reference_path,
                result.overview_path,
            ],
        )
        if args.json:
            expected_output_count = len(
                [
                    path
                    for path in result.dataset_export.expected_output_root.rglob("*")
                    if path.is_file()
                ]
            )
            _print_result(
                build_command_result(
                    command="demo",
                    inputs=[],
                    outputs=outputs,
                    metrics={
                        "artifact_count": len(outputs),
                        "taxon_count": result.dataset.taxon_count,
                        "tree_count": result.dataset.tree_count,
                        "case_count": result.dataset.case_count,
                        "selection_review_case_count": (
                            result.workflow_bundle.selection_review_case_count
                        ),
                        "selection_match_count": (
                            result.workflow_bundle.selection_match_count
                        ),
                        "geiger_selection_match_count": (
                            result.workflow_bundle.geiger_selection_match_count
                        ),
                        "parameter_pass_count": (
                            result.workflow_bundle.parameter_pass_count
                        ),
                        "parameter_row_count": (
                            result.workflow_bundle.parameter_row_count
                        ),
                        "parameter_comparison_row_count": (
                            result.workflow_bundle.parameter_comparison_row_count
                        ),
                        "parameter_closer_to_truth_count_bijux": (
                            result.workflow_bundle.parameter_closer_to_truth_count_bijux
                        ),
                        "parameter_closer_to_truth_count_geiger": (
                            result.workflow_bundle.parameter_closer_to_truth_count_geiger
                        ),
                        "expected_warning_case_count": (
                            result.workflow_bundle.expected_warning_case_count
                        ),
                        "expected_warning_present_count": (
                            result.workflow_bundle.expected_warning_present_count
                        ),
                        "reference_output_count": expected_output_count,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_root)
        return 0

    if args.demo_command == "discrete-mode-recovery-panel":
        result = run_discrete_mode_recovery_panel_demo(args.out)
        outputs = _finalize_outputs(
            args,
            command="demo",
            inputs=[],
            outputs=[
                result.dataset_export.readme_path,
                result.dataset_export.default_tree_path,
                result.dataset_export.simulation_cases_path,
                result.workflow_bundle.workflow_summary_path,
                result.workflow_bundle.recovery_summary_path,
                result.workflow_bundle.rate_recovery_path,
                result.workflow_bundle.rate_comparison_path,
                result.workflow_bundle.model_choice_path,
                result.workflow_bundle.execution_review_path,
                result.workflow_bundle.warning_review_path,
                result.workflow_bundle.geiger_reference_path,
                result.overview_path,
            ],
        )
        if args.json:
            expected_output_count = len(
                [
                    path
                    for path in result.dataset_export.expected_output_root.rglob("*")
                    if path.is_file()
                ]
            )
            _print_result(
                build_command_result(
                    command="demo",
                    inputs=[],
                    outputs=outputs,
                    metrics={
                        "artifact_count": len(outputs),
                        "taxon_count": result.dataset.taxon_count,
                        "tree_count": result.dataset.tree_count,
                        "case_count": result.dataset.case_count,
                        "selection_review_case_count": (
                            result.workflow_bundle.selection_review_case_count
                        ),
                        "selection_match_count": (
                            result.workflow_bundle.selection_match_count
                        ),
                        "geiger_selection_match_count": (
                            result.workflow_bundle.geiger_selection_match_count
                        ),
                        "rate_pass_count": result.workflow_bundle.rate_pass_count,
                        "governed_rate_row_count": (
                            result.workflow_bundle.governed_rate_row_count
                        ),
                        "rate_row_count": result.workflow_bundle.rate_row_count,
                        "governed_rate_comparison_row_count": (
                            result.workflow_bundle.governed_rate_comparison_row_count
                        ),
                        "rate_comparison_row_count": (
                            result.workflow_bundle.rate_comparison_row_count
                        ),
                        "rate_closer_to_truth_count_bijux": (
                            result.workflow_bundle.rate_closer_to_truth_count_bijux
                        ),
                        "rate_closer_to_truth_count_geiger": (
                            result.workflow_bundle.rate_closer_to_truth_count_geiger
                        ),
                        "expected_warning_case_count": (
                            result.workflow_bundle.expected_warning_case_count
                        ),
                        "expected_warning_present_count": (
                            result.workflow_bundle.expected_warning_present_count
                        ),
                        "reference_output_count": expected_output_count,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_root)
        return 0

    if args.demo_command == "known-answer-reference-panel":
        result = run_known_answer_reference_demo(args.out)
        outputs = _finalize_outputs(
            args,
            command="demo",
            inputs=[],
            outputs=[
                result.dataset_export.readme_path,
                result.dataset_export.true_tree_path,
                result.dataset_export.alignment_path,
                result.dataset_export.continuous_traits_path,
                result.dataset_export.ou_traits_path,
                result.dataset_export.discrete_traits_path,
                result.dataset_export.host_traits_path,
                result.dataset_export.geographic_traits_path,
                result.dataset_export.true_parameters_path,
                result.dataset_export.true_continuous_nodes_path,
                result.dataset_export.true_ou_nodes_path,
                result.dataset_export.true_discrete_nodes_path,
                result.dataset_export.true_host_nodes_path,
                result.dataset_export.true_geographic_nodes_path,
                result.dataset_export.true_host_switch_events_path,
                result.dataset_export.true_geographic_transition_events_path,
                result.dataset_export.recovery_thresholds_path,
                result.workflow_bundle.workflow_summary_path,
                result.workflow_bundle.distance_tree_path,
                result.workflow_bundle.tree_recovery_path,
                result.workflow_bundle.parameter_recovery_path,
                result.workflow_bundle.brownian_fit_summary_path,
                result.workflow_bundle.ou_fit_summary_path,
                result.workflow_bundle.continuous_ancestral_summary_path,
                result.workflow_bundle.continuous_ancestral_uncertainty_path,
                result.workflow_bundle.continuous_node_recovery_path,
                result.workflow_bundle.discrete_ancestral_summary_path,
                result.workflow_bundle.discrete_ancestral_probability_path,
                result.workflow_bundle.discrete_node_recovery_path,
                result.workflow_bundle.host_switch_summary_path,
                result.workflow_bundle.host_state_nodes_path,
                result.workflow_bundle.host_switch_branches_path,
                result.workflow_bundle.host_node_recovery_path,
                result.workflow_bundle.host_event_recovery_path,
                result.workflow_bundle.geographic_ancestral_summary_path,
                result.workflow_bundle.geographic_state_probability_path,
                result.workflow_bundle.geographic_transition_summary_path,
                result.workflow_bundle.geographic_node_recovery_path,
                result.workflow_bundle.geographic_event_recovery_path,
                result.workflow_bundle.threshold_evaluation_path,
                result.overview_path,
            ],
        )
        if args.json:
            expected_output_count = len(
                list(result.dataset_export.expected_output_root.glob("*"))
            )
            _print_result(
                build_command_result(
                    command="demo",
                    inputs=[],
                    outputs=outputs,
                    metrics={
                        "artifact_count": len(outputs),
                        "taxon_count": result.dataset.taxon_count,
                        "sequence_length": result.dataset.sequence_length,
                        "distance_method": result.dataset.distance_method,
                        "distance_model": result.dataset.distance_model,
                        "rooted_topology_equal": (
                            result.workflow_bundle.rooted_topology_equal
                        ),
                        "same_unrooted_topology": (
                            result.workflow_bundle.same_unrooted_topology
                        ),
                        "same_taxa_different_rooting": (
                            result.workflow_bundle.same_taxa_different_rooting
                        ),
                        "robinson_foulds_distance": (
                            result.workflow_bundle.robinson_foulds_distance
                        ),
                        "parameter_row_count": (
                            result.workflow_bundle.parameter_row_count
                        ),
                        "threshold_pass_count": (
                            result.workflow_bundle.threshold_pass_count
                        ),
                        "threshold_row_count": (
                            result.workflow_bundle.threshold_row_count
                        ),
                        "continuous_internal_node_mean_absolute_error": (
                            result.workflow_bundle.continuous_internal_node_mean_absolute_error
                        ),
                        "discrete_internal_node_accuracy": (
                            result.workflow_bundle.discrete_internal_node_accuracy
                        ),
                        "host_internal_node_accuracy": (
                            result.workflow_bundle.host_internal_node_accuracy
                        ),
                        "host_event_accuracy": (
                            result.workflow_bundle.host_event_accuracy
                        ),
                        "geographic_internal_node_accuracy": (
                            result.workflow_bundle.geographic_internal_node_accuracy
                        ),
                        "geographic_event_accuracy": (
                            result.workflow_bundle.geographic_event_accuracy
                        ),
                        "reference_output_count": expected_output_count,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_root)
        return 0

    raise NotImplementedError(f"unsupported demo command: {args.demo_command}")
