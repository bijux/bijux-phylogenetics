from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.demo.shared import (
    count_expected_output_entries,
    emit_demo_result,
)
from bijux_phylogenetics.datasets import (
    run_catarrhine_mitogenome_five_locus_panel_demo,
    run_gnathostome_ortholog_protein_benchmark_demo,
    run_influenza_a_ha_reference_demo,
    run_pleistocene_bear_cytb_fragment_demo,
)


def add_sequence_demo_commands(demo_subparsers: Any) -> None:
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


def run_sequence_demo_command(args: Any) -> int | None:
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
        outputs = [
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
        ]
        return emit_demo_result(
            args,
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
                "reference_output_count": count_expected_output_entries(
                    result.dataset_export.expected_output_root
                ),
            },
            data=result,
            output_root=result.output_root,
        )

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
        outputs = [
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
        ]
        return emit_demo_result(
            args,
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
                "reference_output_count": count_expected_output_entries(
                    result.dataset_export.expected_output_root
                ),
            },
            data=result,
            output_root=result.output_root,
        )

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
        outputs = [
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
        ]
        return emit_demo_result(
            args,
            outputs=outputs,
            metrics={
                "artifact_count": len(outputs),
                "sequence_count": result.dataset.sequence_count,
                "degraded_sequence_count": len(result.dataset.degraded_sequence_ids),
                "selected_model": result.workflow_bundle.selected_model,
                "minimum_support": result.workflow_bundle.minimum_support,
                "maximum_support": result.workflow_bundle.maximum_support,
                "removed_column_count": result.workflow_bundle.removed_column_count,
                "cleaned_missing_data_fraction": (
                    result.workflow_bundle.cleaned_missing_data_fraction
                ),
                "reference_output_count": count_expected_output_entries(
                    result.dataset_export.expected_output_root
                ),
            },
            data=result,
            output_root=result.output_root,
        )

    if args.demo_command != "catarrhine-mitogenome-five-locus-panel":
        return None

    result = run_catarrhine_mitogenome_five_locus_panel_demo(
        args.out,
        iqtree_executable=args.iqtree_executable or "iqtree2",
        iqtree_seed=args.iqtree_seed,
        iqtree_threads=args.iqtree_threads,
        bootstrap_replicates=args.bootstrap_replicates,
    )
    outputs = [
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
    ]
    return emit_demo_result(
        args,
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
            "reference_output_count": count_expected_output_entries(
                result.dataset_export.expected_output_root
            ),
        },
        data=result,
        output_root=result.output_root,
    )
