from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import (
    _print_result,
    _write_locus_occupancy_loci_tsv,
    _write_locus_occupancy_matrix_tsv,
    _write_locus_occupancy_taxa_tsv,
)
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.io.fasta import write_fasta_alignment
from bijux_phylogenetics.io.fasta.cleaning import (
    clean_alignment_with_profile,
    compare_alignment_versions,
    compute_pairwise_sequence_identity_matrix,
    list_alignment_filter_profiles,
    trim_alignment,
    write_sequence_identity_matrix,
)
from bijux_phylogenetics.io.fasta.records import summarise_fasta
from bijux_phylogenetics.phylo.alignment.concatenation import (
    concatenate_locus_alignments,
)
from bijux_phylogenetics.phylo.alignment.occupancy import (
    build_locus_occupancy_report,
    filter_locus_occupancy,
    write_locus_partitions,
)
from bijux_phylogenetics.phylo.alignment.partitions import (
    build_partition_summary_report,
    parse_locus_partitions,
    write_partition_summary_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_alignment_matrix_commands(alignment_subparsers: Any) -> None:
    alignment_occupancy = alignment_subparsers.add_parser(
        "occupancy",
        help="Quantify per-taxon and per-locus occupancy across a concatenated multi-locus alignment.",
    )
    alignment_occupancy.add_argument("alignment", type=Path)
    alignment_occupancy.add_argument("partitions", type=Path)
    alignment_occupancy.add_argument("--taxon-coverage-threshold", type=float)
    alignment_occupancy.add_argument("--locus-coverage-threshold", type=float)
    alignment_occupancy.add_argument(
        "--minimum-locus-occupancy",
        type=float,
        default=0.0,
        help=(
            "Require at least this per-taxon per-locus occupancy fraction before "
            "a locus counts as covered for thresholding."
        ),
    )
    alignment_occupancy.add_argument(
        "--taxa-out",
        type=Path,
        help="Write per-taxon locus coverage as TSV.",
    )
    alignment_occupancy.add_argument(
        "--loci-out",
        type=Path,
        help="Write per-locus taxon coverage as TSV.",
    )
    alignment_occupancy.add_argument(
        "--matrix-out",
        type=Path,
        help="Write the taxon-by-locus occupancy matrix as TSV.",
    )
    alignment_occupancy.add_argument(
        "--filtered-alignment-out",
        type=Path,
        help="Write the retained alignment after occupancy filtering.",
    )
    alignment_occupancy.add_argument(
        "--filtered-partitions-out",
        type=Path,
        help="Write the retained partition file after occupancy filtering.",
    )
    alignment_occupancy.add_argument(
        "--json", action="store_true", help="Emit the occupancy report as JSON."
    )
    _add_manifest_argument(alignment_occupancy)

    alignment_concatenate = alignment_subparsers.add_parser(
        "concatenate",
        help="Concatenate aligned per-locus FASTA inputs into one partitioned supermatrix.",
    )
    alignment_concatenate.add_argument("alignments", nargs="+", type=Path)
    alignment_concatenate.add_argument(
        "--locus-name",
        action="append",
        dest="locus_names",
        help="Override one locus name in input order. Repeat once per alignment when stems are not the durable locus names you want in the partition file.",
    )
    alignment_concatenate.add_argument(
        "--data-type",
        action="append",
        dest="data_types",
        help="Override one partition data type in input order. Repeat once per alignment when ambiguous residues make DNA and protein loci indistinguishable by content alone.",
    )
    alignment_concatenate.add_argument("--out", required=True, type=Path)
    alignment_concatenate.add_argument("--partitions-out", required=True, type=Path)
    alignment_concatenate.add_argument("--matrix-out", required=True, type=Path)
    alignment_concatenate.add_argument(
        "--taxa-out",
        type=Path,
        help="Write per-taxon locus coverage as TSV.",
    )
    alignment_concatenate.add_argument(
        "--loci-out",
        type=Path,
        help="Write per-locus taxon coverage as TSV.",
    )
    alignment_concatenate.add_argument(
        "--json", action="store_true", help="Emit the concatenation report as JSON."
    )
    _add_manifest_argument(alignment_concatenate)

    alignment_partition_summary = alignment_subparsers.add_parser(
        "partition-summary",
        help="Validate one partition file against an aligned matrix and summarize each locus as TSV-ready rows.",
    )
    alignment_partition_summary.add_argument("alignment", type=Path)
    alignment_partition_summary.add_argument("partitions", type=Path)
    alignment_partition_summary.add_argument(
        "--out",
        type=Path,
        help="Write the partition summary table as TSV.",
    )
    alignment_partition_summary.add_argument(
        "--json", action="store_true", help="Emit the partition summary report as JSON."
    )
    _add_manifest_argument(alignment_partition_summary)

    alignment_filter = alignment_subparsers.add_parser(
        "filter",
        help="Clean an alignment through one named profile and report what changed.",
    )
    alignment_filter.add_argument("alignment", type=Path)
    alignment_filter.add_argument(
        "--profile",
        required=True,
        choices=tuple(profile.name for profile in list_alignment_filter_profiles()),
    )
    alignment_filter.add_argument("--out", required=True, type=Path)
    alignment_filter.add_argument("--group-table", type=Path)
    alignment_filter.add_argument("--group-columns")
    alignment_filter.add_argument(
        "--json", action="store_true", help="Emit the cleaning report as JSON."
    )
    _add_manifest_argument(alignment_filter)

    alignment_compare = alignment_subparsers.add_parser(
        "compare",
        help="Compare two alignment versions for taxa, sites, missingness, gaps, signal, and composition.",
    )
    alignment_compare.add_argument("left_alignment", type=Path)
    alignment_compare.add_argument("right_alignment", type=Path)
    alignment_compare.add_argument(
        "--json", action="store_true", help="Emit the comparison report as JSON."
    )
    _add_manifest_argument(alignment_compare)

    alignment_trim = alignment_subparsers.add_parser(
        "trim",
        help="Trim all-gap or all-missing sites and optionally drop high-missingness sequences.",
    )
    alignment_trim.add_argument("alignment", type=Path)
    alignment_trim.add_argument("--out", required=True, type=Path)
    alignment_trim.add_argument("--keep-all-gap-sites", action="store_true")
    alignment_trim.add_argument("--keep-all-missing-sites", action="store_true")
    alignment_trim.add_argument("--site-missingness-threshold", type=float)
    alignment_trim.add_argument("--sequence-missingness-threshold", type=float)
    alignment_trim.add_argument(
        "--json", action="store_true", help="Emit the trimming report as JSON."
    )
    _add_manifest_argument(alignment_trim)

    alignment_identity = alignment_subparsers.add_parser(
        "identity-matrix",
        help="Compute a pairwise sequence identity matrix.",
    )
    alignment_identity.add_argument("alignment", type=Path)
    alignment_identity.add_argument("--out", type=Path, help="Write the matrix as TSV.")
    alignment_identity.add_argument(
        "--json", action="store_true", help="Emit the report as JSON."
    )
    _add_manifest_argument(alignment_identity)


def run_alignment_matrix_command(args: Any) -> int | None:
    if args.alignment_command == "concatenate":
        records, partitions, report = concatenate_locus_alignments(
            list(args.alignments),
            locus_names=None if args.locus_names is None else tuple(args.locus_names),
            data_types=None if args.data_types is None else tuple(args.data_types),
            concatenated_alignment_path=args.out,
            concatenated_partition_path=args.partitions_out,
        )
        command_outputs = [
            write_fasta_alignment(args.out, records),
            write_locus_partitions(args.partitions_out, partitions),
            _write_locus_occupancy_matrix_tsv(
                args.matrix_out,
                report.occupancy_report,
            ),
        ]
        if args.taxa_out is not None:
            command_outputs.append(
                _write_locus_occupancy_taxa_tsv(
                    args.taxa_out,
                    report.occupancy_report,
                )
            )
        if args.loci_out is not None:
            command_outputs.append(
                _write_locus_occupancy_loci_tsv(
                    args.loci_out,
                    report.occupancy_report,
                )
            )
        outputs = _finalize_outputs(
            args,
            command="alignment",
            inputs=list(args.alignments),
            outputs=command_outputs,
        )
        _print_result(
            build_command_result(
                command="alignment",
                inputs=list(args.alignments),
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "taxon_count": report.taxon_count,
                    "locus_count": report.locus_count,
                    "alignment_length": report.alignment_length,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.alignment_command == "occupancy":
        report = build_locus_occupancy_report(
            args.alignment,
            args.partitions,
            taxon_coverage_threshold=args.taxon_coverage_threshold,
            locus_coverage_threshold=args.locus_coverage_threshold,
            minimum_locus_occupancy=args.minimum_locus_occupancy,
        )
        command_outputs: list[Path] = []
        if args.taxa_out is not None:
            command_outputs.append(
                _write_locus_occupancy_taxa_tsv(args.taxa_out, report)
            )
        if args.loci_out is not None:
            command_outputs.append(
                _write_locus_occupancy_loci_tsv(args.loci_out, report)
            )
        if args.matrix_out is not None:
            command_outputs.append(
                _write_locus_occupancy_matrix_tsv(args.matrix_out, report)
            )

        filter_report = None
        if (
            args.taxon_coverage_threshold is not None
            or args.locus_coverage_threshold is not None
            or args.filtered_alignment_out is not None
            or args.filtered_partitions_out is not None
        ):
            filtered_records, filtered_partitions, filter_report = (
                filter_locus_occupancy(
                    args.alignment,
                    args.partitions,
                    taxon_coverage_threshold=args.taxon_coverage_threshold,
                    locus_coverage_threshold=args.locus_coverage_threshold,
                    minimum_locus_occupancy=args.minimum_locus_occupancy,
                )
            )
            if args.filtered_alignment_out is not None:
                command_outputs.append(
                    write_fasta_alignment(args.filtered_alignment_out, filtered_records)
                )
            if args.filtered_partitions_out is not None:
                command_outputs.append(
                    write_locus_partitions(
                        args.filtered_partitions_out,
                        filtered_partitions,
                    )
                )

        outputs = _finalize_outputs(
            args,
            command="alignment",
            inputs=[args.alignment, args.partitions],
            outputs=command_outputs,
        )
        warnings = list(
            dict.fromkeys(
                report.warnings
                + ([] if filter_report is None else filter_report.final_report.warnings)
            )
        )
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment, args.partitions],
                outputs=outputs,
                warnings=warnings,
                metrics={
                    "taxon_count": report.taxon_count,
                    "locus_count": report.locus_count,
                    "low_coverage_taxon_count": len(report.low_coverage_taxa),
                    "low_coverage_locus_count": len(report.low_coverage_loci),
                    "minimum_locus_occupancy": report.minimum_locus_occupancy,
                    "filtered_taxon_count": (
                        report.taxon_count
                        if filter_report is None
                        else len(filter_report.retained_taxa)
                    ),
                    "filtered_locus_count": (
                        report.locus_count
                        if filter_report is None
                        else len(filter_report.retained_loci)
                    ),
                },
                data={"report": report, "filter_report": filter_report},
            ),
            json_output=args.json,
        )
        return 0
    if args.alignment_command == "partition-summary":
        alignment_summary = summarise_fasta(args.alignment)
        report = build_partition_summary_report(
            parse_locus_partitions(args.partitions),
            alignment_length=alignment_summary.alignment_length,
        )
        command_outputs: list[Path] = []
        if args.out is not None:
            command_outputs.append(write_partition_summary_table(args.out, report))
        outputs = _finalize_outputs(
            args,
            command="alignment",
            inputs=[args.alignment, args.partitions],
            outputs=command_outputs,
        )
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment, args.partitions],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "partition_count": report.partition_count,
                    "assigned_site_count": report.assigned_site_count,
                    "unassigned_site_count": report.unassigned_site_count,
                    "mixed_data_types": report.mixed_data_types,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.alignment_command == "filter":
        group_columns = (
            None
            if not args.group_columns
            else [
                column.strip()
                for column in args.group_columns.split(",")
                if column.strip()
            ]
        )
        records, report = clean_alignment_with_profile(
            args.alignment,
            profile_name=args.profile,
            group_table_path=args.group_table,
            group_columns=group_columns,
        )
        output_path = write_fasta_alignment(args.out, records)
        filter_inputs = [
            args.alignment,
            *([args.group_table] if args.group_table is not None else []),
        ]
        outputs = _finalize_outputs(
            args,
            command="alignment",
            inputs=filter_inputs,
            outputs=[output_path],
        )
        _print_result(
            build_command_result(
                command="alignment",
                inputs=filter_inputs,
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "profile": report.profile.name,
                    "trimmed_sequence_count": report.trim.trimmed_sequence_count,
                    "trimmed_alignment_length": report.trim.trimmed_alignment_length,
                    "signal_warning_count": len(report.signal_warnings),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.alignment_command == "compare":
        report = compare_alignment_versions(args.left_alignment, args.right_alignment)
        outputs = _finalize_outputs(
            args,
            command="alignment",
            inputs=[args.left_alignment, args.right_alignment],
        )
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.left_alignment, args.right_alignment],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "shared_taxa": len(report.shared_taxa),
                    "left_only_taxa": len(report.left_only_taxa),
                    "right_only_taxa": len(report.right_only_taxa),
                    "left_alignment_length": report.left_alignment_length,
                    "right_alignment_length": report.right_alignment_length,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.alignment_command == "trim":
        records, report = trim_alignment(
            args.alignment,
            remove_all_gap_sites=not args.keep_all_gap_sites,
            remove_all_missing_sites=not args.keep_all_missing_sites,
            site_missingness_threshold=args.site_missingness_threshold,
            sequence_missingness_threshold=args.sequence_missingness_threshold,
        )
        output_path = write_fasta_alignment(args.out, records)
        outputs = _finalize_outputs(
            args,
            command="alignment",
            inputs=[args.alignment],
            outputs=[output_path],
        )
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                metrics={
                    "trimmed_sequence_count": report.trimmed_sequence_count,
                    "trimmed_alignment_length": report.trimmed_alignment_length,
                    "removed_column_count": len(report.removed_columns),
                    "removed_sequence_count": len(report.removed_sequences),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.alignment_command == "identity-matrix":
        report = compute_pairwise_sequence_identity_matrix(args.alignment)
        outputs: list[Path | str] = []
        if args.out is not None:
            outputs.append(write_sequence_identity_matrix(args.out, report))
        outputs = _finalize_outputs(
            args,
            command="alignment",
            inputs=[args.alignment],
            outputs=outputs,
        )
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                metrics={
                    "sequence_count": len(report.identifiers),
                    "pair_count": len(report.pairs),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    return None
