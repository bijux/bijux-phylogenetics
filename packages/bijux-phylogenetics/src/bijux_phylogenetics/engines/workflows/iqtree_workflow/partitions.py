from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.io.fasta import load_fasta_alignment, write_fasta_alignment
from bijux_phylogenetics.io.fasta.records import summarise_fasta
from bijux_phylogenetics.phylo.alignment import AlignmentRecord
from bijux_phylogenetics.phylo.alignment.partitions import (
    build_partition_summary_report,
    normalize_partition_data_type,
    parse_locus_partitions,
    slice_partition_sequence,
    write_locus_partitions,
    write_partition_summary_table,
)
from bijux_phylogenetics.runtime.errors import EngineWorkflowError

from ..models import PreparedIqtreePartitions as _PreparedIqtreePartitions
from ..state import _partition_alignment_file_name, _partition_support_path


def _prepare_iqtree_partitions(
    input_path: Path,
    partition_path: Path,
    *,
    prefix_path: Path,
) -> _PreparedIqtreePartitions:
    records = load_fasta_alignment(input_path)
    alignment_summary = summarise_fasta(input_path)
    partitions = parse_locus_partitions(partition_path)
    summary = build_partition_summary_report(
        partitions,
        alignment_length=alignment_summary.alignment_length,
    )
    summary_path = _partition_support_path(prefix_path, "partition-summary.tsv")
    write_partition_summary_table(summary_path, summary)
    notes = [
        f"validated {summary.partition_count} partitions across {summary.assigned_site_count} assigned sites",
    ]
    output_paths: dict[str, Path] = {
        "partition_summary": summary_path,
    }

    declared_types = {
        normalize_partition_data_type(partition.data_type)
        for partition in partitions
        if partition.data_type is not None
    }
    if len(declared_types) <= 1:
        normalized_partition_path = _partition_support_path(
            prefix_path, "partition-scheme.partitions"
        )
        write_locus_partitions(normalized_partition_path, partitions)
        output_paths["partition_scheme"] = normalized_partition_path
        notes.append(
            "prepared a normalized partition scheme for single-alignment IQ-TREE analysis"
        )
        return _PreparedIqtreePartitions(
            command_args=[
                "-s",
                str(input_path.resolve()),
                "-p",
                str(normalized_partition_path.resolve()),
            ],
            summary=summary,
            output_paths=output_paths,
            notes=notes,
            mixed_data_types=False,
        )

    if any(partition.data_type is None for partition in partitions):
        raise EngineWorkflowError(
            "mixed partition analyses require every partition to declare a data_type"
        )
    unsupported_types = sorted(
        {
            data_type
            for data_type in declared_types
            if data_type is not None and data_type not in {"DNA", "RNA", "PROTEIN"}
        }
    )
    if unsupported_types:
        raise EngineWorkflowError(
            "mixed partition analyses currently support only DNA, RNA, and PROTEIN datatypes; "
            f"got: {', '.join(unsupported_types)}"
        )

    partition_alignment_dir = _partition_support_path(
        prefix_path, "partition-alignments"
    )
    partition_alignment_dir.mkdir(parents=True, exist_ok=True)
    lines = ["#nexus", "begin sets;"]
    for partition in partitions:
        partition_alignment_path = (
            partition_alignment_dir / _partition_alignment_file_name(partition)
        )
        write_fasta_alignment(
            partition_alignment_path,
            [
                AlignmentRecord(
                    identifier=record.identifier,
                    sequence=slice_partition_sequence(record.sequence, partition),
                )
                for record in records
            ],
        )
        output_paths[f"partition_alignment_{partition.name}"] = partition_alignment_path
        lines.append(
            f"    charset {partition.name} = {partition_alignment_path.name}: *;"
        )
    lines.append("end;")
    mixed_partition_path = _partition_support_path(prefix_path, "partition-scheme.nex")
    mixed_partition_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    output_paths["partition_scheme"] = mixed_partition_path
    notes.append(
        "prepared a mixed-datatype NEXUS partition scheme with one extracted alignment per partition"
    )
    return _PreparedIqtreePartitions(
        command_args=["-p", str(mixed_partition_path.resolve())],
        summary=summary,
        output_paths=output_paths,
        notes=notes,
        mixed_data_types=True,
    )
