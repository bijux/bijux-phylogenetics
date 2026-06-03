from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.engines.workflows.state import (
    _ensure_inference_ready_alignment,
)
from bijux_phylogenetics.io.fasta import (
    infer_alignment_alphabet,
    load_fasta_alignment,
)
from bijux_phylogenetics.phylo.alignment import AlignmentAlphabet
from bijux_phylogenetics.phylo.alignment.partitions import (
    LocusPartition,
    build_partition_summary_report,
    normalize_partition_data_type,
    parse_locus_partitions,
    partition_coordinate_text,
)
from bijux_phylogenetics.runtime.errors import EngineWorkflowError

from .models import MrBayesPreparationReport


def _mrbayes_datatype(alphabet: AlignmentAlphabet) -> str:
    if alphabet in {"dna", "rna"}:
        return "dna"
    if alphabet == "protein":
        return "protein"
    raise EngineWorkflowError(
        f"MrBayes preparation requires a recognized alignment alphabet, got {alphabet}"
    )


def _mrbayes_model_commands(
    *,
    alphabet: AlignmentAlphabet,
    model: str,
    rates: str,
    partition_count: int,
) -> list[str]:
    normalized_model = model.lower()
    normalized_rates = rates.lower()
    applyto_prefix = " applyto=(all)" if partition_count > 1 else ""
    if alphabet in {"dna", "rna"}:
        if normalized_model == "jc69":
            base = [f"lset{applyto_prefix} nst=1 rates=equal;"]
        elif normalized_model == "hky":
            base = [f"lset{applyto_prefix} nst=2 rates=equal;"]
        elif normalized_model == "gtr":
            base = [f"lset{applyto_prefix} nst=6 rates=equal;"]
        else:
            raise EngineWorkflowError(f"unsupported nucleotide MrBayes model: {model}")
    else:
        if normalized_model not in {"wag", "jones", "dayhoff", "poisson"}:
            raise EngineWorkflowError(f"unsupported protein MrBayes model: {model}")
        base = [
            f"prset{applyto_prefix} aamodelpr=fixed({normalized_model});",
            f"lset{applyto_prefix} rates=equal;",
        ]

    if normalized_rates == "equal":
        return base
    if normalized_rates in {"gamma", "invgamma", "propinv"}:
        updated: list[str] = []
        for line in base:
            if line.startswith("lset "):
                updated.append(line.replace("rates=equal", f"rates={normalized_rates}"))
            else:
                updated.append(line)
        return updated
    raise EngineWorkflowError(f"unsupported MrBayes rate model: {rates}")


def _mrbayes_partition_datatype(alphabet: AlignmentAlphabet) -> str:
    if alphabet in {"dna", "rna"}:
        return "DNA"
    if alphabet == "protein":
        return "PROTEIN"
    raise EngineWorkflowError(
        f"MrBayes preparation requires a recognized alignment alphabet, got {alphabet}"
    )


def _validate_mrbayes_partitions(
    partitions: tuple[LocusPartition, ...],
    *,
    alphabet: AlignmentAlphabet,
    alignment_length: int,
) -> tuple[list[str], list[str]]:
    summary = build_partition_summary_report(
        partitions, alignment_length=alignment_length
    )
    expected_data_type = _mrbayes_partition_datatype(alphabet)
    declared_data_types = [
        normalize_partition_data_type(partition.data_type)
        for partition in partitions
        if partition.data_type is not None
    ]
    mismatched = sorted(
        {
            data_type
            for data_type in declared_data_types
            if data_type is not None and data_type != expected_data_type
        }
    )
    if mismatched:
        raise EngineWorkflowError(
            "MrBayes preparation requires partition datatypes that match the "
            f"alignment alphabet {expected_data_type}, got {', '.join(mismatched)}"
        )
    return summary.declared_data_types, summary.warnings


def _mrbayes_charset_commands(partitions: tuple[LocusPartition, ...]) -> list[str]:
    return [
        f"charset {partition.name} = {partition_coordinate_text(partition)};"
        for partition in partitions
    ]


def _mrbayes_partition_commands(partitions: tuple[LocusPartition, ...]) -> list[str]:
    if len(partitions) <= 1:
        return []
    partition_name = "loci"
    joined_names = ", ".join(partition.name for partition in partitions)
    return [
        f"partition {partition_name} = {len(partitions)}: {joined_names};",
        f"set partition={partition_name};",
        "prset applyto=(all) ratepr=variable;",
    ]


def prepare_mrbayes_analysis(
    alignment_path: Path,
    nexus_path: Path,
    *,
    partition_path: Path | None = None,
    model: str = "gtr",
    rates: str = "gamma",
    ngen: int = 10000,
    nchains: int = 4,
    samplefreq: int = 100,
    printfreq: int = 100,
    burnin_fraction: float = 0.25,
) -> MrBayesPreparationReport:
    """Prepare a MrBayes NEXUS analysis specification from an aligned FASTA file."""
    _ensure_inference_ready_alignment(alignment_path)
    if not 0.0 <= burnin_fraction < 1.0:
        raise ValueError(
            f"burnin_fraction must be between 0 and 1, got {burnin_fraction}"
        )
    records = load_fasta_alignment(alignment_path)
    alphabet = infer_alignment_alphabet(records)
    datatype = _mrbayes_datatype(alphabet)
    partitions: tuple[LocusPartition, ...] = ()
    partition_data_types: list[str] = []
    partition_warnings: list[str] = []
    if partition_path is not None:
        partitions = parse_locus_partitions(partition_path)
        partition_data_types, partition_warnings = _validate_mrbayes_partitions(
            partitions,
            alphabet=alphabet,
            alignment_length=len(records[0].sequence),
        )
    model_commands = _mrbayes_model_commands(
        alphabet=alphabet,
        model=model,
        rates=rates,
        partition_count=len(partitions),
    )
    matrix_lines = "\n".join(
        f"{record.identifier} {record.sequence}" for record in records
    )
    command_block = "\n".join(
        [
            "begin mrbayes;",
            "  set autoclose=yes nowarn=yes;",
            *[f"  {line}" for line in _mrbayes_charset_commands(partitions)],
            *[f"  {line}" for line in _mrbayes_partition_commands(partitions)],
            *[f"  {line}" for line in model_commands],
            f"  mcmcp ngen={ngen} nchains={nchains} samplefreq={samplefreq} printfreq={printfreq};",
            "  mcmc;",
            f"  sump burninfrac={burnin_fraction:.6f};",
            f"  sumt burninfrac={burnin_fraction:.6f};",
            "end;",
        ]
    )
    nexus_text = "\n".join(
        [
            "#NEXUS",
            "",
            "begin data;",
            f"  dimensions ntax={len(records)} nchar={len(records[0].sequence)};",
            f"  format datatype={datatype} missing=? gap=-;",
            "  matrix",
            matrix_lines,
            "  ;",
            "end;",
            "",
            command_block,
            "",
        ]
    )
    nexus_path.parent.mkdir(parents=True, exist_ok=True)
    nexus_path.write_text(nexus_text, encoding="utf-8")
    return MrBayesPreparationReport(
        alignment_path=alignment_path,
        nexus_path=nexus_path,
        partition_path=partition_path,
        taxon_count=len(records),
        character_count=len(records[0].sequence),
        inferred_alphabet=alphabet,
        partition_count=max(len(partitions), 1),
        partition_names=[partition.name for partition in partitions],
        partition_data_types=partition_data_types,
        partition_warnings=partition_warnings,
        model=model,
        rates=rates,
        ngen=ngen,
        nchains=nchains,
        samplefreq=samplefreq,
        printfreq=printfreq,
        burnin_fraction=burnin_fraction,
    )
