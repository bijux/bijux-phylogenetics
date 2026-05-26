from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.bayesian.beast.execution import (
    run_beast_posterior_inference,
)
from bijux_phylogenetics.bayesian.beast.xml_analysis import (
    prepare_beast_time_tree_analysis,
)
from bijux_phylogenetics.bayesian.mrbayes import (
    prepare_mrbayes_analysis,
    run_mrbayes_posterior_inference,
)
from bijux_phylogenetics.engines import (
    run_alignment_trimming,
    run_bootstrap_support_estimation,
    run_fast_tree_inference,
    run_model_selection,
    run_multiple_sequence_alignment,
)
from bijux_phylogenetics.engines.validation.matrix import (
    ExternalEngineValidationCase,
    build_external_engine_validation_case,
    build_governed_beast_fixture_validation_case,
)
from bijux_phylogenetics.datasets.shared_fixtures import (
    get_shared_beast_posterior_fixture,
)

from .external_engines import (
    real_beast_executable,
    require_alignment_validation_matrix_executables,
    require_bayesian_validation_matrix_executables,
)

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
FIXTURE_GROUPS = ("trees", "alignments", "metadata", "expected")


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def build_real_alignment_validation_cases(
    tmp_path: Path,
) -> list[ExternalEngineValidationCase]:
    executables = require_alignment_validation_matrix_executables()

    mafft_report = run_multiple_sequence_alignment(
        fixture("alignments/example_sequences_raw.fasta"),
        tmp_path / "alignment" / "real-mafft-alignment.fasta",
        executable=executables["mafft"],
        mode="linsi",
    )
    trimal_report = run_alignment_trimming(
        fixture("alignments/example_alignment_trim.fasta"),
        tmp_path / "trim" / "real-trimal-trimmed.fasta",
        executable=executables["trimal"],
        mode="gap-threshold",
        gap_threshold=0.8,
    )
    model_report = run_model_selection(
        fixture("alignments/example_alignment.fasta"),
        out_dir=tmp_path / "model",
        executable=executables["iqtree"],
        prefix="real",
        sequence_type="dna",
    )
    bootstrap_report = run_bootstrap_support_estimation(
        fixture("alignments/example_alignment.fasta"),
        out_dir=tmp_path / "bootstrap",
        model=model_report.selected_model or "GTR+G",
        executable=executables["iqtree"],
        prefix="real",
        sequence_type="dna",
        replicates=1000,
    )
    fasttree_report = run_fast_tree_inference(
        fixture("alignments/example_alignment.fasta"),
        tmp_path / "fasttree" / "real-fasttree.nwk",
        executable=executables["fasttree"],
        sequence_type="dna",
    )

    return [
        build_external_engine_validation_case("mafft alignment", mafft_report),
        build_external_engine_validation_case("trimal trimming", trimal_report),
        build_external_engine_validation_case(
            "iqtree model selection",
            model_report,
            notes=[
                f"selected model: {model_report.selected_model}",
            ],
        ),
        build_external_engine_validation_case(
            "iqtree bootstrap support",
            bootstrap_report,
            notes=[
                f"support value count: {bootstrap_report.iqtree_summary.support_value_count if bootstrap_report.iqtree_summary is not None else 0}",
            ],
        ),
        build_external_engine_validation_case("fasttree inference", fasttree_report),
    ]


def _build_beast_validation_case(tmp_path: Path) -> ExternalEngineValidationCase:
    executable = real_beast_executable()
    if executable is None:
        return build_governed_beast_fixture_validation_case(
            "beast fixture parser acceptance",
            get_shared_beast_posterior_fixture("strict_yule_real_posterior"),
        )
    xml_path = tmp_path / "live-strict-yule.xml"
    prepare_beast_time_tree_analysis(
        fixture("example_alignment.fasta"),
        xml_path,
        clock_model="strict",
        tree_prior="yule",
        chain_length=1000,
        log_every=20,
    )
    report = run_beast_posterior_inference(
        xml_path,
        executable=executable,
        threads=1,
        seed=1,
    )
    return build_external_engine_validation_case(
        "beast posterior inference",
        report,
    )


def build_real_bayesian_validation_cases(
    tmp_path: Path,
) -> list[ExternalEngineValidationCase]:
    executables = require_bayesian_validation_matrix_executables()

    nexus_path = tmp_path / "partitioned-analysis.nex"
    prepare_mrbayes_analysis(
        fixture("alignments/example_multilocus_alignment.fasta"),
        nexus_path,
        partition_path=fixture("alignments/example_multilocus_partitions.txt"),
        model="gtr",
        rates="gamma",
        ngen=20,
        samplefreq=10,
        printfreq=10,
        burnin_fraction=0.25,
    )
    mrbayes_report = run_mrbayes_posterior_inference(
        nexus_path,
        executable=executables["mrbayes"],
        resume=False,
    )

    return [
        build_external_engine_validation_case(
            "mrbayes posterior inference",
            mrbayes_report,
        ),
        _build_beast_validation_case(tmp_path),
    ]
