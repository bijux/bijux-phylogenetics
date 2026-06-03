from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.datasets.shared_fixtures import (
    get_shared_beast_posterior_fixture,
)

from .matrix import (
    ExternalEngineValidationMatrixReport,
    build_external_engine_validation_case,
    build_external_engine_validation_matrix,
    build_governed_beast_fixture_validation_case,
    merge_external_engine_validation_matrices,
)


@dataclass(frozen=True, slots=True)
class AlignmentValidationMatrixInputs:
    raw_sequence_path: Path
    trimming_alignment_path: Path
    inference_alignment_path: Path


@dataclass(frozen=True, slots=True)
class BayesianValidationMatrixInputs:
    mrbayes_alignment_path: Path
    mrbayes_partition_path: Path
    beast_alignment_path: Path


def run_alignment_engine_validation_matrix(
    *,
    inputs: AlignmentValidationMatrixInputs,
    out_dir: Path,
    mafft_executable: str | Path,
    trimal_executable: str | Path,
    iqtree_executable: str | Path,
    fasttree_executable: str | Path,
) -> ExternalEngineValidationMatrixReport:
    """Run the governed real-engine alignment matrix and collect reviewer metadata."""
    from bijux_phylogenetics.engines.workflows import (
        run_alignment_trimming,
        run_bootstrap_support_estimation,
        run_fast_tree_inference,
        run_model_selection,
        run_multiple_sequence_alignment,
    )

    mafft_report = run_multiple_sequence_alignment(
        inputs.raw_sequence_path,
        out_dir / "alignment" / "real-mafft-alignment.fasta",
        executable=mafft_executable,
        mode="linsi",
    )
    trimal_report = run_alignment_trimming(
        inputs.trimming_alignment_path,
        out_dir / "trim" / "real-trimal-trimmed.fasta",
        executable=trimal_executable,
        mode="gap-threshold",
        gap_threshold=0.8,
    )
    model_report = run_model_selection(
        inputs.inference_alignment_path,
        out_dir=out_dir / "model",
        executable=iqtree_executable,
        prefix="real",
        sequence_type="dna",
    )
    bootstrap_report = run_bootstrap_support_estimation(
        inputs.inference_alignment_path,
        out_dir=out_dir / "bootstrap",
        model=model_report.selected_model or "GTR+G",
        executable=iqtree_executable,
        prefix="real",
        sequence_type="dna",
        replicates=1000,
    )
    fasttree_report = run_fast_tree_inference(
        inputs.inference_alignment_path,
        out_dir / "fasttree" / "real-fasttree.nwk",
        executable=fasttree_executable,
        sequence_type="dna",
    )
    return build_external_engine_validation_matrix(
        [
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
                    (
                        "support value count: "
                        f"{bootstrap_report.iqtree_summary.support_value_count if bootstrap_report.iqtree_summary is not None else 0}"
                    ),
                ],
            ),
            build_external_engine_validation_case(
                "fasttree inference",
                fasttree_report,
            ),
        ]
    )


def run_bayesian_engine_validation_matrix(
    *,
    inputs: BayesianValidationMatrixInputs,
    out_dir: Path,
    mrbayes_executable: str | Path,
    beast_executable: str | Path | None = None,
    governed_beast_fixture_id: str = "strict_yule_real_posterior",
) -> ExternalEngineValidationMatrixReport:
    """Run the governed Bayesian matrix with live MrBayes and live-or-governed BEAST."""
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

    nexus_path = out_dir / "partitioned-analysis.nex"
    prepare_mrbayes_analysis(
        inputs.mrbayes_alignment_path,
        nexus_path,
        partition_path=inputs.mrbayes_partition_path,
        model="gtr",
        rates="gamma",
        ngen=20,
        samplefreq=10,
        printfreq=10,
        burnin_fraction=0.25,
    )
    mrbayes_report = run_mrbayes_posterior_inference(
        nexus_path,
        executable=mrbayes_executable,
        resume=False,
    )
    cases = [
        build_external_engine_validation_case(
            "mrbayes posterior inference",
            mrbayes_report,
        )
    ]
    if beast_executable is None:
        cases.append(
            build_governed_beast_fixture_validation_case(
                "beast fixture parser acceptance",
                get_shared_beast_posterior_fixture(governed_beast_fixture_id),
            )
        )
        return build_external_engine_validation_matrix(cases)

    xml_path = out_dir / "live-strict-yule.xml"
    prepare_beast_time_tree_analysis(
        inputs.beast_alignment_path,
        xml_path,
        clock_model="strict",
        tree_prior="yule",
        chain_length=1000,
        log_every=20,
    )
    beast_report = run_beast_posterior_inference(
        xml_path,
        executable=beast_executable,
        threads=1,
        seed=1,
    )
    cases.append(
        build_external_engine_validation_case(
            "beast posterior inference",
            beast_report,
        )
    )
    return build_external_engine_validation_matrix(cases)


def run_external_engine_validation_matrix(
    *,
    alignment_inputs: AlignmentValidationMatrixInputs,
    bayesian_inputs: BayesianValidationMatrixInputs,
    out_dir: Path,
    mafft_executable: str | Path,
    trimal_executable: str | Path,
    iqtree_executable: str | Path,
    fasttree_executable: str | Path,
    mrbayes_executable: str | Path,
    beast_executable: str | Path | None = None,
    governed_beast_fixture_id: str = "strict_yule_real_posterior",
) -> ExternalEngineValidationMatrixReport:
    """Run the full governed external-engine matrix over alignment and Bayesian engines."""
    alignment_matrix = run_alignment_engine_validation_matrix(
        inputs=alignment_inputs,
        out_dir=out_dir / "alignment",
        mafft_executable=mafft_executable,
        trimal_executable=trimal_executable,
        iqtree_executable=iqtree_executable,
        fasttree_executable=fasttree_executable,
    )
    bayesian_matrix = run_bayesian_engine_validation_matrix(
        inputs=bayesian_inputs,
        out_dir=out_dir / "bayesian",
        mrbayes_executable=mrbayes_executable,
        beast_executable=beast_executable,
        governed_beast_fixture_id=governed_beast_fixture_id,
    )
    return merge_external_engine_validation_matrices(
        [alignment_matrix, bayesian_matrix]
    )
