from __future__ import annotations

from pathlib import Path
import subprocess

import pytest

from bijux_phylogenetics.bayesian.beast.execution import (
    run_beast_posterior_inference,
)
from bijux_phylogenetics.bayesian.beast.logs import (
    summarize_beast_log,
)
from bijux_phylogenetics.bayesian.beast.posterior_trees import (
    parse_beast_posterior_tree_samples,
    summarize_beast_posterior_topology_diversity,
    summarize_beast_posterior_trees,
)
from bijux_phylogenetics.bayesian.beast.xml_analysis import (
    prepare_beast_time_tree_analysis,
)
from bijux_phylogenetics.bayesian.mrbayes import (
    parse_mrbayes_consensus_tree,
    parse_mrbayes_mcmc_diagnostics,
    parse_mrbayes_parameter_traces,
    parse_mrbayes_posterior_tree_samples,
    prepare_mrbayes_analysis,
    run_mrbayes_posterior_inference,
)

from ..support.external_engines import (
    real_beast_executable,
    real_mrbayes_executable,
)

pytestmark = [pytest.mark.real_local, pytest.mark.engine_real]

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


def test_prepare_mrbayes_analysis_is_accepted_by_real_mrbayes_on_partitioned_input(
    tmp_path: Path,
) -> None:
    executable = real_mrbayes_executable()
    if executable is None:
        pytest.skip("real MrBayes executable is not available")
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

    report = run_mrbayes_posterior_inference(
        nexus_path, executable=executable, resume=False
    )

    assert report.output_paths["posterior_trees"].exists()
    assert report.output_paths["parameter_traces"].exists()
    assert report.output_paths["mcmc_diagnostics"].exists()
    assert report.output_paths["consensus_tree"].exists()
    assert (
        parse_mrbayes_parameter_traces(
            report.output_paths["parameter_traces"]
        ).row_count
        > 0
    )
    assert (
        parse_mrbayes_posterior_tree_samples(
            report.output_paths["posterior_trees"]
        ).tree_count
        > 0
    )
    assert (
        parse_mrbayes_mcmc_diagnostics(
            report.output_paths["mcmc_diagnostics"]
        ).row_count
        > 0
    )
    assert (
        parse_mrbayes_consensus_tree(report.output_paths["consensus_tree"])[
            1
        ].annotated_node_count
        > 0
    )


def test_run_mrbayes_posterior_inference_resumes_verified_real_outputs(
    tmp_path: Path,
) -> None:
    executable = real_mrbayes_executable()
    if executable is None:
        pytest.skip("real MrBayes executable is not available")

    nexus_path = tmp_path / "analysis.nex"
    prepare_mrbayes_analysis(
        fixture("alignments/example_alignment.fasta"),
        nexus_path,
    )

    first = run_mrbayes_posterior_inference(
        nexus_path,
        executable=executable,
        resume=False,
    )
    resumed = run_mrbayes_posterior_inference(
        nexus_path,
        executable=executable,
        resume=True,
    )

    assert first.output_paths["posterior_trees"].exists()
    assert first.output_paths["parameter_traces"].exists()
    assert first.output_paths["mcmc_diagnostics"].exists()
    assert first.output_paths["consensus_tree"].exists()
    assert resumed.resumed is True


def _run_small_beast_analysis(tmp_path: Path) -> Path:
    executable = real_beast_executable()
    if executable is None:
        pytest.skip("real BEAST executable is not available for integration coverage")

    output_path = tmp_path / "live-strict-yule.xml"
    prepare_beast_time_tree_analysis(
        fixture("example_alignment.fasta"),
        output_path,
        clock_model="strict",
        tree_prior="yule",
        chain_length=1000,
        log_every=20,
    )
    subprocess.run(
        [
            str(executable),
            "-overwrite",
            "-threads",
            "1",
            "-seed",
            "1",
            output_path.name,
        ],
        cwd=tmp_path,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    return output_path


def test_summarize_beast_log_with_real_executable_on_small_alignment(
    tmp_path: Path,
) -> None:
    output_path = _run_small_beast_analysis(tmp_path)
    log_path = output_path.with_suffix(".1.log")
    summary = summarize_beast_log(log_path, burnin_fraction=0.1)

    assert log_path.exists()
    assert summary.kept_row_count >= 40
    assert summary.posterior_parameters == ["posterior"]
    assert "clockRate" in summary.clock_parameters
    assert "birthRate" in summary.tree_parameters


def test_parse_beast_posterior_tree_samples_with_real_executable_on_small_alignment(
    tmp_path: Path,
) -> None:
    output_path = _run_small_beast_analysis(tmp_path)
    report = parse_beast_posterior_tree_samples(
        output_path.with_suffix(".1.trees"),
        burnin_fraction=0.1,
    )

    assert report.kept_tree_count >= 40
    assert report.sampled_states[0] == 100
    assert report.tip_names == ["A", "B", "C", "D"]
    assert report.clades


def test_summarize_beast_posterior_trees_with_real_executable_on_small_alignment(
    tmp_path: Path,
) -> None:
    output_path = _run_small_beast_analysis(tmp_path)
    consensus_tree, report = summarize_beast_posterior_trees(
        output_path.with_suffix(".1.trees"),
        burnin_fraction=0.1,
    )

    assert report.kept_tree_count >= 40
    assert report.annotated_node_count >= 0
    if report.maximum_posterior_probability is not None:
        assert 0.0 < report.maximum_posterior_probability <= 1.0
    assert consensus_tree.tip_names == ["A", "B", "C", "D"]


def test_summarize_beast_posterior_topology_diversity_with_real_executable_on_small_alignment(
    tmp_path: Path,
) -> None:
    output_path = _run_small_beast_analysis(tmp_path)
    report = summarize_beast_posterior_topology_diversity(
        output_path.with_suffix(".1.trees"),
        burnin_fraction=0.1,
    )

    assert report.kept_tree_count >= 40
    assert report.rooted_topology_count >= 1
    assert report.dominant_topology_frequency > 0.0
    assert report.effective_topology_count >= 1.0

@pytest.mark.slow
def test_run_beast_posterior_inference_resumes_verified_real_outputs(
    tmp_path: Path,
) -> None:
    executable = real_beast_executable()
    if executable is None:
        pytest.skip("real BEAST executable is not available for integration coverage")

    output_path = tmp_path / "live-strict-yule.xml"
    prepare_beast_time_tree_analysis(
        fixture("example_alignment.fasta"),
        output_path,
        clock_model="strict",
        tree_prior="yule",
        chain_length=1000,
        log_every=20,
    )

    first = run_beast_posterior_inference(
        output_path,
        executable=executable,
        threads=1,
        seed=1,
    )
    resumed = run_beast_posterior_inference(
        output_path,
        executable=executable,
        threads=1,
        seed=1,
        resume=True,
    )

    assert first.output_paths["posterior_log"].exists()
    assert first.output_paths["posterior_trees"].exists()
    assert resumed.resumed is True
