from __future__ import annotations

import csv
import json
from pathlib import Path

from bijux_phylogenetics.phytools_parity import (
    _load_rows_table,
    list_phytools_parity_cases,
    run_phytools_parity_cases,
    write_phytools_parity_observation_table,
    write_phytools_parity_summary_table,
)
from tests.support.fake_phytools_parity import fake_phytools_rscript


def test_list_phytools_parity_cases_returns_governed_registry() -> None:
    cases = list_phytools_parity_cases()

    assert [case.case_id for case in cases] == [
        "phylosig-lambda-non-ultrametric-strong-signal-twenty-four-taxa",
        "phylosig-lambda-weak-signal-twenty-four-taxa",
        "phylosig-k-strong-signal-twenty-four-taxa",
        "phylosig-k-weak-signal-twenty-four-taxa",
        "fitmk-er-binary-twenty-four-taxa",
        "fitmk-er-multistate-twenty-four-taxa",
        "fitmk-er-binary-missing-twenty-four-taxa",
        "fitmk-er-multistate-missing-twenty-four-taxa",
        "fitmk-sym-multistate-twenty-four-taxa",
        "fitmk-sym-multistate-missing-twenty-four-taxa",
        "fitmk-ard-binary-twenty-four-taxa",
        "fitmk-ard-multistate-twenty-four-taxa",
        "fitmk-ard-binary-missing-twenty-four-taxa",
        "fitmk-ard-multistate-missing-twenty-four-taxa",
        "simmap-er-binary-twenty-four-taxa",
        "simmap-er-multistate-twenty-four-taxa",
        "simmap-er-binary-missing-twenty-four-taxa",
        "simmap-sym-multistate-twenty-four-taxa",
        "simmap-sym-multistate-missing-twenty-four-taxa",
        "simmap-ard-binary-twenty-four-taxa",
        "simmap-ard-multistate-twenty-four-taxa",
        "simmap-ard-binary-missing-twenty-four-taxa",
        "simmap-ard-multistate-missing-twenty-four-taxa",
        "count-simmap-er-binary-twenty-four-taxa",
        "count-simmap-er-multistate-twenty-four-taxa",
        "count-simmap-sym-multistate-twenty-four-taxa",
        "count-simmap-er-binary-missing-twenty-four-taxa",
        "describe-simmap-er-binary-twenty-four-taxa",
        "describe-simmap-er-multistate-twenty-four-taxa",
        "describe-simmap-sym-multistate-twenty-four-taxa",
        "describe-simmap-er-binary-missing-twenty-four-taxa",
        "density-map-er-binary-twenty-four-taxa",
        "density-map-er-binary-missing-twenty-four-taxa",
        "sim-history-binary-no-change-example-tree",
        "sim-history-binary-high-rate-example-tree",
        "sim-history-multistate-no-change-six-taxa",
        "sim-history-multistate-high-rate-six-taxa",
        "rerooting-er-binary-twenty-four-taxa",
        "rerooting-er-multistate-twenty-four-taxa",
        "rerooting-er-binary-missing-twenty-four-taxa",
        "rerooting-sym-multistate-twenty-four-taxa",
        "rerooting-sym-multistate-missing-twenty-four-taxa",
        "fast-anc-strong-signal-twenty-four-taxa",
        "fast-anc-weak-signal-twenty-four-taxa",
        "fast-anc-non-ultrametric-strong-signal-twenty-four-taxa",
        "fast-anc-missing-values-twenty-four-taxa",
        "anc-ml-strong-signal-twenty-four-taxa",
        "anc-ml-weak-signal-twenty-four-taxa",
        "anc-ml-non-ultrametric-strong-signal-twenty-four-taxa",
        "anc-ml-missing-values-twenty-four-taxa",
        "fastbm-example-tree-low-variance",
        "fastbm-example-tree-root-shift-high-variance",
        "fastbm-six-taxa-root-shift",
        "simcorrs-example-tree-low-correlation",
        "simcorrs-example-tree-negative-correlation-root-shift",
        "simcorrs-six-taxa-three-trait",
    ]
    assert cases[0].function_name == "phytools::phylosig(method='lambda')"
    assert cases[1].function_name == "phytools::phylosig(method='lambda')"
    assert cases[2].function_name == "phytools::phylosig(method='K')"
    assert cases[3].function_name == "phytools::phylosig(method='K')"
    assert cases[4].function_name == "phytools::fitMk(model='ER')"
    assert cases[7].function_name == "phytools::fitMk(model='ER')"
    assert cases[8].function_name == "phytools::fitMk(model='SYM')"
    assert cases[9].function_name == "phytools::fitMk(model='SYM')"
    assert cases[10].function_name == "phytools::fitMk(model='ARD')"
    assert cases[13].function_name == "phytools::fitMk(model='ARD')"
    assert cases[14].function_name == "phytools::make.simmap(model='ER')"
    assert cases[16].function_name == "phytools::make.simmap(model='ER')"
    assert cases[17].function_name == "phytools::make.simmap(model='SYM')"
    assert cases[18].function_name == "phytools::make.simmap(model='SYM')"
    assert cases[19].function_name == "phytools::make.simmap(model='ARD')"
    assert cases[22].function_name == "phytools::make.simmap(model='ARD')"
    assert cases[23].function_name == "phytools::countSimmap"
    assert cases[26].function_name == "phytools::countSimmap"
    assert cases[27].function_name == "phytools::describe.simmap"
    assert cases[30].function_name == "phytools::describe.simmap"
    assert cases[31].function_name == "phytools::densityMap"
    assert cases[32].function_name == "phytools::densityMap"
    assert cases[33].function_name == "phytools::sim.history"
    assert cases[36].function_name == "phytools::sim.history"
    assert cases[37].function_name == "phytools::rerootingMethod"
    assert cases[41].function_name == "phytools::rerootingMethod"
    assert cases[42].function_name == "phytools::fastAnc"
    assert cases[45].function_name == "phytools::fastAnc"
    assert cases[46].function_name == "phytools::anc.ML"
    assert cases[49].function_name == "phytools::anc.ML"
    assert cases[50].function_name == "phytools::fastBM"
    assert cases[52].function_name == "phytools::fastBM"
    assert cases[53].function_name == "phytools::sim.corrs"
    assert cases[55].function_name == "phytools::sim.corrs"
    assert (
        cases[0].fixture_id == "phytools_continuous_strong_signal_non_ultrametric_twenty_four_taxa"
    )
    assert cases[2].permutation_count == 199
    assert cases[3].permutation_seed == 17
    assert cases[7].row_field_tolerances == {"rate": 1e-5}
    assert cases[8].row_field_tolerances == {"rate": 1e-4}
    assert cases[10].row_field_tolerances == {"rate": 1e-3}
    assert cases[14].stochastic_map_replicate_count == 128
    assert cases[16].stochastic_map_seed == 17
    assert cases[17].row_field_tolerances == {
        "mean_value": 2.5,
        "lower_95_interval": 5.0,
        "upper_95_interval": 5.0,
        "presence_fraction": 0.25,
    }
    assert cases[20].compare_rows is False
    assert cases[22].compare_rows is False
    assert cases[23].row_field_tolerances == {
        "mean_value": 1.5,
        "lower_95_interval": 3.0,
        "upper_95_interval": 3.0,
        "presence_fraction": 0.2,
    }
    assert cases[25].row_field_tolerances == {
        "mean_value": 2.5,
        "lower_95_interval": 5.0,
        "upper_95_interval": 5.0,
        "presence_fraction": 0.25,
    }
    assert cases[27].row_field_tolerances == {
        "mean_value": 1.5,
        "lower_95_interval": 3.0,
        "upper_95_interval": 3.0,
        "presence_fraction": 0.2,
    }
    assert cases[29].row_field_tolerances == {
        "mean_value": 2.5,
        "lower_95_interval": 5.0,
        "upper_95_interval": 5.0,
        "presence_fraction": 0.25,
    }
    assert cases[31].row_field_tolerances == {
        "mean_posterior_probability": 2e-3,
        "minimum_posterior_probability": 2e-3,
        "maximum_posterior_probability": 2e-3,
        "uncertainty": 2e-3,
    }
    assert cases[33].row_field_tolerances == {
        "mean_value": 0.25,
        "lower_95_interval": 0.35,
        "upper_95_interval": 0.35,
        "presence_fraction": 0.2,
    }
    assert cases[36].row_field_tolerances == {
        "mean_value": 1.5,
        "lower_95_interval": 3.0,
        "upper_95_interval": 3.0,
        "presence_fraction": 0.2,
    }
    assert cases[37].row_field_tolerances == {"probability": 1e-5}
    assert cases[40].row_field_tolerances == {"probability": 5e-5}
    assert cases[41].row_field_tolerances == {"probability": 5e-5}
    assert cases[49].row_field_tolerances == {
        "estimate": 1e-8,
        "standard_error": 5e-8,
        "lower_95_interval": 5e-8,
        "upper_95_interval": 5e-8,
    }
    assert cases[45].row_field_tolerances == {
        "estimate": 1e-8,
        "standard_error": 1e-8,
    }
    assert cases[50].continuous_sigma_squared == 0.25
    assert cases[51].continuous_root_state == 2.5
    assert cases[52].continuous_replicate_count == 512
    assert cases[53].continuous_trait_names == ("trait_alpha", "trait_beta")
    assert cases[54].continuous_root_states == (2.0, -1.0)
    assert cases[55].continuous_replicate_count == 512


def test_run_phytools_parity_cases_passes_against_fake_reference_runner(
    tmp_path: Path,
) -> None:
    rscript = fake_phytools_rscript(tmp_path / "fake-phytools-rscript")

    report = run_phytools_parity_cases(rscript_executable=str(rscript))

    assert report.all_passed is True
    assert report.case_count == 56
    assert report.failed_case_count == 0
    assert report.skipped_case_count == 0
    assert [row.function_name for row in report.summary_rows] == [
        "phytools::anc.ML",
        "phytools::countSimmap",
        "phytools::densityMap",
        "phytools::describe.simmap",
        "phytools::fastAnc",
        "phytools::fastBM",
        "phytools::fitMk(model='ARD')",
        "phytools::fitMk(model='ER')",
        "phytools::fitMk(model='SYM')",
        "phytools::make.simmap(model='ARD')",
        "phytools::make.simmap(model='ER')",
        "phytools::make.simmap(model='SYM')",
        "phytools::phylosig(method='K')",
        "phytools::phylosig(method='lambda')",
        "phytools::rerootingMethod",
        "phytools::sim.corrs",
        "phytools::sim.history",
    ]
    first = report.observations[0]
    assert first.phytools_version == "2.5.2"
    assert first.r_version == "4.6.0"
    assert first.reference_summary is not None


def test_run_phytools_parity_cases_passes_density_map_cases_against_fake_reference_runner(
    tmp_path: Path,
) -> None:
    rscript = fake_phytools_rscript(tmp_path / "fake-phytools-rscript")

    report = run_phytools_parity_cases(
        case_ids=[
            "density-map-er-binary-twenty-four-taxa",
            "density-map-er-binary-missing-twenty-four-taxa",
        ],
        rscript_executable=str(rscript),
    )

    assert report.all_passed is True
    assert report.case_count == 2
    assert report.failed_case_count == 0
    assert [row.function_name for row in report.summary_rows] == [
        "phytools::densityMap"
    ]
    assert all(
        observation.function_name == "phytools::densityMap"
        for observation in report.observations
    )


def test_run_phytools_parity_cases_passes_sim_history_cases_against_fake_reference_runner(
    tmp_path: Path,
) -> None:
    rscript = fake_phytools_rscript(tmp_path / "fake-phytools-rscript")

    report = run_phytools_parity_cases(
        case_ids=[
            "sim-history-binary-no-change-example-tree",
            "sim-history-binary-high-rate-example-tree",
            "sim-history-multistate-no-change-six-taxa",
            "sim-history-multistate-high-rate-six-taxa",
        ],
        rscript_executable=str(rscript),
    )

    assert report.all_passed is True
    assert report.case_count == 4
    assert report.failed_case_count == 0
    assert [row.function_name for row in report.summary_rows] == [
        "phytools::sim.history"
    ]
    assert all(
        observation.function_name == "phytools::sim.history"
        for observation in report.observations
    )


def test_run_phytools_parity_cases_passes_fastbm_cases_against_fake_reference_runner(
    tmp_path: Path,
) -> None:
    rscript = fake_phytools_rscript(tmp_path / "fake-phytools-rscript")

    report = run_phytools_parity_cases(
        case_ids=[
            "fastbm-example-tree-low-variance",
            "fastbm-example-tree-root-shift-high-variance",
            "fastbm-six-taxa-root-shift",
        ],
        rscript_executable=str(rscript),
    )

    assert report.all_passed is True
    assert report.case_count == 3
    assert report.failed_case_count == 0
    assert [row.function_name for row in report.summary_rows] == ["phytools::fastBM"]
    assert all(
        observation.function_name == "phytools::fastBM"
        for observation in report.observations
    )


def test_run_phytools_parity_cases_passes_simcorrs_cases_against_fake_reference_runner(
    tmp_path: Path,
) -> None:
    rscript = fake_phytools_rscript(tmp_path / "fake-phytools-rscript")

    report = run_phytools_parity_cases(
        case_ids=[
            "simcorrs-example-tree-low-correlation",
            "simcorrs-example-tree-negative-correlation-root-shift",
            "simcorrs-six-taxa-three-trait",
        ],
        rscript_executable=str(rscript),
    )

    assert report.all_passed is True
    assert report.case_count == 3
    assert report.failed_case_count == 0
    assert [row.function_name for row in report.summary_rows] == [
        "phytools::sim.corrs"
    ]
    assert all(
        observation.function_name == "phytools::sim.corrs"
        for observation in report.observations
    )


def test_run_phytools_parity_cases_records_failure_artifacts(
    tmp_path: Path,
) -> None:
    rscript = fake_phytools_rscript(
        tmp_path / "fake-phytools-rscript",
        summary_overrides={
            "phylosig-k-strong-signal-twenty-four-taxa": {
                "taxon_count": 24,
                "trait_name": "signal_strong",
                "k": 0.5,
                "p_value": 0.005025125628140704,
                "permutation_count": 199,
                "permutation_seed": 17,
                "simulated_k_minimum": 0.004227447597570447,
                "simulated_k_mean": 0.029614826253456215,
            }
        },
    )

    report = run_phytools_parity_cases(
        case_ids=["phylosig-k-strong-signal-twenty-four-taxa"],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "failures",
    )

    assert report.all_passed is False
    assert report.failed_case_count == 1
    observation = report.observations[0]
    assert observation.status == "failed"
    assert observation.mismatch_reason == "summary_mismatch:k"
    assert observation.reproducible_artifact_root is not None
    assert (observation.reproducible_artifact_root / "reference-summary.json").exists()
    assert (observation.reproducible_artifact_root / "bijux-summary.json").exists()


def test_run_phytools_parity_cases_marks_missing_rscript_as_skipped(
    tmp_path: Path,
) -> None:
    report = run_phytools_parity_cases(
        case_ids=["phylosig-lambda-non-ultrametric-strong-signal-twenty-four-taxa"],
        rscript_executable=str(tmp_path / "missing-rscript"),
        failure_root=tmp_path / "failures",
    )

    assert report.all_passed is False
    assert report.skipped_case_count == 1
    observation = report.observations[0]
    assert observation.status == "skipped"
    assert observation.mismatch_reason == "rscript_unavailable"
    assert observation.reproducible_artifact_root is not None


def test_load_rows_table_preserves_discrete_state_identity_and_booleans(
    tmp_path: Path,
) -> None:
    rows_path = tmp_path / "fitmk-rate-matrix.tsv"
    rows_path.write_text(
        "\n".join(
            [
                "source_state\ttarget_state\ttransition_allowed\tstep_distance\trate",
                "0\t1\tTrue\t1\t70.77642178117888",
                "north\tsouth\tfalse\t2\t0.15",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    rows = _load_rows_table(rows_path)

    assert rows == [
        {
            "source_state": "0",
            "target_state": "1",
            "transition_allowed": True,
            "step_distance": 1,
            "rate": 70.77642178117888,
        },
        {
            "source_state": "north",
            "target_state": "south",
            "transition_allowed": False,
            "step_distance": 2,
            "rate": 0.15,
        },
    ]


def test_load_rows_table_preserves_stochastic_branch_summary_identity(
    tmp_path: Path,
) -> None:
    rows_path = tmp_path / "stochastic-map-summary-rows.tsv"
    rows_path.write_text(
        "\n".join(
            [
                "row_kind\tlabel\tmean_value\tlower_95_interval\tupper_95_interval\tpresence_fraction",
                "branch_state_occupancy\tA|B->A:0\t1.0\t1.0\t1.0\t1.0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    rows = _load_rows_table(rows_path)

    assert rows == [
        {
            "row_kind": "branch_state_occupancy",
            "label": "A|B->A:0",
            "mean_value": 1.0,
            "lower_95_interval": 1.0,
            "upper_95_interval": 1.0,
            "presence_fraction": 1.0,
        }
    ]


def test_load_rows_table_preserves_rerooting_state_identity(
    tmp_path: Path,
) -> None:
    rows_path = tmp_path / "rerooting-method-node-probabilities.tsv"
    rows_path.write_text(
        "\n".join(
            [
                "node\tstate\tprobability",
                "A|B|C\t0\t0.75",
                "A|B|C\tnorth\t0.25",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    rows = _load_rows_table(rows_path)

    assert rows == [
        {
            "node": "A|B|C",
            "state": "0",
            "probability": 0.75,
        },
        {
            "node": "A|B|C",
            "state": "north",
            "probability": 0.25,
        },
    ]


def test_load_rows_table_preserves_stochastic_map_row_identity(
    tmp_path: Path,
) -> None:
    rows_path = tmp_path / "stochastic-map-summary-rows.tsv"
    rows_path.write_text(
        "\n".join(
            [
                "row_kind\tlabel\tmean_value\tlower_95_interval\tupper_95_interval\tpresence_fraction",
                "transition_count\tnorth->south\t1.5\t0\t4\t0.75",
                "state_time\t0\t2.25\t1\t3\t1.0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    rows = _load_rows_table(rows_path)

    assert rows == [
        {
            "row_kind": "transition_count",
            "label": "north->south",
            "mean_value": 1.5,
            "lower_95_interval": 0,
            "upper_95_interval": 4,
            "presence_fraction": 0.75,
        },
        {
            "row_kind": "state_time",
            "label": "0",
            "mean_value": 2.25,
            "lower_95_interval": 1,
            "upper_95_interval": 3,
            "presence_fraction": 1.0,
        },
    ]


def test_write_phytools_parity_tables_writes_summary_and_observations(
    tmp_path: Path,
) -> None:
    rscript = fake_phytools_rscript(tmp_path / "fake-phytools-rscript")
    report = run_phytools_parity_cases(rscript_executable=str(rscript))
    summary_path = tmp_path / "phytools-parity-summary.tsv"
    observation_path = tmp_path / "phytools-parity-observations.tsv"

    write_phytools_parity_summary_table(summary_path, report)
    write_phytools_parity_observation_table(observation_path, report)

    summary_rows = summary_path.read_text(encoding="utf-8").splitlines()
    assert summary_rows[0].startswith("function_name\tcase_count")
    with observation_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows
    input_fixtures = json.loads(rows[0]["input_fixtures"])
    assert isinstance(input_fixtures, list)
    assert rows[0]["phytools_version"] == "2.5.2"
