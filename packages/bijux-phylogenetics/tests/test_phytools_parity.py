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
        "fast-anc-strong-signal-twenty-four-taxa",
        "fast-anc-weak-signal-twenty-four-taxa",
        "fast-anc-non-ultrametric-strong-signal-twenty-four-taxa",
        "fast-anc-missing-values-twenty-four-taxa",
        "anc-ml-strong-signal-twenty-four-taxa",
        "anc-ml-weak-signal-twenty-four-taxa",
        "anc-ml-non-ultrametric-strong-signal-twenty-four-taxa",
        "anc-ml-missing-values-twenty-four-taxa",
    ]
    assert cases[0].function_name == "phytools::phylosig(method='lambda')"
    assert cases[1].function_name == "phytools::phylosig(method='lambda')"
    assert cases[2].function_name == "phytools::phylosig(method='K')"
    assert cases[3].function_name == "phytools::phylosig(method='K')"
    assert cases[4].function_name == "phytools::fitMk(model='ER')"
    assert cases[7].function_name == "phytools::fitMk(model='ER')"
    assert cases[8].function_name == "phytools::fastAnc"
    assert cases[11].function_name == "phytools::fastAnc"
    assert cases[12].function_name == "phytools::anc.ML"
    assert cases[15].function_name == "phytools::anc.ML"
    assert (
        cases[0].fixture_id == "phytools_continuous_strong_signal_non_ultrametric_twenty_four_taxa"
    )
    assert cases[2].permutation_count == 199
    assert cases[3].permutation_seed == 17
    assert cases[7].row_field_tolerances == {"rate": 1e-5}
    assert cases[15].row_field_tolerances == {
        "estimate": 1e-8,
        "standard_error": 5e-8,
        "lower_95_interval": 5e-8,
        "upper_95_interval": 5e-8,
    }
    assert cases[11].row_field_tolerances == {
        "estimate": 1e-8,
        "standard_error": 1e-8,
    }


def test_run_phytools_parity_cases_passes_against_fake_reference_runner(
    tmp_path: Path,
) -> None:
    rscript = fake_phytools_rscript(tmp_path / "fake-phytools-rscript")

    report = run_phytools_parity_cases(rscript_executable=str(rscript))

    assert report.all_passed is True
    assert report.case_count == 16
    assert report.failed_case_count == 0
    assert report.skipped_case_count == 0
    assert [row.function_name for row in report.summary_rows] == [
        "phytools::anc.ML",
        "phytools::fastAnc",
        "phytools::fitMk(model='ER')",
        "phytools::phylosig(method='K')",
        "phytools::phylosig(method='lambda')",
    ]
    first = report.observations[0]
    assert first.phytools_version == "2.5.2"
    assert first.r_version == "4.6.0"
    assert first.reference_summary is not None


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
