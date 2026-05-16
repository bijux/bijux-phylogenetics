from __future__ import annotations

from pathlib import Path

import pytest

import bijux_phylogenetics
from bijux_phylogenetics.comparative.continuous_mode_recovery import (
    ContinuousModeRecoveryScenario,
    run_continuous_mode_recovery,
    write_continuous_mode_recovery_model_choice_table,
    write_continuous_mode_recovery_parameter_table,
    write_continuous_mode_recovery_summary_table,
    write_continuous_mode_recovery_warning_table,
)

_REFERENCE_TREE = (
    "((((A:0.8,B:0.4):0.7,(C:0.6,D:0.3):0.5):1.2,"
    "((E:0.7,F:0.35):0.8,(G:0.5,H:0.25):0.9):1.0):0.6,"
    "((I:0.9,J:0.45):0.6,(K:0.55,L:0.2):1.1):0.85);"
)


@pytest.mark.slow
def test_run_continuous_mode_recovery_validates_strong_and_weak_cases(
    tmp_path: Path,
) -> None:
    tree_path = _write_reference_tree(tmp_path)

    report = run_continuous_mode_recovery(tree_path, _scenarios())

    case_by_id = {case.scenario.case_id: case for case in report.case_reports}
    assert set(case_by_id) == {
        "brownian-sigma-recovery",
        "ou-parameter-recovery",
        "early-burst-rate-recovery",
        "weak-ou-identifiability",
    }
    assert case_by_id["brownian-sigma-recovery"].selected_model == "brownian"
    assert case_by_id["ou-parameter-recovery"].selected_model == ("ornstein-uhlenbeck")
    assert case_by_id["early-burst-rate-recovery"].selected_model == "early-burst"
    assert case_by_id["weak-ou-identifiability"].selected_model == "brownian"
    assert case_by_id["weak-ou-identifiability"].expected_warning_kinds_present is True
    assert all(
        row.within_tolerance
        for case_id in (
            "brownian-sigma-recovery",
            "ou-parameter-recovery",
            "early-burst-rate-recovery",
        )
        for row in case_by_id[case_id].parameter_rows
    )
    assert case_by_id["weak-ou-identifiability"].parameter_rows == []


@pytest.mark.slow
def test_continuous_mode_recovery_writers_emit_review_ledgers(
    tmp_path: Path,
) -> None:
    tree_path = _write_reference_tree(tmp_path)
    report = run_continuous_mode_recovery(tree_path, _scenarios())

    summary_path = write_continuous_mode_recovery_summary_table(
        tmp_path / "summary.tsv",
        report,
    )
    parameter_path = write_continuous_mode_recovery_parameter_table(
        tmp_path / "parameter-recovery.tsv",
        report,
    )
    model_choice_path = write_continuous_mode_recovery_model_choice_table(
        tmp_path / "model-choice.tsv",
        report,
    )
    warning_path = write_continuous_mode_recovery_warning_table(
        tmp_path / "warning-review.tsv",
        report,
    )

    summary_rows = summary_path.read_text(encoding="utf-8").splitlines()
    parameter_rows = parameter_path.read_text(encoding="utf-8").splitlines()
    model_choice_rows = model_choice_path.read_text(encoding="utf-8").splitlines()
    warning_rows = warning_path.read_text(encoding="utf-8").splitlines()

    assert summary_rows[0].startswith("case_id\tlabel\tgenerating_model")
    assert any(
        row.startswith("weak-ou-identifiability\tWeak OU identifiability review\t")
        for row in summary_rows[1:]
    )
    assert parameter_rows[0].startswith("case_id\tgenerating_model\tfitted_model")
    assert any(
        row.startswith(
            "early-burst-rate-recovery\tearly-burst\tearly-burst\trate_change\t"
        )
        for row in parameter_rows[1:]
    )
    assert model_choice_rows[0].startswith("case_id\tgenerating_model")
    assert any(
        row.startswith(
            "ou-parameter-recovery\tornstein-uhlenbeck\tornstein-uhlenbeck\tornstein-uhlenbeck\t"
        )
        for row in model_choice_rows[1:]
    )
    assert warning_rows[0] == "case_id\tfitted_model\tkind\tmessage"
    assert any(
        row.startswith("weak-ou-identifiability\tornstein-uhlenbeck\tflat_likelihood\t")
        for row in warning_rows[1:]
    )


def test_public_runtime_exports_include_continuous_mode_recovery_surface() -> None:
    assert (
        bijux_phylogenetics.run_continuous_mode_recovery is run_continuous_mode_recovery
    )
    assert (
        bijux_phylogenetics.write_continuous_mode_recovery_summary_table
        is write_continuous_mode_recovery_summary_table
    )
    assert (
        bijux_phylogenetics.write_continuous_mode_recovery_parameter_table
        is write_continuous_mode_recovery_parameter_table
    )
    assert (
        bijux_phylogenetics.write_continuous_mode_recovery_model_choice_table
        is write_continuous_mode_recovery_model_choice_table
    )
    assert (
        bijux_phylogenetics.write_continuous_mode_recovery_warning_table
        is write_continuous_mode_recovery_warning_table
    )


def _write_reference_tree(tmp_path: Path) -> Path:
    tree_path = tmp_path / "reference-tree.nwk"
    tree_path.write_text(_REFERENCE_TREE, encoding="utf-8")
    return tree_path


def _scenarios() -> list[ContinuousModeRecoveryScenario]:
    return [
        ContinuousModeRecoveryScenario(
            case_id="brownian-sigma-recovery",
            label="Brownian sigma recovery",
            generating_model="brownian",
            expected_selected_model="brownian",
            root_state=2.5,
            sigma=0.8,
            seed=18,
            parameter_tolerances={"sigma_squared": 0.1},
            notes="Brownian simulation should recover sigma squared and remain preferred by AICc.",
        ),
        ContinuousModeRecoveryScenario(
            case_id="ou-parameter-recovery",
            label="OU parameter recovery",
            generating_model="ornstein-uhlenbeck",
            expected_selected_model="ornstein-uhlenbeck",
            root_state=2.0,
            sigma=0.7,
            alpha=1.0,
            theta=2.8,
            seed=2,
            parameter_tolerances={
                "alpha": 0.6,
                "sigma_squared": 0.25,
                "theta": 0.5,
            },
            notes="OU simulation should recover alpha, sigma squared, and theta while remaining preferred by AICc.",
        ),
        ContinuousModeRecoveryScenario(
            case_id="early-burst-rate-recovery",
            label="Early-burst rate recovery",
            generating_model="early-burst",
            expected_selected_model="early-burst",
            root_state=3.0,
            sigma=0.8,
            rate_change=4.0,
            seed=4,
            early_burst_bounds=(0.0, 10.0),
            parameter_tolerances={"rate_change": 0.25},
            notes="Early-burst simulation should recover a positive rate-change parameter and remain preferred by AICc.",
        ),
        ContinuousModeRecoveryScenario(
            case_id="weak-ou-identifiability",
            label="Weak OU identifiability review",
            generating_model="ornstein-uhlenbeck",
            expected_selected_model="brownian",
            root_state=2.0,
            sigma=0.3,
            alpha=0.05,
            theta=2.1,
            seed=17,
            expected_warning_kinds=["flat_likelihood"],
            notes="Weak OU pull should trigger identifiability warnings and stay close to Brownian support.",
        ),
    ]
