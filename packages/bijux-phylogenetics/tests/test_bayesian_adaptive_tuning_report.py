from __future__ import annotations

import pytest

from bijux_phylogenetics.bayesian.adaptive_tuning import (
    AdaptiveTuningReport,
    AdaptiveTuningWindowRow,
    build_adaptive_tuning_controller,
    build_adaptive_tuning_report,
    build_adaptive_tuning_window_row,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError


def test_build_adaptive_tuning_window_row_records_acceptance_summary() -> None:
    row = build_adaptive_tuning_window_row(
        window_index=1,
        window_start_iteration=1,
        window_end_iteration=4,
        within_burnin=True,
        attempted_count=4,
        accepted_count=1,
        target_acceptance_rate=0.3,
        scale_before_window=2.0,
        scale_after_window=1.0,
        action="decrease",
    )

    assert isinstance(row, AdaptiveTuningWindowRow)
    assert row.acceptance_rate == 0.25
    assert row.action == "decrease"


def test_build_adaptive_tuning_report_records_freeze_point_and_final_scale() -> None:
    controller = build_adaptive_tuning_controller(
        proposal_name="normal-random-walk",
        scale_parameter_name="step-size",
        initial_scale=4.0,
        target_acceptance_rate=0.3,
        burnin_iteration_count=6,
        adaptation_window_size=2,
    )
    report = build_adaptive_tuning_report(
        controller=controller,
        freeze_iteration_index=7,
        burnin_sample_count=3,
        retained_sample_count=2,
        window_rows=[
            build_adaptive_tuning_window_row(
                window_index=1,
                window_start_iteration=1,
                window_end_iteration=2,
                within_burnin=True,
                attempted_count=2,
                accepted_count=0,
                target_acceptance_rate=0.3,
                scale_before_window=4.0,
                scale_after_window=2.0,
                action="decrease",
            ),
            build_adaptive_tuning_window_row(
                window_index=2,
                window_start_iteration=3,
                window_end_iteration=4,
                within_burnin=True,
                attempted_count=2,
                accepted_count=0,
                target_acceptance_rate=0.3,
                scale_before_window=2.0,
                scale_after_window=1.0,
                action="decrease",
            ),
            build_adaptive_tuning_window_row(
                window_index=3,
                window_start_iteration=5,
                window_end_iteration=6,
                within_burnin=True,
                attempted_count=2,
                accepted_count=2,
                target_acceptance_rate=0.3,
                scale_before_window=1.0,
                scale_after_window=2.0,
                action="increase",
            ),
            build_adaptive_tuning_window_row(
                window_index=4,
                window_start_iteration=7,
                window_end_iteration=8,
                within_burnin=False,
                attempted_count=2,
                accepted_count=0,
                target_acceptance_rate=0.3,
                scale_before_window=2.0,
                scale_after_window=2.0,
                action="frozen",
            ),
        ],
    )

    assert isinstance(report, AdaptiveTuningReport)
    assert report.freeze_iteration_index == 7
    assert report.final_scale == 2.0
    assert report.window_rows[-1].action == "frozen"


def test_build_adaptive_tuning_report_rejects_nonfrozen_post_burnin_window() -> None:
    controller = build_adaptive_tuning_controller(
        proposal_name="normal-random-walk",
        scale_parameter_name="step-size",
        initial_scale=4.0,
        target_acceptance_rate=0.3,
        burnin_iteration_count=2,
        adaptation_window_size=2,
    )

    with pytest.raises(
        PhylogeneticsError,
        match="post-burn-in windows to record the frozen action",
    ):
        build_adaptive_tuning_report(
            controller=controller,
            freeze_iteration_index=3,
            burnin_sample_count=1,
            retained_sample_count=1,
            window_rows=[
                build_adaptive_tuning_window_row(
                    window_index=1,
                    window_start_iteration=1,
                    window_end_iteration=2,
                    within_burnin=True,
                    attempted_count=2,
                    accepted_count=0,
                    target_acceptance_rate=0.3,
                    scale_before_window=4.0,
                    scale_after_window=2.0,
                    action="decrease",
                ),
                build_adaptive_tuning_window_row(
                    window_index=2,
                    window_start_iteration=3,
                    window_end_iteration=4,
                    within_burnin=False,
                    attempted_count=2,
                    accepted_count=2,
                    target_acceptance_rate=0.3,
                    scale_before_window=2.0,
                    scale_after_window=4.0,
                    action="increase",
                ),
            ],
        )
