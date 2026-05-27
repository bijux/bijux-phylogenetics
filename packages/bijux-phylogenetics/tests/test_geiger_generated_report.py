from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.parity.geiger.boundary_warning_registry import (
    GeigerBoundaryWarningRow,
)
from bijux_phylogenetics.parity.geiger.generated_report import (
    _load_large_tree_benchmark_summary,
    _load_real_dataset_benchmark_summary,
    build_generated_geiger_parity_report,
    write_generated_geiger_parity_report_json,
    write_generated_geiger_parity_report_markdown,
)
from bijux_phylogenetics.parity.geiger.likelihood_policy import (
    GeigerLikelihoodPolicyRow,
)
from bijux_phylogenetics.parity.geiger.model_confidence import (
    GeigerModelConfidenceRow,
)
from bijux_phylogenetics.parity.geiger.optimizer_triage import (
    GeigerOptimizerTriageRow,
)
from bijux_phylogenetics.parity.geiger.parameterization_registry import (
    GeigerParameterizationRegistryRow,
)
from bijux_phylogenetics.parity.geiger.registry import list_geiger_parity_cases
from bijux_phylogenetics.parity.geiger.runner import (
    GeigerParityObservation,
    GeigerParityReport,
    GeigerParitySummaryRow,
)


def test_build_generated_geiger_parity_report_aggregates_live_and_artifact_evidence(
    monkeypatch,
    tmp_path: Path,
) -> None:
    def _mock_recovery_summary(path: Path) -> dict[str, int | str]:
        path_text = str(path)
        if "macroevolution_recovery_suite" in path_text:
            return {
                "dataset_id": "macroevolution_recovery_suite",
                "component_count": 3,
                "geiger_component_count": 2,
                "case_count": 11,
                "geiger_case_count": 11,
                "max_taxon_count": 24,
                "selection_review_case_count": 6,
                "selection_match_count": 6,
                "geiger_selection_match_count": 5,
                "governed_value_pass_count": 42,
                "governed_value_row_count": 46,
                "governed_comparison_row_count": 23,
                "expected_warning_case_count": 8,
                "expected_warning_present_count": 8,
                "truth_threshold_pass_count": 11,
                "truth_threshold_row_count": 11,
                "sim_char_case_count": 3,
                "requirement_pass_count": 11,
                "requirement_row_count": 11,
            }
        if "discrete" in path_text:
            return {
                "dataset_id": "discrete_mode_recovery_panel",
                "case_count": 4,
                "selection_review_case_count": 2,
                "selection_match_count": 2,
                "geiger_selection_match_count": 2,
                "parameter_pass_count": 0,
                "parameter_row_count": 0,
                "parameter_comparison_row_count": 0,
                "rate_pass_count": 20,
                "governed_rate_row_count": 24,
                "governed_rate_comparison_row_count": 12,
                "expected_warning_case_count": 2,
                "expected_warning_present_count": 2,
            }
        return {
            "dataset_id": "continuous_mode_recovery_panel",
            "case_count": 7,
            "selection_review_case_count": 4,
            "selection_match_count": 4,
            "geiger_selection_match_count": 3,
            "parameter_pass_count": 22,
            "parameter_row_count": 22,
            "parameter_comparison_row_count": 11,
            "rate_pass_count": 0,
            "governed_rate_row_count": 0,
            "governed_rate_comparison_row_count": 0,
            "expected_warning_case_count": 2,
            "expected_warning_present_count": 2,
        }

    monkeypatch.setattr(
        "bijux_phylogenetics.parity.geiger.generated_report._load_recovery_summary",
        _mock_recovery_summary,
    )
    monkeypatch.setattr(
        "bijux_phylogenetics.parity.geiger.generated_report._load_sim_char_summary",
        lambda: {"case_count": 3, "all_passed": True},
    )
    monkeypatch.setattr(
        "bijux_phylogenetics.parity.geiger.generated_report._load_large_tree_benchmark_summary",
        lambda: {
            "case_count": 2,
            "geiger_match_case_count": 2,
            "threshold_pass_case_count": 2,
            "too_slow_case_count": 0,
            "unstable_case_count": 0,
        },
    )
    monkeypatch.setattr(
        "bijux_phylogenetics.parity.geiger.generated_report._load_real_dataset_benchmark_summary",
        lambda: {
            "summary_row_count": 4,
            "model_row_count": 8,
            "alignment_review_row_count": 2,
            "parity_row_count": 10,
            "selection_match_count": 2,
            "unstable_review_count": 1,
        },
    )

    cases = list_geiger_parity_cases()
    summary_rows = [
        GeigerParitySummaryRow(
            function_name=function_name,
            case_count=sum(1 for case in cases if case.function_name == function_name),
            passed_case_count=sum(
                1 for case in cases if case.function_name == function_name
            ),
            failed_case_count=0,
            skipped_case_count=0,
        )
        for function_name in sorted({case.function_name for case in cases})
    ]
    live_report = GeigerParityReport(
        observations=[
            GeigerParityObservation(
                case_id="fitcontinuous-bm-example-tree",
                fixture_id="geiger_continuous_brownian_signal_twenty_four_taxa",
                function_name="geiger::fitContinuous(model='BM')",
                python_function_name="fit_continuous_evolutionary_mode",
                input_fixtures=(),
                model_name="BM",
                optimizer_settings=None,
                tolerance=1e-6,
                r_version="R version 4.4.1",
                geiger_version="2.0.11",
                bijux_version="0.1.0",
                bijux_commit="abc1234",
                status="passed",
                passed=True,
                mismatch_reason=None,
                reproducible_artifact_root=None,
                reference_summary={"log_likelihood": -1.0},
                bijux_summary={"log_likelihood": -1.0},
                reference_rows=None,
                bijux_rows=None,
                reference_error=None,
                bijux_error=None,
            )
        ],
        optimizer_triage_rows=[
            GeigerOptimizerTriageRow(
                case_id="fitcontinuous-lambda-weak-signal-review",
                function_name="geiger::fitContinuous(model='lambda')",
                model_name="lambda",
                status="passed",
                mismatch_reason=None,
                mismatch_type="boundary_solution_review",
                parameter_surface_comparable=True,
                same_likelihood_within_tolerance=True,
                same_parameter_surface_within_tolerance=True,
                reference_log_likelihood=-1.0,
                bijux_log_likelihood=-1.0,
                objective_delta=0.0,
                reference_parameter_name="lambda",
                bijux_parameter_name="lambda",
                reference_parameter_value=0.0,
                bijux_parameter_value=0.0,
                parameter_delta=0.0,
                reference_boundary_detected=True,
                bijux_boundary_detected=True,
                boundary_solution_detected=True,
                reference_trace_row_count=3,
                bijux_trace_row_count=3,
                reference_local_optimum_count=1,
                bijux_local_optimum_count=1,
                reference_optimizer_trace=None,
                bijux_optimizer_trace=None,
            )
        ],
        boundary_warning_rows=[
            GeigerBoundaryWarningRow(
                case_id="fitcontinuous-lambda-weak-signal-review",
                function_name="geiger::fitContinuous(model='lambda')",
                model_name="lambda",
                status="passed",
                scope="fitcontinuous-parameter-boundary-review",
                affected_parameter="lambda",
                lower_bound=0.0,
                upper_bound=1.0,
                reference_parameter_value=0.0,
                bijux_parameter_value=0.0,
                reference_hit_lower_boundary=True,
                reference_hit_upper_boundary=False,
                bijux_hit_lower_boundary=True,
                bijux_hit_upper_boundary=False,
                reference_near_lower_boundary=True,
                reference_near_upper_boundary=False,
                bijux_near_lower_boundary=True,
                bijux_near_upper_boundary=False,
                reference_flat_likelihood_near_boundary=True,
                bijux_flat_likelihood_near_boundary=True,
                reference_boundary_warning_kinds=["lower-boundary-hit"],
                bijux_boundary_warning_kinds=["lower-boundary-hit"],
                reference_boundary_dominates_interpretation=True,
                bijux_boundary_dominates_interpretation=True,
                reference_stable_conclusion_supported=False,
                bijux_stable_conclusion_supported=False,
                affected_parameter_match=True,
                stable_conclusion_supported_match=True,
                boundary_evidence="boundary review",
            )
        ],
        likelihood_policy_rows=[
            GeigerLikelihoodPolicyRow(
                case_id="fitcontinuous-bm-example-tree",
                function_name="geiger::fitContinuous(model='BM')",
                model_name="BM",
                status="passed",
                reference_likelihood_constant_policy="gaussian-full-constant",
                bijux_likelihood_constant_policy="gaussian-full-constant",
                case_level_raw_log_likelihood_comparable=True,
                raw_log_likelihood_match_within_tolerance=True,
                reference_aic_matches_raw_log_likelihood=True,
                bijux_aic_matches_raw_log_likelihood=True,
                relative_aic_comparable=True,
                ranking_permitted=True,
                ranking_guard_outcome="shared-fitcontinuous-policy-ranking-permitted",
                policy_evidence="policy review",
            )
        ],
        model_confidence_rows=[
            GeigerModelConfidenceRow(
                case_id="fitcontinuous-model-comparison-brownian-review",
                function_name="geiger::fitContinuous(model comparison)",
                model_name="brownian-review",
                status="passed",
                candidate_model="brownian",
                weight_basis="AICc",
                delta_threshold=2.0,
                reference_best_model="brownian",
                bijux_best_model="brownian",
                reference_rank=1,
                bijux_rank=1,
                reference_comparable=True,
                bijux_comparable=True,
                reference_delta_aic=0.0,
                bijux_delta_aic=0.0,
                reference_delta_aicc=0.0,
                bijux_delta_aicc=0.0,
                reference_akaike_weight=0.7,
                bijux_akaike_weight=0.7,
                akaike_weight_match_within_tolerance=True,
                reference_within_delta_aic_threshold=True,
                bijux_within_delta_aic_threshold=True,
                within_delta_aic_threshold_match=True,
                reference_within_delta_aicc_threshold=True,
                bijux_within_delta_aicc_threshold=True,
                within_delta_aicc_threshold_match=True,
                reference_selected=True,
                bijux_selected=True,
                selected_match=True,
                reference_uncertainty_class="supported",
                bijux_uncertainty_class="supported",
                uncertainty_class_match=True,
                reference_uncertainty_language="supported",
                bijux_uncertainty_language="supported",
                uncertainty_language_match=True,
                confidence_evidence="confidence review",
            )
        ],
        parameterization_registry_rows=[
            GeigerParameterizationRegistryRow(
                case_id="fitcontinuous-eb-early-burst-rate-recovery",
                function_name="geiger::fitContinuous(model='EB')",
                model_name="EB",
                reference_surface_contract="geiger::fitContinuous(model='EB')",
                bijux_surface_contract="bijux early-burst",
                status="passed",
                canonical_parameter_name="rate_change",
                reference_parameter_name="a",
                bijux_parameter_name="rate_change",
                reference_parameter_value=-1.0,
                bijux_parameter_value=1.0,
                converted_reference_parameter_value=1.0,
                converted_bijux_parameter_value=1.0,
                parameter_match_after_conversion=True,
                parameter_conversion_rule="rate_change=-a",
                reference_parameter_bounds={"lower": -10.0, "upper": 10.0},
                bijux_parameter_bounds={"lower": 0.0, "upper": 50.0},
                converted_reference_parameter_bounds={"lower": -10.0, "upper": 10.0},
                parameter_bounds_match_after_conversion=False,
                bounds_conversion_rule="sign-flip-only",
                root_state_parameterization_reference="ancestral-state-mean",
                root_state_parameterization_bijux="ancestral-state-mean",
                root_state_match_within_tolerance=True,
                variance_parameterization_reference="sigma-squared",
                variance_parameterization_bijux="sigma-squared",
                variance_match_within_tolerance=True,
                likelihood_constants_policy_reference="gaussian-full-constant",
                likelihood_constants_policy_bijux="gaussian-full-constant",
                likelihood_constants_comparison_policy="direct",
                log_likelihood_match_within_tolerance=True,
                expected_divergence=True,
                expected_divergence_kind="continuous-early-burst-sign-and-bound-convention",
                expected_divergence_evidence="expected",
            )
        ],
        summary_rows=summary_rows,
        case_count=len(cases),
        passed_case_count=len(cases),
        failed_case_count=0,
        skipped_case_count=0,
        all_passed=True,
        limitations=["live test limitation"],
    )

    report = build_generated_geiger_parity_report(parity_report=live_report)

    assert report.goal_start == 251
    assert report.goal_end == 289
    assert report.live_case_count == len(cases)
    assert report.r_version == "R version 4.4.1"
    assert report.geiger_version == "2.0.11"
    assert "geiger::fitContinuous(model='BM')" in report.covered_models
    assert any(row.goal_id == 278 for row in report.excluded_models)
    assert any(
        row.mismatch_type == "boundary_solution_review"
        for row in report.optimizer_mismatch_categories
    )
    assert any(
        row.panel_id == "continuous_mode_recovery_panel"
        and row.case_count == 7
        and row.geiger_selection_match_count == 3
        for row in report.simulation_recovery_rows
    )
    assert any(
        row.panel_id == "macroevolution_recovery_suite"
        and row.case_count == 11
        and row.governed_value_pass_count == 42
        for row in report.simulation_recovery_rows
    )
    assert any(
        row.benchmark_id == "real_dataset_macroevolution" and row.parity_row_count == 10
        for row in report.benchmark_rows
    )
    assert any(
        row.goal_id == 289 and row.status == "artifact-backed"
        for row in report.goal_coverage_rows
    )
    assert any(
        row.goal_id == 269
        and any("nondeterministic plateau lambda values" in note for note in row.notes)
        for row in report.goal_coverage_rows
    )
    assert any(
        row.surface == "geiger::fitDiscrete(model='ER', transform='lambda')"
        and any("nondeterministic plateau lambda value" in note for note in row.notes)
        for row in report.tolerance_rules
    )

    markdown_path = write_generated_geiger_parity_report_markdown(
        tmp_path / "geiger-parity-report.md",
        report,
    )
    json_path = write_generated_geiger_parity_report_json(
        tmp_path / "geiger-parity-report.json",
        report,
    )

    markdown_text = markdown_path.read_text(encoding="utf-8")
    json_payload = json.loads(json_path.read_text(encoding="utf-8"))

    assert "## Goal Coverage" in markdown_text
    assert "## Simulation Recovery" in markdown_text
    assert json_payload["goal_start"] == 251
    assert json_payload["sim_char_all_passed"] is True


def test_load_large_tree_benchmark_summary_reads_packaged_artifact() -> None:
    summary = _load_large_tree_benchmark_summary()

    assert summary == {
        "case_count": 2,
        "geiger_match_case_count": 2,
        "threshold_pass_case_count": 2,
        "too_slow_case_count": 0,
        "unstable_case_count": 0,
    }


def test_load_real_dataset_benchmark_summary_reads_packaged_artifact() -> None:
    summary = _load_real_dataset_benchmark_summary()

    assert summary == {
        "summary_row_count": 4,
        "model_row_count": 8,
        "alignment_review_row_count": 2,
        "parity_row_count": 10,
        "selection_match_count": 2,
        "unstable_review_count": 1,
    }
