from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory

from bijux_phylogenetics.simulation import write_discrete_trait_table

from .comparison import (
    _build_parameter_comparison_rows,
    _build_parameter_rows,
    _build_rate_comparison_rows,
    _build_rate_rows,
    _selection_match,
)
from .fitting import (
    _build_bijux_fit_snapshot,
    _build_bijux_model_choice_rows,
    _build_execution_row,
    _build_geiger_fit_snapshot,
    _build_geiger_model_choice_rows,
    _build_model_comparison_execution_row,
    _geiger_payload_for_case,
    _simulate_case,
)
from .models import (
    DiscreteModeRecoveryCaseReport,
    DiscreteModeRecoveryReport,
    DiscreteModeRecoveryScenario,
)


def run_discrete_mode_recovery(
    tree_path: Path,
    scenarios: list[DiscreteModeRecoveryScenario],
    *,
    artifacts_root: Path | None = None,
) -> DiscreteModeRecoveryReport:
    """Simulate, refit, and compare governed discrete Mk recovery cases."""
    if artifacts_root is not None:
        artifacts_root.mkdir(parents=True, exist_ok=True)
        case_reports = [
            _run_case(
                default_tree_path=tree_path,
                scenario=scenario,
                working_root=artifacts_root,
                persist_traits=True,
            )
            for scenario in scenarios
        ]
        return DiscreteModeRecoveryReport(
            default_tree_path=tree_path,
            case_reports=case_reports,
        )
    with TemporaryDirectory(prefix="discrete-mode-recovery-") as temporary_root:
        temporary_path = Path(temporary_root)
        case_reports = [
            _run_case(
                default_tree_path=tree_path,
                scenario=scenario,
                working_root=temporary_path,
                persist_traits=False,
            )
            for scenario in scenarios
        ]
    return DiscreteModeRecoveryReport(
        default_tree_path=tree_path,
        case_reports=case_reports,
    )


def _run_case(
    *,
    default_tree_path: Path,
    scenario: DiscreteModeRecoveryScenario,
    working_root: Path,
    persist_traits: bool,
) -> DiscreteModeRecoveryCaseReport:
    case_tree_path = (
        default_tree_path if scenario.tree_path is None else scenario.tree_path
    )
    simulation = _simulate_case(case_tree_path, scenario)
    traits_path = write_discrete_trait_table(
        working_root / f"{scenario.case_id}-traits.tsv",
        simulation,
    )
    bijux_fit = _build_bijux_fit_snapshot(case_tree_path, traits_path, scenario)
    bijux_model_rows, bijux_selected_model = _build_bijux_model_choice_rows(
        case_tree_path,
        traits_path,
        scenario,
    )
    geiger_payload = _geiger_payload_for_case(scenario.case_id)
    geiger_fit = _build_geiger_fit_snapshot(
        scenario,
        geiger_payload["fit_summary"],
        geiger_payload["fit_rows"],
    )
    geiger_model_rows, geiger_selected_model = _build_geiger_model_choice_rows(
        scenario,
        geiger_payload["comparison_summary"],
        geiger_payload["comparison_rows"],
    )
    parameter_rows = _build_parameter_rows(
        scenario=scenario,
        bijux_fit=bijux_fit,
        geiger_fit=geiger_fit,
    )
    parameter_comparison_rows = _build_parameter_comparison_rows(parameter_rows)
    rate_rows = _build_rate_rows(
        scenario=scenario,
        bijux_fit=bijux_fit,
        geiger_fit=geiger_fit,
    )
    rate_comparison_rows = _build_rate_comparison_rows(rate_rows)
    warning_rows = [*bijux_fit.warning_rows, *geiger_fit.warning_rows]
    observed_warning_kinds = {
        row.kind for row in warning_rows if row.recovery_engine == "bijux"
    }
    overparameterized_review_matches_expectation = (
        bijux_fit.overparameterized == scenario.expected_overparameterized
        and geiger_fit.overparameterized == scenario.expected_overparameterized
    )
    execution_rows = [
        _build_execution_row(scenario.case_id, bijux_fit, operation="fit"),
        _build_execution_row(scenario.case_id, geiger_fit, operation="fit"),
        _build_model_comparison_execution_row(
            scenario.case_id,
            "bijux",
            bijux_selected_model,
            bijux_model_rows,
        ),
        _build_model_comparison_execution_row(
            scenario.case_id,
            "geiger",
            geiger_selected_model,
            geiger_model_rows,
        ),
    ]
    return DiscreteModeRecoveryCaseReport(
        scenario=scenario,
        tree_path=case_tree_path,
        traits_path=traits_path if persist_traits else None,
        simulation=simulation,
        parameter_rows=parameter_rows,
        parameter_comparison_rows=parameter_comparison_rows,
        rate_rows=rate_rows,
        rate_comparison_rows=rate_comparison_rows,
        model_choice_rows=[*bijux_model_rows, *geiger_model_rows],
        execution_rows=execution_rows,
        warning_rows=warning_rows,
        selected_model=bijux_selected_model,
        geiger_selected_model=geiger_selected_model,
        selection_matches_expectation=_selection_match(
            scenario.expected_selected_model,
            bijux_selected_model,
        ),
        geiger_selection_matches_expectation=_selection_match(
            scenario.expected_selected_model,
            geiger_selected_model,
        ),
        overparameterized_review_matches_expectation=overparameterized_review_matches_expectation,
        expected_warning_kinds_present=all(
            kind in observed_warning_kinds for kind in scenario.expected_warning_kinds
        ),
    )
