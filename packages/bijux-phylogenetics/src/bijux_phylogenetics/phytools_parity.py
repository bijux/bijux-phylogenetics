from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
import json
from importlib import metadata
import math
import os
from pathlib import Path
import subprocess
import tempfile

from bijux_phylogenetics.ancestral.continuous import (
    reconstruct_continuous_ancestral_states,
)
from bijux_phylogenetics.comparative.discrete_mk import fit_discrete_mk_model
from bijux_phylogenetics.comparative.signal import (
    compute_blombergs_k,
    compute_phylogenetic_signal_test,
    estimate_pagels_lambda,
)
from bijux_phylogenetics.shared_phytools_comparative_fixtures import (
    get_shared_phytools_comparative_fixture,
)


@dataclass(frozen=True, slots=True)
class PhytoolsParityCase:
    """One governed live `phytools` parity case."""

    case_id: str
    fixture_id: str
    function_name: str
    python_function_name: str
    operation: str
    input_fixtures: tuple[Path, ...]
    tolerance: float
    trait_name: str
    taxon_column: str | None = None
    discrete_model: str | None = None
    permutation_count: int | None = None
    permutation_seed: int | None = None
    field_tolerances: dict[str, float] | None = None
    row_field_tolerances: dict[str, float] | None = None


@dataclass(frozen=True, slots=True)
class PhytoolsParityObservation:
    """One live parity comparison between Bijux and `phytools`."""

    case_id: str
    fixture_id: str
    function_name: str
    python_function_name: str
    input_fixtures: tuple[Path, ...]
    tolerance: float
    r_version: str | None
    phytools_version: str | None
    bijux_version: str
    bijux_commit: str | None
    status: str
    passed: bool
    mismatch_reason: str | None
    reproducible_artifact_root: Path | None
    reference_summary: dict[str, object] | None
    bijux_summary: dict[str, object] | None
    reference_rows: list[dict[str, object]] | None
    bijux_rows: list[dict[str, object]] | None
    reference_error: dict[str, object] | None
    bijux_error: dict[str, object] | None


@dataclass(frozen=True, slots=True)
class PhytoolsParitySummaryRow:
    """One function-level summary across governed `phytools` parity cases."""

    function_name: str
    case_count: int
    passed_case_count: int
    failed_case_count: int
    skipped_case_count: int


@dataclass(slots=True)
class PhytoolsParityReport:
    """Aggregate report for governed live `phytools` parity cases."""

    observations: list[PhytoolsParityObservation]
    summary_rows: list[PhytoolsParitySummaryRow]
    case_count: int
    passed_case_count: int
    failed_case_count: int
    skipped_case_count: int
    all_passed: bool
    limitations: list[str]


def _package_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _repository_root() -> Path:
    return Path(__file__).resolve().parents[4]


def _phytools_runner_path() -> Path:
    return (
        Path(__file__).resolve().parent
        / "resources"
        / "reference"
        / "phytools_parity_runner.R"
    )


def _failure_root() -> Path:
    return _repository_root() / "artifacts" / "phytools-parity-failures"


def _reference_environment() -> dict[str, str]:
    environment = dict(os.environ)
    r_library = _repository_root() / "artifacts" / "r-lib"
    if "R_LIBS_USER" not in environment and r_library.is_dir():
        environment["R_LIBS_USER"] = str(r_library)
    return environment


def _bijux_version() -> str:
    try:
        return metadata.version("bijux-phylogenetics")
    except metadata.PackageNotFoundError:
        return "0.1.0"


def _bijux_commit() -> str | None:
    result = subprocess.run(
        ["git", "rev-parse", "--short", "HEAD"],
        capture_output=True,
        check=False,
        cwd=_repository_root(),
        text=True,
    )
    if result.returncode != 0:
        return None
    commit = result.stdout.strip()
    return commit or None


def list_phytools_parity_cases() -> list[PhytoolsParityCase]:
    """Return the governed live `phytools` parity cases."""
    strong_signal_fixture = get_shared_phytools_comparative_fixture(
        "phytools_continuous_strong_signal_twenty_four_taxa"
    )
    nonultrametric_signal_fixture = get_shared_phytools_comparative_fixture(
        "phytools_continuous_strong_signal_non_ultrametric_twenty_four_taxa"
    )
    weak_signal_fixture = get_shared_phytools_comparative_fixture(
        "phytools_continuous_weak_signal_twenty_four_taxa"
    )
    missing_signal_fixture = get_shared_phytools_comparative_fixture(
        "phytools_continuous_missing_values_twenty_four_taxa"
    )
    binary_discrete_fixture = get_shared_phytools_comparative_fixture(
        "phytools_discrete_binary_twenty_four_taxa"
    )
    multistate_discrete_fixture = get_shared_phytools_comparative_fixture(
        "phytools_discrete_multistate_twenty_four_taxa"
    )
    binary_discrete_missing_fixture = get_shared_phytools_comparative_fixture(
        "phytools_discrete_binary_missing_twenty_four_taxa"
    )
    multistate_discrete_missing_fixture = get_shared_phytools_comparative_fixture(
        "phytools_discrete_multistate_missing_twenty_four_taxa"
    )
    return [
        PhytoolsParityCase(
            case_id="phylosig-lambda-non-ultrametric-strong-signal-twenty-four-taxa",
            fixture_id=nonultrametric_signal_fixture.fixture_id,
            function_name="phytools::phylosig(method='lambda')",
            python_function_name="estimate_pagels_lambda",
            operation="phylogenetic-signal-lambda",
            input_fixtures=(
                nonultrametric_signal_fixture.tree_path,
                nonultrametric_signal_fixture.traits_path,
            ),
            tolerance=5e-4,
            trait_name=nonultrametric_signal_fixture.trait_name,
            taxon_column=nonultrametric_signal_fixture.taxon_column,
        ),
        PhytoolsParityCase(
            case_id="phylosig-lambda-weak-signal-twenty-four-taxa",
            fixture_id=weak_signal_fixture.fixture_id,
            function_name="phytools::phylosig(method='lambda')",
            python_function_name="estimate_pagels_lambda",
            operation="phylogenetic-signal-lambda",
            input_fixtures=(
                weak_signal_fixture.tree_path,
                weak_signal_fixture.traits_path,
            ),
            tolerance=1e-3,
            trait_name=weak_signal_fixture.trait_name,
            taxon_column=weak_signal_fixture.taxon_column,
        ),
        PhytoolsParityCase(
            case_id="phylosig-k-strong-signal-twenty-four-taxa",
            fixture_id=strong_signal_fixture.fixture_id,
            function_name="phytools::phylosig(method='K')",
            python_function_name="compute_phylogenetic_signal_test",
            operation="phylogenetic-signal-k",
            input_fixtures=(
                strong_signal_fixture.tree_path,
                strong_signal_fixture.traits_path,
            ),
            tolerance=1e-6,
            trait_name=strong_signal_fixture.trait_name,
            taxon_column=strong_signal_fixture.taxon_column,
            permutation_count=199,
            permutation_seed=17,
            field_tolerances={
                "p_value": 0.03,
                "simulated_k_minimum": 0.01,
                "simulated_k_mean": 0.01,
            },
        ),
        PhytoolsParityCase(
            case_id="phylosig-k-weak-signal-twenty-four-taxa",
            fixture_id=weak_signal_fixture.fixture_id,
            function_name="phytools::phylosig(method='K')",
            python_function_name="compute_phylogenetic_signal_test",
            operation="phylogenetic-signal-k",
            input_fixtures=(
                weak_signal_fixture.tree_path,
                weak_signal_fixture.traits_path,
            ),
            tolerance=1e-6,
            trait_name=weak_signal_fixture.trait_name,
            taxon_column=weak_signal_fixture.taxon_column,
            permutation_count=199,
            permutation_seed=17,
            field_tolerances={
                "p_value": 0.03,
                "simulated_k_minimum": 0.01,
                "simulated_k_mean": 0.01,
            },
        ),
        PhytoolsParityCase(
            case_id="fitmk-er-binary-twenty-four-taxa",
            fixture_id=binary_discrete_fixture.fixture_id,
            function_name="phytools::fitMk(model='ER')",
            python_function_name="fit_discrete_mk_model",
            operation="discrete-fit-mk",
            input_fixtures=(
                binary_discrete_fixture.tree_path,
                binary_discrete_fixture.traits_path,
            ),
            tolerance=1e-6,
            trait_name=binary_discrete_fixture.trait_name,
            taxon_column=binary_discrete_fixture.taxon_column,
            discrete_model="equal-rates",
            row_field_tolerances={"rate": 1e-5},
        ),
        PhytoolsParityCase(
            case_id="fitmk-er-multistate-twenty-four-taxa",
            fixture_id=multistate_discrete_fixture.fixture_id,
            function_name="phytools::fitMk(model='ER')",
            python_function_name="fit_discrete_mk_model",
            operation="discrete-fit-mk",
            input_fixtures=(
                multistate_discrete_fixture.tree_path,
                multistate_discrete_fixture.traits_path,
            ),
            tolerance=1e-6,
            trait_name=multistate_discrete_fixture.trait_name,
            taxon_column=multistate_discrete_fixture.taxon_column,
            discrete_model="equal-rates",
            row_field_tolerances={"rate": 1e-5},
        ),
        PhytoolsParityCase(
            case_id="fitmk-er-binary-missing-twenty-four-taxa",
            fixture_id=binary_discrete_missing_fixture.fixture_id,
            function_name="phytools::fitMk(model='ER')",
            python_function_name="fit_discrete_mk_model",
            operation="discrete-fit-mk",
            input_fixtures=(
                binary_discrete_missing_fixture.tree_path,
                binary_discrete_missing_fixture.traits_path,
            ),
            tolerance=1e-6,
            trait_name=binary_discrete_missing_fixture.trait_name,
            taxon_column=binary_discrete_missing_fixture.taxon_column,
            discrete_model="equal-rates",
            row_field_tolerances={"rate": 1e-5},
        ),
        PhytoolsParityCase(
            case_id="fitmk-er-multistate-missing-twenty-four-taxa",
            fixture_id=multistate_discrete_missing_fixture.fixture_id,
            function_name="phytools::fitMk(model='ER')",
            python_function_name="fit_discrete_mk_model",
            operation="discrete-fit-mk",
            input_fixtures=(
                multistate_discrete_missing_fixture.tree_path,
                multistate_discrete_missing_fixture.traits_path,
            ),
            tolerance=1e-6,
            trait_name=multistate_discrete_missing_fixture.trait_name,
            taxon_column=multistate_discrete_missing_fixture.taxon_column,
            discrete_model="equal-rates",
            row_field_tolerances={"rate": 1e-5},
        ),
        PhytoolsParityCase(
            case_id="fitmk-sym-multistate-twenty-four-taxa",
            fixture_id=multistate_discrete_fixture.fixture_id,
            function_name="phytools::fitMk(model='SYM')",
            python_function_name="fit_discrete_mk_model",
            operation="discrete-fit-mk",
            input_fixtures=(
                multistate_discrete_fixture.tree_path,
                multistate_discrete_fixture.traits_path,
            ),
            tolerance=1e-6,
            trait_name=multistate_discrete_fixture.trait_name,
            taxon_column=multistate_discrete_fixture.taxon_column,
            discrete_model="symmetric",
            field_tolerances={
                "log_likelihood": 2e-4,
                "aic": 2e-4,
                "aicc": 2e-4,
            },
            row_field_tolerances={"rate": 1e-4},
        ),
        PhytoolsParityCase(
            case_id="fitmk-sym-multistate-missing-twenty-four-taxa",
            fixture_id=multistate_discrete_missing_fixture.fixture_id,
            function_name="phytools::fitMk(model='SYM')",
            python_function_name="fit_discrete_mk_model",
            operation="discrete-fit-mk",
            input_fixtures=(
                multistate_discrete_missing_fixture.tree_path,
                multistate_discrete_missing_fixture.traits_path,
            ),
            tolerance=1e-6,
            trait_name=multistate_discrete_missing_fixture.trait_name,
            taxon_column=multistate_discrete_missing_fixture.taxon_column,
            discrete_model="symmetric",
            field_tolerances={
                "log_likelihood": 2e-4,
                "aic": 2e-4,
                "aicc": 2e-4,
            },
            row_field_tolerances={"rate": 1e-4},
        ),
        PhytoolsParityCase(
            case_id="fast-anc-strong-signal-twenty-four-taxa",
            fixture_id=strong_signal_fixture.fixture_id,
            function_name="phytools::fastAnc",
            python_function_name="reconstruct_continuous_ancestral_states",
            operation="continuous-ancestral-fast-anc",
            input_fixtures=(
                strong_signal_fixture.tree_path,
                strong_signal_fixture.traits_path,
            ),
            tolerance=1e-8,
            trait_name=strong_signal_fixture.trait_name,
            taxon_column=strong_signal_fixture.taxon_column,
            row_field_tolerances={
                "estimate": 1e-8,
                "standard_error": 1e-8,
            },
        ),
        PhytoolsParityCase(
            case_id="fast-anc-weak-signal-twenty-four-taxa",
            fixture_id=weak_signal_fixture.fixture_id,
            function_name="phytools::fastAnc",
            python_function_name="reconstruct_continuous_ancestral_states",
            operation="continuous-ancestral-fast-anc",
            input_fixtures=(
                weak_signal_fixture.tree_path,
                weak_signal_fixture.traits_path,
            ),
            tolerance=1e-8,
            trait_name=weak_signal_fixture.trait_name,
            taxon_column=weak_signal_fixture.taxon_column,
            row_field_tolerances={
                "estimate": 1e-8,
                "standard_error": 1e-8,
            },
        ),
        PhytoolsParityCase(
            case_id="fast-anc-non-ultrametric-strong-signal-twenty-four-taxa",
            fixture_id=nonultrametric_signal_fixture.fixture_id,
            function_name="phytools::fastAnc",
            python_function_name="reconstruct_continuous_ancestral_states",
            operation="continuous-ancestral-fast-anc",
            input_fixtures=(
                nonultrametric_signal_fixture.tree_path,
                nonultrametric_signal_fixture.traits_path,
            ),
            tolerance=1e-8,
            trait_name=nonultrametric_signal_fixture.trait_name,
            taxon_column=nonultrametric_signal_fixture.taxon_column,
            row_field_tolerances={
                "estimate": 1e-8,
                "standard_error": 1e-8,
            },
        ),
        PhytoolsParityCase(
            case_id="fast-anc-missing-values-twenty-four-taxa",
            fixture_id=missing_signal_fixture.fixture_id,
            function_name="phytools::fastAnc",
            python_function_name="reconstruct_continuous_ancestral_states",
            operation="continuous-ancestral-fast-anc",
            input_fixtures=(
                missing_signal_fixture.tree_path,
                missing_signal_fixture.traits_path,
            ),
            tolerance=1e-8,
            trait_name=missing_signal_fixture.trait_name,
            taxon_column=missing_signal_fixture.taxon_column,
            row_field_tolerances={
                "estimate": 1e-8,
                "standard_error": 1e-8,
            },
        ),
        PhytoolsParityCase(
            case_id="anc-ml-strong-signal-twenty-four-taxa",
            fixture_id=strong_signal_fixture.fixture_id,
            function_name="phytools::anc.ML",
            python_function_name="reconstruct_continuous_ancestral_states",
            operation="continuous-ancestral-anc-ml",
            input_fixtures=(
                strong_signal_fixture.tree_path,
                strong_signal_fixture.traits_path,
            ),
            tolerance=1e-8,
            trait_name=strong_signal_fixture.trait_name,
            taxon_column=strong_signal_fixture.taxon_column,
            field_tolerances={
                "sigma_squared": 5e-8,
            },
            row_field_tolerances={
                "estimate": 1e-8,
                "standard_error": 5e-8,
                "lower_95_interval": 5e-8,
                "upper_95_interval": 5e-8,
            },
        ),
        PhytoolsParityCase(
            case_id="anc-ml-weak-signal-twenty-four-taxa",
            fixture_id=weak_signal_fixture.fixture_id,
            function_name="phytools::anc.ML",
            python_function_name="reconstruct_continuous_ancestral_states",
            operation="continuous-ancestral-anc-ml",
            input_fixtures=(
                weak_signal_fixture.tree_path,
                weak_signal_fixture.traits_path,
            ),
            tolerance=1e-8,
            trait_name=weak_signal_fixture.trait_name,
            taxon_column=weak_signal_fixture.taxon_column,
            field_tolerances={
                "sigma_squared": 5e-8,
            },
            row_field_tolerances={
                "estimate": 1e-8,
                "standard_error": 5e-8,
                "lower_95_interval": 5e-8,
                "upper_95_interval": 5e-8,
            },
        ),
        PhytoolsParityCase(
            case_id="anc-ml-non-ultrametric-strong-signal-twenty-four-taxa",
            fixture_id=nonultrametric_signal_fixture.fixture_id,
            function_name="phytools::anc.ML",
            python_function_name="reconstruct_continuous_ancestral_states",
            operation="continuous-ancestral-anc-ml",
            input_fixtures=(
                nonultrametric_signal_fixture.tree_path,
                nonultrametric_signal_fixture.traits_path,
            ),
            tolerance=1e-8,
            trait_name=nonultrametric_signal_fixture.trait_name,
            taxon_column=nonultrametric_signal_fixture.taxon_column,
            field_tolerances={
                "sigma_squared": 5e-8,
            },
            row_field_tolerances={
                "estimate": 1e-8,
                "standard_error": 5e-8,
                "lower_95_interval": 5e-8,
                "upper_95_interval": 5e-8,
            },
        ),
        PhytoolsParityCase(
            case_id="anc-ml-missing-values-twenty-four-taxa",
            fixture_id=missing_signal_fixture.fixture_id,
            function_name="phytools::anc.ML",
            python_function_name="reconstruct_continuous_ancestral_states",
            operation="continuous-ancestral-anc-ml",
            input_fixtures=(
                missing_signal_fixture.tree_path,
                missing_signal_fixture.traits_path,
            ),
            tolerance=1e-8,
            trait_name=missing_signal_fixture.trait_name,
            taxon_column=missing_signal_fixture.taxon_column,
            field_tolerances={
                "sigma_squared": 5e-8,
            },
            row_field_tolerances={
                "estimate": 1e-8,
                "standard_error": 5e-8,
                "lower_95_interval": 5e-8,
                "upper_95_interval": 5e-8,
            },
        ),
    ]


def _selected_cases(case_ids: list[str] | None) -> list[PhytoolsParityCase]:
    registry = {case.case_id: case for case in list_phytools_parity_cases()}
    if case_ids is None:
        return list(registry.values())
    missing = [case_id for case_id in case_ids if case_id not in registry]
    if missing:
        missing_text = ", ".join(sorted(missing))
        raise ValueError(f"unknown phytools parity case id(s): {missing_text}")
    return [registry[case_id] for case_id in case_ids]


def _write_case_file(path: Path, case: PhytoolsParityCase) -> Path:
    payload = {
        "case_id": case.case_id,
        "fixture_id": case.fixture_id,
        "function_name": case.function_name,
        "operation": case.operation,
        "input_fixtures": [str(path) for path in case.input_fixtures],
        "trait_name": case.trait_name,
        "taxon_column": case.taxon_column,
        "discrete_model": case.discrete_model,
        "tolerance": case.tolerance,
        "permutation_count": case.permutation_count,
        "permutation_seed": case.permutation_seed,
    }
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return path


def _build_bijux_case_payload(
    case: PhytoolsParityCase,
) -> tuple[dict[str, object], list[dict[str, object]] | None]:
    tree_path, traits_path = case.input_fixtures
    if case.operation == "phylogenetic-signal-lambda":
        report = estimate_pagels_lambda(
            tree_path,
            traits_path,
            trait=case.trait_name,
            taxon_column=case.taxon_column,
        )
        return (
            {
                "taxon_count": report.taxon_count,
                "trait_name": report.trait,
                "lambda_value": report.lambda_value,
                "log_likelihood": report.log_likelihood,
                "null_log_likelihood": report.null_log_likelihood,
                "brownian_log_likelihood": report.brownian_log_likelihood,
                "tree_is_ultrametric": report.input_audit.tree_is_ultrametric,
                "pruned_missing_value_taxa": list(
                    report.input_audit.pruned_missing_value_taxa
                ),
                "warning_count": len(report.input_audit.warnings),
            },
            None,
        )
    if case.operation == "phylogenetic-signal-k":
        signal_test = compute_phylogenetic_signal_test(
            tree_path,
            traits_path,
            trait=case.trait_name,
            taxon_column=case.taxon_column,
            permutations=case.permutation_count or 199,
            seed=case.permutation_seed or 1,
        )
        report = compute_blombergs_k(
            tree_path,
            traits_path,
            trait=case.trait_name,
            taxon_column=case.taxon_column,
        )
        return (
            {
                "taxon_count": report.taxon_count,
                "trait_name": report.trait,
                "k": report.k,
                "p_value": signal_test.p_value,
                "permutation_count": signal_test.permutations,
                "permutation_seed": signal_test.seed,
                "null_distribution_count": len(signal_test.permutation_rows),
                "simulated_k_minimum": signal_test.null_distribution_minimum,
                "simulated_k_mean": signal_test.null_distribution_mean,
                "simulated_k_maximum": signal_test.null_distribution_maximum,
                "generalized_mean": report.generalized_mean,
                "observed_mean_square": report.observed_mean_square,
                "phylogenetic_mean_square": report.phylogenetic_mean_square,
                "expected_mean_square_ratio": report.expected_mean_square_ratio,
                "tree_is_ultrametric": report.input_audit.tree_is_ultrametric,
                "pruned_missing_value_taxa": list(
                    report.input_audit.pruned_missing_value_taxa
                ),
                "warning_count": len(report.input_audit.warnings),
            },
            None,
        )
    if case.operation == "discrete-fit-mk":
        report = fit_discrete_mk_model(
            tree_path,
            traits_path,
            trait=case.trait_name,
            taxon_column=case.taxon_column,
            model=case.discrete_model or "equal-rates",
        )
        rows = sorted(
            [
                {
                    "source_state": row.source_state,
                    "target_state": row.target_state,
                    "transition_allowed": row.transition_allowed,
                    "step_distance": row.step_distance,
                    "rate": row.rate,
                }
                for row in report.transition_rate_rows
            ],
            key=lambda row: (str(row["source_state"]), str(row["target_state"])),
        )
        return (
            {
                "taxon_count": report.taxon_count,
                "trait_name": report.trait,
                "excluded_taxon_count": len(report.input_audit.pruned_missing_value_taxa),
                "excluded_taxa": list(report.input_audit.pruned_missing_value_taxa),
                "model": report.model,
                "state_count": len(report.input_audit.observed_states),
                "parameter_count": report.parameter_count,
                "log_likelihood": report.log_likelihood,
                "aic": report.aic,
                "aicc": report.aicc,
                "overparameterized": report.overparameterized,
                "baseline_model": (
                    None
                    if report.baseline_comparison is None
                    else report.baseline_comparison.baseline_model
                ),
                "preferred_model_by_aic": (
                    None
                    if report.baseline_comparison is None
                    else report.baseline_comparison.preferred_model_by_aic
                ),
            },
            rows,
        )
    if case.operation == "continuous-ancestral-fast-anc":
        report = reconstruct_continuous_ancestral_states(
            tree_path,
            traits_path,
            trait=case.trait_name,
            taxon_column=case.taxon_column,
            model="brownian",
            estimator="fast-anc",
        )
        rows = sorted(
            [
                {
                    "node": estimate.node,
                    "estimate": estimate.estimate,
                    "standard_error": estimate.standard_error,
                    "lower_95_interval": estimate.lower_95_interval,
                    "upper_95_interval": estimate.upper_95_interval,
                }
                for estimate in report.estimates
                if not estimate.is_tip
            ],
            key=lambda row: str(row["node"]),
        )
        return (
            {
                "taxon_count": report.taxon_count,
                "trait_name": report.trait,
                "internal_node_count": len(rows),
                "excluded_taxon_count": len(report.dropped_missing_taxa)
                + len(report.dropped_non_numeric_taxa),
                "excluded_taxa": sorted(
                    report.dropped_missing_taxa + report.dropped_non_numeric_taxa
                ),
                "tree_is_ultrametric": (
                    None
                    if report.brownian_fit_diagnostics is None
                    else report.brownian_fit_diagnostics.tree_is_ultrametric
                ),
                "covariance_condition_number": (
                    None
                    if report.brownian_fit_diagnostics is None
                    else report.brownian_fit_diagnostics.covariance_condition_number
                ),
                "log_likelihood": (
                    None
                    if report.brownian_fit_diagnostics is None
                    else report.brownian_fit_diagnostics.log_likelihood
                ),
                "warning_count": len(report.warnings),
            },
            rows,
        )
    if case.operation == "continuous-ancestral-anc-ml":
        report = reconstruct_continuous_ancestral_states(
            tree_path,
            traits_path,
            trait=case.trait_name,
            taxon_column=case.taxon_column,
            model="brownian",
            estimator="anc-ml",
        )
        rows = sorted(
            [
                {
                    "node": estimate.node,
                    "estimate": estimate.estimate,
                    "standard_error": estimate.standard_error,
                    "lower_95_interval": estimate.lower_95_interval,
                    "upper_95_interval": estimate.upper_95_interval,
                }
                for estimate in report.estimates
                if not estimate.is_tip
            ],
            key=lambda row: str(row["node"]),
        )
        return (
            {
                "taxon_count": report.taxon_count,
                "trait_name": report.trait,
                "internal_node_count": len(rows),
                "excluded_taxon_count": len(report.dropped_missing_taxa)
                + len(report.dropped_non_numeric_taxa),
                "excluded_taxa": sorted(
                    report.dropped_missing_taxa + report.dropped_non_numeric_taxa
                ),
                "tree_is_ultrametric": (
                    None
                    if report.brownian_fit_diagnostics is None
                    else report.brownian_fit_diagnostics.tree_is_ultrametric
                ),
                "sigma_squared": (
                    None
                    if report.brownian_fit_diagnostics is None
                    else report.brownian_fit_diagnostics.residual_sigma_squared
                ),
                "log_likelihood": (
                    None
                    if report.brownian_fit_diagnostics is None
                    else report.brownian_fit_diagnostics.log_likelihood
                ),
                "warning_count": len(report.warnings),
            },
            rows,
        )
    raise ValueError(f"unsupported phytools parity operation: {case.operation}")


def _load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _load_rows_table(path: Path) -> list[dict[str, object]]:
    string_identity_fields = {
        "source_state",
        "target_state",
        "node",
    }
    boolean_fields = {
        "transition_allowed",
    }
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows: list[dict[str, object]] = []
        for row in reader:
            parsed: dict[str, object] = {}
            for key, value in row.items():
                if value is None or value == "":
                    parsed[key] = ""
                    continue
                if key in string_identity_fields:
                    parsed[key] = value
                    continue
                if key in boolean_fields:
                    lowered = value.lower()
                    if lowered == "true":
                        parsed[key] = True
                        continue
                    if lowered == "false":
                        parsed[key] = False
                        continue
                try:
                    parsed[key] = int(value)
                    continue
                except ValueError:
                    pass
                try:
                    parsed[key] = float(value)
                    continue
                except ValueError:
                    parsed[key] = value
            rows.append(parsed)
        return rows


def _isclose(left: object, right: object, *, tolerance: float) -> bool:
    if isinstance(left, bool) or isinstance(right, bool):
        return left == right
    if isinstance(left, (int, float)) and isinstance(right, (int, float)):
        return math.isclose(
            float(left),
            float(right),
            rel_tol=tolerance,
            abs_tol=tolerance,
        )
    return left == right


def _field_tolerance(case: PhytoolsParityCase, key: str) -> float:
    if case.field_tolerances and key in case.field_tolerances:
        return case.field_tolerances[key]
    return case.tolerance


def _row_field_tolerance(case: PhytoolsParityCase, key: str) -> float:
    if case.row_field_tolerances and key in case.row_field_tolerances:
        return case.row_field_tolerances[key]
    return case.tolerance


def _mismatch_reason(
    case: PhytoolsParityCase,
    *,
    reference_summary: dict[str, object] | None,
    bijux_summary: dict[str, object] | None,
) -> str | None:
    if reference_summary is None or bijux_summary is None:
        return "summary_missing"
    if case.operation == "phylogenetic-signal-lambda":
        compare_keys = ("taxon_count", "trait_name", "lambda_value", "log_likelihood")
    elif case.operation == "phylogenetic-signal-k":
        compare_keys = (
            "taxon_count",
            "trait_name",
            "k",
            "p_value",
            "permutation_count",
            "permutation_seed",
            "simulated_k_minimum",
            "simulated_k_mean",
        )
    elif case.operation == "discrete-fit-mk":
        compare_keys = (
            "taxon_count",
            "trait_name",
            "excluded_taxon_count",
            "excluded_taxa",
            "model",
            "state_count",
            "parameter_count",
            "log_likelihood",
            "aic",
            "aicc",
            "overparameterized",
            "baseline_model",
            "preferred_model_by_aic",
        )
    elif case.operation == "continuous-ancestral-fast-anc":
        compare_keys = (
            "taxon_count",
            "trait_name",
            "internal_node_count",
            "excluded_taxon_count",
            "excluded_taxa",
            "tree_is_ultrametric",
        )
    elif case.operation == "continuous-ancestral-anc-ml":
        compare_keys = (
            "taxon_count",
            "trait_name",
            "internal_node_count",
            "excluded_taxon_count",
            "excluded_taxa",
            "tree_is_ultrametric",
            "sigma_squared",
            "log_likelihood",
        )
    else:
        return "unsupported_operation"
    for key in compare_keys:
        if key not in reference_summary or key not in bijux_summary:
            return f"summary_field_missing:{key}"
        if not _isclose(
            reference_summary[key],
            bijux_summary[key],
            tolerance=_field_tolerance(case, key),
        ):
            return f"summary_mismatch:{key}"
    return None


def _row_mismatch_reason(
    case: PhytoolsParityCase,
    *,
    reference_rows: list[dict[str, object]] | None,
    bijux_rows: list[dict[str, object]] | None,
) -> str | None:
    if case.operation not in {
        "discrete-fit-mk",
        "continuous-ancestral-fast-anc",
        "continuous-ancestral-anc-ml",
    }:
        return None
    if reference_rows is None or bijux_rows is None:
        return "rows_missing"
    if case.operation == "discrete-fit-mk":
        reference_rows = sorted(
            reference_rows,
            key=lambda row: (
                str(row.get("source_state", "")),
                str(row.get("target_state", "")),
            ),
        )
        bijux_rows = sorted(
            bijux_rows,
            key=lambda row: (
                str(row.get("source_state", "")),
                str(row.get("target_state", "")),
            ),
        )
    else:
        reference_rows = sorted(reference_rows, key=lambda row: str(row.get("node", "")))
        bijux_rows = sorted(bijux_rows, key=lambda row: str(row.get("node", "")))
    if len(reference_rows) != len(bijux_rows):
        return "row_count_mismatch"
    if case.operation == "discrete-fit-mk":
        compare_keys = (
            "source_state",
            "target_state",
            "transition_allowed",
            "step_distance",
            "rate",
        )
    else:
        compare_keys = (
            ("node", "estimate", "standard_error")
            if case.operation == "continuous-ancestral-fast-anc"
            else (
                "node",
                "estimate",
                "standard_error",
                "lower_95_interval",
                "upper_95_interval",
            )
        )
    for reference_row, bijux_row in zip(reference_rows, bijux_rows, strict=True):
        for key in compare_keys:
            if key not in reference_row or key not in bijux_row:
                return f"row_field_missing:{key}"
            if not _isclose(
                reference_row[key],
                bijux_row[key],
                tolerance=_row_field_tolerance(case, key),
            ):
                return f"row_mismatch:{key}"
    return None


def _persist_failure_bundle(
    *,
    failure_root: Path,
    case: PhytoolsParityCase,
    case_file: Path,
    execution_root: Path,
    execution_payload: dict[str, object] | None,
    reference_summary: dict[str, object] | None,
    bijux_summary: dict[str, object] | None,
    reference_rows: list[dict[str, object]] | None,
    bijux_rows: list[dict[str, object]] | None,
    reference_error: dict[str, object] | None,
    bijux_error: dict[str, object] | None,
    mismatch_reason: str,
) -> Path:
    artifact_root = failure_root / case.case_id
    if artifact_root.exists():
        for child in artifact_root.iterdir():
            if child.is_file():
                child.unlink()
    artifact_root.mkdir(parents=True, exist_ok=True)
    (artifact_root / "case.json").write_text(
        case_file.read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    if execution_root.exists():
        for source in execution_root.iterdir():
            if source.is_file():
                (artifact_root / source.name).write_text(
                    source.read_text(encoding="utf-8"),
                    encoding="utf-8",
                )
    if execution_payload is not None:
        (artifact_root / "reference-execution-copy.json").write_text(
            json.dumps(execution_payload, indent=2),
            encoding="utf-8",
        )
    if reference_summary is not None:
        (artifact_root / "reference-summary-copy.json").write_text(
            json.dumps(reference_summary, indent=2),
            encoding="utf-8",
        )
    if bijux_summary is not None:
        (artifact_root / "bijux-summary.json").write_text(
            json.dumps(bijux_summary, indent=2),
            encoding="utf-8",
        )
    if reference_rows is not None:
        (artifact_root / "reference-rows.json").write_text(
            json.dumps(reference_rows, indent=2),
            encoding="utf-8",
        )
    if bijux_rows is not None:
        (artifact_root / "bijux-rows.json").write_text(
            json.dumps(bijux_rows, indent=2),
            encoding="utf-8",
        )
    if reference_error is not None:
        (artifact_root / "reference-error.json").write_text(
            json.dumps(reference_error, indent=2),
            encoding="utf-8",
        )
    if bijux_error is not None:
        (artifact_root / "bijux-error.json").write_text(
            json.dumps(bijux_error, indent=2),
            encoding="utf-8",
        )
    (artifact_root / "mismatch-reason.txt").write_text(
        mismatch_reason,
        encoding="utf-8",
    )
    return artifact_root


def _summary_rows(
    observations: list[PhytoolsParityObservation],
) -> list[PhytoolsParitySummaryRow]:
    function_names = sorted({item.function_name for item in observations})
    rows: list[PhytoolsParitySummaryRow] = []
    for function_name in function_names:
        matching = [
            observation
            for observation in observations
            if observation.function_name == function_name
        ]
        rows.append(
            PhytoolsParitySummaryRow(
                function_name=function_name,
                case_count=len(matching),
                passed_case_count=sum(
                    1 for observation in matching if observation.status == "passed"
                ),
                failed_case_count=sum(
                    1 for observation in matching if observation.status == "failed"
                ),
                skipped_case_count=sum(
                    1 for observation in matching if observation.status == "skipped"
                ),
            )
        )
    return rows


def run_phytools_parity_cases(
    *,
    case_ids: list[str] | None = None,
    rscript_executable: str = "Rscript",
    failure_root: Path | None = None,
) -> PhytoolsParityReport:
    """Run governed live `phytools` parity cases through the checked-in R runner."""
    selected = _selected_cases(case_ids)
    observations: list[PhytoolsParityObservation] = []
    active_failure_root = _failure_root() if failure_root is None else failure_root
    bijux_version = _bijux_version()
    bijux_commit = _bijux_commit()
    for case in selected:
        with tempfile.TemporaryDirectory(
            prefix=f"bijux-phytools-parity-{case.case_id}-"
        ) as tmpdir:
            working_root = Path(tmpdir)
            case_file = _write_case_file(working_root / "case.json", case)
            execution_root = working_root / "reference"
            execution_root.mkdir(parents=True, exist_ok=True)
            bijux_summary: dict[str, object] | None = None
            bijux_rows: list[dict[str, object]] | None = None
            bijux_error: dict[str, object] | None = None
            try:
                bijux_summary, bijux_rows = _build_bijux_case_payload(case)
            except Exception as error:
                bijux_error = {
                    "error_type": type(error).__name__,
                    "message": str(error),
                }
            execution_payload: dict[str, object] | None = None
            reference_summary: dict[str, object] | None = None
            reference_rows: list[dict[str, object]] | None = None
            reference_error: dict[str, object] | None = None
            status = "failed"
            mismatch_reason: str | None = None
            artifact_root: Path | None = None
            r_version: str | None = None
            phytools_version: str | None = None
            process_stdout = ""
            process_stderr = ""
            try:
                process = subprocess.run(
                    [
                        rscript_executable,
                        str(_phytools_runner_path()),
                        str(case_file),
                        str(execution_root),
                    ],
                    capture_output=True,
                    check=False,
                    cwd=_repository_root(),
                    env=_reference_environment(),
                    text=True,
                )
                process_stdout = process.stdout
                process_stderr = process.stderr
            except FileNotFoundError:
                process = None
                status = "skipped"
                mismatch_reason = "rscript_unavailable"
            if process is not None and process.returncode == 0:
                execution_path = execution_root / "reference-execution.json"
                summary_path = execution_root / "reference-summary.json"
                if not execution_path.exists():
                    mismatch_reason = "reference_execution_failed"
                else:
                    execution_payload = _load_json(execution_path)
                    r_version = execution_payload.get("r_version")  # type: ignore[assignment]
                    phytools_version = execution_payload.get("phytools_version")  # type: ignore[assignment]
                    execution_status = execution_payload.get("status")
                    if execution_status == "unavailable":
                        status = "skipped"
                        mismatch_reason = str(
                            execution_payload.get(
                                "mismatch_reason",
                                "phytools_package_unavailable",
                            )
                        )
                    elif execution_status != "ok":
                        reference_error = {
                            "error_type": str(
                                execution_payload.get(
                                    "error_type",
                                    execution_payload.get(
                                        "mismatch_reason",
                                        "reference_execution_failed",
                                    ),
                                )
                            ),
                            "message": str(execution_payload.get("message", "")),
                        }
                        mismatch_reason = str(
                            execution_payload.get(
                                "mismatch_reason",
                                "reference_execution_failed",
                            )
                        )
                    elif not summary_path.exists():
                        mismatch_reason = "reference_summary_missing"
                    else:
                        reference_summary = _load_json(summary_path)
                        mismatch_reason = _mismatch_reason(
                            case,
                            reference_summary=reference_summary,
                            bijux_summary=bijux_summary,
                        )
                        if mismatch_reason is None and case.operation in {
                            "discrete-fit-mk",
                            "continuous-ancestral-fast-anc",
                            "continuous-ancestral-anc-ml",
                        }:
                            rows_path = execution_root / (
                                "fitmk-rate-matrix.tsv"
                                if case.operation == "discrete-fit-mk"
                                else (
                                    "fast-anc-node-estimates.tsv"
                                    if case.operation == "continuous-ancestral-fast-anc"
                                    else "anc-ml-node-estimates.tsv"
                                )
                            )
                            if not rows_path.exists():
                                mismatch_reason = "reference_rows_missing"
                            else:
                                reference_rows = _load_rows_table(rows_path)
                                mismatch_reason = _row_mismatch_reason(
                                    case,
                                    reference_rows=reference_rows,
                                    bijux_rows=bijux_rows,
                                )
                        if mismatch_reason is None:
                            status = "passed"
            elif process is not None and process.returncode != 0:
                mismatch_reason = "reference_execution_failed"
            if status != "passed":
                artifact_root = _persist_failure_bundle(
                    failure_root=active_failure_root,
                    case=case,
                    case_file=case_file,
                    execution_root=execution_root,
                    execution_payload=execution_payload,
                    reference_summary=reference_summary,
                    bijux_summary=bijux_summary,
                    reference_rows=reference_rows,
                    bijux_rows=bijux_rows,
                    reference_error=reference_error,
                    bijux_error=bijux_error,
                    mismatch_reason=mismatch_reason or "reference_execution_failed",
                )
                if process_stdout:
                    (artifact_root / "reference-stdout.txt").write_text(
                        process_stdout,
                        encoding="utf-8",
                    )
                if process_stderr:
                    (artifact_root / "reference-stderr.txt").write_text(
                        process_stderr,
                        encoding="utf-8",
                    )
            observations.append(
                PhytoolsParityObservation(
                    case_id=case.case_id,
                    fixture_id=case.fixture_id,
                    function_name=case.function_name,
                    python_function_name=case.python_function_name,
                    input_fixtures=case.input_fixtures,
                    tolerance=case.tolerance,
                    r_version=r_version,
                    phytools_version=phytools_version,
                    bijux_version=bijux_version,
                    bijux_commit=bijux_commit,
                    status=status,
                    passed=status == "passed",
                    mismatch_reason=mismatch_reason,
                    reproducible_artifact_root=artifact_root,
                    reference_summary=reference_summary,
                    bijux_summary=bijux_summary,
                    reference_rows=reference_rows,
                    bijux_rows=bijux_rows,
                    reference_error=reference_error,
                    bijux_error=bijux_error,
                )
            )
    case_count = len(observations)
    passed_case_count = sum(1 for item in observations if item.status == "passed")
    failed_case_count = sum(1 for item in observations if item.status == "failed")
    skipped_case_count = sum(1 for item in observations if item.status == "skipped")
    return PhytoolsParityReport(
        observations=observations,
        summary_rows=_summary_rows(observations),
        case_count=case_count,
        passed_case_count=passed_case_count,
        failed_case_count=failed_case_count,
        skipped_case_count=skipped_case_count,
        all_passed=case_count > 0
        and passed_case_count == case_count
        and failed_case_count == 0
        and skipped_case_count == 0,
        limitations=[
            "The governed live `phytools` parity registry is intentionally narrow until later rounds expand the comparative fixture surface.",
            "This harness requires Rscript plus the `phytools` and `jsonlite` R packages for live reference execution.",
        ],
    )


def write_phytools_parity_summary_table(
    path: Path,
    report: PhytoolsParityReport,
) -> Path:
    """Write one row per governed `phytools` function summary."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "function_name",
                "case_count",
                "passed_case_count",
                "failed_case_count",
                "skipped_case_count",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.summary_rows:
            writer.writerow(asdict(row))
    return path


def write_phytools_parity_observation_table(
    path: Path,
    report: PhytoolsParityReport,
) -> Path:
    """Write one row per governed `phytools` parity observation."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "case_id",
                "fixture_id",
                "function_name",
                "python_function_name",
                "input_fixtures",
                "tolerance",
                "r_version",
                "phytools_version",
                "bijux_version",
                "bijux_commit",
                "status",
                "passed",
                "mismatch_reason",
                "reproducible_artifact_root",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for observation in report.observations:
            writer.writerow(
                {
                    "case_id": observation.case_id,
                    "fixture_id": observation.fixture_id,
                    "function_name": observation.function_name,
                    "python_function_name": observation.python_function_name,
                    "input_fixtures": json.dumps(
                        [str(path) for path in observation.input_fixtures]
                    ),
                    "tolerance": format(observation.tolerance, ".12g"),
                    "r_version": observation.r_version or "",
                    "phytools_version": observation.phytools_version or "",
                    "bijux_version": observation.bijux_version,
                    "bijux_commit": observation.bijux_commit or "",
                    "status": observation.status,
                    "passed": str(observation.passed).lower(),
                    "mismatch_reason": observation.mismatch_reason or "",
                    "reproducible_artifact_root": ""
                    if observation.reproducible_artifact_root is None
                    else str(observation.reproducible_artifact_root),
                }
            )
    return path
