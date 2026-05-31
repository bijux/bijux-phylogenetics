from __future__ import annotations

import csv
from dataclasses import dataclass
import json
from pathlib import Path
import re
import tempfile

from Bio import Phylo

from bijux_phylogenetics.comparative import validate_comparative_reference_examples
from bijux_phylogenetics.comparative.evolutionary_modes import (
    fit_continuous_evolutionary_mode,
)
from bijux_phylogenetics.comparative.pgls import build_pgls_model_matrix, run_pgls
from bijux_phylogenetics.comparative.signal import (
    compute_blombergs_k,
    estimate_pagels_lambda,
)
from bijux_phylogenetics.compare.topology import (
    compare_branch_score_distance,
    compare_robinson_foulds,
)
from bijux_phylogenetics.diagnostics.validation.structure import _load_tree
from bijux_phylogenetics.io.biopython import tree_from_biophylo
from bijux_phylogenetics.io.newick import dumps_newick, write_newick
from bijux_phylogenetics.io.trees import detect_tree_format
from bijux_phylogenetics.phylo.topology.clades import informative_unrooted_splits
from bijux_phylogenetics.trees import (
    compute_clade_frequency_table,
    compute_consensus_tree,
)


@dataclass(frozen=True, slots=True)
class ReferenceParityObservation:
    """One parity comparison between Bijux and a checked-in external reference."""

    suite: str
    method: str
    case: str
    input_fixtures: list[Path]
    reference_tool: str
    reference_version: str
    reference_source: str
    tolerance: float
    tolerance_reason: str
    expected_failure_mode: str
    taxon_overlap_policy: str | None
    shared_taxa: list[str]
    left_only_taxa: list[str]
    right_only_taxa: list[str]
    expected_output: dict[str, object]
    observed_output: dict[str, object]
    mismatch_kind: str | None
    passed: bool


@dataclass(frozen=True, slots=True)
class ReferenceParitySummaryRow:
    """One method-level summary row across multiple parity observations."""

    suite: str
    method: str
    case_count: int
    passed_case_count: int
    failed_case_count: int
    reference_tools: list[str]


@dataclass(slots=True)
class ReferenceParityReport:
    """Integrated reference-parity report across core numerical methods."""

    suite: str
    observations: list[ReferenceParityObservation]
    summary_rows: list[ReferenceParitySummaryRow]
    reference_tools: dict[str, str]
    covered_methods: list[str]
    case_count: int
    passed_case_count: int
    failed_case_count: int
    all_passed: bool
    limitations: list[str]


def _fixtures_root() -> Path:
    return _package_root() / "tests" / "fixtures"


def _package_root() -> Path:
    return Path(__file__).resolve().parents[3]


def _repository_root() -> Path:
    return Path(__file__).resolve().parents[5]


def _expected_root() -> Path:
    return _fixtures_root() / "expected"


def _load_fixture_document(name: str) -> dict[str, object]:
    return json.loads((_expected_root() / name).read_text(encoding="utf-8"))


def _resolve_fixture_path(relative_path: str) -> Path:
    return _package_root() / relative_path


def _resolve_repository_path(relative_path: str) -> Path:
    return _repository_root() / relative_path


def _resolve_input_path(relative_path: str) -> Path:
    package_candidate = _resolve_fixture_path(relative_path)
    if package_candidate.exists():
        return package_candidate
    repository_candidate = _resolve_repository_path(relative_path)
    if repository_candidate.exists():
        return repository_candidate
    return package_candidate


def _strip_nexus_rooted_flag(tree_text: str) -> str:
    """Remove a leading Nexus rootedness marker before plain Newick comparison."""
    return re.sub(r"^\s*\[\&[RU]\]\s*", "", tree_text, count=1)


def _numeric_outputs_match(
    expected: dict[str, object],
    observed: dict[str, object],
    *,
    tolerance: float,
) -> bool:
    if set(expected) != set(observed):
        return False
    for key, expected_value in expected.items():
        observed_value = observed[key]
        if isinstance(expected_value, (int, float)):
            if not isinstance(observed_value, (int, float)):
                return False
            if abs(float(observed_value) - float(expected_value)) > tolerance:
                return False
            continue
        if observed_value != expected_value:
            return False
    return True


def _classify_comparative_mismatch(method: str) -> str:
    if method in {
        "pgls",
        "pagels-lambda",
        "ornstein-uhlenbeck-trait-model",
    }:
        return "model_assumption"
    return "numerical_tolerance"


def _build_comparative_observations() -> tuple[
    list[ReferenceParityObservation], dict[str, str]
]:
    fixture = _load_fixture_document("reference_parity_core.json")
    comparative_cases = [
        entry
        for entry in fixture["observations"]
        if entry["method"]
        in {
            "pgls",
            "pagels-lambda",
            "brownian-trait-model",
            "ornstein-uhlenbeck-trait-model",
            "phylogenetic-independent-contrasts",
            "blombergs-k",
        }
    ]
    reference_report = validate_comparative_reference_examples()
    observed_by_case = {
        observation.case: observation for observation in reference_report.observations
    }
    observations: list[ReferenceParityObservation] = []
    reference_tools: dict[str, str] = {}
    for entry in comparative_cases:
        reference_tools[str(entry["reference_tool"])] = str(entry["reference_version"])
        if str(entry["method"]) == "pgls":
            observations.append(_build_pgls_observation(entry, suite="core"))
            continue
        observed = observed_by_case[str(entry["case"])]
        expected_output = {
            key: float(value) for key, value in dict(entry["expected_output"]).items()
        }
        observed_output = {
            key: float(value) for key, value in observed.observed_parameters.items()
        }
        passed = _numeric_outputs_match(
            expected_output,
            observed_output,
            tolerance=float(entry["tolerance"]),
        )
        observations.append(
            ReferenceParityObservation(
                suite="core",
                method=str(entry["method"]),
                case=str(entry["case"]),
                input_fixtures=[
                    _resolve_input_path(path) for path in entry["input_fixtures"]
                ],
                reference_tool=str(entry["reference_tool"]),
                reference_version=str(entry["reference_version"]),
                reference_source=str(entry["reference_source"]),
                tolerance=float(entry["tolerance"]),
                tolerance_reason=str(entry["tolerance_reason"]),
                expected_failure_mode=str(entry["failure_mode"]),
                taxon_overlap_policy=None,
                shared_taxa=[],
                left_only_taxa=[],
                right_only_taxa=[],
                expected_output=expected_output,
                observed_output=observed_output,
                mismatch_kind=None
                if passed
                else _classify_comparative_mismatch(str(entry["method"])),
                passed=passed,
            )
        )
    return observations, reference_tools


def _format_unrooted_split(split: frozenset[str], universe: set[str]) -> str:
    left = sorted(split)
    right = sorted(universe - set(split))
    return f"{'|'.join(left)}||{'|'.join(right)}"


def _consensus_splits(path: Path) -> list[str]:
    tree = _load_tree(path)
    shared_taxa = set(tree.tip_names)
    return sorted(
        _format_unrooted_split(split, shared_taxa)
        for split in informative_unrooted_splits(tree, shared_taxa)
    )


def _clade_frequency_observed_output(path: Path) -> dict[str, float]:
    with tempfile.TemporaryDirectory(prefix="bijux-reference-parity-") as tmp_dir:
        normalized_path = _normalize_tree_set_path(path, Path(tmp_dir))
        report = compute_clade_frequency_table(normalized_path)
    return {row.clade: row.frequency for row in report.clade_frequencies}


def _normalize_tree_set_path(path: Path, work_dir: Path) -> Path:
    source_format = detect_tree_format(path)
    if source_format == "newick":
        return path

    trees = [
        tree_from_biophylo(tree, source_format=source_format)
        for tree in Phylo.parse(path, source_format)
    ]
    normalized_path = work_dir / f"{path.stem}.normalized.nwk"
    normalized_path.write_text(
        "".join(f"{dumps_newick(tree)}\n" for tree in trees),
        encoding="utf-8",
    )
    return normalized_path


def _build_tree_observation(
    entry: dict[str, object],
    *,
    suite: str,
) -> ReferenceParityObservation:
    input_paths = [_resolve_input_path(path) for path in entry["input_fixtures"]]
    method = str(entry["method"])
    expected_output: dict[str, object] = dict(entry["expected_output"])
    tolerance = float(entry["tolerance"])
    expected_failure_mode = str(entry["failure_mode"])
    taxon_overlap_policy = (
        str(entry["taxon_overlap_policy"])
        if entry.get("taxon_overlap_policy") is not None
        else None
    )
    mismatch_kind: str | None = None
    shared_taxa: list[str] = []
    left_only_taxa: list[str] = []
    right_only_taxa: list[str] = []

    if method == "robinson-foulds-distance":
        report = compare_robinson_foulds(
            input_paths[0],
            input_paths[1],
            rf_mode=str(entry.get("rf_mode", "unrooted")),
            taxon_overlap_policy=taxon_overlap_policy or "prune-to-shared",
        )
        shared_taxa = report.shared_taxa
        left_only_taxa = report.left_only_taxa
        right_only_taxa = report.right_only_taxa
        taxon_overlap_policy = report.taxon_overlap_policy
        observed_output: dict[str, object] = {
            "robinson_foulds_distance": report.robinson_foulds_distance,
            "normalized_robinson_foulds": report.normalized_robinson_foulds,
        }
        passed = _numeric_outputs_match(
            expected_output,
            observed_output,
            tolerance=tolerance,
        )
        mismatch_kind = None if passed else "topology"
    elif method == "branch-score-distance":
        report = compare_branch_score_distance(
            input_paths[0],
            input_paths[1],
            taxon_overlap_policy=taxon_overlap_policy or "prune-to-shared",
        )
        shared_taxa = report.shared_taxa
        left_only_taxa = report.left_only_taxa
        right_only_taxa = report.right_only_taxa
        taxon_overlap_policy = report.taxon_overlap_policy
        observed_output = {
            "branch_score_distance": report.branch_score_distance,
        }
        passed = _numeric_outputs_match(
            expected_output,
            observed_output,
            tolerance=tolerance,
        )
        mismatch_kind = None if passed else "branch_length"
    elif method == "posterior-clade-frequencies":
        observed_output = _clade_frequency_observed_output(input_paths[0])
        passed = _numeric_outputs_match(
            expected_output,
            observed_output,
            tolerance=tolerance,
        )
        mismatch_kind = None if passed else "numerical_tolerance"
    elif method == "consensus-tree-generation":
        with tempfile.TemporaryDirectory(prefix="bijux-reference-parity-") as tmp_dir:
            normalized_input_path = _normalize_tree_set_path(
                input_paths[0],
                Path(tmp_dir),
            )
            tree, _consensus = compute_consensus_tree(normalized_input_path)
            observed_path = Path(tmp_dir) / "observed-consensus.nwk"
            reference_path = Path(tmp_dir) / "reference-consensus.nwk"
            write_newick(observed_path, tree)
            reference_path.write_text(
                f"{_strip_nexus_rooted_flag(str(expected_output['reference_consensus_newick']))}\n",
                encoding="utf-8",
            )
            rf_report = compare_robinson_foulds(
                observed_path,
                reference_path,
                rf_mode="unrooted",
            )
            observed_output = {
                "consensus_splits": _consensus_splits(observed_path),
                "observed_consensus_newick": observed_path.read_text(
                    encoding="utf-8"
                ).strip(),
                "unrooted_robinson_foulds": rf_report.robinson_foulds_distance,
            }
        passed = (
            observed_output["consensus_splits"] == expected_output["consensus_splits"]
            and observed_output["unrooted_robinson_foulds"]
            == expected_output["unrooted_robinson_foulds"]
        )
        mismatch_kind = None if passed else "topology"
    else:  # pragma: no cover - guarded by checked-in fixtures
        raise ValueError(f"unsupported parity method '{method}'")

    return ReferenceParityObservation(
        suite=suite,
        method=method,
        case=str(entry["case"]),
        input_fixtures=input_paths,
        reference_tool=str(entry["reference_tool"]),
        reference_version=str(entry["reference_version"]),
        reference_source=str(entry["reference_source"]),
        tolerance=tolerance,
        tolerance_reason=str(entry["tolerance_reason"]),
        expected_failure_mode=expected_failure_mode,
        taxon_overlap_policy=taxon_overlap_policy,
        shared_taxa=shared_taxa,
        left_only_taxa=left_only_taxa,
        right_only_taxa=right_only_taxa,
        expected_output=expected_output,
        observed_output=observed_output,
        mismatch_kind=mismatch_kind,
        passed=passed,
    )


def _extended_primate_pgls_estimated_lambda_observation(
    input_paths: list[Path],
) -> dict[str, float]:
    report = run_pgls(
        input_paths[0],
        input_paths[1],
        response="longevity",
        predictors=["social_group_size"],
        taxon_column="species",
        lambda_value="estimate",
    )
    coefficients = {
        coefficient.name: coefficient.estimate for coefficient in report.coefficients
    }
    return {
        "intercept": coefficients["intercept"],
        "social_group_size": coefficients["social_group_size"],
        "log_likelihood": report.log_likelihood,
        "lambda_value": report.lambda_value,
    }


def _flatten_pgls_observed_output(
    *,
    report,
    model_matrix=None,
) -> dict[str, float]:
    observed_output: dict[str, float] = {
        "log_likelihood": report.log_likelihood,
        "aic": report.aic,
        "lambda_value": report.lambda_value,
    }
    for coefficient in report.coefficients:
        observed_output[f"coefficient.{coefficient.name}.estimate"] = (
            coefficient.estimate
        )
        observed_output[f"coefficient.{coefficient.name}.standard_error"] = (
            coefficient.standard_error
        )
        observed_output[f"coefficient.{coefficient.name}.p_value"] = coefficient.p_value
    if model_matrix is not None:
        for row in model_matrix.rows:
            for column, value in row.encoded_values.items():
                observed_output[f"model_matrix.{row.taxon}.{column}"] = value
    return observed_output


def _build_pgls_observation(
    entry: dict[str, object],
    *,
    suite: str,
) -> ReferenceParityObservation:
    input_paths = [_resolve_input_path(path) for path in entry["input_fixtures"]]
    formula = None if entry.get("formula") is None else str(entry.get("formula"))
    predictors = (
        None
        if entry.get("predictors") is None
        else [str(value) for value in list(entry["predictors"])]
    )
    lambda_value_entry = entry.get("lambda_value", "estimate")
    lambda_value: float | str
    if lambda_value_entry == "estimate":
        lambda_value = "estimate"
    else:
        lambda_value = float(lambda_value_entry)
    report = run_pgls(
        input_paths[0],
        input_paths[1],
        response=None if entry.get("response") is None else str(entry["response"]),
        predictors=predictors,
        formula=formula,
        taxon_column=(
            None if entry.get("taxon_column") is None else str(entry["taxon_column"])
        ),
        lambda_value=lambda_value,
    )
    model_matrix = None
    if bool(entry.get("include_model_matrix", False)):
        model_matrix = build_pgls_model_matrix(
            input_paths[0],
            input_paths[1],
            response=None if entry.get("response") is None else str(entry["response"]),
            predictors=predictors,
            formula=formula,
            taxon_column=(
                None
                if entry.get("taxon_column") is None
                else str(entry["taxon_column"])
            ),
        )
    expected_output = {
        key: float(value) for key, value in dict(entry["expected_output"]).items()
    }
    observed_output = _flatten_pgls_observed_output(
        report=report,
        model_matrix=model_matrix,
    )
    passed = _numeric_outputs_match(
        expected_output,
        observed_output,
        tolerance=float(entry["tolerance"]),
    )
    return ReferenceParityObservation(
        suite=suite,
        method="pgls",
        case=str(entry["case"]),
        input_fixtures=input_paths,
        reference_tool=str(entry["reference_tool"]),
        reference_version=str(entry["reference_version"]),
        reference_source=str(entry["reference_source"]),
        tolerance=float(entry["tolerance"]),
        tolerance_reason=str(entry["tolerance_reason"]),
        expected_failure_mode=str(entry["failure_mode"]),
        taxon_overlap_policy=None,
        shared_taxa=[],
        left_only_taxa=[],
        right_only_taxa=[],
        expected_output=expected_output,
        observed_output=observed_output,
        mismatch_kind=None if passed else _classify_comparative_mismatch("pgls"),
        passed=passed,
    )


def _extended_primate_pagels_lambda_observation(
    input_paths: list[Path],
) -> dict[str, float]:
    report = estimate_pagels_lambda(
        input_paths[0],
        input_paths[1],
        trait="longevity",
        taxon_column="species",
    )
    return {
        "lambda_value": report.lambda_value,
        "log_likelihood": report.log_likelihood,
    }


def _extended_signal_pagels_lambda_observation(
    input_paths: list[Path],
    *,
    trait: str,
    taxon_column: str,
) -> dict[str, float]:
    report = estimate_pagels_lambda(
        input_paths[0],
        input_paths[1],
        trait=trait,
        taxon_column=taxon_column,
    )
    return {
        "lambda_value": report.lambda_value,
        "log_likelihood": report.log_likelihood,
    }


def _extended_signal_blombergs_k_observation(
    input_paths: list[Path],
    *,
    trait: str,
    taxon_column: str,
) -> dict[str, float]:
    report = compute_blombergs_k(
        input_paths[0],
        input_paths[1],
        trait=trait,
        taxon_column=taxon_column,
    )
    return {"k": report.k}


def _extended_primate_brownian_mode_observation(
    input_paths: list[Path],
) -> dict[str, float]:
    report = fit_continuous_evolutionary_mode(
        input_paths[0],
        input_paths[1],
        trait="longevity",
        taxon_column="species",
        mode="brownian",
    )
    return {
        "aic": report.aic,
        "log_likelihood": report.log_likelihood,
        "rate": report.rate,
        "root_state": report.root_state,
    }


def _extended_primate_ou_mode_observation(
    input_paths: list[Path],
) -> dict[str, float]:
    report = fit_continuous_evolutionary_mode(
        input_paths[0],
        input_paths[1],
        trait="longevity",
        taxon_column="species",
        mode="ornstein-uhlenbeck",
    )
    return {
        "aic": report.aic,
        "alpha": report.parameter_value or 0.0,
        "log_likelihood": report.log_likelihood,
        "rate": report.rate,
        "root_state": report.root_state,
    }


def _build_extended_comparative_observations() -> tuple[
    list[ReferenceParityObservation], dict[str, str]
]:
    fixture = _load_fixture_document("reference_parity_extended_comparative.json")
    observed_builders = {
        "pagel-lambda-non-ultrametric-strong-signal-twenty-four-taxa": (
            lambda input_paths: _extended_signal_pagels_lambda_observation(
                input_paths,
                trait="signal_strong",
                taxon_column="taxon",
            )
        ),
        "pagel-lambda-weak-signal-twenty-four-taxa": (
            lambda input_paths: _extended_signal_pagels_lambda_observation(
                input_paths,
                trait="signal_weak",
                taxon_column="taxon",
            )
        ),
        "blomberg-k-strong-signal-twenty-four-taxa": (
            lambda input_paths: _extended_signal_blombergs_k_observation(
                input_paths,
                trait="signal_strong",
                taxon_column="taxon",
            )
        ),
        "blomberg-k-weak-signal-twenty-four-taxa": (
            lambda input_paths: _extended_signal_blombergs_k_observation(
                input_paths,
                trait="signal_weak",
                taxon_column="taxon",
            )
        ),
        "pgls-primate-longevity-estimated-lambda": (
            _extended_primate_pgls_estimated_lambda_observation
        ),
        "pagel-lambda-primate-longevity": (_extended_primate_pagels_lambda_observation),
        "brownian-primate-longevity-intercept-only": (
            _extended_primate_brownian_mode_observation
        ),
        "ornstein-uhlenbeck-primate-longevity-intercept-only": (
            _extended_primate_ou_mode_observation
        ),
    }
    observations: list[ReferenceParityObservation] = []
    for entry in fixture["observations"]:
        case = str(entry["case"])
        if str(entry["method"]) == "pgls":
            observations.append(_build_pgls_observation(entry, suite="extended"))
            continue
        input_paths = [_resolve_input_path(path) for path in entry["input_fixtures"]]
        expected_output = {
            key: float(value) for key, value in dict(entry["expected_output"]).items()
        }
        observed_output = observed_builders[case](input_paths)
        passed = _numeric_outputs_match(
            expected_output,
            observed_output,
            tolerance=float(entry["tolerance"]),
        )
        observations.append(
            ReferenceParityObservation(
                suite="extended",
                method=str(entry["method"]),
                case=case,
                input_fixtures=input_paths,
                reference_tool=str(entry["reference_tool"]),
                reference_version=str(entry["reference_version"]),
                reference_source=str(entry["reference_source"]),
                tolerance=float(entry["tolerance"]),
                tolerance_reason=str(entry["tolerance_reason"]),
                expected_failure_mode=str(entry["failure_mode"]),
                taxon_overlap_policy=None,
                shared_taxa=[],
                left_only_taxa=[],
                right_only_taxa=[],
                expected_output=expected_output,
                observed_output=observed_output,
                mismatch_kind=None if passed else str(entry["failure_mode"]),
                passed=passed,
            )
        )
    return observations, {
        str(tool): str(version)
        for tool, version in dict(fixture["source_packages"]).items()
    }


def _build_tree_observations(
    fixture_name: str,
    *,
    suite: str,
    methods: set[str],
) -> tuple[list[ReferenceParityObservation], dict[str, str]]:
    fixture = _load_fixture_document(fixture_name)
    observations = [
        _build_tree_observation(entry, suite=suite)
        for entry in fixture["observations"]
        if str(entry["method"]) in methods
    ]
    return observations, {
        str(tool): str(version)
        for tool, version in dict(fixture["source_packages"]).items()
    }


def _build_summary_rows(
    observations: list[ReferenceParityObservation],
) -> list[ReferenceParitySummaryRow]:
    methods = sorted({observation.method for observation in observations})
    rows: list[ReferenceParitySummaryRow] = []
    for method in methods:
        selected = [item for item in observations if item.method == method]
        suites = sorted({item.suite for item in selected})
        rows.append(
            ReferenceParitySummaryRow(
                suite=suites[0] if len(suites) == 1 else "mixed",
                method=method,
                case_count=len(selected),
                passed_case_count=sum(1 for item in selected if item.passed),
                failed_case_count=sum(1 for item in selected if not item.passed),
                reference_tools=sorted(
                    {
                        f"{item.reference_tool} {item.reference_version}"
                        for item in selected
                    }
                ),
            )
        )
    return rows


def validate_reference_parity_examples(
    *,
    include_extended: bool = False,
) -> ReferenceParityReport:
    """Validate core numerical methods against checked-in external reference outputs."""
    observations, comparative_tools = _build_comparative_observations()
    tree_observations, tree_tools = _build_tree_observations(
        "reference_parity_core.json",
        suite="core",
        methods={
            "robinson-foulds-distance",
            "branch-score-distance",
            "posterior-clade-frequencies",
            "consensus-tree-generation",
        },
    )
    observations.extend(tree_observations)
    reference_tools = {**comparative_tools, **tree_tools}
    limitations = [
        "The default parity suite uses small checked-in fixtures so it can run routinely in CI.",
        "Consensus-tree parity is evaluated on unrooted majority-rule topology against DendroPy rather than on branch-length aggregation semantics.",
    ]
    suite = "core"
    if include_extended:
        extended_comparative_observations, extended_comparative_tools = (
            _build_extended_comparative_observations()
        )
        observations.extend(extended_comparative_observations)
        reference_tools.update(extended_comparative_tools)
        extended_observations, extended_tools = _build_tree_observations(
            "reference_parity_extended.json",
            suite="extended",
            methods={
                "posterior-clade-frequencies",
                "consensus-tree-generation",
            },
        )
        observations.extend(extended_observations)
        reference_tools.update(extended_tools)
        suite = "core+extended"
        limitations.append(
            "The extended suite adds governed primate comparative parity checks and a larger posterior tree-set parity check for optional validation lanes."
        )
    summary_rows = _build_summary_rows(observations)
    case_count = len(observations)
    passed_case_count = sum(1 for observation in observations if observation.passed)
    failed_case_count = case_count - passed_case_count
    return ReferenceParityReport(
        suite=suite,
        observations=observations,
        summary_rows=summary_rows,
        reference_tools=reference_tools,
        covered_methods=sorted({observation.method for observation in observations}),
        case_count=case_count,
        passed_case_count=passed_case_count,
        failed_case_count=failed_case_count,
        all_passed=failed_case_count == 0,
        limitations=limitations,
    )


def write_reference_parity_summary_table(
    path: Path,
    report: ReferenceParityReport,
) -> Path:
    """Write one row per method in the reference parity report."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "suite",
                "method",
                "case_count",
                "passed_case_count",
                "failed_case_count",
                "reference_tools",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.summary_rows:
            writer.writerow(
                {
                    "suite": row.suite,
                    "method": row.method,
                    "case_count": row.case_count,
                    "passed_case_count": row.passed_case_count,
                    "failed_case_count": row.failed_case_count,
                    "reference_tools": ",".join(row.reference_tools),
                }
            )
    return path


def write_reference_parity_observation_table(
    path: Path,
    report: ReferenceParityReport,
) -> Path:
    """Write one row per checked reference parity observation."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "suite",
                "method",
                "case",
                "input_fixtures",
                "reference_tool",
                "reference_version",
                "reference_source",
                "tolerance",
                "tolerance_reason",
                "expected_failure_mode",
                "taxon_overlap_policy",
                "shared_taxa",
                "left_only_taxa",
                "right_only_taxa",
                "passed",
                "mismatch_kind",
                "expected_output",
                "observed_output",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for observation in report.observations:
            writer.writerow(
                {
                    "suite": observation.suite,
                    "method": observation.method,
                    "case": observation.case,
                    "input_fixtures": ",".join(
                        str(path) for path in observation.input_fixtures
                    ),
                    "reference_tool": observation.reference_tool,
                    "reference_version": observation.reference_version,
                    "reference_source": observation.reference_source,
                    "tolerance": format(observation.tolerance, ".12g"),
                    "tolerance_reason": observation.tolerance_reason,
                    "expected_failure_mode": observation.expected_failure_mode,
                    "taxon_overlap_policy": observation.taxon_overlap_policy or "",
                    "shared_taxa": "|".join(observation.shared_taxa),
                    "left_only_taxa": "|".join(observation.left_only_taxa),
                    "right_only_taxa": "|".join(observation.right_only_taxa),
                    "passed": str(observation.passed).lower(),
                    "mismatch_kind": observation.mismatch_kind or "",
                    "expected_output": json.dumps(
                        observation.expected_output, sort_keys=True
                    ),
                    "observed_output": json.dumps(
                        observation.observed_output, sort_keys=True
                    ),
                }
            )
    return path
