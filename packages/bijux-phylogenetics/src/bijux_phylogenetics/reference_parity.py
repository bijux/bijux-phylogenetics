from __future__ import annotations

import csv
from dataclasses import dataclass
import json
from pathlib import Path
import tempfile

from bijux_phylogenetics.compare.topology import (
    _unrooted_splits,
    compare_branch_score_distance,
    compare_robinson_foulds,
)
from bijux_phylogenetics.comparative import validate_comparative_reference_examples
from bijux_phylogenetics.diagnostics.validation import _load_tree
from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.tree_set import (
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
    return Path(__file__).resolve().parents[2] / "tests" / "fixtures"


def _expected_root() -> Path:
    return _fixtures_root() / "expected"


def _load_fixture_document(name: str) -> dict[str, object]:
    return json.loads((_expected_root() / name).read_text(encoding="utf-8"))


def _resolve_fixture_path(relative_path: str) -> Path:
    return Path(__file__).resolve().parents[2] / relative_path


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


def _build_comparative_observations() -> tuple[list[ReferenceParityObservation], dict[str, str]]:
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
        reference_tools[str(entry["reference_tool"])] = str(entry["reference_version"])
        observations.append(
            ReferenceParityObservation(
                suite="core",
                method=str(entry["method"]),
                case=str(entry["case"]),
                input_fixtures=[
                    _resolve_fixture_path(path) for path in entry["input_fixtures"]
                ],
                reference_tool=str(entry["reference_tool"]),
                reference_version=str(entry["reference_version"]),
                reference_source=str(entry["reference_source"]),
                tolerance=float(entry["tolerance"]),
                tolerance_reason=str(entry["tolerance_reason"]),
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
        for split in _unrooted_splits(tree, shared_taxa)
    )


def _build_tree_observation(
    entry: dict[str, object],
    *,
    suite: str,
) -> ReferenceParityObservation:
    input_paths = [_resolve_fixture_path(path) for path in entry["input_fixtures"]]
    method = str(entry["method"])
    expected_output = dict(entry["expected_output"])
    tolerance = float(entry["tolerance"])
    mismatch_kind: str | None = None

    if method == "robinson-foulds-distance":
        report = compare_robinson_foulds(
            input_paths[0],
            input_paths[1],
            rf_mode=str(entry.get("rf_mode", "unrooted")),
        )
        observed_output = {
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
        report = compare_branch_score_distance(input_paths[0], input_paths[1])
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
        report = compute_clade_frequency_table(input_paths[0])
        observed_output = {
            row.clade: row.frequency for row in report.clade_frequencies
        }
        passed = _numeric_outputs_match(
            expected_output,
            observed_output,
            tolerance=tolerance,
        )
        mismatch_kind = None if passed else "numerical_tolerance"
    elif method == "consensus-tree-generation":
        tree, _consensus = compute_consensus_tree(input_paths[0])
        with tempfile.TemporaryDirectory(prefix="bijux-reference-parity-") as tmp_dir:
            observed_path = Path(tmp_dir) / "observed-consensus.nwk"
            reference_path = Path(tmp_dir) / "reference-consensus.nwk"
            write_newick(observed_path, tree)
            reference_path.write_text(
                f"{expected_output['reference_consensus_newick']}\n",
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
        expected_output=expected_output,
        observed_output=observed_output,
        mismatch_kind=mismatch_kind,
        passed=passed,
    )


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
                    {f"{item.reference_tool} {item.reference_version}" for item in selected}
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
            "The extended suite adds a larger posterior tree-set parity check and is intended for optional validation lanes."
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
