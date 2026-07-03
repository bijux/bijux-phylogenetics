"""Audit classification coverage for long-running pytest cases."""

from __future__ import annotations

import argparse
import ast
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import re
from typing import Literal

DEFAULT_JSON_OUT = Path("artifacts/root/runtime-budget-report.json")
DEFAULT_MINIMUM_SECONDS = 100.0
_DURATION_PATTERN = re.compile(
    r"^(?P<duration>\d+\.\d+)s (?P<phase>\w+)\s+(?P<nodeid>packages/[^ ]+::[^\n]+)$",
    re.MULTILINE,
)
_CLASSIFICATION_CHOICES = ("split", "frozen", "slow", "optimize", "delete")


@dataclass(frozen=True, slots=True)
class RuntimeBudgetPolicy:
    """One durable classification rule for a family of long-running tests."""

    nodeid_prefix: str
    classification: Literal["split", "frozen", "slow", "optimize", "delete"]
    rationale: str


@dataclass(frozen=True, slots=True)
class RuntimeBudgetObservation:
    """One parsed pytest duration entry above the governed runtime threshold."""

    nodeid: str
    phase: str
    duration_seconds: float


@dataclass(frozen=True, slots=True)
class RuntimeBudgetIssue:
    """One classification or marker-governance issue."""

    code: str
    nodeid: str
    message: str


_RUNTIME_BUDGET_POLICIES: tuple[RuntimeBudgetPolicy, ...] = (
    RuntimeBudgetPolicy(
        nodeid_prefix=(
            "packages/bijux-phylogenetics/tests/test_geiger_parity_cli.py::"
        ),
        classification="split",
        rationale=(
            "the live geiger parity CLI lane exercises multiple responsibilities in one run and should remain outside routine maintainer checks until the workflow is separated into narrower assertions"
        ),
    ),
    RuntimeBudgetPolicy(
        nodeid_prefix=(
            "packages/bijux-phylogenetics/tests/test_phytools_parity_cli.py::"
        ),
        classification="split",
        rationale=(
            "the live phytools parity CLI lane couples external execution, payload rendering, and table governance, so it stays as an explicitly split slow surface"
        ),
    ),
    RuntimeBudgetPolicy(
        nodeid_prefix=(
            "packages/bijux-phylogenetics/tests/parity/geiger/test_discrete_runner.py::"
        ),
        classification="slow",
        rationale=(
            "the discrete geiger parity corpus is scientifically useful reference evidence, but its cross-language likelihood fits are inherently slow and belong on the explicit slow lane"
        ),
    ),
    RuntimeBudgetPolicy(
        nodeid_prefix=(
            "packages/bijux-phylogenetics/tests/test_comparative_model_figure_reference_validation.py::"
        ),
        classification="frozen",
        rationale=(
            "publication-facing comparative figure validation is governed reference evidence and should be preserved as a frozen audit surface rather than treated as a routine local check"
        ),
    ),
    RuntimeBudgetPolicy(
        nodeid_prefix=(
            "packages/bijux-phylogenetics/tests/test_comparative_model_figure_package.py::"
        ),
        classification="frozen",
        rationale=(
            "the comparative figure package build is a governed artifact review surface whose cost is justified by frozen reference coverage, not by routine repetition"
        ),
    ),
    RuntimeBudgetPolicy(
        nodeid_prefix=(
            "packages/bijux-phylogenetics/tests/test_avian_reproductive_trait_dataset.py::"
        ),
        classification="frozen",
        rationale=(
            "the avian reproductive workflow bundle tests materialize governed dataset outputs and are best preserved as explicit frozen evidence checks"
        ),
    ),
    RuntimeBudgetPolicy(
        nodeid_prefix=(
            "packages/bijux-phylogenetics/tests/test_discrete_mode_recovery.py::"
        ),
        classification="frozen",
        rationale=(
            "the discrete mode recovery workflow protects benchmark ledgers and governed review output, so it stays as a frozen slow lane instead of a routine package check"
        ),
    ),
    RuntimeBudgetPolicy(
        nodeid_prefix=(
            "packages/bijux-phylogenetics/tests/test_large_fasta_streaming.py::"
        ),
        classification="slow",
        rationale=(
            "the large FASTA streaming resource check is intentionally stress-oriented and should remain available as an explicit slow scalability surface"
        ),
    ),
    RuntimeBudgetPolicy(
        nodeid_prefix=(
            "packages/bijux-phylogenetics/tests/test_macroevolution_recovery_suite_dataset.py::"
        ),
        classification="frozen",
        rationale=(
            "the macroevolution recovery suite materializes governed workflow bundles and should be preserved as frozen evidence rather than rerun casually"
        ),
    ),
    RuntimeBudgetPolicy(
        nodeid_prefix=(
            "packages/bijux-phylogenetics/tests/test_discrete_mk.py::test_fit_discrete_mk_model_"
        ),
        classification="optimize",
        rationale=(
            "the owned discrete Mk boundary and parity fits expose important native behavior, but their current runtime should be reduced instead of accepted as the steady-state budget"
        ),
    ),
    RuntimeBudgetPolicy(
        nodeid_prefix=(
            "packages/bijux-phylogenetics/tests/test_real_dataset_macroevolution_benchmark.py::"
        ),
        classification="frozen",
        rationale=(
            "the real-dataset macroevolution benchmark governs stored reference evidence and should stay as an explicit frozen review surface"
        ),
    ),
    RuntimeBudgetPolicy(
        nodeid_prefix=(
            "packages/bijux-phylogenetics/tests/test_large_dataset_stress.py::"
        ),
        classification="slow",
        rationale=(
            "the large-dataset stress suite exists specifically to protect heavy resource envelopes and belongs on a deliberate slow lane"
        ),
    ),
    RuntimeBudgetPolicy(
        nodeid_prefix=(
            "packages/bijux-phylogenetics/tests/test_shared_phytools_comparative_fixtures.py::"
        ),
        classification="slow",
        rationale=(
            "the shared phytools comparative fixture catalog provides governed scientific parity evidence and remains an explicit slow validation lane"
        ),
    ),
    RuntimeBudgetPolicy(
        nodeid_prefix=(
            "packages/bijux-phylogenetics/tests/test_shared_geiger_continuous_fixtures.py::"
        ),
        classification="slow",
        rationale=(
            "the shared geiger continuous fixture catalog captures parity-sensitive signal-strength boundaries and remains intentionally slow"
        ),
    ),
    RuntimeBudgetPolicy(
        nodeid_prefix=(
            "packages/bijux-phylogenetics/tests/test_fasta_to_tree_real_workflows.py::"
        ),
        classification="frozen",
        rationale=(
            "the real workflow goldens are governed end-to-end evidence and should be preserved as frozen slow workflows rather than routine local checks"
        ),
    ),
    RuntimeBudgetPolicy(
        nodeid_prefix=(
            "packages/bijux-phylogenetics/tests/test_primate_comparative_dataset.py::"
        ),
        classification="frozen",
        rationale=(
            "the primate comparative dataset demo and bundle tests materialize governed outputs and belong on a frozen evidence lane"
        ),
    ),
    RuntimeBudgetPolicy(
        nodeid_prefix=(
            "packages/bijux-phylogenetics/tests/test_native_maximum_likelihood_accuracy_benchmark.py::"
        ),
        classification="optimize",
        rationale=(
            "native maximum-likelihood accuracy evidence is essential, but the owned benchmark should be optimized rather than accepted as a long-term 100+ second routine"
        ),
    ),
    RuntimeBudgetPolicy(
        nodeid_prefix=(
            "packages/bijux-phylogenetics/tests/test_reference_parity_extended.py::"
        ),
        classification="frozen",
        rationale=(
            "the extended reference parity corpus is governed scientific validation evidence and should be handled as a frozen audit surface"
        ),
    ),
    RuntimeBudgetPolicy(
        nodeid_prefix=(
            "packages/bijux-phylogenetics/tests/test_runtime.py::test_cli_evidence_book_"
        ),
        classification="split",
        rationale=(
            "the evidence-book CLI build surfaces combine broad repository work and should stay split from routine package checks until the heaviest subflows are isolated"
        ),
    ),
    RuntimeBudgetPolicy(
        nodeid_prefix=(
            "packages/bijux-phylogenetics/tests/test_discrete_mode_recovery_dataset.py::"
        ),
        classification="frozen",
        rationale=(
            "the discrete mode recovery dataset demos and bundle outputs are governed reference artifacts and belong on a frozen review lane"
        ),
    ),
    RuntimeBudgetPolicy(
        nodeid_prefix=(
            "packages/bijux-phylogenetics/tests/test_continuous_mode_recovery.py::"
        ),
        classification="frozen",
        rationale=(
            "the continuous mode recovery workflow preserves governed benchmark evidence and is best kept as an explicit frozen audit surface"
        ),
    ),
)


def _is_pytest_slow_marker(node: ast.AST) -> bool:
    return (
        isinstance(node, ast.Attribute)
        and node.attr == "slow"
        and isinstance(node.value, ast.Attribute)
        and node.value.attr == "mark"
        and isinstance(node.value.value, ast.Name)
        and node.value.value.id == "pytest"
    )


def _contains_slow_marker(node: ast.AST) -> bool:
    if _is_pytest_slow_marker(node):
        return True
    if isinstance(node, (ast.List, ast.Tuple, ast.Set)):
        return any(_contains_slow_marker(child) for child in node.elts)
    return False


def _module_has_slow_marker(module: ast.Module) -> bool:
    for node in module.body:
        if not isinstance(node, ast.Assign):
            continue
        if not any(
            isinstance(target, ast.Name) and target.id == "pytestmark"
            for target in node.targets
        ):
            continue
        if _contains_slow_marker(node.value):
            return True
    return False


def _function_has_slow_marker(module: ast.Module, function_name: str) -> bool:
    for node in module.body:
        if not isinstance(node, ast.FunctionDef) or node.name != function_name:
            continue
        return any(_is_pytest_slow_marker(decorator) for decorator in node.decorator_list)
    return False


def _is_slow_marked(repo_root: Path, nodeid: str) -> bool:
    module_path_text, function_name = nodeid.split("::", 1)
    module_path = repo_root / module_path_text
    if not module_path.exists():
        return False
    module = ast.parse(module_path.read_text(encoding="utf-8"))
    return _module_has_slow_marker(module) or _function_has_slow_marker(
        module, function_name
    )


def _matching_policy(nodeid: str) -> RuntimeBudgetPolicy | None:
    matching = [
        policy for policy in _RUNTIME_BUDGET_POLICIES if nodeid.startswith(policy.nodeid_prefix)
    ]
    if not matching:
        return None
    return max(matching, key=lambda policy: len(policy.nodeid_prefix))


def parse_runtime_budget_observations(
    log_path: Path,
    *,
    minimum_seconds: float = DEFAULT_MINIMUM_SECONDS,
) -> list[RuntimeBudgetObservation]:
    """Parse the long-running pytest duration rows from one captured log."""
    observations = [
        RuntimeBudgetObservation(
            nodeid=match.group("nodeid"),
            phase=match.group("phase"),
            duration_seconds=float(match.group("duration")),
        )
        for match in _DURATION_PATTERN.finditer(log_path.read_text(encoding="utf-8"))
        if match.group("phase") == "call"
        and float(match.group("duration")) >= minimum_seconds
    ]
    observations.sort(
        key=lambda observation: (-observation.duration_seconds, observation.nodeid)
    )
    return observations


def build_runtime_budget_report(
    repo_root: Path,
    log_path: Path,
    *,
    minimum_seconds: float = DEFAULT_MINIMUM_SECONDS,
) -> dict[str, object]:
    """Build the runtime-budget classification report for one pytest duration log."""
    repo_root = repo_root.resolve()
    log_path = log_path.resolve()
    observations = parse_runtime_budget_observations(
        log_path,
        minimum_seconds=minimum_seconds,
    )
    issues: list[RuntimeBudgetIssue] = []
    observed_rows: list[dict[str, object]] = []
    counts_by_classification = {
        classification: 0 for classification in _CLASSIFICATION_CHOICES
    }

    for observation in observations:
        policy = _matching_policy(observation.nodeid)
        slow_marked = _is_slow_marked(repo_root, observation.nodeid)
        classification = None if policy is None else policy.classification
        rationale = "" if policy is None else policy.rationale
        if policy is None:
            issues.append(
                RuntimeBudgetIssue(
                    code="missing-runtime-budget-classification",
                    nodeid=observation.nodeid,
                    message=(
                        "long-running test exceeds the governed runtime threshold but has no durable budget classification"
                    ),
                )
            )
        else:
            counts_by_classification[policy.classification] += 1
            if policy.classification != "delete" and not slow_marked:
                issues.append(
                    RuntimeBudgetIssue(
                        code="long-running-test-missing-slow-marker",
                        nodeid=observation.nodeid,
                        message=(
                            f"{policy.classification} runtime budget classification still requires pytest.mark.slow so the case stays out of default package checks"
                        ),
                    )
                )
        observed_rows.append(
            {
                "nodeid": observation.nodeid,
                "phase": observation.phase,
                "duration_seconds": observation.duration_seconds,
                "classification": classification,
                "slow_marked": slow_marked,
                "rationale": rationale,
            }
        )

    return {
        "schema_version": 1,
        "minimum_seconds": minimum_seconds,
        "log_path": str(log_path),
        "policy_count": len(_RUNTIME_BUDGET_POLICIES),
        "observation_count": len(observations),
        "classified_observation_count": sum(
            1 for row in observed_rows if row["classification"] is not None
        ),
        "counts_by_classification": counts_by_classification,
        "issue_count": len(issues),
        "issues": [asdict(issue) for issue in issues],
        "observations": observed_rows,
    }


def _write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def check_runtime_budgets(
    repo_root: Path,
    log_path: Path,
    *,
    minimum_seconds: float = DEFAULT_MINIMUM_SECONDS,
    json_out: Path | None = None,
) -> dict[str, object]:
    """Raise when the governed long-running tests are not fully classified."""
    payload = build_runtime_budget_report(
        repo_root,
        log_path,
        minimum_seconds=minimum_seconds,
    )
    if json_out is not None:
        _write_json(json_out, payload)
    if payload["issue_count"]:
        raise SystemExit("runtime budget audit failed")
    return payload


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the runtime-budget audit."""
    parser = argparse.ArgumentParser(
        description="Audit classification coverage for long-running pytest cases."
    )
    parser.add_argument("command", choices=("report", "check"))
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--log-path", required=True)
    parser.add_argument("--minimum-seconds", type=float, default=DEFAULT_MINIMUM_SECONDS)
    parser.add_argument("--json-out", default=str(DEFAULT_JSON_OUT))
    return parser.parse_args()


def main() -> int:
    """Run the runtime-budget CLI entry point."""
    args = parse_args()
    repo_root = Path(args.repo_root).resolve()
    log_path = Path(args.log_path).resolve()
    json_out = Path(args.json_out)
    if not json_out.is_absolute():
        json_out = repo_root / json_out
    if args.command == "check":
        payload = check_runtime_budgets(
            repo_root,
            log_path,
            minimum_seconds=args.minimum_seconds,
            json_out=json_out,
        )
    else:
        payload = build_runtime_budget_report(
            repo_root,
            log_path,
            minimum_seconds=args.minimum_seconds,
        )
        _write_json(json_out, payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
