from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
import json
from importlib import metadata
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

from bijux_phylogenetics.diagnostics.validation import inspect_tree_path
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree


@dataclass(frozen=True, slots=True)
class ApeParityCase:
    """One governed live `ape` parity case."""

    case_id: str
    function_name: str
    python_function_name: str
    operation: str
    input_fixture: Path
    tolerance: float


@dataclass(frozen=True, slots=True)
class ApeParityObservation:
    """One live parity comparison between Bijux and `ape`."""

    case_id: str
    function_name: str
    python_function_name: str
    input_fixture: Path
    tolerance: float
    r_version: str | None
    ape_version: str | None
    bijux_version: str
    bijux_commit: str | None
    status: str
    passed: bool
    mismatch_reason: str | None
    reproducible_artifact_root: Path | None
    reference_summary: dict[str, object] | None
    bijux_summary: dict[str, object] | None


@dataclass(frozen=True, slots=True)
class ApeParitySummaryRow:
    """One function-level summary across governed `ape` parity cases."""

    function_name: str
    case_count: int
    passed_case_count: int
    failed_case_count: int
    skipped_case_count: int


@dataclass(slots=True)
class ApeParityReport:
    """Aggregate report for governed live `ape` parity cases."""

    observations: list[ApeParityObservation]
    summary_rows: list[ApeParitySummaryRow]
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


def _fixtures_root() -> Path:
    return _package_root() / "tests" / "fixtures"


def _ape_runner_path() -> Path:
    return (
        Path(__file__).resolve().parent
        / "resources"
        / "reference"
        / "ape_parity_runner.R"
    )


def _failure_root() -> Path:
    return _repository_root() / "artifacts" / "ape-parity-failures"


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


def list_ape_parity_cases(fixtures_root: Path | None = None) -> list[ApeParityCase]:
    """Return the governed live `ape` parity cases."""
    root = _fixtures_root() if fixtures_root is None else fixtures_root
    return [
        ApeParityCase(
            case_id="read-tree-example-rooted",
            function_name="ape::read.tree",
            python_function_name="load_tree+inspect_tree_path",
            operation="read-tree-summary",
            input_fixture=root / "trees" / "example_tree.nwk",
            tolerance=0.0,
        ),
        ApeParityCase(
            case_id="read-tree-example-unrooted",
            function_name="ape::read.tree",
            python_function_name="load_tree+inspect_tree_path",
            operation="read-tree-summary",
            input_fixture=root / "trees" / "example_tree_unrooted.nwk",
            tolerance=0.0,
        ),
    ]


def _build_case_lookup(fixtures_root: Path | None = None) -> dict[str, ApeParityCase]:
    return {case.case_id: case for case in list_ape_parity_cases(fixtures_root)}


def _selected_cases(
    *,
    case_ids: list[str] | None,
    fixtures_root: Path | None = None,
) -> list[ApeParityCase]:
    cases = _build_case_lookup(fixtures_root)
    if case_ids is None:
        return list(cases.values())
    selected: list[ApeParityCase] = []
    for case_id in case_ids:
        try:
            selected.append(cases[case_id])
        except KeyError as error:
            supported = ", ".join(sorted(cases))
            raise ValueError(
                f"unsupported ape parity case '{case_id}'; expected one of: {supported}"
            ) from error
    return selected


def _write_case_file(path: Path, case: ApeParityCase) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            {
                "case_id": case.case_id,
                "function_name": case.function_name,
                "operation": case.operation,
                "input_fixture": str(case.input_fixture),
                "tolerance": case.tolerance,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return path


def _build_bijux_tree_summary(input_fixture: Path) -> tuple[dict[str, object], list[dict[str, object]], str]:
    tree = load_tree(input_fixture)
    inspection = inspect_tree_path(input_fixture)
    summary = {
        "tip_count": inspection.tip_count,
        "internal_node_count": inspection.internal_node_count,
        "edge_count": inspection.edge_count,
        "rooted": inspection.rooted,
        "tip_labels": tree.tip_names,
        "branch_length_count": sum(
            1 for branch_length in tree.branch_lengths() if branch_length is not None
        ),
    }
    tips = [
        {"position": index, "label": label}
        for index, label in enumerate(tree.tip_names, start=1)
    ]
    return summary, tips, dumps_newick(tree)


def _load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _load_tip_table(path: Path) -> list[dict[str, object]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return [
        {"position": int(row["position"]), "label": row["label"]}
        for row in rows
    ]


def _compare_scalar(expected: object, observed: object, *, tolerance: float) -> bool:
    if isinstance(expected, (int, float)) and isinstance(observed, (int, float)):
        return abs(float(expected) - float(observed)) <= tolerance
    return expected == observed


def _compare_json(expected: object, observed: object, *, tolerance: float) -> bool:
    if isinstance(expected, dict) and isinstance(observed, dict):
        if set(expected) != set(observed):
            return False
        return all(
            _compare_json(expected[key], observed[key], tolerance=tolerance)
            for key in expected
        )
    if isinstance(expected, list) and isinstance(observed, list):
        if len(expected) != len(observed):
            return False
        return all(
            _compare_json(left, right, tolerance=tolerance)
            for left, right in zip(expected, observed, strict=True)
        )
    return _compare_scalar(expected, observed, tolerance=tolerance)


def _canonical_newick(path: Path) -> str:
    return dumps_newick(load_tree(path))


def _copy_if_exists(source: Path, destination: Path) -> None:
    if source.exists():
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)


def _write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )


def _write_tip_table(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["position", "label"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _persist_failure_bundle(
    *,
    failure_root: Path,
    case: ApeParityCase,
    case_file: Path,
    execution_root: Path,
    execution_payload: dict[str, object] | None,
    reference_summary: dict[str, object] | None,
    bijux_summary: dict[str, object] | None,
    bijux_tips: list[dict[str, object]],
    bijux_newick: str,
    mismatch_reason: str,
) -> Path:
    artifact_root = failure_root / case.case_id
    artifact_root.mkdir(parents=True, exist_ok=True)
    _copy_if_exists(case_file, artifact_root / "case.json")
    _copy_if_exists(
        execution_root / "reference-execution.json",
        artifact_root / "reference-execution.json",
    )
    _copy_if_exists(
        execution_root / "summary.json",
        artifact_root / "reference-summary.json",
    )
    _copy_if_exists(
        execution_root / "tips.tsv",
        artifact_root / "reference-tips.tsv",
    )
    _copy_if_exists(
        execution_root / "normalized-tree.nwk",
        artifact_root / "reference-normalized-tree.nwk",
    )
    if execution_payload is not None:
        _write_json(artifact_root / "reference-execution.observed.json", execution_payload)
    if reference_summary is not None:
        _write_json(artifact_root / "reference-summary.observed.json", reference_summary)
    if bijux_summary is not None:
        _write_json(artifact_root / "bijux-summary.json", bijux_summary)
    _write_tip_table(artifact_root / "bijux-tips.tsv", bijux_tips)
    (artifact_root / "bijux-normalized-tree.nwk").write_text(
        f"{bijux_newick}\n", encoding="utf-8"
    )
    _write_json(
        artifact_root / "comparison.json",
        {
            "case_id": case.case_id,
            "function_name": case.function_name,
            "mismatch_reason": mismatch_reason,
        },
    )
    return artifact_root


def _summary_rows(observations: list[ApeParityObservation]) -> list[ApeParitySummaryRow]:
    rows: list[ApeParitySummaryRow] = []
    for function_name in sorted({item.function_name for item in observations}):
        selected = [item for item in observations if item.function_name == function_name]
        rows.append(
            ApeParitySummaryRow(
                function_name=function_name,
                case_count=len(selected),
                passed_case_count=sum(1 for item in selected if item.status == "passed"),
                failed_case_count=sum(1 for item in selected if item.status == "failed"),
                skipped_case_count=sum(1 for item in selected if item.status == "skipped"),
            )
        )
    return rows


def run_ape_parity_cases(
    *,
    case_ids: list[str] | None = None,
    rscript_executable: str = "Rscript",
    failure_root: Path | None = None,
    fixtures_root: Path | None = None,
) -> ApeParityReport:
    """Run governed live `ape` parity cases through the checked-in R runner."""
    selected = _selected_cases(case_ids=case_ids, fixtures_root=fixtures_root)
    observations: list[ApeParityObservation] = []
    active_failure_root = _failure_root() if failure_root is None else failure_root
    bijux_version = _bijux_version()
    bijux_commit = _bijux_commit()
    for case in selected:
        with tempfile.TemporaryDirectory(prefix=f"bijux-ape-parity-{case.case_id}-") as tmpdir:
            working_root = Path(tmpdir)
            case_file = _write_case_file(working_root / "case.json", case)
            execution_root = working_root / "reference"
            execution_root.mkdir(parents=True, exist_ok=True)
            bijux_summary, bijux_tips, bijux_newick = _build_bijux_tree_summary(
                case.input_fixture
            )
            execution_payload: dict[str, object] | None = None
            reference_summary: dict[str, object] | None = None
            reference_tips: list[dict[str, object]] | None = None
            status = "failed"
            mismatch_reason: str | None = None
            artifact_root: Path | None = None
            r_version: str | None = None
            ape_version: str | None = None
            process_stdout = ""
            process_stderr = ""
            try:
                process = subprocess.run(
                    [
                        rscript_executable,
                        str(_ape_runner_path()),
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
            if process is None:
                pass
            elif process.returncode != 0:
                mismatch_reason = "reference_execution_failed"
            else:
                execution_path = execution_root / "reference-execution.json"
                if not execution_path.exists():
                    mismatch_reason = "reference_execution_failed"
                else:
                    execution_payload = _load_json(execution_path)
                    r_version = execution_payload.get("r_version")  # type: ignore[assignment]
                    ape_version = execution_payload.get("ape_version")  # type: ignore[assignment]
                    execution_status = execution_payload.get("status")
                    if execution_status == "unavailable":
                        status = "skipped"
                        mismatch_reason = str(
                            execution_payload.get(
                                "mismatch_reason", "ape_package_unavailable"
                            )
                        )
                    elif execution_status != "ok":
                        mismatch_reason = str(
                            execution_payload.get(
                                "mismatch_reason", "reference_execution_failed"
                            )
                        )
                    else:
                        reference_summary = _load_json(execution_root / "summary.json")
                        reference_tips = _load_tip_table(execution_root / "tips.tsv")
                        reference_newick = _canonical_newick(
                            execution_root / "normalized-tree.nwk"
                        )
                        if not _compare_json(
                            reference_summary, bijux_summary, tolerance=case.tolerance
                        ):
                            mismatch_reason = "summary_mismatch"
                        elif reference_tips != bijux_tips:
                            mismatch_reason = "tip_table_mismatch"
                        elif reference_newick != bijux_newick:
                            mismatch_reason = "normalized_tree_mismatch"
                        else:
                            status = "passed"
            if status != "passed":
                artifact_root = _persist_failure_bundle(
                    failure_root=active_failure_root,
                    case=case,
                    case_file=case_file,
                    execution_root=execution_root,
                    execution_payload=execution_payload,
                    reference_summary=reference_summary,
                    bijux_summary=bijux_summary,
                    bijux_tips=bijux_tips,
                    bijux_newick=bijux_newick,
                    mismatch_reason=mismatch_reason or "reference_execution_failed",
                )
                if process_stdout:
                    (artifact_root / "reference-stdout.txt").write_text(
                        process_stdout, encoding="utf-8"
                    )
                if process_stderr:
                    (artifact_root / "reference-stderr.txt").write_text(
                        process_stderr, encoding="utf-8"
                    )
            observations.append(
                ApeParityObservation(
                    case_id=case.case_id,
                    function_name=case.function_name,
                    python_function_name=case.python_function_name,
                    input_fixture=case.input_fixture,
                    tolerance=case.tolerance,
                    r_version=r_version,
                    ape_version=ape_version,
                    bijux_version=bijux_version,
                    bijux_commit=bijux_commit,
                    status=status,
                    passed=status == "passed",
                    mismatch_reason=mismatch_reason,
                    reproducible_artifact_root=artifact_root,
                    reference_summary=reference_summary,
                    bijux_summary=bijux_summary,
                )
            )
    case_count = len(observations)
    passed_case_count = sum(1 for item in observations if item.status == "passed")
    failed_case_count = sum(1 for item in observations if item.status == "failed")
    skipped_case_count = sum(1 for item in observations if item.status == "skipped")
    return ApeParityReport(
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
            "The governed live `ape` parity registry is intentionally narrow until later rounds expand the shared fixture surface.",
            "This harness requires Rscript plus the `ape` and `jsonlite` R packages for live reference execution.",
        ],
    )


def write_ape_parity_summary_table(path: Path, report: ApeParityReport) -> Path:
    """Write one row per governed `ape` function summary."""
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


def write_ape_parity_observation_table(path: Path, report: ApeParityReport) -> Path:
    """Write one row per governed `ape` parity observation."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "case_id",
                "function_name",
                "python_function_name",
                "input_fixture",
                "tolerance",
                "r_version",
                "ape_version",
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
                    "function_name": observation.function_name,
                    "python_function_name": observation.python_function_name,
                    "input_fixture": str(observation.input_fixture),
                    "tolerance": format(observation.tolerance, ".12g"),
                    "r_version": observation.r_version or "",
                    "ape_version": observation.ape_version or "",
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
