from __future__ import annotations

from pathlib import Path
import shutil
import stat
import subprocess

from bijux_phylogenetics_dev.quality.frozen_test_report import (
    build_output_manifest,
    build_frozen_test_report,
    build_scientific_drift_records,
    diff_output_manifests,
    diff_test_records,
    parse_junit_report,
)


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _git() -> str:
    resolved = shutil.which("git")
    if resolved is None:
        raise AssertionError("git executable not available for test")
    return resolved


def _run(command: list[str], *, cwd: Path) -> None:
    subprocess.run(command, cwd=cwd, check=True)


def _synthetic_repo(tmp_path: Path) -> tuple[Path, str]:
    repo_root = tmp_path / "repo"
    repo_root.mkdir()
    _run([_git(), "init"], cwd=repo_root)
    _run([_git(), "config", "user.name", "Bijux Tests"], cwd=repo_root)
    _run([_git(), "config", "user.email", "tests@bijux.io"], cwd=repo_root)
    _write(repo_root / "packages" / "demo-package" / "placeholder.txt", "demo\n")
    _write(repo_root / "makes" / "packages" / "demo-package.mk", "test-all:\n\t@true\n")
    _write(repo_root / "science.txt", "frozen\n")
    _run([_git(), "add", "."], cwd=repo_root)
    _run([_git(), "commit", "-m", "seed frozen fixture"], cwd=repo_root)
    frozen_ref = (
        subprocess.run(
            [_git(), "rev-parse", "HEAD"],
            cwd=repo_root,
            check=True,
            capture_output=True,
            text=True,
        )
        .stdout.strip()
    )
    _write(repo_root / "science.txt", "current\n")
    _run([_git(), "add", "science.txt"], cwd=repo_root)
    _run([_git(), "commit", "-m", "change scientific output"], cwd=repo_root)
    return repo_root, frozen_ref


def _fake_make_executable(tmp_path: Path) -> str:
    script = tmp_path / "fake-make"
    script.write_text(
        """#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import sys


def _arg_value(args: list[str], prefix: str) -> str:
    for arg in args:
        if arg.startswith(prefix):
            return arg.split("=", 1)[1]
    raise SystemExit(f"missing {prefix}")


args = sys.argv[1:]
repo_root = Path(_arg_value(args, "MONOREPO_ROOT="))
artifacts_root = Path(_arg_value(args, "PROJECT_ARTIFACTS_DIR="))
science_value = (repo_root / "science.txt").read_text(encoding="utf-8").strip()
junit_xml = artifacts_root / "test" / "junit.xml"
junit_xml.parent.mkdir(parents=True, exist_ok=True)
outcome_tag = ""
suite_attrs = 'tests="1"'
if science_value == "current":
    outcome_tag = '<failure message="scientific drift">observed output changed</failure>'
    suite_attrs = 'tests="1" failures="1"'
junit_xml.write_text(
    f'<testsuites><testsuite name="pytest" {suite_attrs}><testcase file="tests/test_demo.py" name="test_science">{outcome_tag}</testcase></testsuite></testsuites>\\n',
    encoding="utf-8",
)
output_path = artifacts_root / "tree-report.json"
output_path.parent.mkdir(parents=True, exist_ok=True)
output_path.write_text(f'{{"science": "{science_value}"}}\\n', encoding="utf-8")
print(f"ran {science_value}")
if science_value == "current":
    raise SystemExit(1)
""",
        encoding="utf-8",
    )
    script.chmod(script.stat().st_mode | stat.S_IXUSR)
    return script.as_posix()


def test_parse_junit_report_normalizes_pytest_testcase_outcomes(tmp_path: Path) -> None:
    junit_xml = tmp_path / "junit.xml"
    junit_xml.write_text(
        """
<testsuites>
  <testsuite name="pytest" tests="3" failures="1" skipped="1">
    <testcase classname="pkg.module" file="tests/test_alpha.py" name="test_pass" />
    <testcase classname="pkg.module" file="tests/test_alpha.py" name="test_fail">
      <failure message="assertion failed">expected parity</failure>
    </testcase>
    <testcase classname="pkg.module" name="test_skip">
      <skipped message="external engine unavailable" />
    </testcase>
  </testsuite>
</testsuites>
""".strip()
        + "\n",
        encoding="utf-8",
    )

    records = parse_junit_report(junit_xml)

    assert records["tests/test_alpha.py::test_pass"].outcome == "passed"
    assert records["tests/test_alpha.py::test_fail"].outcome == "failed"
    assert "expected parity" in records["tests/test_alpha.py::test_fail"].message
    assert records["pkg.module::test_skip"].outcome == "skipped"


def test_build_output_manifest_ignores_harness_files_and_hashes_governed_outputs(
    tmp_path: Path,
) -> None:
    _write(tmp_path / "tree-report.json", '{"likelihood": -12.3}\n')
    _write(tmp_path / "benchmarks" / "bench.json", "{}\n")
    _write(tmp_path / "htmlcov" / "index.html", "<html></html>\n")
    _write(tmp_path / "nested" / "support.tsv", "clade\tsupport\nA\t0.8\n")
    _write(tmp_path / "junit.xml", "<testsuites />\n")
    _write(tmp_path / ".coverage", "coverage-data\n")

    manifest = build_output_manifest(tmp_path)

    assert set(manifest) == {"nested/support.tsv", "tree-report.json"}
    assert manifest["tree-report.json"].size_bytes > 0
    assert len(manifest["tree-report.json"].sha256) == 64


def test_diff_helpers_detect_changed_tests_outputs_and_scientific_drift(
    tmp_path: Path,
) -> None:
    current_junit = tmp_path / "current-junit.xml"
    frozen_junit = tmp_path / "frozen-junit.xml"
    current_junit.write_text(
        """
<testsuites>
  <testsuite name="pytest" tests="2" failures="1">
    <testcase file="tests/test_beta.py" name="test_same" />
    <testcase file="tests/test_beta.py" name="test_changed">
      <failure message="likelihood drift">current result changed</failure>
    </testcase>
  </testsuite>
</testsuites>
""".strip()
        + "\n",
        encoding="utf-8",
    )
    frozen_junit.write_text(
        """
<testsuites>
  <testsuite name="pytest" tests="2">
    <testcase file="tests/test_beta.py" name="test_same" />
    <testcase file="tests/test_beta.py" name="test_changed" />
  </testsuite>
</testsuites>
""".strip()
        + "\n",
        encoding="utf-8",
    )

    current_outputs = tmp_path / "current"
    frozen_outputs = tmp_path / "frozen"
    _write(current_outputs / "report.json", '{"score": 2}\n')
    _write(current_outputs / "only-current.tsv", "present\n")
    _write(frozen_outputs / "report.json", '{"score": 1}\n')
    _write(frozen_outputs / "only-frozen.tsv", "present\n")

    test_differences = diff_test_records(
        parse_junit_report(current_junit),
        parse_junit_report(frozen_junit),
    )
    output_differences = diff_output_manifests(
        build_output_manifest(current_outputs),
        build_output_manifest(frozen_outputs),
    )
    drift_records = build_scientific_drift_records(
        test_differences, output_differences
    )

    assert [difference.nodeid for difference in test_differences] == [
        "tests/test_beta.py::test_changed"
    ]
    assert {difference.relative_path for difference in output_differences} == {
        "only-current.tsv",
        "only-frozen.tsv",
        "report.json",
    }
    assert {record.kind for record in drift_records} == {"output", "test"}
    assert any(record.identifier == "report.json" for record in drift_records)


def test_build_frozen_test_report_compares_current_run_against_frozen_ref(
    tmp_path: Path,
) -> None:
    repo_root, frozen_ref = _synthetic_repo(tmp_path)
    artifacts_root = tmp_path / "artifacts"
    _write(artifacts_root / "current" / "stale.txt", "stale\n")
    _write(artifacts_root / "frozen" / "stale.txt", "stale\n")
    payload = build_frozen_test_report(
        repo_root,
        package="demo-package",
        frozen_ref=frozen_ref,
        artifacts_root=artifacts_root,
        full_suite_worker_count=6,
        git_executable=_git(),
        make_executable=_fake_make_executable(tmp_path),
    )

    assert payload["current_ref"]["input_ref"] == "HEAD"
    assert payload["frozen_ref"]["resolved_ref"] == frozen_ref
    assert payload["current_run_succeeded"] is False
    assert payload["frozen_run_succeeded"] is True
    assert payload["changed_test_count"] == 1
    assert payload["changed_output_count"] == 1
    assert payload["scientific_drift_count"] == 2
    assert payload["current_run"]["full_suite_worker_count"] == 6
    assert payload["frozen_run"]["full_suite_worker_count"] == 6
    assert "FULL_SUITE_PYTEST_WORKERS=6" in payload["current_run"]["command"]
    assert "FULL_SUITE_PYTEST_WORKERS=6" in payload["frozen_run"]["command"]
    assert payload["changed_outputs"][0]["relative_path"] == "tree-report.json"
    assert all(
        difference["relative_path"] != "stale.txt"
        for difference in payload["changed_outputs"]
    )
