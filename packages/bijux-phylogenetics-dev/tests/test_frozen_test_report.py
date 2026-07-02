from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics_dev.quality.frozen_test_report import (
    build_output_manifest,
    build_scientific_drift_records,
    diff_output_manifests,
    diff_test_records,
    parse_junit_report,
)


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


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
