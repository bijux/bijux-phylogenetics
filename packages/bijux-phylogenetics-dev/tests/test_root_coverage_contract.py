from __future__ import annotations

from configparser import ConfigParser
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]


def _coverage_config() -> ConfigParser:
    parser = ConfigParser()
    parser.read(REPO_ROOT / "configs" / "coveragerc.ini", encoding="utf-8")
    return parser


def _package_import_roots() -> set[str]:
    import_roots: set[str] = set()
    for source_dir in (REPO_ROOT / "packages").glob("*/src"):
        import_roots.update(
            child.name for child in source_dir.iterdir() if child.is_dir()
        )
    return import_roots


def test_root_coverage_configuration_matches_shared_python_baseline() -> None:
    coverage_config = _coverage_config()
    run_section = coverage_config["run"]
    report_section = coverage_config["report"]
    html_section = coverage_config["html"]

    assert run_section["branch"] == "True"
    assert run_section["parallel"] == "False"
    assert run_section["data_file"] == "artifacts/test/.coverage"
    assert {
        line.strip() for line in run_section["source"].splitlines() if line.strip()
    } == _package_import_roots()
    assert {
        line.strip() for line in run_section["omit"].splitlines() if line.strip()
    } == {
        "*/.tox/*",
        "*/.venv/*",
        "*/site-packages/*",
        "*/__main__.py",
        "tests/*",
    }

    assert report_section["show_missing"] == "True"
    assert report_section["skip_covered"] == "False"
    assert report_section["precision"] == "1"
    assert html_section["directory"] == "artifacts/test/htmlcov"
