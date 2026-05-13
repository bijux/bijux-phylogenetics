from __future__ import annotations

from configparser import ConfigParser
from pathlib import Path
import tomllib
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[3]
CONFIGS_ROOT = REPO_ROOT / "configs"
SHARED_PYTEST_MARKERS = {
    "api: HTTP API tests (manual, not for CI)",
    "e2e: end-to-end tests",
    "evaluation: evaluation benchmarks (deterministic, no regressions)",
    "gpu: requires CUDA",
    "integration: integration tests",
    "live: live provider integration",
    "real: real local model tests (slow, manual, not for CI)",
    "real_local: requires local models or hardware",
    "regression: regression tests",
    "slow: mark test as slow",
    "smoke: smoke tests",
    "unit: unit tests",
    "windows: mark tests for Windows-only",
}
REPOSITORY_PYTEST_MARKERS = {
    "engine_contract: fake-executable and parser contract tests for external engines",
    "engine_real: real executable integration tests for external engines",
    "scientific_validation: slower workflow checks against governed reference outputs",
    "stress_heavy: governed heavy stress tier for large-dataset resource checks",
    "stress_small: governed small stress tier for routine large-dataset resource checks",
}


def _config_parser(path: Path) -> ConfigParser:
    parser = ConfigParser()
    parser.read(path, encoding="utf-8")
    return parser


def _ruff_config() -> dict[str, Any]:
    with (CONFIGS_ROOT / "ruff.toml").open("rb") as handle:
        return tomllib.load(handle)


def _resolve_config_relative_path(config_path: Path, configured_path: str) -> Path:
    return (config_path.parent / configured_path).resolve()


def _package_roots(kind: str) -> set[str]:
    return {
        path.relative_to(REPO_ROOT).as_posix()
        for path in (REPO_ROOT / "packages").glob(f"*/{kind}")
    }


def _package_import_roots() -> set[str]:
    import_roots: set[str] = set()
    for source_dir in (REPO_ROOT / "packages").glob("*/src"):
        import_roots.update(
            child.name for child in source_dir.iterdir() if child.is_dir()
        )
    return import_roots


def test_root_pytest_configuration_matches_shared_python_baseline() -> None:
    pytest_ini = CONFIGS_ROOT / "pytest.ini"
    pytest_config = _config_parser(pytest_ini)["pytest"]

    assert pytest_config["minversion"] == "8.0"
    assert pytest_config["python_files"] == "test_*.py"
    assert pytest_config["python_classes"] == "Test*"
    assert pytest_config["python_functions"] == "test_*"
    assert pytest_config["asyncio_mode"] == "auto"
    assert pytest_config["cache_dir"] == "../artifacts/root/pytest-cache"
    assert (
        _resolve_config_relative_path(
            pytest_ini,
            pytest_config["cache_dir"],
        )
        == REPO_ROOT / "artifacts" / "root" / "pytest-cache"
    )
    assert pytest_config["timeout"] == "120"
    assert pytest_config["timeout_method"] == "thread"
    assert pytest_config["timeout_func_only"] == "true"
    assert pytest_config["xfail_strict"] == "true"

    assert {
        line.strip()
        for line in pytest_config["norecursedirs"].splitlines()
        if line.strip()
    } == {
        ".venv",
        ".tox",
        ".pytest_cache",
        ".ruff_cache",
        ".mypy_cache",
        ".hypothesis",
        ".benchmarks",
        "build",
        "dist",
        "htmlcov",
        "docs",
        "artifacts",
        "node_modules",
        "site",
    }
    assert [
        line.strip() for line in pytest_config["addopts"].splitlines() if line.strip()
    ] == [
        "-ra",
        "--import-mode=importlib",
        "--strict-markers",
        "--tb=short",
    ]
    configured_markers = {
        line.strip() for line in pytest_config["markers"].splitlines() if line.strip()
    }
    assert SHARED_PYTEST_MARKERS.issubset(configured_markers)
    assert configured_markers - SHARED_PYTEST_MARKERS == REPOSITORY_PYTEST_MARKERS
    assert [
        line.strip()
        for line in pytest_config["filterwarnings"].splitlines()
        if line.strip()
    ] == [
        "ignore:Not saving anything, no benchmarks have been run!",
        "ignore:jsonschema\\.exceptions\\.RefResolutionError is deprecated:DeprecationWarning",
        "ignore:jsonschema\\.exceptions\\.RefResolutionError is deprecated:DeprecationWarning:schemathesis.generation.coverage",
        "ignore:.*forkpty.*:DeprecationWarning",
        "ignore:datetime\\.datetime\\.utcnow\\(\\) is deprecated:DeprecationWarning",
        "ignore:'asyncio\\.iscoroutinefunction' is deprecated:DeprecationWarning",
        "ignore:'asyncio\\.get_event_loop_policy' is deprecated:DeprecationWarning",
        "ignore:'HTTP_422_UNPROCESSABLE_ENTITY' is deprecated:DeprecationWarning:anyio._backends._asyncio",
    ]


def test_root_ruff_configuration_matches_shared_python_baseline() -> None:
    ruff_config = _ruff_config()

    assert ruff_config["target-version"] == "py311"
    assert ruff_config["line-length"] == 88
    assert ruff_config["respect-gitignore"] is True
    assert ruff_config["cache-dir"] == "../artifacts/root/ruff-cache"
    assert (
        _resolve_config_relative_path(
            CONFIGS_ROOT / "ruff.toml",
            ruff_config["cache-dir"],
        )
        == REPO_ROOT / "artifacts" / "root" / "ruff-cache"
    )
    assert set(ruff_config["src"]) == _package_roots("src") | _package_roots("tests")
    assert ruff_config["exclude"] == [
        ".git",
        ".hg",
        ".mypy_cache",
        ".pytest_cache",
        ".ruff_cache",
        ".tox",
        ".venv",
        "artifacts",
        "build",
        "dist",
        "docs/report",
        "htmlcov",
        "__pycache__",
        "migrations",
        "node_modules",
        "*.egg-info",
        "site",
    ]

    lint = ruff_config["lint"]
    assert lint["select"] == [
        "E",
        "F",
        "I",
        "B",
        "UP",
        "SIM",
        "C4",
        "PIE",
        "RET",
        "ISC",
    ]
    assert lint["ignore"] == ["E501", "E203"]
    assert lint["per-file-ignores"] == {"__init__.py": ["F401"]}
    assert lint["isort"]["force-sort-within-sections"] is True
    assert set(lint["isort"]["known-first-party"]) == _package_import_roots() | {
        "tests"
    }
    assert lint["mccabe"]["max-complexity"] == 10


def test_root_mypy_configuration_matches_shared_python_baseline() -> None:
    mypy_config = _config_parser(REPO_ROOT / "configs" / "mypy.ini")
    root_mypy = mypy_config["mypy"]

    assert root_mypy["python_version"] == "3.11"
    assert root_mypy["strict"] == "true"
    assert root_mypy["pretty"] == "true"
    assert root_mypy["show_error_codes"] == "true"
    assert root_mypy["warn_unreachable"] == "true"
    assert root_mypy["warn_unused_configs"] == "true"
    assert root_mypy["warn_unused_ignores"] == "true"
    assert root_mypy["namespace_packages"] == "true"
    assert root_mypy["plugins"] == "pydantic.mypy"
    assert root_mypy["exclude"] == (
        "^(\\.venv|build|dist|docs|htmlcov|\\.mypy_cache|\\.pytest_cache|"
        "\\.ruff_cache|\\.tox|__pycache__|migrations|\\.egg-info|node_modules|"
        "artifacts|site)/"
    )

    configured_files = {
        entry.strip() for entry in root_mypy["files"].split(",") if entry.strip()
    }
    assert configured_files == _package_roots("src") | _package_roots("tests")

    configured_paths = {
        entry.strip() for entry in root_mypy["mypy_path"].split(":") if entry.strip()
    }
    assert configured_paths == _package_roots("src")
