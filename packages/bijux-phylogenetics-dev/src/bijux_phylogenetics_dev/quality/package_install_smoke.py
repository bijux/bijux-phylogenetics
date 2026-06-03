"""Run installability smoke tests against built runtime distribution artifacts."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import os
from pathlib import Path
import shutil
import sys
from typing import Any

from ..trusted_process import run_text

DEFAULT_ARTIFACTS_ROOT = Path("artifacts/root/package-install-smoke")
_ARTIFACT_KIND_PATTERNS = {
    "wheel": "bijux_phylogenetics-*.whl",
    "sdist": "bijux_phylogenetics-*.tar.gz",
}
_RUNTIME_SOURCE_ROOT = (
    Path("packages") / "bijux-phylogenetics" / "src" / "bijux_phylogenetics"
)
_RESOURCE_SENTINELS: dict[str, tuple[str, ...]] = {
    "example_alignment": ("examples", "alignments", "example_alignment.fasta"),
    "example_tree": ("examples", "trees", "example_tree.nwk"),
    "primate_tree": ("datasets", "mammals", "primate_comparative", "tree.nwk"),
    "primate_traits": ("datasets", "mammals", "primate_comparative", "traits.csv"),
    "rabies_workflow_config": (
        "datasets",
        "pathogens",
        "rabies_cross_host_geography_panel",
        "workflow-config.json",
    ),
}
_EXPECTED_EXAMPLE_INPUTS = ("alignment", "alt_tree", "metadata", "traits", "tree")
_TREE_PACKAGE_OUTPUTS = (
    "tree-report.html",
    "tree-image.svg",
    "support-table.tsv",
    "clade-table.tsv",
    "branch-stats.tsv",
    "tree-report.manifest.json",
)

JsonObject = dict[str, Any]


@dataclass(frozen=True)
class SmokeIssue:
    """Describe one installability smoke-test failure."""

    code: str
    path: str
    message: str


def _artifact_kind_choices(value: str) -> tuple[str, ...]:
    if value == "both":
        return ("wheel", "sdist")
    if value in _ARTIFACT_KIND_PATTERNS:
        return (value,)
    raise ValueError(f"unsupported artifact kind: {value}")


def _as_dict(value: object) -> JsonObject:
    return value if isinstance(value, dict) else {}


def _as_list(value: object) -> list[object]:
    return value if isinstance(value, list) else []


def _as_str(value: object) -> str:
    return value if isinstance(value, str) else ""


def _write_json(path: Path, payload: JsonObject) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _smoke_env() -> dict[str, str]:
    env = dict(os.environ)
    for key in ("PYTHONPATH", "PYTHONHOME", "VIRTUAL_ENV", "__PYVENV_LAUNCHER__"):
        env.pop(key, None)
    env["PIP_DISABLE_PIP_VERSION_CHECK"] = "1"
    return env


def _venv_executable(venv_root: Path, name: str) -> Path:
    bin_dir = "Scripts" if os.name == "nt" else "bin"
    suffix = ".exe" if os.name == "nt" else ""
    return venv_root / bin_dir / f"{name}{suffix}"


def _select_newest_artifact(dist_dir: Path, pattern: str, artifact_kind: str) -> Path:
    matches = sorted(dist_dir.glob(pattern))
    if not matches:
        raise ValueError(
            f"expected at least one {artifact_kind} artifact for pattern {pattern} in {dist_dir}"
        )
    return max(matches, key=lambda path: (path.stat().st_mtime_ns, path.name))


def select_artifact_paths(dist_dir: Path, artifact_kind: str) -> list[tuple[str, Path]]:
    """Select the built runtime artifacts to smoke-test."""
    selections: list[tuple[str, Path]] = []
    for kind in _artifact_kind_choices(artifact_kind):
        pattern = _ARTIFACT_KIND_PATTERNS[kind]
        selections.append((kind, _select_newest_artifact(dist_dir, pattern, kind)))
    return selections


def _probe_installed_resources(
    venv_python: Path, cwd: Path, env: dict[str, str]
) -> JsonObject:
    sentinel_json = json.dumps(_RESOURCE_SENTINELS, sort_keys=True)
    probe_script = f"""
import importlib.resources
import json
from pathlib import Path
import bijux_phylogenetics

sentinels = json.loads({sentinel_json!r})
package_root = Path(bijux_phylogenetics.__file__).resolve().parent
resource_root = importlib.resources.files("bijux_phylogenetics").joinpath("resources")
resource_paths = {{}}
for name, parts in sentinels.items():
    candidate = resource_root
    for part in parts:
        candidate = candidate.joinpath(part)
    resource_paths[name] = str(Path(str(candidate)).resolve())
missing_resources = [
    name for name, path in resource_paths.items() if not Path(path).exists()
]
print(
    json.dumps(
        {{
            "package_root": str(package_root),
            "resource_root": str(Path(str(resource_root)).resolve()),
            "resource_paths": resource_paths,
            "missing_resources": missing_resources,
        }},
        sort_keys=True,
    )
)
"""
    completed = run_text(
        [venv_python, "-c", probe_script],
        cwd=cwd,
        check=True,
        capture_output=True,
        env=env,
    )
    return _as_dict(json.loads(completed.stdout))


def _probe_copied_example_inputs(
    venv_python: Path,
    *,
    cwd: Path,
    env: dict[str, str],
    destination: Path,
) -> JsonObject:
    probe_script = """
import json
from pathlib import Path
import sys

from bijux_phylogenetics.core import copy_example_inputs

copied = copy_example_inputs(Path(sys.argv[1]))
print(json.dumps({name: str(path.resolve()) for name, path in copied.items()}, sort_keys=True))
"""
    completed = run_text(
        [venv_python, "-c", probe_script, str(destination)],
        cwd=cwd,
        check=True,
        capture_output=True,
        env=env,
    )
    copied_paths = _as_dict(json.loads(completed.stdout))
    return {
        "destination": destination.as_posix(),
        "copied_paths": copied_paths,
    }


def validate_resource_probe(report: JsonObject, repo_root: Path) -> list[SmokeIssue]:
    """Validate that packaged resources resolve from the installed artifact."""
    issues: list[SmokeIssue] = []
    source_root = (repo_root / _RUNTIME_SOURCE_ROOT).resolve()
    package_root_text = _as_str(report.get("package_root"))
    package_root = Path(package_root_text) if package_root_text else None
    if package_root is None or not package_root.exists():
        issues.append(
            SmokeIssue(
                code="missing-package-root",
                path=package_root_text,
                message="installed package probe did not resolve a package root",
            )
        )
    elif package_root == source_root or source_root in package_root.parents:
        issues.append(
            SmokeIssue(
                code="source-tree-package-root",
                path=package_root.as_posix(),
                message="installed smoke test imported the runtime from the repository source tree",
            )
        )
    resource_paths = {
        name: _as_str(path)
        for name, path in _as_dict(report.get("resource_paths")).items()
    }
    for sentinel_name in sorted(_RESOURCE_SENTINELS):
        resource_path_text = resource_paths.get(sentinel_name, "")
        if not resource_path_text:
            issues.append(
                SmokeIssue(
                    code="missing-resource-sentinel",
                    path=sentinel_name,
                    message="installed resource probe omitted a required packaged resource",
                )
            )
            continue
        resource_path = Path(resource_path_text)
        if not resource_path.exists():
            issues.append(
                SmokeIssue(
                    code="missing-installed-resource",
                    path=resource_path_text,
                    message="installed package does not contain one required packaged resource",
                )
            )
        elif resource_path == source_root or source_root in resource_path.parents:
            issues.append(
                SmokeIssue(
                    code="source-tree-resource-path",
                    path=resource_path.as_posix(),
                    message="resource probe resolved a repository source-tree path instead of an installed package resource",
                )
            )
    for missing in _as_list(report.get("missing_resources")):
        issues.append(
            SmokeIssue(
                code="reported-missing-installed-resource",
                path=_as_str(missing),
                message="resource probe reported a missing packaged resource",
            )
        )
    return issues


def validate_example_input_probe(
    report: JsonObject, repo_root: Path
) -> list[SmokeIssue]:
    """Validate that installed example-copy helpers create usable writable inputs."""
    issues: list[SmokeIssue] = []
    source_root = (repo_root / _RUNTIME_SOURCE_ROOT).resolve()
    copied_paths = {
        name: _as_str(path)
        for name, path in _as_dict(report.get("copied_paths")).items()
    }
    for input_name in _EXPECTED_EXAMPLE_INPUTS:
        path_text = copied_paths.get(input_name, "")
        if not path_text:
            issues.append(
                SmokeIssue(
                    code="missing-example-input",
                    path=input_name,
                    message="installed example-input copy omitted one required example file",
                )
            )
            continue
        copied_path = Path(path_text)
        if not copied_path.exists():
            issues.append(
                SmokeIssue(
                    code="missing-copied-example-input",
                    path=path_text,
                    message="installed example-input copy did not materialize one required file",
                )
            )
        elif copied_path == source_root or source_root in copied_path.parents:
            issues.append(
                SmokeIssue(
                    code="source-tree-example-input",
                    path=path_text,
                    message="installed example-input copy resolved back into the repository source tree",
                )
            )
    return issues


def _run_command(
    command: list[str | os.PathLike[str]],
    *,
    cwd: Path,
    env: dict[str, str],
) -> tuple[int, str, str]:
    completed = run_text(
        command,
        cwd=cwd,
        check=False,
        capture_output=True,
        env=env,
    )
    return completed.returncode, completed.stdout, completed.stderr


def _parse_json_payload(
    name: str,
    stdout: str,
    *,
    output_path: Path,
    issues: list[SmokeIssue],
) -> JsonObject:
    _write_text(output_path, stdout)
    try:
        return _as_dict(json.loads(stdout))
    except json.JSONDecodeError as error:
        issues.append(
            SmokeIssue(
                code="invalid-json-output",
                path=output_path.as_posix(),
                message=f"{name} did not emit valid JSON: {error}",
            )
        )
        return {}


def validate_alignment_payload(payload: JsonObject) -> list[SmokeIssue]:
    """Validate the installed FASTA validation smoke payload."""
    metrics = _as_dict(payload.get("metrics"))
    issues: list[SmokeIssue] = []
    if _as_str(payload.get("status")) != "ok":
        issues.append(
            SmokeIssue(
                code="alignment-validate-status",
                path="alignment-validate.json",
                message="installed FASTA validation did not report status=ok",
            )
        )
    if metrics.get("sequence_count") != 4:
        issues.append(
            SmokeIssue(
                code="alignment-sequence-count",
                path="alignment-validate.json",
                message="installed FASTA validation did not preserve the packaged example sequence count",
            )
        )
    if metrics.get("detected_type") != "dna" or metrics.get("selected_type") != "dna":
        issues.append(
            SmokeIssue(
                code="alignment-sequence-type",
                path="alignment-validate.json",
                message="installed FASTA validation did not preserve the packaged example DNA classification",
            )
        )
    return issues


def validate_tree_package_payload(
    payload: JsonObject, out_dir: Path
) -> list[SmokeIssue]:
    """Validate the installed tree report package smoke payload."""
    metrics = _as_dict(payload.get("metrics"))
    issues: list[SmokeIssue] = []
    if _as_str(payload.get("status")) != "ok":
        issues.append(
            SmokeIssue(
                code="tree-package-status",
                path="tree-package.json",
                message="installed tree report package did not report status=ok",
            )
        )
    if metrics.get("tip_count") != 4:
        issues.append(
            SmokeIssue(
                code="tree-package-tip-count",
                path="tree-package.json",
                message="installed tree report package did not preserve the packaged example tip count",
            )
        )
    output_names = (
        {path.name for path in out_dir.iterdir()} if out_dir.is_dir() else set()
    )
    missing_outputs = [
        filename for filename in _TREE_PACKAGE_OUTPUTS if filename not in output_names
    ]
    for filename in missing_outputs:
        issues.append(
            SmokeIssue(
                code="tree-package-output-missing",
                path=(out_dir / filename).as_posix(),
                message="installed tree report package did not write one required artifact",
            )
        )
    return issues


def validate_pgls_payload(payload: JsonObject) -> list[SmokeIssue]:
    """Validate the installed comparative PGLS smoke payload."""
    metrics = _as_dict(payload.get("metrics"))
    issues: list[SmokeIssue] = []
    if _as_str(payload.get("status")) != "ok":
        issues.append(
            SmokeIssue(
                code="comparative-pgls-status",
                path="comparative-pgls.json",
                message="installed comparative PGLS did not report status=ok",
            )
        )
    if metrics.get("taxon_count") != 75:
        issues.append(
            SmokeIssue(
                code="comparative-pgls-taxon-count",
                path="comparative-pgls.json",
                message="installed comparative PGLS did not preserve the packaged primate taxon count",
            )
        )
    if metrics.get("predictor_count") != 1 or metrics.get("coefficient_count") != 2:
        issues.append(
            SmokeIssue(
                code="comparative-pgls-shape",
                path="comparative-pgls.json",
                message="installed comparative PGLS did not preserve the expected model shape for the packaged primate dataset",
            )
        )
    return issues


def run_artifact_install_smoke(
    *,
    artifact_kind: str,
    artifact_path: Path,
    repo_root: Path,
    artifacts_root: Path,
    build_python: Path,
) -> JsonObject:
    """Run the governed installability smoke checks for one built artifact."""
    smoke_root = artifacts_root / artifact_kind
    if smoke_root.exists():
        shutil.rmtree(smoke_root)
    smoke_root.mkdir(parents=True, exist_ok=True)
    env = _smoke_env()
    venv_root = smoke_root / "venv"
    run_text(
        [build_python, "-m", "venv", venv_root],
        cwd=repo_root,
        check=True,
        env=env,
    )
    venv_python = _venv_executable(venv_root, "python")
    cli_path = _venv_executable(venv_root, "bijux-phylogenetics")
    run_text(
        [venv_python, "-m", "pip", "install", artifact_path],
        cwd=smoke_root,
        check=True,
        env=env,
    )

    issues: list[SmokeIssue] = []
    version_code, version_stdout, version_stderr = _run_command(
        [str(cli_path), "--version"],
        cwd=smoke_root,
        env=env,
    )
    _write_text(smoke_root / "version.txt", version_stdout)
    if version_stderr:
        _write_text(smoke_root / "version.stderr.txt", version_stderr)
    if version_code != 0 or not version_stdout.strip():
        issues.append(
            SmokeIssue(
                code="cli-version-failed",
                path="version.txt",
                message="installed CLI version command did not succeed cleanly",
            )
        )
    help_code, help_stdout, help_stderr = _run_command(
        [str(cli_path), "--help"],
        cwd=smoke_root,
        env=env,
    )
    _write_text(smoke_root / "help.txt", help_stdout)
    if help_stderr:
        _write_text(smoke_root / "help.stderr.txt", help_stderr)
    if help_code != 0 or "usage:" not in help_stdout.lower():
        issues.append(
            SmokeIssue(
                code="cli-help-failed",
                path="help.txt",
                message="installed CLI help did not return one usage surface",
            )
        )

    resource_probe = _probe_installed_resources(venv_python, smoke_root, env)
    _write_json(smoke_root / "resource-probe.json", resource_probe)
    issues.extend(validate_resource_probe(resource_probe, repo_root))
    resource_paths = {
        name: Path(path)
        for name, path in _as_dict(resource_probe.get("resource_paths")).items()
        if isinstance(path, str) and path
    }
    example_input_probe = _probe_copied_example_inputs(
        venv_python,
        cwd=smoke_root,
        env=env,
        destination=smoke_root / "example-inputs",
    )
    _write_json(smoke_root / "example-input-probe.json", example_input_probe)
    issues.extend(validate_example_input_probe(example_input_probe, repo_root))
    copied_inputs = {
        name: Path(path)
        for name, path in _as_dict(example_input_probe.get("copied_paths")).items()
        if isinstance(path, str) and path
    }

    alignment_code, alignment_stdout, alignment_stderr = _run_command(
        [
            str(cli_path),
            "alignment",
            "validate-input",
            str(copied_inputs["alignment"]),
            "--json",
        ],
        cwd=smoke_root,
        env=env,
    )
    if alignment_stderr:
        _write_text(smoke_root / "alignment-validate.stderr.txt", alignment_stderr)
    alignment_payload = _parse_json_payload(
        "alignment validate-input",
        alignment_stdout,
        output_path=smoke_root / "alignment-validate.json",
        issues=issues,
    )
    if alignment_code != 0:
        issues.append(
            SmokeIssue(
                code="alignment-validate-failed",
                path="alignment-validate.json",
                message="installed FASTA validation command exited non-zero",
            )
        )
    issues.extend(validate_alignment_payload(alignment_payload))

    tree_package_dir = smoke_root / "tree-package"
    tree_code, tree_stdout, tree_stderr = _run_command(
        [
            str(cli_path),
            "report",
            "tree-package",
            str(copied_inputs["tree"]),
            "--out-dir",
            str(tree_package_dir),
            "--json",
        ],
        cwd=smoke_root,
        env=env,
    )
    if tree_stderr:
        _write_text(smoke_root / "tree-package.stderr.txt", tree_stderr)
    tree_payload = _parse_json_payload(
        "report tree-package",
        tree_stdout,
        output_path=smoke_root / "tree-package.json",
        issues=issues,
    )
    if tree_code != 0:
        issues.append(
            SmokeIssue(
                code="tree-package-failed",
                path="tree-package.json",
                message="installed tree report package command exited non-zero",
            )
        )
    issues.extend(validate_tree_package_payload(tree_payload, tree_package_dir))

    pgls_code, pgls_stdout, pgls_stderr = _run_command(
        [
            str(cli_path),
            "comparative",
            "pgls",
            str(resource_paths["primate_tree"]),
            str(resource_paths["primate_traits"]),
            "--response",
            "longevity",
            "--predictors",
            "social_group_size",
            "--taxon-column",
            "species",
            "--json",
        ],
        cwd=smoke_root,
        env=env,
    )
    if pgls_stderr:
        _write_text(smoke_root / "comparative-pgls.stderr.txt", pgls_stderr)
    pgls_payload = _parse_json_payload(
        "comparative pgls",
        pgls_stdout,
        output_path=smoke_root / "comparative-pgls.json",
        issues=issues,
    )
    if pgls_code != 0:
        issues.append(
            SmokeIssue(
                code="comparative-pgls-failed",
                path="comparative-pgls.json",
                message="installed comparative PGLS command exited non-zero",
            )
        )
    issues.extend(validate_pgls_payload(pgls_payload))

    report = {
        "schema_version": 1,
        "artifact_kind": artifact_kind,
        "artifact_path": artifact_path.as_posix(),
        "artifacts_root": smoke_root.as_posix(),
        "package_root": _as_str(resource_probe.get("package_root")),
        "resource_check_count": len(_RESOURCE_SENTINELS),
        "copied_example_input_count": len(copied_inputs),
        "command_count": 5,
        "issue_count": len(issues),
        "all_passed": not issues,
        "issues": [asdict(issue) for issue in issues],
        "resource_probe": resource_probe,
        "example_input_probe": example_input_probe,
        "commands": {
            "version": {
                "returncode": version_code,
                "stdout_path": (smoke_root / "version.txt").as_posix(),
            },
            "help": {
                "returncode": help_code,
                "stdout_path": (smoke_root / "help.txt").as_posix(),
            },
            "copy_example_inputs": {
                "output_path": (smoke_root / "example-input-probe.json").as_posix(),
            },
            "alignment_validate": {
                "returncode": alignment_code,
                "output_path": (smoke_root / "alignment-validate.json").as_posix(),
            },
            "tree_package": {
                "returncode": tree_code,
                "output_path": (smoke_root / "tree-package.json").as_posix(),
                "output_dir": tree_package_dir.as_posix(),
            },
            "comparative_pgls": {
                "returncode": pgls_code,
                "output_path": (smoke_root / "comparative-pgls.json").as_posix(),
            },
        },
    }
    _write_json(smoke_root / "install-smoke-report.json", report)
    return report


def run_install_smoke_suite(
    *,
    dist_dir: Path,
    repo_root: Path,
    artifacts_root: Path,
    build_python: Path,
    artifact_kind: str,
) -> JsonObject:
    """Run the governed installability smoke suite for one or both artifacts."""
    reports = [
        run_artifact_install_smoke(
            artifact_kind=kind,
            artifact_path=path,
            repo_root=repo_root,
            artifacts_root=artifacts_root,
            build_python=build_python,
        )
        for kind, path in select_artifact_paths(dist_dir, artifact_kind)
    ]
    suite_report = {
        "schema_version": 1,
        "artifact_kind": artifact_kind,
        "dist_dir": dist_dir.as_posix(),
        "repo_root": repo_root.as_posix(),
        "artifacts_root": artifacts_root.as_posix(),
        "artifact_count": len(reports),
        "issue_count": sum(int(report["issue_count"]) for report in reports),
        "all_passed": all(bool(report["all_passed"]) for report in reports),
        "artifacts": reports,
    }
    _write_json(artifacts_root / "install-smoke-report.json", suite_report)
    return suite_report


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Run installability smoke tests against built runtime artifacts."
    )
    parser.add_argument("--dist-dir", required=True, type=Path)
    parser.add_argument("--repo-root", default=Path.cwd(), type=Path)
    parser.add_argument("--artifacts-root", default=DEFAULT_ARTIFACTS_ROOT, type=Path)
    parser.add_argument(
        "--build-python",
        default=Path(sys.executable).resolve(),
        type=Path,
        help="Python interpreter used to create the clean smoke-test virtual environments.",
    )
    parser.add_argument(
        "--artifact-kind",
        choices=("wheel", "sdist", "both"),
        default="both",
        help="Which built artifact shape to install into a clean virtual environment.",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Emit the suite report as JSON to stdout.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """Run the installability smoke-test entry point."""
    args = parse_args(argv)
    report = run_install_smoke_suite(
        dist_dir=args.dist_dir.resolve(),
        repo_root=args.repo_root.resolve(),
        artifacts_root=args.artifacts_root.resolve(),
        build_python=args.build_python.resolve(),
        artifact_kind=args.artifact_kind,
    )
    if args.json:
        print(json.dumps(report, indent=2, sort_keys=True))
    else:
        print(args.artifacts_root.resolve() / "install-smoke-report.json")
    return 0 if report["all_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
