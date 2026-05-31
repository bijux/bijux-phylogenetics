"""Audit repository package ownership, exports, and cross-package boundaries."""

from __future__ import annotations

import argparse
import ast
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import tomllib
from typing import Any

from .policies import PACKAGE_BOUNDARIES_POLICY_PATH

DEFAULT_JSON_OUT = Path("artifacts/root/package-boundaries.json")


@dataclass(frozen=True)
class BoundaryIssue:
    """Describe one package-boundary contract drift."""

    code: str
    path: str
    message: str


@dataclass(frozen=True)
class PackageRolePolicy:
    """Define one governed repository package role."""

    package_name: str
    role: str
    package_dir: str
    module_root: str
    allowed_repo_import_roots: tuple[str, ...]
    owned_module_prefixes: tuple[str, ...]
    required_install_dependencies: tuple[str, ...]


@dataclass(frozen=True)
class TargetPackageRole:
    """Define one governed external-facing target package role."""

    package_name: str
    role: str
    target_module_root: str
    required_runtime_dependency: str


@dataclass(frozen=True)
class RuntimeEvidenceCompatibilityContract:
    """Define the minimal runtime API promised to evidence consumers."""

    runtime_version_spec: str
    supported_api_modules: tuple[str, ...]
    supported_api_locators: tuple[str, ...]
    notes: str


@dataclass(frozen=True)
class PackageBoundaryPolicy:
    """Collect the repository package-boundary policy surfaces."""

    known_repo_module_roots: tuple[str, ...]
    forbidden_runtime_top_level_exports: tuple[str, ...]
    alias_allowed_local_files: tuple[str, ...]
    package_roles: dict[str, PackageRolePolicy]
    target_package_roles: dict[str, TargetPackageRole]
    runtime_evidence_compatibility: RuntimeEvidenceCompatibilityContract


def _load_toml(path: Path) -> dict[str, Any]:
    with path.open("rb") as handle:
        return tomllib.load(handle)


def _as_dict(value: object) -> dict[str, Any]:
    return value if isinstance(value, dict) else {}


def _as_str_tuple(value: object) -> tuple[str, ...]:
    if not isinstance(value, list):
        return ()
    return tuple(entry for entry in value if isinstance(entry, str))


def load_package_boundary_policy(repo_root: Path) -> PackageBoundaryPolicy:
    """Load the repository package-boundary policy from the governed TOML file."""
    payload = _load_toml(repo_root / PACKAGE_BOUNDARIES_POLICY_PATH)
    tool = _as_dict(payload.get("tool"))
    workspace = _as_dict(tool.get("bijux_phylogenetics"))
    section = _as_dict(workspace.get("package_boundaries"))
    package_roles = {
        package_name: PackageRolePolicy(
            package_name=package_name,
            role=str(values["role"]),
            package_dir=str(values["package_dir"]),
            module_root=str(values["module_root"]),
            allowed_repo_import_roots=_as_str_tuple(
                values.get("allowed_repo_import_roots")
            ),
            owned_module_prefixes=_as_str_tuple(values.get("owned_module_prefixes")),
            required_install_dependencies=_as_str_tuple(
                values.get("required_install_dependencies")
            ),
        )
        for package_name, values in _as_dict(section.get("package_roles")).items()
        if isinstance(values, dict)
    }
    target_package_roles = {
        package_name: TargetPackageRole(
            package_name=package_name,
            role=str(values["role"]),
            target_module_root=str(values["target_module_root"]),
            required_runtime_dependency=str(values["required_runtime_dependency"]),
        )
        for package_name, values in _as_dict(
            section.get("target_package_roles")
        ).items()
        if isinstance(values, dict)
    }
    compatibility = _as_dict(section.get("runtime_evidence_compatibility"))
    return PackageBoundaryPolicy(
        known_repo_module_roots=_as_str_tuple(section.get("known_repo_module_roots")),
        forbidden_runtime_top_level_exports=_as_str_tuple(
            section.get("forbidden_runtime_top_level_exports")
        ),
        alias_allowed_local_files=_as_str_tuple(
            section.get("alias_allowed_local_files")
        ),
        package_roles=package_roles,
        target_package_roles=target_package_roles,
        runtime_evidence_compatibility=RuntimeEvidenceCompatibilityContract(
            runtime_version_spec=str(compatibility["runtime_version_spec"]),
            supported_api_modules=_as_str_tuple(
                compatibility.get("supported_api_modules")
            ),
            supported_api_locators=_as_str_tuple(
                compatibility.get("supported_api_locators")
            ),
            notes=str(compatibility.get("notes", "")),
        ),
    )


def _runtime_init_exports(path: Path) -> list[str]:
    module = ast.parse(path.read_text(encoding="utf-8"))
    for node in module.body:
        if not isinstance(node, ast.Assign):
            continue
        if not any(
            isinstance(target, ast.Name) and target.id == "__all__"
            for target in node.targets
        ):
            continue
        if not isinstance(node.value, ast.List):
            continue
        return [
            element.value
            for element in node.value.elts
            if isinstance(element, ast.Constant) and isinstance(element.value, str)
        ]
    return []


def _package_dependencies(repo_root: Path, package_dir: str) -> tuple[str, ...]:
    pyproject = _load_toml(repo_root / package_dir / "pyproject.toml")
    project = _as_dict(pyproject.get("project"))
    dependencies = project.get("dependencies", [])
    return tuple(entry for entry in dependencies if isinstance(entry, str))


def _iter_python_files(root: Path) -> list[Path]:
    return sorted(path for path in root.rglob("*.py") if path.is_file())


def _is_generated_python_cache(path: Path) -> bool:
    return "__pycache__" in path.parts or path.suffix in {".pyc", ".pyo"}


def _repo_import_roots(path: Path, known_roots: tuple[str, ...]) -> list[str]:
    module = ast.parse(path.read_text(encoding="utf-8"))
    imports: list[str] = []
    for node in ast.walk(module):
        if isinstance(node, ast.Import):
            for alias in node.names:
                root = alias.name.split(".", 1)[0]
                if root in known_roots:
                    imports.append(root)
        elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            root = node.module.split(".", 1)[0]
            if root in known_roots:
                imports.append(root)
    return sorted(set(imports))


def _string_literal_values(node: ast.AST) -> set[str] | None:
    if not isinstance(node, (ast.List, ast.Tuple)):
        return None
    return {
        element.value
        for element in node.elts
        if isinstance(element, ast.Constant) and isinstance(element.value, str)
    }


def _relative_import_module_path(
    module_path: Path,
    *,
    module: str | None,
    level: int,
) -> Path | None:
    anchor = module_path.parent
    for _ in range(max(level - 1, 0)):
        anchor = anchor.parent
    relative_parts = module.split(".") if module else []
    candidate_root = anchor.joinpath(*relative_parts)
    package_candidate = candidate_root / "__init__.py"
    if package_candidate.is_file():
        return package_candidate
    module_candidate = candidate_root.with_suffix(".py")
    return module_candidate if module_candidate.is_file() else None


def _resolve_imported_symbol(
    module_path: Path,
    module: ast.Module,
    symbol_name: str,
) -> tuple[Path, str] | None:
    for node in module.body:
        if not isinstance(node, ast.ImportFrom):
            continue
        for alias in node.names:
            local_name = alias.asname or alias.name
            if local_name != symbol_name:
                continue
            imported_module_path = _relative_import_module_path(
                module_path,
                module=node.module or alias.name,
                level=node.level,
            )
            if imported_module_path is None:
                return None
            return imported_module_path, alias.name
    return None


def _resolve_export_module_paths(
    module_path: Path,
    module: ast.Module,
    symbol_name: str,
) -> list[Path] | None:
    for node in module.body:
        if not isinstance(node, ast.Assign):
            continue
        if not any(
            isinstance(target, ast.Name) and target.id == symbol_name
            for target in node.targets
        ):
            continue
        if not isinstance(node.value, (ast.List, ast.Tuple)):
            return None
        module_paths: list[Path] = []
        for element in node.value.elts:
            if not isinstance(element, ast.Name):
                return None
            imported_symbol = _resolve_imported_symbol(module_path, module, element.id)
            if imported_symbol is None:
                return None
            imported_module_path, _ = imported_symbol
            module_paths.append(imported_module_path)
        return module_paths
    return None


def _resolve_symbol_exports(
    module_path: Path,
    symbol_name: str,
    *,
    seen: set[tuple[Path, str]],
) -> set[str] | None:
    key = (module_path, symbol_name)
    if key in seen:
        return None
    seen.add(key)
    module = ast.parse(module_path.read_text(encoding="utf-8"))
    for node in module.body:
        if not isinstance(node, ast.Assign):
            continue
        if not any(
            isinstance(target, ast.Name) and target.id == symbol_name
            for target in node.targets
        ):
            continue
        literal_values = _string_literal_values(node.value)
        if literal_values is not None:
            return literal_values
        if (
            isinstance(node.value, ast.Call)
            and isinstance(node.value.func, ast.Name)
            and node.value.func.id in {"list", "tuple"}
            and len(node.value.args) == 1
        ):
            argument = node.value.args[0]
            if isinstance(argument, ast.Name):
                imported_symbol = _resolve_imported_symbol(
                    module_path, module, argument.id
                )
                if imported_symbol is None:
                    return _resolve_symbol_exports(module_path, argument.id, seen=seen)
                imported_module_path, imported_name = imported_symbol
                return _resolve_symbol_exports(
                    imported_module_path, imported_name, seen=seen
                )
            if (
                isinstance(argument, ast.GeneratorExp)
                and len(argument.generators) == 2
                and isinstance(argument.generators[0].target, ast.Name)
                and isinstance(argument.generators[0].iter, ast.Name)
                and isinstance(argument.generators[1].iter, ast.Attribute)
                and isinstance(argument.generators[1].iter.value, ast.Name)
                and argument.generators[1].iter.attr == "__all__"
                and argument.generators[0].target.id
                == argument.generators[1].iter.value.id
            ):
                export_module_paths = _resolve_export_module_paths(
                    module_path, module, argument.generators[0].iter.id
                )
                if export_module_paths is None:
                    return None
                exports: set[str] = set()
                for export_module_path in export_module_paths:
                    exports.update(_module_exports(export_module_path))
                return exports
    imported_symbol = _resolve_imported_symbol(module_path, module, symbol_name)
    if imported_symbol is None:
        return None
    imported_module_path, imported_name = imported_symbol
    return _resolve_symbol_exports(imported_module_path, imported_name, seen=seen)


def _module_exports(module_path: Path, *, module_name: str | None = None) -> set[str]:
    del module_name
    resolved_exports = _resolve_symbol_exports(module_path, "__all__", seen=set())
    if resolved_exports is not None:
        return resolved_exports
    exports: set[str] = set()
    module = ast.parse(module_path.read_text(encoding="utf-8"))
    for node in module.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            exports.add(node.name)
    return exports


def _tuple_constant_strings(module_path: Path, constant_name: str) -> tuple[str, ...]:
    module = ast.parse(module_path.read_text(encoding="utf-8"))
    for node in module.body:
        if not isinstance(node, ast.Assign):
            continue
        if not any(
            isinstance(target, ast.Name) and target.id == constant_name
            for target in node.targets
        ):
            continue
        if not isinstance(node.value, (ast.Tuple, ast.List)):
            continue
        values: list[str] = []
        for element in node.value.elts:
            if isinstance(element, ast.Constant) and isinstance(element.value, str):
                values.append(element.value)
        return tuple(values)
    return ()


def _module_source_path(
    repo_root: Path,
    policy: PackageRolePolicy,
    module_name: str,
) -> Path | None:
    prefix = policy.module_root.split(".")
    parts = module_name.split(".")
    if parts[: len(prefix)] != prefix:
        return None
    relative = parts[len(prefix) :]
    module_root = repo_root / policy.package_dir / "src" / Path(*prefix)
    if not relative:
        candidate = module_root / "__init__.py"
        return candidate if candidate.is_file() else None
    package_candidate = module_root / Path(*relative) / "__init__.py"
    if package_candidate.is_file():
        return package_candidate
    module_candidate = module_root / Path(*relative[:-1]) / f"{relative[-1]}.py"
    return module_candidate if module_candidate.is_file() else None


def build_package_boundary_report(repo_root: Path) -> dict[str, Any]:
    """Build the package-boundary ownership and export report."""
    repo_root = repo_root.resolve()
    policy = load_package_boundary_policy(repo_root)
    issues: list[BoundaryIssue] = []

    runtime_policy = next(
        role_policy
        for role_policy in policy.package_roles.values()
        if role_policy.role == "runtime"
    )
    alias_policy = next(
        role_policy
        for role_policy in policy.package_roles.values()
        if role_policy.role == "compatibility-alias"
    )
    runtime_init = (
        repo_root
        / runtime_policy.package_dir
        / "src"
        / runtime_policy.module_root
        / "__init__.py"
    )
    runtime_exports = sorted(_runtime_init_exports(runtime_init))
    forbidden_exports_present = sorted(
        export
        for export in policy.forbidden_runtime_top_level_exports
        if export in runtime_exports
    )
    for export in forbidden_exports_present:
        issues.append(
            BoundaryIssue(
                code="forbidden-runtime-top-level-export",
                path=runtime_init.relative_to(repo_root).as_posix(),
                message=f"runtime top-level export {export} leaks evidence-only helper surface",
            )
        )
    evidence_contract_module = (
        repo_root
        / runtime_policy.package_dir
        / "src"
        / runtime_policy.module_root
        / "comparative"
        / "evidence_contract.py"
    )
    contract_modules: tuple[str, ...] = ()
    contract_locators: tuple[str, ...] = ()
    if not evidence_contract_module.is_file():
        issues.append(
            BoundaryIssue(
                code="missing-runtime-evidence-contract-module",
                path=evidence_contract_module.relative_to(repo_root).as_posix(),
                message="runtime package must expose a governed evidence compatibility contract module",
            )
        )
    else:
        contract_modules = _tuple_constant_strings(
            evidence_contract_module, "SUPPORTED_EVIDENCE_API_MODULES"
        )
        contract_locators = _tuple_constant_strings(
            evidence_contract_module, "SUPPORTED_EVIDENCE_API_LOCATORS"
        )
        missing_contract_modules = tuple(
            module_name
            for module_name in policy.runtime_evidence_compatibility.supported_api_modules
            if module_name not in contract_modules
        )
        if missing_contract_modules:
            issues.append(
                BoundaryIssue(
                    code="runtime-evidence-module-contract-drift",
                    path=evidence_contract_module.relative_to(repo_root).as_posix(),
                    message="runtime evidence API module contract is missing one or more governed compatibility modules",
                )
            )
        missing_contract_locators = tuple(
            locator
            for locator in policy.runtime_evidence_compatibility.supported_api_locators
            if locator not in contract_locators
        )
        if missing_contract_locators:
            issues.append(
                BoundaryIssue(
                    code="runtime-evidence-locator-contract-drift",
                    path=evidence_contract_module.relative_to(repo_root).as_posix(),
                    message="runtime evidence API locator contract is missing one or more governed compatibility locators",
                )
            )

    alias_root = repo_root / alias_policy.package_dir / "src" / alias_policy.module_root
    actual_alias_files = sorted(
        path.relative_to(alias_root).as_posix()
        for path in alias_root.rglob("*")
        if path.is_file() and not _is_generated_python_cache(path)
    )
    allowed_alias_files = sorted(policy.alias_allowed_local_files)
    unexpected_alias_files = sorted(
        set(actual_alias_files).difference(allowed_alias_files)
    )
    missing_alias_files = sorted(
        set(allowed_alias_files).difference(actual_alias_files)
    )
    if unexpected_alias_files or missing_alias_files:
        issues.append(
            BoundaryIssue(
                code="alias-local-surface-drift",
                path=alias_root.relative_to(repo_root).as_posix(),
                message=(
                    "compatibility alias package local surface drifted: "
                    f"missing={missing_alias_files or ['<none>']}, "
                    f"unexpected={unexpected_alias_files or ['<none>']}"
                ),
            )
        )

    package_reports: list[dict[str, Any]] = []
    for package_name, role_policy in sorted(policy.package_roles.items()):
        package_dir = repo_root / role_policy.package_dir
        src_root = package_dir / "src" / role_policy.module_root
        dependencies = _package_dependencies(repo_root, role_policy.package_dir)
        dependency_ok = dependencies == role_policy.required_install_dependencies
        if not dependency_ok:
            issues.append(
                BoundaryIssue(
                    code="dependency-role-drift",
                    path=f"{role_policy.package_dir}/pyproject.toml",
                    message=(
                        f"{package_name} dependencies {dependencies!r} do not match "
                        f"the governed role contract {role_policy.required_install_dependencies!r}"
                    ),
                )
            )
        disallowed_imports: list[dict[str, Any]] = []
        for path in _iter_python_files(src_root):
            imported_roots = _repo_import_roots(path, policy.known_repo_module_roots)
            forbidden = sorted(
                root
                for root in imported_roots
                if root not in role_policy.allowed_repo_import_roots
            )
            if forbidden:
                disallowed_imports.append(
                    {
                        "path": path.relative_to(repo_root).as_posix(),
                        "forbidden_repo_import_roots": forbidden,
                    }
                )
                issues.append(
                    BoundaryIssue(
                        code="repo-import-boundary-drift",
                        path=path.relative_to(repo_root).as_posix(),
                        message=(
                            f"{package_name} imports forbidden repository package roots: "
                            + ", ".join(forbidden)
                        ),
                    )
                )
        package_reports.append(
            {
                "package_name": package_name,
                "role": role_policy.role,
                "package_dir": role_policy.package_dir,
                "module_root": role_policy.module_root,
                "owned_module_prefixes": list(role_policy.owned_module_prefixes),
                "allowed_repo_import_roots": list(
                    role_policy.allowed_repo_import_roots
                ),
                "required_install_dependencies": list(
                    role_policy.required_install_dependencies
                ),
                "dependencies": list(dependencies),
                "dependency_ok": dependency_ok,
                "disallowed_imports": disallowed_imports,
            }
        )

    compatibility = policy.runtime_evidence_compatibility
    compatibility_results: list[dict[str, Any]] = []
    for locator in compatibility.supported_api_locators:
        module_name, separator, attribute_name = locator.partition(":")
        valid = True
        issue_code = ""
        message = ""
        if not separator or not attribute_name:
            valid = False
            issue_code = "invalid-compatibility-locator"
            message = "supported locator must use module:attribute syntax"
        elif module_name not in compatibility.supported_api_modules:
            valid = False
            issue_code = "unsupported-compatibility-module"
            message = f"{module_name} is not an allowed supported API module for evidence consumers"
        else:
            module_path = _module_source_path(repo_root, runtime_policy, module_name)
            if module_path is None:
                valid = False
                issue_code = "missing-compatibility-module"
                message = f"cannot resolve supported API module {module_name}"
            else:
                exports = _module_exports(module_path, module_name=module_name)
                if attribute_name not in exports:
                    valid = False
                    issue_code = "missing-compatibility-export"
                    message = f"{attribute_name} is not exported by supported API module {module_name}"
        compatibility_results.append(
            {
                "locator": locator,
                "valid": valid,
                "issue_code": issue_code,
                "message": message,
            }
        )
        if not valid:
            issues.append(
                BoundaryIssue(
                    code=issue_code,
                    path=PACKAGE_BOUNDARIES_POLICY_PATH.as_posix(),
                    message=message,
                )
            )

    target_role_reports = [
        {
            "package_name": target.package_name,
            "role": target.role,
            "target_module_root": target.target_module_root,
            "required_runtime_dependency": target.required_runtime_dependency,
            "package_exists": (repo_root / "packages" / target.package_name).is_dir(),
        }
        for target in policy.target_package_roles.values()
    ]

    return {
        "schema_version": 1,
        "runtime_public_api": {
            "top_level_export_count": len(runtime_exports),
            "forbidden_top_level_exports": list(
                policy.forbidden_runtime_top_level_exports
            ),
            "forbidden_top_level_exports_present": forbidden_exports_present,
        },
        "alias_surface": {
            "allowed_local_files": allowed_alias_files,
            "actual_local_files": actual_alias_files,
            "missing_local_files": missing_alias_files,
            "unexpected_local_files": unexpected_alias_files,
        },
        "package_roles": package_reports,
        "target_package_roles": target_role_reports,
        "runtime_evidence_compatibility": {
            "contract_module_path": evidence_contract_module.relative_to(
                repo_root
            ).as_posix(),
            "contract_modules": list(contract_modules),
            "contract_locators": list(contract_locators),
            "runtime_version_spec": compatibility.runtime_version_spec,
            "supported_api_modules": list(compatibility.supported_api_modules),
            "supported_api_locators": list(compatibility.supported_api_locators),
            "notes": compatibility.notes,
            "locator_results": compatibility_results,
        },
        "issue_count": len(issues),
        "issues": [asdict(issue) for issue in issues],
    }


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def check_package_boundaries(
    repo_root: Path,
    *,
    json_out: Path | None = None,
) -> dict[str, Any]:
    """Raise when the package-boundary report contains governance issues."""
    payload = build_package_boundary_report(repo_root.resolve())
    if json_out is not None:
        _write_json(json_out, payload)
    if payload["issue_count"]:
        raise SystemExit("package boundary audit failed")
    return payload


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the package-boundary audit."""
    parser = argparse.ArgumentParser(
        description="Audit package ownership, exports, and cross-package boundaries."
    )
    parser.add_argument("command", choices=("report", "check"))
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--json-out", default=str(DEFAULT_JSON_OUT))
    return parser.parse_args()


def main() -> int:
    """Run the package-boundary CLI entry point."""
    args = parse_args()
    repo_root = Path(args.repo_root).resolve()
    json_out = Path(args.json_out)
    if not json_out.is_absolute():
        json_out = repo_root / json_out
    if args.command == "check":
        payload = check_package_boundaries(repo_root, json_out=json_out)
    else:
        payload = build_package_boundary_report(repo_root)
        _write_json(json_out, payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
