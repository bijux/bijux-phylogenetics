"""Publication guard helpers."""

from __future__ import annotations

import argparse
from pathlib import Path

from bijux_phylogenetics_dev.quality.config_ssot import check_config_ssot
from bijux_phylogenetics_dev.quality.package_boundaries import check_package_boundaries
from bijux_phylogenetics_dev.quality.package_bundles import check_package_bundles
from bijux_phylogenetics_dev.quality.publish_readiness import check_publish_readiness

from .version_resolver import resolve_version


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate that a package version is safe to publish."
    )
    parser.add_argument("--pyproject", required=True, help="Path to pyproject.toml")
    parser.add_argument("--package-name", required=True, help="Package slug/tag prefix")
    parser.add_argument(
        "--dist-dir",
        help="Optional dist directory whose artifact versions must match the resolved version",
    )
    parser.add_argument(
        "--allow-prerelease",
        action="store_true",
        help="Allow prerelease versions such as .devN or rcN",
    )
    parser.add_argument(
        "--allow-local-version",
        action="store_true",
        help="Allow local version segments such as +dirty",
    )
    parser.add_argument(
        "--repo-root",
        default="",
        help="Optional repository root for repo-level publish guards.",
    )
    parser.add_argument(
        "--require-config-ssot",
        action="store_true",
        help="Require a clean repository config SSOT audit before publish checks pass.",
    )
    parser.add_argument(
        "--require-package-bundles",
        action="store_true",
        help="Require a clean publish-artifact bundle audit before publish checks pass.",
    )
    parser.add_argument(
        "--require-package-boundaries",
        action="store_true",
        help="Require a clean package boundary audit before publish checks pass.",
    )
    parser.add_argument(
        "--require-publish-readiness",
        action="store_true",
        help="Require a clean repository publish-readiness report before publish checks pass.",
    )
    return parser.parse_args()


def _artifact_version(path: Path) -> str:
    """Parse the version embedded in an sdist or wheel filename."""
    if path.name.endswith(".whl"):
        parts = path.name[:-4].split("-")
        if len(parts) < 2:
            raise ValueError(f"unrecognized wheel filename: {path.name}")
        return parts[1]
    if path.name.endswith(".tar.gz"):
        stem = path.name[: -len(".tar.gz")]
        if "-" not in stem:
            raise ValueError(f"unrecognized sdist filename: {path.name}")
        return stem.rsplit("-", 1)[1]
    raise ValueError(f"unsupported artifact extension: {path.name}")


def artifact_versions(dist_dir: Path) -> dict[str, str]:
    """Collect resolved artifact versions from a dist directory."""
    versions: dict[str, str] = {}
    for path in sorted(dist_dir.glob("*.whl")) + sorted(dist_dir.glob("*.tar.gz")):
        versions[path.name] = _artifact_version(path)
    return versions


def assert_publishable_version(
    version: str,
    *,
    allow_prerelease: bool = False,
    allow_local_version: bool = False,
) -> None:
    """Reject prerelease and local-only versions unless explicitly allowed."""
    lowered = version.lower()
    prerelease_markers = (".dev", "a", "b", "rc")
    if not allow_prerelease and any(marker in lowered for marker in prerelease_markers):
        raise ValueError(
            f"{version} is a prerelease version; create the release tag or set "
            "PUBLISH_ALLOW_PRERELEASE=1 for an intentional prerelease publish"
        )
    if not allow_local_version and ("+" in version or "dirty" in lowered):
        raise ValueError(
            f"{version} includes a local build marker; clean the checkout or set "
            "PUBLISH_ALLOW_LOCAL_VERSION=1 for an intentional local publish"
        )


def assert_artifacts_match_version(dist_dir: Path, version: str) -> None:
    """Ensure all publish artifacts align with the resolved version."""
    versions = artifact_versions(dist_dir)
    if not versions:
        raise ValueError(f"no artifacts found under {dist_dir}")
    mismatched = {
        name: artifact_version
        for name, artifact_version in versions.items()
        if artifact_version != version
    }
    if mismatched:
        details = ", ".join(
            f"{name}={artifact_version}"
            for name, artifact_version in sorted(mismatched.items())
        )
        raise ValueError(
            f"artifact versions do not match resolved version {version}: {details}"
        )


def assert_publishable_repository(
    *,
    repo_root: Path | None = None,
    require_config_ssot: bool = False,
    require_package_boundaries: bool = False,
    require_package_bundles: bool = False,
    require_publish_readiness: bool = False,
) -> None:
    """Reject repository states that are not safe to publish."""
    if (
        not require_config_ssot
        and not require_package_boundaries
        and not require_package_bundles
        and not require_publish_readiness
    ):
        return
    if repo_root is None:
        raise ValueError(
            "repo_root is required when repository-level publish guards are enabled"
        )
    if require_config_ssot:
        check_config_ssot(repo_root)
    if require_package_boundaries:
        check_package_boundaries(repo_root)
    if require_package_bundles:
        check_package_bundles(repo_root)
    if require_publish_readiness:
        check_publish_readiness(repo_root)


def main() -> int:
    """Run the command-line entry point."""
    args = parse_args()
    version = resolve_version(Path(args.pyproject), args.package_name)
    if version == "0.0.0":
        raise SystemExit("unable to resolve package version")
    assert_publishable_version(
        version,
        allow_prerelease=args.allow_prerelease,
        allow_local_version=args.allow_local_version,
    )
    repo_root = Path(args.repo_root).resolve() if args.repo_root else None
    assert_publishable_repository(
        repo_root=repo_root,
        require_config_ssot=args.require_config_ssot,
        require_package_boundaries=args.require_package_boundaries,
        require_package_bundles=args.require_package_bundles,
        require_publish_readiness=args.require_publish_readiness,
    )
    if args.dist_dir:
        assert_artifacts_match_version(Path(args.dist_dir), version)
    print(version)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
