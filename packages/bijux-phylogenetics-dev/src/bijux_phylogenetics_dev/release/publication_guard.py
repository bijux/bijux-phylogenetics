"""Publication guard helpers."""

from __future__ import annotations

import argparse
from pathlib import Path

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
    if args.dist_dir:
        assert_artifacts_match_version(Path(args.dist_dir), version)
    print(version)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
