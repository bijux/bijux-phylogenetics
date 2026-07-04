"""Synchronize shared badge blocks into managed README surfaces."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import re
import tomllib
from typing import Any, cast

REPO_ROOT = Path(__file__).resolve().parents[5]
BADGE_SOURCE_PATH = REPO_ROOT / "docs" / "badges.md"
START_MARKER = "<!-- bijux-phylogenetics-badges:generated:start -->"
END_MARKER = "<!-- bijux-phylogenetics-badges:generated:end -->"
BADGE_BLOCK_RE = re.compile(
    r"<!-- bijux-phylogenetics-badges:(?P<name>[a-z0-9-]+):start -->\n"
    r"(?P<body>.*?)\n"
    r"<!-- bijux-phylogenetics-badges:(?P=name):end -->",
    re.DOTALL,
)
TOKEN_RE = re.compile(r"{{\s*(?P<name>[a-z0-9_]+)\s*}}")
BADGE_GROUPS: tuple[str, ...] = (
    "family-pypi-badge",
    "family-ghcr-badge",
    "family-docs-badge",
)


@dataclass(frozen=True)
class PackageBadgeRecord:
    """Rendered badge metadata for one public package."""

    package_slug: str
    distribution_name: str
    docs_url: str
    package_pypi_url: str
    package_ghcr_url: str
    pypi_badge_label: str
    docs_badge_label: str
    docs_badge_alt: str


@dataclass(frozen=True)
class BadgeTarget:
    """A README-like surface that consumes generated badge content."""

    path: Path
    kind: str
    package_slug: str | None = None


def _workspace_metadata() -> dict[str, Any]:
    """Load repository workspace metadata from the root `pyproject.toml`."""
    with (REPO_ROOT / "pyproject.toml").open("rb") as handle:
        data = tomllib.load(handle)
    return cast(dict[str, Any], data["tool"]["bijux_phylogenetics"])


def _package_project(package_slug: str) -> dict[str, Any]:
    """Load one package's project table from its `pyproject.toml`."""
    pyproject_path = REPO_ROOT / "packages" / package_slug / "pyproject.toml"
    with pyproject_path.open("rb") as handle:
        data = tomllib.load(handle)
    return cast(dict[str, Any], data["project"])


def _public_package_slugs() -> tuple[str, ...]:
    """Return publishable package slugs, excluding the maintainer package."""
    workspace = _workspace_metadata()
    docs_package = cast(str, workspace["docs_package"])
    packages = cast(list[str], workspace["packages"])
    return tuple(
        package_slug for package_slug in packages if package_slug != docs_package
    )


def _shield_text(value: str) -> str:
    """Escape a label so it renders correctly in Shields badge URLs."""
    return value.replace("-", "--").replace(" ", "%20")


def _docs_badge_alt(distribution_name: str) -> str:
    """Build the docs badge alt text for one distribution."""
    return f"{distribution_name} docs"


def _short_family_label(distribution_name: str) -> str:
    """Shorten the family prefix so badge labels stay compact."""
    if distribution_name.startswith("bijux-phylogenetics-"):
        return distribution_name.removeprefix("bijux-phylogenetics-")
    return distribution_name


def _package_record(package_slug: str) -> PackageBadgeRecord:
    """Assemble the badge metadata for one public package."""
    project = _package_project(package_slug)
    distribution_name = str(project["name"])
    docs_url = str(project.get("urls", {}).get("Documentation", ""))
    label_text = _short_family_label(distribution_name)
    return PackageBadgeRecord(
        package_slug=package_slug,
        distribution_name=distribution_name,
        docs_url=docs_url,
        package_pypi_url=f"https://pypi.org/project/{distribution_name}/",
        package_ghcr_url=(
            "https://github.com/bijux/bijux-phylogenetics/pkgs/container/"
            f"bijux-phylogenetics%2F{distribution_name}"
        ),
        pypi_badge_label=_shield_text(label_text),
        docs_badge_label=_shield_text(label_text),
        docs_badge_alt=_docs_badge_alt(distribution_name),
    )


def public_package_records() -> tuple[PackageBadgeRecord, ...]:
    """Return badge metadata for every public package in workspace order."""
    return tuple(
        _package_record(package_slug) for package_slug in _public_package_slugs()
    )


def load_badge_catalog() -> dict[str, str]:
    """Load named badge templates from the repository badge catalog."""
    text = BADGE_SOURCE_PATH.read_text(encoding="utf-8")
    catalog = {
        match.group("name"): match.group("body").strip()
        for match in BADGE_BLOCK_RE.finditer(text)
    }
    if not catalog:
        raise ValueError(f"No badge blocks found in {BADGE_SOURCE_PATH}")
    return catalog


def _render_template(template: str, context: dict[str, str]) -> str:
    """Render one badge template with the provided placeholder values."""

    def replace(match: re.Match[str]) -> str:
        key = match.group("name")
        try:
            return context[key]
        except KeyError as exc:
            raise KeyError(f"Missing badge token {key!r}") from exc

    return TOKEN_RE.sub(replace, template)


def _record_context(record: PackageBadgeRecord) -> dict[str, str]:
    """Convert a badge record into the template context mapping."""
    return {
        "distribution_name": record.distribution_name,
        "docs_badge_alt": record.docs_badge_alt,
        "docs_badge_label": record.docs_badge_label,
        "docs_url": record.docs_url,
        "package_ghcr_url": record.package_ghcr_url,
        "package_pypi_url": record.package_pypi_url,
        "pypi_badge_label": record.pypi_badge_label,
    }


def _render_badge_group(template: str, records: tuple[PackageBadgeRecord, ...]) -> str:
    """Render one badge family line for every record in order."""
    badges = [_render_template(template, _record_context(record)) for record in records]
    if not badges:
        return ""
    return "\n".join(badges)


def _prioritize_record(
    records: tuple[PackageBadgeRecord, ...], current: PackageBadgeRecord
) -> tuple[PackageBadgeRecord, ...]:
    """Move the current package to the front of the rendered badge order."""
    return (current,) + tuple(
        record for record in records if record.package_slug != current.package_slug
    )


def _render_badge_groups(
    catalog: dict[str, str], records: tuple[PackageBadgeRecord, ...]
) -> list[str]:
    """Render the configured badge families in their canonical order."""
    return [
        _render_badge_group(catalog[template_name], records)
        for template_name in BADGE_GROUPS
    ]


def _render_repository_badges(
    catalog: dict[str, str], records: tuple[PackageBadgeRecord, ...]
) -> str:
    """Render the repository-level badge block."""
    sections = [
        _render_template(
            catalog["repository-summary"],
            {"public_package_count": str(len(records))},
        ),
        *_render_badge_groups(catalog, records),
    ]
    return "\n\n".join(section for section in sections if section)


def _render_package_badges(
    catalog: dict[str, str],
    record: PackageBadgeRecord,
    records: tuple[PackageBadgeRecord, ...],
) -> str:
    """Render the package-level badge block with the current package first."""
    ordered_records = _prioritize_record(records, record)
    sections = [
        _render_template(catalog["package-summary"], _record_context(record)),
        *_render_badge_groups(catalog, ordered_records),
    ]
    return "\n\n".join(section for section in sections if section)


def _render_maintainer_badges(catalog: dict[str, str]) -> str:
    """Render the maintainer package badge block."""
    return catalog["maintainer-summary"]


def iter_badge_targets() -> tuple[BadgeTarget, ...]:
    """Enumerate documentation surfaces that consume managed badge blocks."""
    workspace = _workspace_metadata()
    docs_package = cast(str, workspace["docs_package"])
    package_targets = tuple(
        BadgeTarget(
            path=REPO_ROOT / "packages" / package_slug / "README.md",
            kind="package",
            package_slug=package_slug,
        )
        for package_slug in _public_package_slugs()
    )
    return (
        BadgeTarget(path=REPO_ROOT / "README.md", kind="repository"),
        BadgeTarget(path=REPO_ROOT / "docs" / "index.md", kind="repository"),
        *package_targets,
        BadgeTarget(
            path=REPO_ROOT / "packages" / docs_package / "README.md",
            kind="maintainer",
        ),
    )


def render_badge_block(target: BadgeTarget) -> str:
    """Render the managed badge block for one documentation surface."""
    catalog = load_badge_catalog()
    records = public_package_records()
    if target.kind == "repository":
        return _render_repository_badges(catalog, records)
    if target.kind == "maintainer":
        return _render_maintainer_badges(catalog)
    if target.kind == "package" and target.package_slug is not None:
        record = next(
            record for record in records if record.package_slug == target.package_slug
        )
        return _render_package_badges(catalog, record, records)
    raise ValueError(f"Unsupported badge target: {target}")


def _managed_block(rendered_badges: str) -> str:
    """Wrap rendered badges in the managed-section markers."""
    return f"{START_MARKER}\n{rendered_badges.rstrip()}\n{END_MARKER}"


def render_target_text(target: BadgeTarget, current_text: str) -> str:
    """Replace a surface badge section with the generated managed block."""
    managed_block = _managed_block(render_badge_block(target))
    if START_MARKER in current_text and END_MARKER in current_text:
        start = current_text.index(START_MARKER)
        end = current_text.index(END_MARKER) + len(END_MARKER)
        updated = current_text[:start] + managed_block + current_text[end:]
        return updated if updated.endswith("\n") else updated + "\n"

    lines = current_text.splitlines()
    start_index = next(
        (index for index, line in enumerate(lines) if line.startswith("[![")),
        None,
    )
    if start_index is None:
        raise ValueError(f"Unable to locate badge block in {target.path}")

    end_index = start_index
    while end_index < len(lines):
        line = lines[end_index]
        if line.startswith("[![") or not line.strip():
            end_index += 1
            continue
        break

    replacement = managed_block.splitlines()
    replacement.append("")
    updated_lines = lines[:start_index] + replacement + lines[end_index:]
    return "\n".join(updated_lines) + "\n"


def synchronize_badges(*, check: bool) -> list[Path]:
    """Sync managed badge blocks and return the surfaces that changed."""
    changed: list[Path] = []
    for target in iter_badge_targets():
        current_text = target.path.read_text(encoding="utf-8")
        expected_text = render_target_text(target, current_text)
        if expected_text == current_text:
            continue
        changed.append(target.path)
        if not check:
            target.path.write_text(expected_text, encoding="utf-8")
    return changed


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser for badge synchronization."""
    parser = argparse.ArgumentParser(
        description="Synchronize shared badge templates into managed README surfaces."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser(
        "sync", help="Render badge blocks into managed README surfaces."
    )
    subparsers.add_parser(
        "check",
        help="Fail when managed README surfaces differ from the generated badge contract.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the badge synchronization entry point."""
    parser = build_parser()
    args = parser.parse_args(argv)
    check = args.command == "check"
    changed = synchronize_badges(check=check)
    if check and changed:
        for path in changed:
            print(path.relative_to(REPO_ROOT))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
