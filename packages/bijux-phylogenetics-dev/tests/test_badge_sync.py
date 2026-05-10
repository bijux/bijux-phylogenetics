from __future__ import annotations

from pathlib import Path
import re

from bijux_phylogenetics_dev.docs.badge_sync import (
    BadgeTarget,
    load_badge_catalog,
    render_badge_block,
    synchronize_badges,
)

GENERATED_BLOCK_RE = re.compile(
    r"<!-- bijux-phylogenetics-badges:generated:start -->.*?<!-- bijux-phylogenetics-badges:generated:end -->",
    re.DOTALL,
)


def test_badge_catalog_exposes_expected_templates() -> None:
    catalog = load_badge_catalog()
    assert set(catalog) == {
        "family-docs-badge",
        "family-ghcr-badge",
        "family-pypi-badge",
        "maintainer-summary",
        "package-summary",
        "repository-summary",
    }


def test_repository_badge_block_uses_the_shared_workflow_badge_syntax() -> None:
    rendered = render_badge_block(
        BadgeTarget(path=Path("README.md"), kind="repository")
    )

    assert (
        "https://github.com/bijux/bijux-phylogenetics/actions/workflows/verify.yml/badge.svg?branch=main"
        in rendered
    )
    assert (
        "https://img.shields.io/badge/release-pypi%20workflow-2563EB?logo=githubactions&logoColor=white"
        in rendered
    )
    assert (
        "https://img.shields.io/badge/release-ghcr%20workflow-2563EB?logo=githubactions&logoColor=white"
        in rendered
    )
    assert (
        "https://img.shields.io/badge/release-github%20workflow-2563EB?logo=githubactions&logoColor=white"
        in rendered
    )
    assert "workflows/release-pypi/badge.svg" not in rendered
    assert rendered.count("https://img.shields.io/pypi/v/") == 2
    assert rendered.count("/pkgs/container/") == 2


def test_package_badge_block_prioritizes_the_public_distribution() -> None:
    rendered = render_badge_block(
        BadgeTarget(
            path=Path("packages/bijux-phylogenetics/README.md"),
            kind="package",
            package_slug="bijux-phylogenetics",
        )
    )

    assert (
        "\n[![bijux-phylogenetics](https://img.shields.io/pypi/v/bijux-phylogenetics"
        in rendered
    )
    assert (
        "\n[![bijux-phylogenetics](https://img.shields.io/badge/bijux--phylogenetics-ghcr"
        in rendered
    )
    assert (
        "\n[![bijux-phylogenetics docs](https://img.shields.io/badge/docs-bijux--phylogenetics"
        in rendered
    )


def test_alias_package_badge_block_prioritizes_the_alias_distribution() -> None:
    rendered = render_badge_block(
        BadgeTarget(
            path=Path("packages/phylogenetic/README.md"),
            kind="package",
            package_slug="phylogenetic",
        )
    )

    assert "\n[![phylogenetic](https://img.shields.io/pypi/v/phylogenetic" in rendered
    assert (
        "\n[![phylogenetic](https://img.shields.io/badge/phylogenetic-ghcr" in rendered
    )
    assert (
        "\n[![phylogenetic docs](https://img.shields.io/badge/docs-phylogenetic"
        in rendered
    )


def test_badge_surfaces_are_synchronized() -> None:
    assert synchronize_badges(check=True) == []


def test_managed_surfaces_only_use_generated_badges() -> None:
    targets = [
        Path("README.md"),
        Path("docs/index.md"),
        Path("packages/bijux-phylogenetics/README.md"),
        Path("packages/phylogenetic/README.md"),
        Path("packages/bijux-phylogenetics-dev/README.md"),
    ]
    for path in targets:
        text = path.read_text(encoding="utf-8")
        stripped = GENERATED_BLOCK_RE.sub("", text)
        assert "[![" not in stripped, (
            f"{path} contains inline badges outside the generated block"
        )
