"""Repository-owned governance policy resources for bijux-phylogenetics."""

from __future__ import annotations

from pathlib import Path

POLICY_ROOT = Path(
    "packages/bijux-phylogenetics-dev/src/bijux_phylogenetics_dev/quality/policies"
)
CONFIG_SSOT_POLICY_PATH = POLICY_ROOT / "config_ssot.toml"
EXECUTION_SURFACES_POLICY_PATH = POLICY_ROOT / "execution_surfaces.toml"
PACKAGE_BOUNDARIES_POLICY_PATH = POLICY_ROOT / "package_boundaries.toml"
PUBLICATION_READINESS_POLICY_PATH = POLICY_ROOT / "publication_readiness.toml"
