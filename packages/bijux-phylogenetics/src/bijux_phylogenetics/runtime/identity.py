from __future__ import annotations

from dataclasses import dataclass

PACKAGE_NAME = "bijux-phylogenetics"
IMPORT_NAME = "bijux_phylogenetics"
PRODUCT_NAME = "Bijux Phylogenetics"
CLI_NAME = "bijux-phylogenetics"
UMBRELLA_COMMAND = "bijux phylogenetics"
CLI_ALIASES = ("bijux phylo",)


@dataclass(frozen=True, slots=True)
class PackageIdentity:
    """Canonical product naming contract for the runtime package."""

    package_name: str = PACKAGE_NAME
    import_name: str = IMPORT_NAME
    product_name: str = PRODUCT_NAME
    cli_name: str = CLI_NAME
    umbrella_command: str = UMBRELLA_COMMAND
    cli_aliases: tuple[str, ...] = CLI_ALIASES


IDENTITY = PackageIdentity()
