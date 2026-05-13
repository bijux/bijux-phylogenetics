"""Release support checks for repository maintenance."""

from .version_resolver import resolve_version

__all__ = [
    "artifact_versions",
    "assert_artifacts_match_version",
    "assert_publishable_version",
    "resolve_version",
]


def __getattr__(name: str) -> object:
    """Load publication guard helpers lazily to keep `python -m` warning-free."""
    if name in {
        "artifact_versions",
        "assert_artifacts_match_version",
        "assert_publishable_version",
    }:
        from . import publication_guard

        return getattr(publication_guard, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
