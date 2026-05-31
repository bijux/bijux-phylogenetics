from __future__ import annotations

from typing import TypeVar

from bijux_phylogenetics.runtime.errors import PhylogeneticsError

_T = TypeVar("_T")


def require_present(
    value: _T | None,
    *,
    owner_name: str,
    field_name: str,
) -> _T:
    """Return one required Bayesian value or raise one stable domain error."""
    if value is None:
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be present",
            code="bayesian_required_field_missing",
            details={"owner_name": owner_name, "field_name": field_name},
        )
    return value
