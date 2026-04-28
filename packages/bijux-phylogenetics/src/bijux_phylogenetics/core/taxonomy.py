from __future__ import annotations

from dataclasses import dataclass


@dataclass(slots=True)
class TaxonomyReference:
    """Placeholder taxonomy reference used by future adapters."""

    authority: str
    identifier: str

