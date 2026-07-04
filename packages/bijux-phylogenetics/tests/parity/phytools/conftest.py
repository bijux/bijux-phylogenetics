from __future__ import annotations

import pytest


@pytest.fixture(autouse=True)
def _enable_live_parity_access(enable_live_parity_access: None) -> None:
    """Enable live parity access for phytools parity runner coverage."""
