from __future__ import annotations

import pytest

from bijux_phylogenetics.parity.live_access import (
    RUN_LIVE_PARITY_ENV_VAR,
    is_live_parity_enabled,
    require_live_parity_access,
)
from bijux_phylogenetics.runtime.errors import (
    ExternalEngineAccessError,
    PhylogeneticsError,
)


def test_live_parity_is_disabled_by_default(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.delenv(RUN_LIVE_PARITY_ENV_VAR, raising=False)

    assert is_live_parity_enabled() is False


@pytest.mark.parametrize("raw_value", ["1", "true", "TRUE", " yes ", "on"])
def test_live_parity_truthy_values_enable_access(
    raw_value: str,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv(RUN_LIVE_PARITY_ENV_VAR, raw_value)

    assert is_live_parity_enabled() is True


@pytest.mark.parametrize("raw_value", ["0", "false", "FALSE", " no ", "off", ""])
def test_live_parity_falsey_values_disable_access(
    raw_value: str,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv(RUN_LIVE_PARITY_ENV_VAR, raw_value)

    assert is_live_parity_enabled() is False


def test_live_parity_rejects_invalid_environment_values(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv(RUN_LIVE_PARITY_ENV_VAR, "sometimes")

    with pytest.raises(
        PhylogeneticsError,
        match="RUN_LIVE_PARITY must be one of",
    ) as error:
        is_live_parity_enabled()

    assert error.value.code == "live_parity_opt_in_invalid"
    assert error.value.details["environment_variable"] == RUN_LIVE_PARITY_ENV_VAR
    assert error.value.details["received_value"] == "sometimes"


def test_live_parity_access_reports_opt_in_requirement(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.delenv(RUN_LIVE_PARITY_ENV_VAR, raising=False)

    with pytest.raises(
        ExternalEngineAccessError,
        match="set RUN_LIVE_PARITY=1 to continue",
    ) as error:
        require_live_parity_access(
            surface="ape live parity runner",
            reference_source="ape-live",
            executable="Rscript",
        )

    assert error.value.code == "external_engine_access_forbidden"
    assert error.value.details["surface"] == "ape live parity runner"
    assert error.value.details["reference_source"] == "ape-live"
    assert error.value.details["executable"] == "Rscript"
    assert error.value.details["environment_variable"] == RUN_LIVE_PARITY_ENV_VAR
