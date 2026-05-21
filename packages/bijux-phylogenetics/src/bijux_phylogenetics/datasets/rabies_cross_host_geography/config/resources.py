from __future__ import annotations

from pathlib import Path

from ..models import _DATASET_ID, _WORKFLOW_CONFIG_NAME


def _resource_root() -> Path:
    return (
        Path(__file__).resolve().parent.parent.parent.parent
        / "resources"
        / "datasets"
        / "pathogens"
        / _DATASET_ID
    )


def _default_workflow_config_path() -> Path:
    return _resource_root() / _WORKFLOW_CONFIG_NAME
