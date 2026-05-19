from __future__ import annotations

from pathlib import Path

_DATASET_ID = "rabies_method_sensitivity_panel"
_DATASET_LABEL = "Rabies method-sensitivity panel"
_SEQUENCE_TYPE = "dna"
_WORKFLOW_PREFIX = "rabies-method-sensitivity-panel"
_SOURCE_ACCESSIONS = (
    "MG458305",
    "MG458304",
    "PV641713",
    "PX845689",
    "OQ693985",
    "PX845683",
    "PX845681",
    "PX845678",
    "PX845676",
)


def _resource_root() -> Path:
    return (
        Path(__file__).resolve().parent.parent.parent
        / "resources"
        / "datasets"
        / "pathogens"
        / _DATASET_ID
    )
