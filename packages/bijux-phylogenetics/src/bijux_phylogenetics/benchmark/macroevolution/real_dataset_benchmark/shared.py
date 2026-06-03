from __future__ import annotations

import math

from .contracts import RealDatasetMacroevolutionModelRow

CONTINUOUS_SURFACE_ID = "seed-mass-native-model-table"
DISCRETE_SURFACE_ID = "lifeform-native-model-table"
CONTINUOUS_REVIEW_SURFACE_ID = "seed-mass-alignment-review"
DISCRETE_REVIEW_SURFACE_ID = "lifeform-alignment-review"
CONTINUOUS_MODES = (
    "brownian",
    "white-noise",
    "pagel-lambda",
    "ornstein-uhlenbeck",
    "early-burst",
)
DISCRETE_MODELS = ("equal-rates", "symmetric", "all-rates-different")
REMOVED_TREE_TAXON = "Triglochin_maritimum"
EXTRA_TRAIT_TAXON = "unmatched_review_taxon"
CONTINUOUS_MISSING_VALUE_TAXON = "Juncus_maritimus"
DISCRETE_MISSING_VALUE_TAXON = "Juncus_gerardii"
PROVENANCE_CITATION = (
    "Vandelook, Filip; Janssens, Steven B.; Matthies, Diethart (2018). "
    "Data from: Ecological niche and phylogeny explain distribution of seed mass "
    "in the Central European flora. Dryad."
)
PROVENANCE_DOI = "10.5061/dryad.0st06f0"


def akaike_weights(aicc_by_model: dict[str, float]) -> dict[str, float]:
    finite_rows = {
        model: value for model, value in aicc_by_model.items() if math.isfinite(value)
    }
    minimum = min(finite_rows.values())
    unnormalized = {
        model: math.exp(-0.5 * (value - minimum))
        for model, value in finite_rows.items()
    }
    total = sum(unnormalized.values())
    return {model: weight / total for model, weight in unnormalized.items()}


def selected_model_from_aicc(aicc_by_model: dict[str, float]) -> str:
    return min(aicc_by_model.items(), key=lambda item: item[1])[0]


def apply_akaike_weights_from_report(
    rows: list[RealDatasetMacroevolutionModelRow],
    *,
    surface_id: str,
    engine: str,
) -> None:
    surface_rows = [row for row in rows if row.surface_id == surface_id]
    if engine == "bijux":
        weights = akaike_weights({row.model: row.bijux_aicc for row in surface_rows})
        for row in surface_rows:
            row.bijux_akaike_weight = weights[row.model]
        return
    weights = akaike_weights({row.model: row.geiger_aicc for row in surface_rows})
    for row in surface_rows:
        row.geiger_akaike_weight = weights[row.model]


def apply_geiger_akaike_weights(
    rows: list[RealDatasetMacroevolutionModelRow],
    *,
    surface_id: str,
) -> None:
    apply_akaike_weights_from_report(rows, surface_id=surface_id, engine="geiger")


def geiger_selected_weight(payload: dict[str, object]) -> float:
    weights = akaike_weights(
        {row["model"]: float(row["aicc"]) for row in payload["comparison_rows"]}
    )
    return weights[str(payload["selected_model"])]


def optional_float(value: object) -> float | None:
    if value in (None, "", []):
        return None
    return float(value)


def optional_str(value: object) -> str | None:
    if value in (None, "", []):
        return None
    return str(value)


def format_float(value: float) -> str:
    return f"{value:.12f}"


def format_optional_float(value: float | None) -> str:
    return "" if value is None else format_float(value)
