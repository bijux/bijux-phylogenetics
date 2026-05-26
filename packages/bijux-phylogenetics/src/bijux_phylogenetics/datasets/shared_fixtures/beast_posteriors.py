from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path

from ._paths import governed_fixture_root


@dataclass(frozen=True, slots=True)
class SharedBeastPosteriorBurninReference:
    """Expected retained-versus-burned counts for one governed burn-in fraction."""

    burnin_fraction: float
    burnin_row_count: int
    burnin_tree_count: int
    kept_row_count: int
    kept_tree_count: int


@dataclass(frozen=True, slots=True)
class SharedBeastPosteriorParameterReference:
    """Expected posterior summary statistics for one governed BEAST parameter."""

    parameter: str
    effective_sample_size: float
    mean: float
    median: float
    hpd_95_lower: float
    hpd_95_upper: float


@dataclass(frozen=True, slots=True)
class SharedBeastPosteriorConsensusReference:
    """Expected majority-rule consensus summary for the governed posterior corpus."""

    burnin_fraction: float
    annotated_node_count: int
    minimum_posterior_probability: float
    maximum_posterior_probability: float
    newick: str


@dataclass(frozen=True, slots=True)
class SharedBeastPosteriorMccReference:
    """Expected maximum clade credibility summary for the governed posterior corpus."""

    burnin_fraction: float
    selected_tree_index: int
    clade_credibility_score: float
    newick: str


@dataclass(frozen=True, slots=True)
class SharedBeastPosteriorReference:
    """Typed governed reference expectations for one BEAST posterior fixture."""

    burnin_reference: dict[float, SharedBeastPosteriorBurninReference]
    parameter_reference: dict[str, SharedBeastPosteriorParameterReference]
    consensus_reference: SharedBeastPosteriorConsensusReference
    mcc_reference: SharedBeastPosteriorMccReference


@dataclass(frozen=True, slots=True)
class SharedBeastPosteriorFixture:
    """One governed real-artifact BEAST posterior corpus."""

    fixture_id: str
    analysis_xml_relative_path: str
    posterior_log_relative_path: str
    posterior_trees_relative_path: str
    consensus_tree_relative_path: str
    mcc_tree_relative_path: str
    reference_json_relative_path: str
    beast_version: str
    recommended_burnin_fraction: float
    posterior_row_count: int
    posterior_tree_count: int
    shared_taxa: tuple[str, ...]
    feature_tags: tuple[str, ...]
    notes: str

    @property
    def analysis_xml_path(self) -> Path:
        return _fixtures_root() / self.analysis_xml_relative_path

    @property
    def posterior_log_path(self) -> Path:
        return _fixtures_root() / self.posterior_log_relative_path

    @property
    def posterior_trees_path(self) -> Path:
        return _fixtures_root() / self.posterior_trees_relative_path

    @property
    def consensus_tree_path(self) -> Path:
        return _fixtures_root() / self.consensus_tree_relative_path

    @property
    def mcc_tree_path(self) -> Path:
        return _fixtures_root() / self.mcc_tree_relative_path

    @property
    def reference_json_path(self) -> Path:
        return _fixtures_root() / self.reference_json_relative_path

    def load_reference(self) -> SharedBeastPosteriorReference:
        return _load_reference(self.reference_json_path)


def _fixtures_root() -> Path:
    return governed_fixture_root()


def _catalog_path() -> Path:
    return _fixtures_root() / "metadata" / "shared_beast_posterior_fixture_catalog.json"


def _load_catalog() -> dict[str, object]:
    return json.loads(_catalog_path().read_text(encoding="utf-8"))


def _load_reference(path: Path) -> SharedBeastPosteriorReference:
    payload = json.loads(path.read_text(encoding="utf-8"))
    consensus_payload = payload["consensus_reference_0.1"]
    mcc_payload = payload["mcc_reference_0.1"]
    return SharedBeastPosteriorReference(
        burnin_reference={
            float(fraction): SharedBeastPosteriorBurninReference(
                burnin_fraction=float(fraction),
                burnin_row_count=counts["burnin_row_count"],
                burnin_tree_count=counts["burnin_tree_count"],
                kept_row_count=counts["kept_row_count"],
                kept_tree_count=counts["kept_tree_count"],
            )
            for fraction, counts in payload["burnin_reference"].items()
        },
        parameter_reference={
            parameter: SharedBeastPosteriorParameterReference(
                parameter=parameter,
                effective_sample_size=summary["effective_sample_size"],
                mean=summary["mean"],
                median=summary["median"],
                hpd_95_lower=summary["hpd_95_lower"],
                hpd_95_upper=summary["hpd_95_upper"],
            )
            for parameter, summary in payload["parameter_reference_0.1"].items()
        },
        consensus_reference=SharedBeastPosteriorConsensusReference(
            burnin_fraction=0.1,
            annotated_node_count=consensus_payload["annotated_node_count"],
            minimum_posterior_probability=consensus_payload[
                "minimum_posterior_probability"
            ],
            maximum_posterior_probability=consensus_payload[
                "maximum_posterior_probability"
            ],
            newick=consensus_payload["newick"],
        ),
        mcc_reference=SharedBeastPosteriorMccReference(
            burnin_fraction=0.1,
            selected_tree_index=mcc_payload["selected_tree_index"],
            clade_credibility_score=mcc_payload["clade_credibility_score"],
            newick=mcc_payload["newick"],
        ),
    )


def list_shared_beast_posterior_fixtures() -> list[SharedBeastPosteriorFixture]:
    """Return the governed real-artifact BEAST posterior corpus."""
    catalog = _load_catalog()
    return [
        SharedBeastPosteriorFixture(
            fixture_id=entry["fixture_id"],
            analysis_xml_relative_path=entry["analysis_xml_relative_path"],
            posterior_log_relative_path=entry["posterior_log_relative_path"],
            posterior_trees_relative_path=entry["posterior_trees_relative_path"],
            consensus_tree_relative_path=entry["consensus_tree_relative_path"],
            mcc_tree_relative_path=entry["mcc_tree_relative_path"],
            reference_json_relative_path=entry["reference_json_relative_path"],
            beast_version=entry["beast_version"],
            recommended_burnin_fraction=entry["recommended_burnin_fraction"],
            posterior_row_count=entry["posterior_row_count"],
            posterior_tree_count=entry["posterior_tree_count"],
            shared_taxa=tuple(entry["shared_taxa"]),
            feature_tags=tuple(entry["feature_tags"]),
            notes=entry["notes"],
        )
        for entry in catalog["fixtures"]
    ]


def get_shared_beast_posterior_fixture(fixture_id: str) -> SharedBeastPosteriorFixture:
    """Resolve one governed BEAST posterior corpus by durable fixture id."""
    for fixture in list_shared_beast_posterior_fixtures():
        if fixture.fixture_id == fixture_id:
            return fixture
    supported = ", ".join(
        sorted(fixture.fixture_id for fixture in list_shared_beast_posterior_fixtures())
    )
    raise ValueError(
        f"unsupported shared BEAST posterior fixture '{fixture_id}'; expected one of: {supported}"
    )
