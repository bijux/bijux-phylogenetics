from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
import json
from pathlib import Path

from bijux_phylogenetics.runtime.errors import PhylogeneticsError

from .metropolis_hastings import MetropolisHastingsRunReport
from .state import (
    BayesianPhylogeneticState,
    BayesianTreeState,
    deserialize_bayesian_phylogenetic_state,
    serialize_bayesian_phylogenetic_state,
)

_POSTERIOR_TREE_SAMPLE_ARCHIVE_KIND = "bayesian-posterior-tree-sample-archive"


@dataclass(frozen=True, slots=True)
class BayesianPosteriorTreeSample:
    """One posterior tree sample with native Bayesian state metadata."""

    sample_index: int
    iteration_index: int | None
    model_id: str
    state: BayesianPhylogeneticState

    @property
    def tree(self) -> BayesianTreeState:
        return self.state.tree

    @property
    def posterior_log_score(self) -> float:
        return self.state.posterior_log_score


@dataclass(frozen=True, slots=True)
class BayesianPosteriorTreeSampleArchive:
    """One native Bijux posterior tree sample archive."""

    archive_kind: str
    sample_count: int
    samples: list[BayesianPosteriorTreeSample]


def infer_bayesian_model_id(
    *,
    state: BayesianPhylogeneticState,
) -> str:
    """Infer one stable model identifier from one sampled Bayesian state."""
    if not isinstance(state, BayesianPhylogeneticState):
        raise PhylogeneticsError(
            "bayesian posterior tree sample model-id inference requires one BayesianPhylogeneticState",
            code="bayesian_posterior_tree_sample_state_type_invalid",
        )
    categorical_parameters = state.model_parameters.categorical_parameters
    if not categorical_parameters:
        raise PhylogeneticsError(
            "bayesian posterior tree sample model-id inference requires at least one categorical model parameter or one explicit model_id override",
            code="bayesian_posterior_tree_sample_model_id_missing",
        )
    return "|".join(
        f"{parameter_name}={parameter_value}"
        for parameter_name, parameter_value in sorted(categorical_parameters.items())
    )


def build_bayesian_posterior_tree_sample(
    *,
    sample_index: int,
    state: BayesianPhylogeneticState,
    iteration_index: int | None = None,
    model_id: str | None = None,
) -> BayesianPosteriorTreeSample:
    """Build one validated posterior tree sample from one Bayesian state."""
    validated_sample_index = _validate_nonnegative_integer(
        value=sample_index,
        field_name="sample_index",
        owner_name="bayesian posterior tree sample",
    )
    validated_iteration_index = (
        _validate_nonnegative_integer(
            value=iteration_index,
            field_name="iteration_index",
            owner_name="bayesian posterior tree sample",
        )
        if iteration_index is not None
        else None
    )
    if not isinstance(state, BayesianPhylogeneticState):
        raise PhylogeneticsError(
            "bayesian posterior tree sample requires one BayesianPhylogeneticState",
            code="bayesian_posterior_tree_sample_state_type_invalid",
        )
    validated_model_id = _validate_nonblank_string(
        value=model_id
        if model_id is not None
        else infer_bayesian_model_id(state=state),
        field_name="model_id",
        owner_name="bayesian posterior tree sample",
    )
    return BayesianPosteriorTreeSample(
        sample_index=validated_sample_index,
        iteration_index=validated_iteration_index,
        model_id=validated_model_id,
        state=state,
    )


def build_bayesian_posterior_tree_sample_archive(
    *,
    samples: Sequence[BayesianPosteriorTreeSample],
) -> BayesianPosteriorTreeSampleArchive:
    """Build one validated native posterior tree sample archive."""
    validated_samples = list(samples)
    if not validated_samples:
        raise PhylogeneticsError(
            "bayesian posterior tree sample archive requires at least one sample",
            code="bayesian_posterior_tree_sample_archive_empty",
        )
    for sample in validated_samples:
        if not isinstance(sample, BayesianPosteriorTreeSample):
            raise PhylogeneticsError(
                "bayesian posterior tree sample archive requires every sample to be one BayesianPosteriorTreeSample",
                code="bayesian_posterior_tree_sample_archive_sample_type_invalid",
            )
    expected_sample_indices = list(range(len(validated_samples)))
    observed_sample_indices = [sample.sample_index for sample in validated_samples]
    if observed_sample_indices != expected_sample_indices:
        raise PhylogeneticsError(
            "bayesian posterior tree sample archive requires zero-based contiguous sample indices",
            code="bayesian_posterior_tree_sample_archive_sample_indices_invalid",
            details={
                "expected_sample_indices": expected_sample_indices,
                "observed_sample_indices": observed_sample_indices,
            },
        )
    return BayesianPosteriorTreeSampleArchive(
        archive_kind=_POSTERIOR_TREE_SAMPLE_ARCHIVE_KIND,
        sample_count=len(validated_samples),
        samples=validated_samples,
    )


def build_metropolis_hastings_posterior_tree_sample_archive(
    *,
    chain_report: MetropolisHastingsRunReport,
    model_id: str | None = None,
) -> BayesianPosteriorTreeSampleArchive:
    """Build one native posterior tree sample archive from one MH chain report."""
    if not isinstance(chain_report, MetropolisHastingsRunReport):
        raise PhylogeneticsError(
            "metropolis-hastings posterior tree sample archive builder requires one MetropolisHastingsRunReport",
            code="bayesian_posterior_tree_sample_archive_chain_report_type_invalid",
        )
    return build_bayesian_posterior_tree_sample_archive(
        samples=[
            build_bayesian_posterior_tree_sample(
                sample_index=sample_index,
                iteration_index=sample_index * chain_report.sample_every,
                model_id=model_id,
                state=state,
            )
            for sample_index, state in enumerate(chain_report.sampled_states)
        ]
    )


def write_bayesian_posterior_tree_sample_archive(
    path: Path | str,
    archive: BayesianPosteriorTreeSampleArchive,
) -> Path:
    """Write one native posterior tree sample archive as JSON."""
    if not isinstance(archive, BayesianPosteriorTreeSampleArchive):
        raise PhylogeneticsError(
            "bayesian posterior tree sample archive writer requires one BayesianPosteriorTreeSampleArchive",
            code="bayesian_posterior_tree_sample_archive_write_type_invalid",
        )
    archive_path = Path(path)
    archive_path.parent.mkdir(parents=True, exist_ok=True)
    archive_path.write_text(
        json.dumps(
            {
                "archive_kind": archive.archive_kind,
                "sample_count": archive.sample_count,
                "samples": [
                    {
                        "sample_index": sample.sample_index,
                        "iteration_index": sample.iteration_index,
                        "model_id": sample.model_id,
                        "state": serialize_bayesian_phylogenetic_state(sample.state),
                    }
                    for sample in archive.samples
                ],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return archive_path


def load_bayesian_posterior_tree_sample_archive(
    path: Path | str,
) -> BayesianPosteriorTreeSampleArchive:
    """Load one native posterior tree sample archive from JSON."""
    archive_path = Path(path)
    payload = json.loads(archive_path.read_text(encoding="utf-8"))
    if not isinstance(payload, Mapping):
        raise PhylogeneticsError(
            "bayesian posterior tree sample archive json must decode to one mapping",
            code="bayesian_posterior_tree_sample_archive_payload_type_invalid",
        )
    archive_kind = _require_string(payload, key="archive_kind")
    if archive_kind != _POSTERIOR_TREE_SAMPLE_ARCHIVE_KIND:
        raise PhylogeneticsError(
            "bayesian posterior tree sample archive loader requires archive_kind 'bayesian-posterior-tree-sample-archive'",
            code="bayesian_posterior_tree_sample_archive_kind_invalid",
            details={"archive_kind": archive_kind},
        )
    sample_payloads = _require_list(payload, key="samples")
    archive = build_bayesian_posterior_tree_sample_archive(
        samples=[
            build_bayesian_posterior_tree_sample(
                sample_index=_require_integer(sample_payload, key="sample_index"),
                iteration_index=_optional_integer(
                    sample_payload.get("iteration_index")
                ),
                model_id=_require_string(sample_payload, key="model_id"),
                state=deserialize_bayesian_phylogenetic_state(
                    _require_mapping(sample_payload, key="state")
                ),
            )
            for sample_payload in (
                _require_mapping(
                    item,
                    owner_name="bayesian posterior tree sample archive sample payload",
                )
                for item in sample_payloads
            )
        ]
    )
    serialized_sample_count = _require_integer(payload, key="sample_count")
    if serialized_sample_count != archive.sample_count:
        raise PhylogeneticsError(
            "bayesian posterior tree sample archive sample_count must match the loaded sample count",
            code="bayesian_posterior_tree_sample_archive_count_mismatch",
            details={
                "serialized_sample_count": serialized_sample_count,
                "loaded_sample_count": archive.sample_count,
            },
        )
    return archive


def _validate_nonnegative_integer(
    *,
    value: object,
    field_name: str,
    owner_name: str,
) -> int:
    if not isinstance(value, int):
        raise PhylogeneticsError(
            f"{owner_name} requires {field_name} to be one integer",
            code="bayesian_posterior_tree_sample_integer_invalid",
            details={"field_name": field_name, "value": value},
        )
    if value < 0:
        raise PhylogeneticsError(
            f"{owner_name} requires {field_name} to be nonnegative",
            code="bayesian_posterior_tree_sample_integer_negative",
            details={"field_name": field_name, "value": value},
        )
    return value


def _validate_nonblank_string(
    *,
    value: object,
    field_name: str,
    owner_name: str,
) -> str:
    if not isinstance(value, str):
        raise PhylogeneticsError(
            f"{owner_name} requires {field_name} to be one string",
            code="bayesian_posterior_tree_sample_string_invalid",
            details={"field_name": field_name, "value": value},
        )
    if value.strip() == "":
        raise PhylogeneticsError(
            f"{owner_name} requires {field_name} to be nonblank",
            code="bayesian_posterior_tree_sample_string_blank",
            details={"field_name": field_name},
        )
    return value


def _require_integer(payload: Mapping[str, object], *, key: str) -> int:
    if key not in payload:
        raise PhylogeneticsError(
            f"bayesian posterior tree sample archive payload requires key '{key}'",
            code="bayesian_posterior_tree_sample_archive_key_missing",
            details={"key": key},
        )
    return _validate_nonnegative_integer(
        value=payload[key],
        field_name=key,
        owner_name="bayesian posterior tree sample archive payload",
    )


def _optional_integer(value: object) -> int | None:
    if value is None:
        return None
    return _validate_nonnegative_integer(
        value=value,
        field_name="iteration_index",
        owner_name="bayesian posterior tree sample archive payload",
    )


def _require_string(payload: Mapping[str, object], *, key: str) -> str:
    if key not in payload:
        raise PhylogeneticsError(
            f"bayesian posterior tree sample archive payload requires key '{key}'",
            code="bayesian_posterior_tree_sample_archive_key_missing",
            details={"key": key},
        )
    return _validate_nonblank_string(
        value=payload[key],
        field_name=key,
        owner_name="bayesian posterior tree sample archive payload",
    )


def _require_list(payload: Mapping[str, object], *, key: str) -> list[object]:
    value = payload.get(key)
    if not isinstance(value, list):
        raise PhylogeneticsError(
            f"bayesian posterior tree sample archive payload requires key '{key}' to be one list",
            code="bayesian_posterior_tree_sample_archive_list_invalid",
            details={"key": key},
        )
    return value


def _require_mapping(
    payload: Mapping[str, object],
    *,
    key: str | None = None,
    owner_name: str | None = None,
) -> Mapping[str, object]:
    value: object = payload if key is None else payload.get(key)
    if not isinstance(value, Mapping):
        subject = owner_name or f"bayesian posterior tree sample archive key '{key}'"
        raise PhylogeneticsError(
            f"{subject} requires one mapping",
            code="bayesian_posterior_tree_sample_archive_mapping_invalid",
            details={"key": key},
        )
    return value
