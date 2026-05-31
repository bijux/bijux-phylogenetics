from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

from bijux_phylogenetics.engines.common import build_file_checksums
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

from .metropolis_hastings import MetropolisHastingsRunReport
from .state import BayesianPhylogeneticState, serialize_bayesian_phylogenetic_state

BAYESIAN_BURNIN_POLICY_NAMES = ("none", "discarded-samples", "adaptive-tuning")


@dataclass(frozen=True, slots=True)
class BayesianRunPriorRow:
    """One prior configuration row recorded in one Bayesian run manifest."""

    prior_name: str
    family: str
    parameter_values: dict[str, float]


@dataclass(frozen=True, slots=True)
class BayesianRunBurninPolicy:
    """One explicit burn-in policy recorded in one Bayesian run manifest."""

    policy_name: str
    discarded_iteration_count: int
    discarded_sample_count: int
    freeze_iteration_index: int | None = None


@dataclass(frozen=True, slots=True)
class BayesianRunManifest:
    """One reproducibility manifest for one Bayesian run."""

    manifest_kind: str
    run_kind: str
    model_name: str
    model_configuration: dict[str, object]
    prior_rows: list[BayesianRunPriorRow]
    proposal_schedule: dict[str, object]
    seed: int
    chain_count: int
    retained_sample_count: int
    burnin_policy: BayesianRunBurninPolicy
    execution_configuration: dict[str, object]
    retained_sample_ids: list[str]
    input_paths: list[str]
    input_checksums: dict[str, str]
    output_paths: list[str]
    output_checksums: dict[str, str]


def build_bayesian_run_burnin_policy(
    *,
    policy_name: str,
    discarded_iteration_count: int = 0,
    discarded_sample_count: int = 0,
    freeze_iteration_index: int | None = None,
) -> BayesianRunBurninPolicy:
    """Build one validated Bayesian burn-in policy record."""
    validated_policy_name = _validate_policy_name(policy_name)
    validated_discarded_iteration_count = _validate_nonnegative_integer(
        value=discarded_iteration_count,
        field_name="discarded_iteration_count",
        owner_name="bayesian run burn-in policy",
    )
    validated_discarded_sample_count = _validate_nonnegative_integer(
        value=discarded_sample_count,
        field_name="discarded_sample_count",
        owner_name="bayesian run burn-in policy",
    )
    validated_freeze_iteration_index = _validate_optional_nonnegative_integer(
        value=freeze_iteration_index,
        field_name="freeze_iteration_index",
        owner_name="bayesian run burn-in policy",
    )
    if validated_policy_name == "none":
        if validated_discarded_iteration_count != 0:
            raise PhylogeneticsError(
                "bayesian run burn-in policy 'none' requires discarded_iteration_count to equal zero",
                code="bayesian_run_manifest_burnin_discarded_iterations_unexpected",
            )
        if validated_discarded_sample_count != 0:
            raise PhylogeneticsError(
                "bayesian run burn-in policy 'none' requires discarded_sample_count to equal zero",
                code="bayesian_run_manifest_burnin_discarded_samples_unexpected",
            )
        if validated_freeze_iteration_index is not None:
            raise PhylogeneticsError(
                "bayesian run burn-in policy 'none' does not accept freeze_iteration_index",
                code="bayesian_run_manifest_burnin_freeze_iteration_unexpected",
            )
    return BayesianRunBurninPolicy(
        policy_name=validated_policy_name,
        discarded_iteration_count=validated_discarded_iteration_count,
        discarded_sample_count=validated_discarded_sample_count,
        freeze_iteration_index=validated_freeze_iteration_index,
    )


def list_metropolis_hastings_retained_sample_ids(
    *,
    chain_report: MetropolisHastingsRunReport,
) -> list[str]:
    """List canonical retained sample identifiers for one Metropolis-Hastings chain."""
    if not isinstance(chain_report, MetropolisHastingsRunReport):
        raise PhylogeneticsError(
            "bayesian retained sample identifiers require one MetropolisHastingsRunReport",
            code="bayesian_run_manifest_chain_report_type_invalid",
        )
    return [
        _fingerprint_state(sampled_state) for sampled_state in chain_report.sampled_states
    ]


def build_bayesian_run_manifest(
    *,
    run_kind: str,
    model_name: str,
    model_configuration: Mapping[str, object],
    prior_rows: Sequence[BayesianRunPriorRow],
    proposal_schedule: Mapping[str, object],
    seed: int,
    chain_count: int,
    burnin_policy: BayesianRunBurninPolicy,
    execution_configuration: Mapping[str, object],
    retained_sample_ids: Sequence[str],
    input_paths: Sequence[Path | str],
    output_paths: Sequence[Path | str],
) -> BayesianRunManifest:
    """Build one validated Bayesian run manifest."""
    validated_run_kind = _validate_nonblank_string(
        value=run_kind,
        field_name="run_kind",
        owner_name="bayesian run manifest",
    )
    validated_model_name = _validate_nonblank_string(
        value=model_name,
        field_name="model_name",
        owner_name="bayesian run manifest",
    )
    validated_prior_rows = _validate_prior_rows(prior_rows)
    if not isinstance(burnin_policy, BayesianRunBurninPolicy):
        raise PhylogeneticsError(
            "bayesian run manifest requires one BayesianRunBurninPolicy",
            code="bayesian_run_manifest_burnin_policy_type_invalid",
        )
    validated_seed = _validate_nonnegative_integer(
        value=seed,
        field_name="seed",
        owner_name="bayesian run manifest",
    )
    validated_chain_count = _validate_positive_integer(
        value=chain_count,
        field_name="chain_count",
        owner_name="bayesian run manifest",
    )
    validated_retained_sample_ids = _validate_retained_sample_ids(retained_sample_ids)
    validated_input_paths = _validate_existing_file_paths(
        paths=input_paths,
        field_name="input_paths",
        owner_name="bayesian run manifest",
        require_nonempty=False,
    )
    validated_output_paths = _validate_existing_file_paths(
        paths=output_paths,
        field_name="output_paths",
        owner_name="bayesian run manifest",
        require_nonempty=True,
    )
    return BayesianRunManifest(
        manifest_kind="bayesian-run",
        run_kind=validated_run_kind,
        model_name=validated_model_name,
        model_configuration=_json_ready_mapping(
            value=model_configuration,
            field_name="model_configuration",
            owner_name="bayesian run manifest",
        ),
        prior_rows=list(validated_prior_rows),
        proposal_schedule=_json_ready_mapping(
            value=proposal_schedule,
            field_name="proposal_schedule",
            owner_name="bayesian run manifest",
        ),
        seed=validated_seed,
        chain_count=validated_chain_count,
        retained_sample_count=len(validated_retained_sample_ids),
        burnin_policy=burnin_policy,
        execution_configuration=_json_ready_mapping(
            value=execution_configuration,
            field_name="execution_configuration",
            owner_name="bayesian run manifest",
        ),
        retained_sample_ids=list(validated_retained_sample_ids),
        input_paths=[str(path) for path in validated_input_paths],
        input_checksums=build_file_checksums(list(validated_input_paths)),
        output_paths=[str(path) for path in validated_output_paths],
        output_checksums=build_file_checksums(list(validated_output_paths)),
    )


def write_bayesian_run_manifest(
    path: Path,
    manifest: BayesianRunManifest,
) -> Path:
    """Write one Bayesian run manifest to one JSON file."""
    if not isinstance(manifest, BayesianRunManifest):
        raise PhylogeneticsError(
            "bayesian run manifest writer requires one BayesianRunManifest",
            code="bayesian_run_manifest_write_type_invalid",
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(asdict(manifest), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def load_bayesian_run_manifest(path: Path | str) -> BayesianRunManifest:
    """Load one Bayesian run manifest from one JSON file."""
    manifest_path = Path(path)
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    if not isinstance(payload, Mapping):
        raise PhylogeneticsError(
            "bayesian run manifest json must decode to one mapping",
            code="bayesian_run_manifest_payload_type_invalid",
        )
    return _load_manifest_payload(payload)


def _load_manifest_payload(payload: Mapping[str, object]) -> BayesianRunManifest:
    manifest_kind = _require_string(payload, key="manifest_kind")
    if manifest_kind != "bayesian-run":
        raise PhylogeneticsError(
            "bayesian run manifest loader requires manifest_kind 'bayesian-run'",
            code="bayesian_run_manifest_kind_invalid",
            details={"manifest_kind": manifest_kind},
        )
    prior_payloads = _require_list(payload, key="prior_rows")
    burnin_payload = _require_mapping(payload, key="burnin_policy")
    model_configuration = _require_mapping(payload, key="model_configuration")
    proposal_schedule = _require_mapping(payload, key="proposal_schedule")
    execution_configuration = _require_mapping(payload, key="execution_configuration")
    retained_sample_ids = _require_string_list(payload, key="retained_sample_ids")
    input_paths = _require_string_list(payload, key="input_paths")
    output_paths = _require_string_list(payload, key="output_paths")
    input_checksums = _require_string_mapping(payload, key="input_checksums")
    output_checksums = _require_string_mapping(payload, key="output_checksums")
    retained_sample_count = _require_integer(payload, key="retained_sample_count")
    if retained_sample_count != len(retained_sample_ids):
        raise PhylogeneticsError(
            "bayesian run manifest retained_sample_count must match retained_sample_ids length",
            code="bayesian_run_manifest_retained_sample_count_mismatch",
            details={
                "retained_sample_count": retained_sample_count,
                "retained_sample_id_count": len(retained_sample_ids),
            },
        )
    return BayesianRunManifest(
        manifest_kind=manifest_kind,
        run_kind=_require_string(payload, key="run_kind"),
        model_name=_require_string(payload, key="model_name"),
        model_configuration=dict(model_configuration),
        prior_rows=[
            BayesianRunPriorRow(
                prior_name=_require_string(prior_payload, key="prior_name"),
                family=_require_string(prior_payload, key="family"),
                parameter_values=_require_float_mapping(
                    prior_payload,
                    key="parameter_values",
                ),
            )
            for prior_payload in (
                _require_mapping(item, owner_name="bayesian run prior row payload")
                for item in prior_payloads
            )
        ],
        proposal_schedule=dict(proposal_schedule),
        seed=_require_integer(payload, key="seed"),
        chain_count=_require_integer(payload, key="chain_count"),
        retained_sample_count=retained_sample_count,
        burnin_policy=build_bayesian_run_burnin_policy(
            policy_name=_require_string(burnin_payload, key="policy_name"),
            discarded_iteration_count=_require_integer(
                burnin_payload,
                key="discarded_iteration_count",
            ),
            discarded_sample_count=_require_integer(
                burnin_payload,
                key="discarded_sample_count",
            ),
            freeze_iteration_index=_optional_integer(
                burnin_payload.get("freeze_iteration_index")
            ),
        ),
        execution_configuration=dict(execution_configuration),
        retained_sample_ids=retained_sample_ids,
        input_paths=input_paths,
        input_checksums=input_checksums,
        output_paths=output_paths,
        output_checksums=output_checksums,
    )


def _fingerprint_state(state: BayesianPhylogeneticState) -> str:
    return _fingerprint_payload(serialize_bayesian_phylogenetic_state(state))


def _fingerprint_payload(payload: object) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True).encode("utf-8")
    ).hexdigest()


def _validate_policy_name(value: str) -> str:
    validated_value = _validate_nonblank_string(
        value=value,
        field_name="policy_name",
        owner_name="bayesian run burn-in policy",
    )
    if validated_value not in BAYESIAN_BURNIN_POLICY_NAMES:
        raise PhylogeneticsError(
            "bayesian run burn-in policy name is unsupported",
            code="bayesian_run_manifest_burnin_policy_name_invalid",
            details={
                "policy_name": validated_value,
                "allowed_policy_names": list(BAYESIAN_BURNIN_POLICY_NAMES),
            },
        )
    return validated_value


def _validate_prior_rows(
    prior_rows: Sequence[BayesianRunPriorRow],
) -> tuple[BayesianRunPriorRow, ...]:
    validated_rows = tuple(prior_rows)
    if not validated_rows:
        raise PhylogeneticsError(
            "bayesian run manifest requires at least one prior row",
            code="bayesian_run_manifest_prior_rows_empty",
        )
    if any(not isinstance(row, BayesianRunPriorRow) for row in validated_rows):
        raise PhylogeneticsError(
            "bayesian run manifest requires every prior row to be one BayesianRunPriorRow",
            code="bayesian_run_manifest_prior_row_type_invalid",
        )
    return validated_rows


def _validate_retained_sample_ids(
    retained_sample_ids: Sequence[str],
) -> tuple[str, ...]:
    validated_values = tuple(retained_sample_ids)
    if not validated_values:
        raise PhylogeneticsError(
            "bayesian run manifest requires at least one retained sample identifier",
            code="bayesian_run_manifest_retained_sample_ids_empty",
        )
    if any(not isinstance(value, str) or not value.strip() for value in validated_values):
        raise PhylogeneticsError(
            "bayesian run manifest requires every retained sample identifier to be one nonblank string",
            code="bayesian_run_manifest_retained_sample_id_invalid",
        )
    return validated_values


def _validate_existing_file_paths(
    *,
    paths: Sequence[Path | str],
    field_name: str,
    owner_name: str,
    require_nonempty: bool,
) -> tuple[Path, ...]:
    validated_paths = tuple(Path(path) for path in paths)
    if require_nonempty and not validated_paths:
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to contain at least one file path",
            code="bayesian_run_manifest_file_paths_empty",
            details={"field_name": field_name},
        )
    for path in validated_paths:
        if not path.exists() or not path.is_file():
            raise PhylogeneticsError(
                f"{owner_name} requires every '{field_name}' entry to be one existing file",
                code="bayesian_run_manifest_file_path_missing",
                details={"field_name": field_name, "path": str(path)},
            )
    return validated_paths


def _json_ready_mapping(
    *,
    value: Mapping[str, object],
    field_name: str,
    owner_name: str,
) -> dict[str, object]:
    if not isinstance(value, Mapping):
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be one mapping",
            code="bayesian_run_manifest_mapping_type_invalid",
            details={"field_name": field_name},
        )
    try:
        normalized_value = json.loads(json.dumps(dict(value), sort_keys=True))
    except TypeError as error:
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be JSON-serializable",
            code="bayesian_run_manifest_mapping_not_json_ready",
            details={"field_name": field_name},
        ) from error
    if not isinstance(normalized_value, dict):
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to serialize to one mapping",
            code="bayesian_run_manifest_mapping_json_shape_invalid",
            details={"field_name": field_name},
        )
    return normalized_value


def _validate_nonblank_string(
    *,
    value: str,
    field_name: str,
    owner_name: str,
) -> str:
    if not isinstance(value, str) or not value.strip():
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be one nonblank string",
            code="bayesian_run_manifest_string_invalid",
            details={"field_name": field_name},
        )
    return value.strip()


def _validate_positive_integer(
    *,
    value: int,
    field_name: str,
    owner_name: str,
) -> int:
    if not isinstance(value, int) or value <= 0:
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be one positive integer",
            code="bayesian_run_manifest_positive_integer_invalid",
            details={"field_name": field_name, "value": value},
        )
    return value


def _validate_nonnegative_integer(
    *,
    value: int,
    field_name: str,
    owner_name: str,
) -> int:
    if not isinstance(value, int) or value < 0:
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be one nonnegative integer",
            code="bayesian_run_manifest_nonnegative_integer_invalid",
            details={"field_name": field_name, "value": value},
        )
    return value


def _validate_optional_nonnegative_integer(
    *,
    value: int | None,
    field_name: str,
    owner_name: str,
) -> int | None:
    if value is None:
        return None
    return _validate_nonnegative_integer(
        value=value,
        field_name=field_name,
        owner_name=owner_name,
    )


def _require_mapping(
    payload: Mapping[str, object],
    *,
    key: str,
    owner_name: str | None = None,
) -> Mapping[str, object]:
    value = payload.get(key)
    if not isinstance(value, Mapping):
        subject = owner_name or f"'{key}'"
        raise PhylogeneticsError(
            f"bayesian run manifest loader requires {subject} to be one mapping",
            code="bayesian_run_manifest_load_mapping_invalid",
            details={"key": key},
        )
    return value


def _require_list(payload: Mapping[str, object], *, key: str) -> list[object]:
    value = payload.get(key)
    if not isinstance(value, list):
        raise PhylogeneticsError(
            f"bayesian run manifest loader requires '{key}' to be one list",
            code="bayesian_run_manifest_load_list_invalid",
            details={"key": key},
        )
    return value


def _require_string(payload: Mapping[str, object], *, key: str) -> str:
    value = payload.get(key)
    if not isinstance(value, str) or not value.strip():
        raise PhylogeneticsError(
            f"bayesian run manifest loader requires '{key}' to be one nonblank string",
            code="bayesian_run_manifest_load_string_invalid",
            details={"key": key},
        )
    return value


def _require_integer(payload: Mapping[str, object], *, key: str) -> int:
    value = payload.get(key)
    if not isinstance(value, int):
        raise PhylogeneticsError(
            f"bayesian run manifest loader requires '{key}' to be one integer",
            code="bayesian_run_manifest_load_integer_invalid",
            details={"key": key},
        )
    return value


def _optional_integer(value: object) -> int | None:
    if value is None:
        return None
    if not isinstance(value, int):
        raise PhylogeneticsError(
            "bayesian run manifest loader requires optional integer fields to be integers when present",
            code="bayesian_run_manifest_load_optional_integer_invalid",
        )
    return value


def _require_string_list(payload: Mapping[str, object], *, key: str) -> list[str]:
    raw_values = _require_list(payload, key=key)
    if any(not isinstance(value, str) for value in raw_values):
        raise PhylogeneticsError(
            f"bayesian run manifest loader requires '{key}' to contain only strings",
            code="bayesian_run_manifest_load_string_list_invalid",
            details={"key": key},
        )
    return list(raw_values)


def _require_string_mapping(
    payload: Mapping[str, object],
    *,
    key: str,
) -> dict[str, str]:
    value = _require_mapping(payload, key=key)
    if any(
        not isinstance(mapping_key, str) or not isinstance(mapping_value, str)
        for mapping_key, mapping_value in value.items()
    ):
        raise PhylogeneticsError(
            f"bayesian run manifest loader requires '{key}' to contain only string keys and values",
            code="bayesian_run_manifest_load_string_mapping_invalid",
            details={"key": key},
        )
    return dict(value)


def _require_float_mapping(
    payload: Mapping[str, object],
    *,
    key: str,
) -> dict[str, float]:
    value = _require_mapping(payload, key=key)
    if any(
        not isinstance(mapping_key, str)
        or not isinstance(mapping_value, (int, float))
        for mapping_key, mapping_value in value.items()
    ):
        raise PhylogeneticsError(
            f"bayesian run manifest loader requires '{key}' to contain only string keys and numeric values",
            code="bayesian_run_manifest_load_float_mapping_invalid",
            details={"key": key},
        )
    return {mapping_key: float(mapping_value) for mapping_key, mapping_value in value.items()}


__all__ = [
    "BAYESIAN_BURNIN_POLICY_NAMES",
    "BayesianRunBurninPolicy",
    "BayesianRunManifest",
    "BayesianRunPriorRow",
    "build_bayesian_run_burnin_policy",
    "build_bayesian_run_manifest",
    "list_metropolis_hastings_retained_sample_ids",
    "load_bayesian_run_manifest",
    "write_bayesian_run_manifest",
]
