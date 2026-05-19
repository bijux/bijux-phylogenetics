from __future__ import annotations

import csv
from dataclasses import replace
import json
from pathlib import Path
import re
import shutil
import subprocess  # nosec B404 - parity helpers invoke repository-owned reference commands
import tempfile

from bijux_phylogenetics.ancestral.common import load_discrete_dataset, node_signature
from bijux_phylogenetics.ancestral.discrete import (
    reconstruct_discrete_ancestral_states,
)
from bijux_phylogenetics.compare.structural_parity import (
    compare_tree_sets_structurally,
    compare_tree_structurally,
)
from bijux_phylogenetics.core._node_identity import build_ape_internal_node_map
from bijux_phylogenetics.core.pruning import (
    drop_tree_taxa,
    prune_tree_to_requested_taxa,
)
from bijux_phylogenetics.core.topology import (
    extract_tree_clade_by_node_id,
    root_tree_on_outgroup,
    unroot_tree,
)
from bijux_phylogenetics.core.tree import TreeNode
from bijux_phylogenetics.distance import build_tree_from_imported_distance_matrix
from bijux_phylogenetics.io.newick import dumps_newick, load_newick_tree_set
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.trees import (
    compute_consensus_tree,
    compute_strict_consensus_tree,
)
from ..registry import ApeParityCase, _selected_cases, _write_case_file
from .models import (
    ApeParityObservation,
    ApeParityReport,
    ApeParitySummaryRow as ApeParitySummaryRow,
)
from .reporting import (
    build_ape_parity_report,
    write_ape_parity_observation_table as write_ape_parity_observation_table,
    write_ape_parity_summary_table as write_ape_parity_summary_table,
)
from .comparative_payloads import (
    _build_bijux_brownian_covariance_rows,
    _build_bijux_branching_time_rows,
    _build_bijux_continuous_ancestral_rows,
    _build_bijux_diversification_gamma_rows,
    _build_bijux_independent_contrast_rows,
    _build_bijux_tree_node_depth_rows,
    _build_bijux_tree_simulation_envelope_rows,
    _build_bijux_tree_ultrametric_rows,
)
from .sequence_payloads import (
    _build_bijux_base_frequency_summary,
    _build_bijux_distance_rows,
    _build_bijux_dnabin_rows,
    _build_bijux_neighbor_joining_structure,
    _build_bijux_segregating_site_rows,
    _build_bijux_translation_rows,
)
from .tree_payloads import (
    _build_bijux_consensus_rows,
    _build_bijux_drop_tip_structure,
    _build_bijux_extract_clade_structure,
    _build_bijux_keep_tip_structure,
    _build_bijux_monophyly_summary,
    _build_bijux_mrca_summary,
    _build_bijux_prop_clades_rows,
    _build_bijux_root_outgroup_structure,
    _build_bijux_topology_distance_rows,
    _build_bijux_tree_set_structure,
    _build_bijux_tree_structure,
    _build_bijux_tree_tip_distance_rows,
    _build_bijux_unroot_structure,
    _materialize_reference_input,
    _sort_parity_rows,
)
from .runtime import (
    ape_runner_path as _ape_runner_path,
    bijux_commit as _bijux_commit,
    bijux_version as _bijux_version,
    failure_root as _failure_root,
    reference_environment as _reference_environment,
    repository_root as _repository_root,
)





























def _build_bijux_discrete_ancestral_rows(
    input_fixture: Path,
    *,
    trait_table_path: Path,
    trait_name: str,
    trait_taxon_column: str,
    ancestral_model: str,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    dataset = load_discrete_dataset(
        input_fixture,
        trait_table_path,
        trait=trait_name,
        taxon_column=trait_taxon_column,
    )
    report = reconstruct_discrete_ancestral_states(
        input_fixture,
        trait_table_path,
        trait=trait_name,
        taxon_column=trait_taxon_column,
        model=ancestral_model,
    )
    internal_node_map = {
        node_signature(node): node_id
        for node_id, node in build_ape_internal_node_map(dataset.tree).items()
    }
    rows = sorted(
        [
            {
                "node_id": internal_node_map[estimate.node],
                "node": estimate.node,
                "state": _coerce_table_cell(state),
                "posterior_probability": probability,
                "most_likely_state": _coerce_table_cell(estimate.most_likely_state),
                "max_posterior_probability": estimate.confidence,
            }
            for estimate in report.estimates
            if not estimate.is_tip
            for state, probability in sorted(estimate.state_probabilities.items())
        ],
        key=lambda row: (int(row["node_id"]), str(row["state"])),
    )
    transition_rows = [
        {
            "source_state": row.source_state,
            "target_state": row.target_state,
            "transition_allowed": row.transition_allowed,
            "step_distance": row.step_distance,
            "rate": row.rate,
        }
        for row in report.transition_rate_rows
    ]
    return {
        "trait": report.trait,
        "taxon_count": report.taxon_count,
        "excluded_taxon_count": len(report.dropped_missing_taxa),
        "dropped_missing_taxa": report.dropped_missing_taxa,
        "internal_node_count": len(
            [estimate for estimate in report.estimates if not estimate.is_tip]
        ),
        "model": report.model,
        "state_count": len(report.observed_states),
        "state_labels": report.observed_states,
        "log_likelihood": report.log_likelihood,
        "parameter_count": report.parameter_count,
        "aic": report.aic,
        "overparameterized": report.overparameterized,
        "baseline_model": (
            None
            if report.baseline_comparison is None
            else report.baseline_comparison.baseline_model
        ),
        "baseline_delta_aic": (
            None
            if report.baseline_comparison is None
            else report.baseline_comparison.delta_aic
        ),
        "preferred_model_by_aic": (
            None
            if report.baseline_comparison is None
            else report.baseline_comparison.preferred_model_by_aic
        ),
        "transition_rate_rows": transition_rows,
    }, rows








def _load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _normalize_newick_label(label: str) -> str:
    if len(label) >= 2 and label.startswith("'") and label.endswith("'"):
        return label[1:-1].replace("''", "'")
    return label


def _normalize_expected_label(
    label: str,
    *,
    expected_tip_labels: set[str] | None,
) -> str:
    normalized = _normalize_newick_label(label)
    if (
        expected_tip_labels
        and normalized not in expected_tip_labels
        and normalized.replace("_", " ") in expected_tip_labels
    ):
        return normalized.replace("_", " ")
    return normalized


def _normalize_joined_labels(
    value: str,
    *,
    expected_tip_labels: set[str] | None,
) -> str:
    if value == "":
        return value
    labels = [
        _normalize_expected_label(label, expected_tip_labels=expected_tip_labels)
        for label in value.split("|")
    ]
    return "|".join(sorted(labels))


def _normalize_reference_summary(summary: dict[str, object]) -> dict[str, object]:
    normalized = dict(summary)
    tip_labels = normalized.get("tip_labels")
    if isinstance(tip_labels, list):
        expected_tip_labels = {
            _normalize_newick_label(str(label)) for label in tip_labels
        }
        normalized["tip_labels"] = [
            _normalize_expected_label(
                str(label),
                expected_tip_labels=expected_tip_labels,
            )
            for label in tip_labels
        ]
    return normalized


def _summary_rooted_flag(summary: dict[str, object]) -> bool:
    rooted = summary.get("rooted")
    if isinstance(rooted, bool):
        return rooted
    raise ValueError("reference summary must include a boolean rooted flag")


def _optional_payload_string(payload: dict[str, object], key: str) -> str | None:
    value = payload.get(key)
    return value if isinstance(value, str) else None


def _coerce_table_cell(value: str) -> object:
    if value == "":
        return ""
    if value.lower() in {"true", "false"}:
        return value.lower() == "true"
    if re.fullmatch(r"-?\d+", value):
        return int(value)
    if re.fullmatch(r"-?(?:\d+\.\d*|\d*\.\d+)(?:[eE][+-]?\d+)?", value):
        return float(value)
    return value



def _load_rows_table(
    path: Path,
    *,
    expected_tip_labels: set[str] | None = None,
    sort_rows: bool = False,
) -> list[dict[str, object]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    normalized_rows = [
        {
            key: _normalize_expected_label(
                value, expected_tip_labels=expected_tip_labels
            )
            if key.endswith("label")
            else _normalize_joined_labels(
                value,
                expected_tip_labels=expected_tip_labels,
            )
            if key in {"clade_id", "taxa", "descendant_taxa", "shared_taxa"}
            else _coerce_table_cell(value)
            for key, value in row.items()
        }
        for row in rows
    ]
    if sort_rows:
        return _sort_parity_rows(normalized_rows)
    return normalized_rows


def _compare_scalar(expected: object, observed: object, *, tolerance: float) -> bool:
    if isinstance(expected, (int, float)) and isinstance(observed, (int, float)):
        return abs(float(expected) - float(observed)) <= tolerance
    return expected == observed


def _compare_json(expected: object, observed: object, *, tolerance: float) -> bool:
    if isinstance(expected, dict) and isinstance(observed, dict):
        if set(expected) != set(observed):
            return False
        return all(
            _compare_json(expected[key], observed[key], tolerance=tolerance)
            for key in expected
        )
    if isinstance(expected, list) and isinstance(observed, list):
        if len(expected) != len(observed):
            return False
        return all(
            _compare_json(left, right, tolerance=tolerance)
            for left, right in zip(expected, observed, strict=True)
        )
    return _compare_scalar(expected, observed, tolerance=tolerance)


def _normalize_tree_labels(
    node: TreeNode,
    *,
    expected_tip_labels: set[str] | None,
) -> None:
    if node.name is not None:
        normalized = _normalize_newick_label(node.name)
        if (
            expected_tip_labels
            and normalized not in expected_tip_labels
            and normalized.replace("_", " ") in expected_tip_labels
        ):
            normalized = normalized.replace("_", " ")
        node.name = normalized
    for child in node.children:
        _normalize_tree_labels(child, expected_tip_labels=expected_tip_labels)


def _clear_branch_lengths(node: TreeNode) -> None:
    node.branch_length = None
    for child in node.children:
        _clear_branch_lengths(child)


def _canonical_newick(
    path: Path,
    *,
    expected_tip_labels: set[str] | None = None,
) -> str:
    tree = load_tree(path)
    _normalize_tree_labels(tree.root, expected_tip_labels=expected_tip_labels)
    return dumps_newick(tree)


def _copy_if_exists(source: Path, destination: Path) -> None:
    if source.exists():
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)


def _write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )


def _write_rows_table(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def _persist_failure_bundle(
    *,
    failure_root: Path,
    case: ApeParityCase,
    case_file: Path,
    execution_root: Path,
    execution_payload: dict[str, object] | None,
    reference_summary: dict[str, object] | None,
    bijux_summary: dict[str, object] | None,
    reference_error: dict[str, object] | None,
    bijux_error: dict[str, object] | None,
    reference_rows: list[dict[str, object]] | None,
    bijux_rows: list[dict[str, object]] | None,
    bijux_normalized_text: str | None,
    mismatch_reason: str,
) -> Path:
    artifact_root = failure_root / case.case_id
    artifact_root.mkdir(parents=True, exist_ok=True)
    _copy_if_exists(case_file, artifact_root / "case.json")
    _copy_if_exists(
        execution_root.parent / "bijux-reference-input.nwk",
        artifact_root / "bijux-reference-input.nwk",
    )
    _copy_if_exists(
        execution_root / "reference-execution.json",
        artifact_root / "reference-execution.json",
    )
    if execution_payload is not None:
        outputs = execution_payload.get("outputs")
        if isinstance(outputs, dict):
            for path_string in outputs.values():
                if isinstance(path_string, str):
                    source = Path(path_string)
                    _copy_if_exists(source, artifact_root / f"reference-{source.name}")
    if execution_payload is not None:
        _write_json(
            artifact_root / "reference-execution.observed.json", execution_payload
        )
    if reference_summary is not None:
        _write_json(
            artifact_root / "reference-summary.observed.json", reference_summary
        )
    if reference_rows is not None:
        _write_json(artifact_root / "reference-rows.observed.json", reference_rows)
        _write_rows_table(artifact_root / "reference-rows.observed.tsv", reference_rows)
    if bijux_summary is not None:
        _write_json(artifact_root / "bijux-summary.json", bijux_summary)
    if reference_error is not None:
        _write_json(artifact_root / "reference-error.observed.json", reference_error)
    if bijux_error is not None:
        _write_json(artifact_root / "bijux-error.json", bijux_error)
    if bijux_rows is not None:
        _write_json(artifact_root / "bijux-rows.json", bijux_rows)
        _write_rows_table(artifact_root / "bijux-rows.tsv", bijux_rows)
    if bijux_normalized_text is not None:
        (artifact_root / "bijux-normalized.txt").write_text(
            f"{bijux_normalized_text}\n",
            encoding="utf-8",
        )
    _write_json(
        artifact_root / "comparison.json",
        {
            "case_id": case.case_id,
            "function_name": case.function_name,
            "mismatch_reason": mismatch_reason,
        },
    )
    return artifact_root


def _build_bijux_case_payload(
    case: ApeParityCase,
) -> tuple[dict[str, object], list[dict[str, object]] | None, str | None]:
    if case.operation in {"read-tree-structure", "write-tree-structure"}:
        summary, rows, normalized_text = _build_bijux_tree_structure(case.input_fixture)
        return summary, rows, normalized_text
    if case.operation == "root-tree-outgroup":
        summary, rows, normalized_text = _build_bijux_root_outgroup_structure(
            case.input_fixture,
            outgroup_taxa=case.outgroup_taxa,
        )
        return summary, rows, normalized_text
    if case.operation == "unroot-tree":
        summary, rows, normalized_text = _build_bijux_unroot_structure(
            case.input_fixture,
        )
        return summary, rows, normalized_text
    if case.operation == "drop-tree-taxa":
        summary, rows, normalized_text = _build_bijux_drop_tip_structure(
            case.input_fixture,
            excluded_taxa=case.excluded_taxa,
        )
        return summary, rows, normalized_text
    if case.operation == "keep-tree-taxa":
        summary, rows, normalized_text = _build_bijux_keep_tip_structure(
            case.input_fixture,
            requested_taxa=case.requested_taxa,
        )
        return summary, rows, normalized_text
    if case.operation == "extract-tree-clade":
        if case.node_id is None:
            raise ValueError(
                f"ape parity case '{case.case_id}' is missing an extraction node id"
            )
        summary, rows, normalized_text = _build_bijux_extract_clade_structure(
            case.input_fixture,
            node_id=case.node_id,
        )
        return summary, rows, normalized_text
    if case.operation == "get-tree-mrca":
        if not case.mrca_taxa:
            raise ValueError(f"ape parity case '{case.case_id}' is missing MRCA taxa")
        summary = _build_bijux_mrca_summary(
            case.input_fixture,
            mrca_taxa=case.mrca_taxa,
        )
        return summary, None, None
    if case.operation == "assess-tree-monophyly":
        if case.monophyly_reroot is None:
            raise ValueError(
                f"ape parity case '{case.case_id}' is missing a monophyly reroot policy"
            )
        summary = _build_bijux_monophyly_summary(
            case.input_fixture,
            requested_taxa=case.requested_taxa,
            reroot=case.monophyly_reroot,
        )
        return summary, None, None
    if case.operation in {"read-tree-set-structure", "write-tree-set-structure"}:
        summary, rows, normalized_text = _build_bijux_tree_set_structure(
            case.input_fixture
        )
        return summary, rows, normalized_text
    if case.operation == "tree-consensus":
        if case.consensus_method is None:
            raise ValueError(
                f"ape parity case '{case.case_id}' is missing a consensus method"
            )
        summary, rows, normalized_text = _build_bijux_consensus_rows(
            case.input_fixture,
            consensus_method=case.consensus_method,
        )
        return summary, rows, normalized_text
    if case.operation == "tree-clade-support":
        if case.reference_tree_path is None:
            raise ValueError(
                f"ape parity case '{case.case_id}' is missing a reference tree path"
            )
        summary, rows = _build_bijux_prop_clades_rows(
            case.reference_tree_path,
            case.input_fixture,
        )
        return summary, rows, None
    if case.operation == "tree-tip-distance":
        summary, rows = _build_bijux_tree_tip_distance_rows(case.input_fixture)
        return summary, rows, None
    if case.operation == "distance-matrix-neighbor-joining":
        return _build_bijux_neighbor_joining_structure(case.input_fixture)
    if case.operation == "tree-topology-distance":
        if case.rf_mode is None:
            raise ValueError(
                f"ape parity case '{case.case_id}' is missing a topology rf mode"
            )
        summary, rows = _build_bijux_topology_distance_rows(
            case.input_fixture,
            rf_mode=case.rf_mode,
        )
        return summary, rows, None
    if case.operation == "tree-brownian-covariance":
        summary, rows = _build_bijux_brownian_covariance_rows(case.input_fixture)
        return summary, rows, None
    if case.operation == "tree-continuous-ancestral-states":
        if (
            case.trait_table_path is None
            or case.trait_name is None
            or case.trait_taxon_column is None
        ):
            raise ValueError(
                f"ape parity case '{case.case_id}' is missing a trait table path, trait name, or taxon column"
            )
        summary, rows = _build_bijux_continuous_ancestral_rows(
            case.input_fixture,
            trait_table_path=case.trait_table_path,
            trait_name=case.trait_name,
            trait_taxon_column=case.trait_taxon_column,
        )
        return summary, rows, None
    if case.operation == "tree-discrete-ancestral-states":
        if (
            case.trait_table_path is None
            or case.trait_name is None
            or case.trait_taxon_column is None
        ):
            raise ValueError(
                f"ape parity case '{case.case_id}' is missing a trait table path, trait name, or taxon column"
            )
        summary, rows = _build_bijux_discrete_ancestral_rows(
            case.input_fixture,
            trait_table_path=case.trait_table_path,
            trait_name=case.trait_name,
            trait_taxon_column=case.trait_taxon_column,
            ancestral_model=case.ancestral_model or "equal-rates",
        )
        return summary, rows, None
    if case.operation == "tree-independent-contrasts":
        if case.trait_table_path is None or case.trait_name is None:
            raise ValueError(
                f"ape parity case '{case.case_id}' is missing a trait table path or trait name"
            )
        summary, rows = _build_bijux_independent_contrast_rows(
            case.input_fixture,
            trait_table_path=case.trait_table_path,
            trait_name=case.trait_name,
        )
        return summary, rows, None
    if case.operation == "tree-node-depth":
        summary, rows = _build_bijux_tree_node_depth_rows(case.input_fixture)
        return summary, rows, None
    if case.operation == "tree-branching-times":
        summary, rows = _build_bijux_branching_time_rows(case.input_fixture)
        return summary, rows, None
    if case.operation == "tree-diversification-gamma-statistic":
        summary, rows = _build_bijux_diversification_gamma_rows(case.input_fixture)
        return summary, rows, None
    if case.operation == "tree-simulation-envelope":
        summary, rows = _build_bijux_tree_simulation_envelope_rows(case.fixture_id)
        return summary, rows, None
    if case.operation == "tree-ultrametricity":
        if case.ultrametric_option is None:
            raise ValueError(
                f"ape parity case '{case.case_id}' is missing an ultrametric option"
            )
        summary, rows = _build_bijux_tree_ultrametric_rows(
            case.input_fixture,
            tolerance=case.tolerance,
            option=case.ultrametric_option,
        )
        return summary, rows, None
    if case.operation == "dna-dnabin-structure":
        summary, rows = _build_bijux_dnabin_rows(case.input_fixture)
        return summary, rows, None
    if case.operation == "dna-base-frequency":
        summary, rows = _build_bijux_base_frequency_summary(case.input_fixture)
        return summary, rows, None
    if case.operation == "dna-segregating-sites":
        summary, rows = _build_bijux_segregating_site_rows(case.input_fixture)
        return summary, rows, None
    if case.operation == "dna-distance":
        if case.pairwise_deletion is None:
            raise ValueError(
                f"ape parity case '{case.case_id}' is missing pairwise deletion policy"
            )
        if case.distance_model is None:
            raise ValueError(
                f"ape parity case '{case.case_id}' is missing a distance model"
            )
        summary, rows = _build_bijux_distance_rows(
            case.input_fixture,
            pairwise_deletion=case.pairwise_deletion,
            distance_model=case.distance_model,
        )
        return summary, rows, None
    if case.operation == "dna-translation":
        if case.genetic_code_id is None:
            raise ValueError(
                f"ape parity case '{case.case_id}' is missing a genetic code id"
            )
        summary, rows = _build_bijux_translation_rows(
            case.input_fixture,
            genetic_code_id=case.genetic_code_id,
        )
        return summary, rows, None
    raise ValueError(f"unsupported ape parity operation '{case.operation}'")


def _load_reference_case_payload(
    case: ApeParityCase,
    execution_root: Path,
) -> tuple[dict[str, object], list[dict[str, object]] | None, str | None]:
    if case.operation in {
        "read-tree-structure",
        "write-tree-structure",
        "root-tree-outgroup",
        "unroot-tree",
        "drop-tree-taxa",
        "keep-tree-taxa",
        "extract-tree-clade",
    }:
        summary = _normalize_reference_summary(
            _load_json(execution_root / "summary.json")
        )
        expected_tip_labels = {str(label) for label in summary.get("tip_labels", [])}
        rows = _load_rows_table(
            execution_root / "clades.tsv",
            expected_tip_labels=expected_tip_labels,
            sort_rows=True,
        )
        normalized_text = _canonical_newick(
            execution_root / "normalized-tree.nwk",
            expected_tip_labels=expected_tip_labels,
        )
        return summary, rows, normalized_text
    if case.operation == "get-tree-mrca":
        summary = _normalize_reference_summary(
            _load_json(execution_root / "summary.json")
        )
        return summary, None, None
    if case.operation == "assess-tree-monophyly":
        summary = _normalize_reference_summary(
            _load_json(execution_root / "summary.json")
        )
        return summary, None, None
    if case.operation in {"read-tree-set-structure", "write-tree-set-structure"}:
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "clades.tsv", sort_rows=True)
        return summary, rows, None
    if case.operation == "tree-consensus":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "clade-frequencies.tsv")
        normalized_text = _canonical_newick(execution_root / "normalized-tree.nwk")
        return summary, rows, normalized_text
    if case.operation == "tree-clade-support":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "support-table.tsv")
        return summary, rows, None
    if case.operation == "tree-tip-distance":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "tip-distance-long.tsv")
        return summary, rows, None
    if case.operation == "distance-matrix-neighbor-joining":
        summary = _normalize_reference_summary(
            _load_json(execution_root / "summary.json")
        )
        normalized_text = _canonical_newick(execution_root / "normalized-tree.nwk")
        return summary, None, normalized_text
    if case.operation == "tree-topology-distance":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "split-table.tsv")
        return summary, rows, None
    if case.operation == "tree-brownian-covariance":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "covariance-long.tsv")
        return summary, rows, None
    if case.operation == "tree-continuous-ancestral-states":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "continuous-ancestral.tsv")
        return summary, rows, None
    if case.operation == "tree-discrete-ancestral-states":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "discrete-ancestral.tsv")
        return summary, rows, None
    if case.operation == "tree-independent-contrasts":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "independent-contrasts.tsv")
        return summary, rows, None
    if case.operation == "tree-node-depth":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "node-depths.tsv")
        return summary, rows, None
    if case.operation == "tree-branching-times":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "branching-times.tsv")
        return summary, rows, None
    if case.operation == "tree-diversification-gamma-statistic":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "gamma-statistic.tsv")
        return summary, rows, None
    if case.operation == "tree-simulation-envelope":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "simulation-envelope.tsv")
        return summary, rows, None
    if case.operation == "tree-ultrametricity":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "ultrametric-diagnostics.tsv")
        return summary, rows, None
    if case.operation == "dna-dnabin-structure":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "dnabin.tsv")
        return summary, rows, None
    if case.operation == "dna-base-frequency":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "base-frequency.tsv")
        return summary, rows, None
    if case.operation == "dna-segregating-sites":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "segregating-sites.tsv")
        return summary, rows, None
    if case.operation == "dna-distance":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "distance-matrix.tsv")
        return summary, rows, None
    if case.operation == "dna-translation":
        summary = _load_json(execution_root / "summary.json")
        rows = _load_rows_table(execution_root / "translation.tsv")
        return summary, rows, None
    raise ValueError(f"unsupported ape parity operation '{case.operation}'")


def _tree_structure_mismatch_reason(
    case: ApeParityCase,
    execution_root: Path,
) -> str | None:
    reference_tree = load_tree(execution_root / "normalized-tree.nwk")
    expected_tree = load_tree(case.input_fixture)
    _normalize_tree_labels(
        reference_tree.root,
        expected_tip_labels=set(expected_tree.tip_names),
    )
    report = compare_tree_structurally(
        expected_tree,
        reference_tree,
        tolerance=case.tolerance,
        compare_internal_labels=True,
    )
    return report.mismatch_reason


def _tree_set_structure_mismatch_reason(
    case: ApeParityCase,
    execution_root: Path,
) -> str | None:
    reference_tree_set = load_newick_tree_set(
        execution_root / "normalized-tree-set.nwk"
    )
    expected_tree_set = load_newick_tree_set(case.input_fixture)
    expected_tip_labels = {
        tip_name for tree in expected_tree_set for tip_name in tree.tip_names
    }
    for tree in reference_tree_set:
        _normalize_tree_labels(tree.root, expected_tip_labels=expected_tip_labels)
    report = compare_tree_sets_structurally(
        expected_tree_set,
        reference_tree_set,
        tolerance=case.tolerance,
        compare_internal_labels=True,
    )
    return report.mismatch_reason


def _root_tree_outgroup_mismatch_reason(
    case: ApeParityCase,
    execution_root: Path,
) -> str | None:
    reference_tree = load_tree(execution_root / "normalized-tree.nwk")
    # Canonical Newick does not preserve rootedness metadata, but ape::root
    # produced this record explicitly as a rooted output for these governed cases.
    reference_tree.rooted = True
    expected_tree, _report = root_tree_on_outgroup(
        case.input_fixture,
        outgroup_taxa=list(case.outgroup_taxa),
    )
    _normalize_tree_labels(
        reference_tree.root,
        expected_tip_labels=set(expected_tree.tip_names),
    )
    report = compare_tree_structurally(
        expected_tree,
        reference_tree,
        tolerance=case.tolerance,
        compare_internal_labels=True,
    )
    return report.mismatch_reason


def _unroot_tree_mismatch_reason(
    case: ApeParityCase,
    execution_root: Path,
) -> str | None:
    reference_tree = load_tree(execution_root / "normalized-tree.nwk")
    reference_tree.rooted = False
    expected_tree, _report = unroot_tree(case.input_fixture)
    _normalize_tree_labels(
        reference_tree.root,
        expected_tip_labels=set(expected_tree.tip_names),
    )
    report = compare_tree_structurally(
        expected_tree,
        reference_tree,
        tolerance=case.tolerance,
        compare_internal_labels=True,
    )
    return report.mismatch_reason


def _drop_tip_tree_mismatch_reason(
    case: ApeParityCase,
    execution_root: Path,
) -> str | None:
    reference_summary = _normalize_reference_summary(
        _load_json(execution_root / "summary.json")
    )
    reference_tree = load_tree(execution_root / "normalized-tree.nwk")
    reference_tree.rooted = _summary_rooted_flag(reference_summary)
    expected_tree, _report = drop_tree_taxa(
        case.input_fixture, list(case.excluded_taxa)
    )
    _normalize_tree_labels(
        reference_tree.root,
        expected_tip_labels=set(expected_tree.tip_names),
    )
    report = compare_tree_structurally(
        expected_tree,
        reference_tree,
        tolerance=case.tolerance,
        compare_internal_labels=True,
    )
    return report.mismatch_reason


def _keep_tip_tree_mismatch_reason(
    case: ApeParityCase,
    execution_root: Path,
) -> str | None:
    reference_summary = _normalize_reference_summary(
        _load_json(execution_root / "summary.json")
    )
    reference_tree = load_tree(execution_root / "normalized-tree.nwk")
    reference_tree.rooted = _summary_rooted_flag(reference_summary)
    expected_tree, _report = prune_tree_to_requested_taxa(
        case.input_fixture,
        list(case.requested_taxa),
    )
    _normalize_tree_labels(
        reference_tree.root,
        expected_tip_labels=set(expected_tree.tip_names),
    )
    report = compare_tree_structurally(
        expected_tree,
        reference_tree,
        tolerance=case.tolerance,
        compare_internal_labels=True,
    )
    return report.mismatch_reason


def _extract_clade_tree_mismatch_reason(
    case: ApeParityCase,
    execution_root: Path,
) -> str | None:
    reference_summary = _normalize_reference_summary(
        _load_json(execution_root / "summary.json")
    )
    reference_tree = load_tree(execution_root / "normalized-tree.nwk")
    reference_tree.rooted = _summary_rooted_flag(reference_summary)
    if case.node_id is None:
        raise ValueError(
            f"ape parity case '{case.case_id}' is missing an extraction node id"
        )
    expected_tree, _report = extract_tree_clade_by_node_id(
        case.input_fixture,
        node_id=case.node_id,
    )
    _normalize_tree_labels(
        reference_tree.root,
        expected_tip_labels=set(expected_tree.tip_names),
    )
    report = compare_tree_structurally(
        expected_tree,
        reference_tree,
        tolerance=case.tolerance,
        compare_internal_labels=True,
    )
    return report.mismatch_reason


def _consensus_tree_mismatch_reason(
    case: ApeParityCase,
    execution_root: Path,
) -> str | None:
    if case.consensus_method == "strict":
        expected_tree, _report = compute_strict_consensus_tree(case.input_fixture)
    elif case.consensus_method == "majority-rule":
        expected_tree, _report = compute_consensus_tree(case.input_fixture)
    else:
        raise ValueError(
            f"ape parity case '{case.case_id}' has unsupported consensus method "
            f"{case.consensus_method!r}"
        )
    expected_tree.rooted = False
    _clear_branch_lengths(expected_tree.root)
    reference_tree = load_tree(execution_root / "normalized-tree.nwk")
    reference_tree.rooted = False
    _clear_branch_lengths(reference_tree.root)
    _normalize_tree_labels(
        reference_tree.root,
        expected_tip_labels=set(expected_tree.tip_names),
    )
    report = compare_tree_structurally(
        expected_tree,
        reference_tree,
        tolerance=case.tolerance,
        compare_internal_labels=False,
    )
    return report.mismatch_reason


def _neighbor_joining_tree_mismatch_reason(
    case: ApeParityCase,
    execution_root: Path,
) -> str | None:
    reference_tree = load_tree(execution_root / "normalized-tree.nwk")
    reference_tree.rooted = False
    expected_tree, _report = build_tree_from_imported_distance_matrix(
        case.input_fixture,
        method="neighbor-joining",
    )
    _normalize_tree_labels(
        reference_tree.root,
        expected_tip_labels=set(expected_tree.tip_names),
    )
    report = compare_tree_structurally(
        expected_tree,
        reference_tree,
        tolerance=case.tolerance,
        compare_internal_labels=False,
    )
    return report.mismatch_reason


def _supports_ard_rate_multiset_equivalence(
    *,
    reference_summary: dict[str, object] | None,
    bijux_summary: dict[str, object] | None,
) -> bool:
    if reference_summary is None or bijux_summary is None:
        return False
    return (
        reference_summary.get("model") == "all-rates-different"
        and bijux_summary.get("model") == "all-rates-different"
        and reference_summary.get("overparameterized") is True
        and bijux_summary.get("overparameterized") is True
    )


def _group_transition_rate_rows(
    rows: list[dict[str, object]],
) -> dict[tuple[bool, int], list[float]]:
    grouped: dict[tuple[bool, int], list[float]] = {}
    for row in rows:
        grouped.setdefault(
            (
                bool(row["transition_allowed"]),
                int(row["step_distance"]),
            ),
            [],
        ).append(float(row["rate"]))
    for values in grouped.values():
        values.sort()
    return grouped


def _transition_rate_rows_match(
    *,
    reference_rows: list[dict[str, object]],
    bijux_rows: list[dict[str, object]],
    reference_summary: dict[str, object] | None,
    bijux_summary: dict[str, object] | None,
    tolerance: float,
) -> bool:
    if _compare_json(reference_rows, bijux_rows, tolerance=tolerance):
        return True
    if not _supports_ard_rate_multiset_equivalence(
        reference_summary=reference_summary,
        bijux_summary=bijux_summary,
    ):
        return False
    return _compare_json(
        _group_transition_rate_rows(reference_rows),
        _group_transition_rate_rows(bijux_rows),
        tolerance=tolerance,
    )


def run_ape_parity_cases(
    *,
    case_ids: list[str] | None = None,
    rscript_executable: str = "Rscript",
    failure_root: Path | None = None,
    fixtures_root: Path | None = None,
) -> ApeParityReport:
    """Run governed live `ape` parity cases through the checked-in R runner."""
    selected = _selected_cases(case_ids=case_ids, fixtures_root=fixtures_root)
    observations: list[ApeParityObservation] = []
    active_failure_root = _failure_root() if failure_root is None else failure_root
    bijux_version = _bijux_version()
    bijux_commit = _bijux_commit()
    for case in selected:
        with tempfile.TemporaryDirectory(
            prefix=f"bijux-ape-parity-{case.case_id}-"
        ) as tmpdir:
            working_root = Path(tmpdir)
            reference_input_path = _materialize_reference_input(case, working_root)
            reference_case = replace(case, input_fixture=reference_input_path)
            case_file = _write_case_file(working_root / "case.json", reference_case)
            execution_root = working_root / "reference"
            execution_root.mkdir(parents=True, exist_ok=True)
            bijux_summary: dict[str, object] | None = None
            bijux_rows: list[dict[str, object]] | None = None
            bijux_normalized_text: str | None = None
            bijux_error: dict[str, object] | None = None
            try:
                (
                    bijux_summary,
                    bijux_rows,
                    bijux_normalized_text,
                ) = _build_bijux_case_payload(case)
            except Exception as error:
                bijux_error = {
                    "error_type": type(error).__name__,
                    "message": str(error),
                }
            execution_payload: dict[str, object] | None = None
            reference_summary: dict[str, object] | None = None
            reference_error: dict[str, object] | None = None
            reference_rows: list[dict[str, object]] | None = None
            reference_normalized_text: str | None = None
            status = "failed"
            mismatch_reason: str | None = None
            artifact_root: Path | None = None
            r_version: str | None = None
            ape_version: str | None = None
            process_stdout = ""
            process_stderr = ""
            try:
                # Repository-owned R parity runner.
                process = subprocess.run(  # nosec
                    [
                        rscript_executable,
                        str(_ape_runner_path()),
                        str(case_file),
                        str(execution_root),
                    ],
                    capture_output=True,
                    check=False,
                    cwd=_repository_root(),
                    env=_reference_environment(),
                    text=True,
                )
                process_stdout = process.stdout
                process_stderr = process.stderr
            except FileNotFoundError:
                process = None
                status = "skipped"
                mismatch_reason = "rscript_unavailable"
            if process is None:
                pass
            elif process.returncode != 0:
                mismatch_reason = "reference_execution_failed"
            else:
                execution_path = execution_root / "reference-execution.json"
                if not execution_path.exists():
                    mismatch_reason = "reference_execution_failed"
                else:
                    execution_payload = _load_json(execution_path)
                    r_version = _optional_payload_string(execution_payload, "r_version")
                    ape_version = _optional_payload_string(
                        execution_payload, "ape_version"
                    )
                    execution_status = execution_payload.get("status")
                    if execution_status == "unavailable":
                        status = "skipped"
                        mismatch_reason = str(
                            execution_payload.get(
                                "mismatch_reason", "ape_package_unavailable"
                            )
                        )
                    elif execution_status != "ok":
                        reference_error = {
                            "error_type": str(
                                execution_payload.get(
                                    "error_type",
                                    execution_payload.get(
                                        "mismatch_reason", "reference_execution_failed"
                                    ),
                                )
                            ),
                            "message": str(execution_payload.get("message", "")),
                        }
                        mismatch_reason = str(
                            execution_payload.get(
                                "mismatch_reason", "reference_execution_failed"
                            )
                        )
                    else:
                        (
                            reference_summary,
                            reference_rows,
                            reference_normalized_text,
                        ) = _load_reference_case_payload(case, execution_root)
                        if case.operation in {
                            "read-tree-structure",
                            "write-tree-structure",
                        }:
                            mismatch_reason = _tree_structure_mismatch_reason(
                                case,
                                execution_root,
                            )
                        elif case.operation == "root-tree-outgroup":
                            mismatch_reason = _root_tree_outgroup_mismatch_reason(
                                case,
                                execution_root,
                            )
                        elif case.operation == "unroot-tree":
                            mismatch_reason = _unroot_tree_mismatch_reason(
                                case,
                                execution_root,
                            )
                        elif case.operation == "drop-tree-taxa":
                            mismatch_reason = _drop_tip_tree_mismatch_reason(
                                case,
                                execution_root,
                            )
                            if mismatch_reason is None and not _compare_json(
                                reference_summary,
                                bijux_summary,
                                tolerance=case.tolerance,
                            ):
                                mismatch_reason = "summary_mismatch"
                        elif case.operation == "keep-tree-taxa":
                            mismatch_reason = _keep_tip_tree_mismatch_reason(
                                case,
                                execution_root,
                            )
                            if mismatch_reason is None and not _compare_json(
                                reference_summary,
                                bijux_summary,
                                tolerance=case.tolerance,
                            ):
                                mismatch_reason = "summary_mismatch"
                        elif case.operation == "extract-tree-clade":
                            mismatch_reason = _extract_clade_tree_mismatch_reason(
                                case,
                                execution_root,
                            )
                            if mismatch_reason is None and not _compare_json(
                                reference_summary,
                                bijux_summary,
                                tolerance=case.tolerance,
                            ):
                                mismatch_reason = "summary_mismatch"
                        elif (
                            case.operation == "get-tree-mrca"
                            or case.operation == "assess-tree-monophyly"
                        ):
                            if not _compare_json(
                                reference_summary,
                                bijux_summary,
                                tolerance=case.tolerance,
                            ):
                                mismatch_reason = "summary_mismatch"
                            else:
                                status = "passed"
                        elif case.operation in {
                            "read-tree-set-structure",
                            "write-tree-set-structure",
                        }:
                            mismatch_reason = _tree_set_structure_mismatch_reason(
                                case,
                                execution_root,
                            )
                        elif case.operation == "tree-consensus":
                            mismatch_reason = _consensus_tree_mismatch_reason(
                                case,
                                execution_root,
                            )
                            if mismatch_reason is None and not _compare_json(
                                reference_summary,
                                bijux_summary,
                                tolerance=case.tolerance,
                            ):
                                mismatch_reason = "summary_mismatch"
                            elif mismatch_reason is None and not _compare_json(
                                reference_rows,
                                bijux_rows,
                                tolerance=case.tolerance,
                            ):
                                mismatch_reason = "rows_mismatch"
                        elif case.operation == "distance-matrix-neighbor-joining":
                            mismatch_reason = _neighbor_joining_tree_mismatch_reason(
                                case,
                                execution_root,
                            )
                            if mismatch_reason is None and not _compare_json(
                                reference_summary,
                                bijux_summary,
                                tolerance=case.tolerance,
                            ):
                                mismatch_reason = "summary_mismatch"
                        elif case.operation == "tree-discrete-ancestral-states":
                            reference_transition_rows = (
                                []
                                if reference_summary is None
                                else reference_summary.get("transition_rate_rows", [])
                            )
                            bijux_transition_rows = (
                                []
                                if bijux_summary is None
                                else bijux_summary.get("transition_rate_rows", [])
                            )
                            reference_summary_without_transition_rows = (
                                {}
                                if reference_summary is None
                                else {
                                    key: value
                                    for key, value in reference_summary.items()
                                    if key != "transition_rate_rows"
                                }
                            )
                            bijux_summary_without_transition_rows = (
                                {}
                                if bijux_summary is None
                                else {
                                    key: value
                                    for key, value in bijux_summary.items()
                                    if key != "transition_rate_rows"
                                }
                            )
                            if not _compare_json(
                                reference_summary_without_transition_rows,
                                bijux_summary_without_transition_rows,
                                tolerance=case.tolerance,
                            ):
                                mismatch_reason = "summary_mismatch"
                            elif not _transition_rate_rows_match(
                                reference_rows=reference_transition_rows,
                                bijux_rows=bijux_transition_rows,
                                reference_summary=reference_summary,
                                bijux_summary=bijux_summary,
                                tolerance=(
                                    case.transition_rate_tolerance
                                    if case.transition_rate_tolerance is not None
                                    else case.tolerance
                                ),
                            ):
                                mismatch_reason = "transition_rate_rows_mismatch"
                        elif not _compare_json(
                            reference_summary, bijux_summary, tolerance=case.tolerance
                        ):
                            mismatch_reason = "summary_mismatch"
                        elif not _compare_json(
                            reference_rows,
                            bijux_rows,
                            tolerance=case.tolerance,
                        ):
                            mismatch_reason = "rows_mismatch"
                        elif reference_normalized_text != bijux_normalized_text:
                            mismatch_reason = "normalized_text_mismatch"
                        else:
                            status = "passed"
                        if mismatch_reason is None:
                            status = "passed"
            if case.expected_status == "parse-error":
                if bijux_error is None:
                    mismatch_reason = "bijux_expected_parse_error_missing"
                elif reference_error is None:
                    mismatch_reason = "reference_expected_parse_error_missing"
                elif not bijux_error.get("message") or not reference_error.get(
                    "message"
                ):
                    mismatch_reason = "parse_error_message_missing"
                else:
                    status = "passed"
                    mismatch_reason = None
            if case.expected_status == "rooting-error":
                if bijux_error is None:
                    mismatch_reason = "bijux_expected_rooting_error_missing"
                elif reference_error is None:
                    mismatch_reason = "reference_expected_rooting_error_missing"
                elif not bijux_error.get("message") or not reference_error.get(
                    "message"
                ):
                    mismatch_reason = "rooting_error_message_missing"
                else:
                    status = "passed"
                    mismatch_reason = None
            if case.expected_status == "clade-extraction-error":
                if bijux_error is None:
                    mismatch_reason = "bijux_expected_clade_extraction_error_missing"
                elif reference_error is None:
                    mismatch_reason = (
                        "reference_expected_clade_extraction_error_missing"
                    )
                elif not bijux_error.get("message") or not reference_error.get(
                    "message"
                ):
                    mismatch_reason = "clade_extraction_error_message_missing"
                else:
                    status = "passed"
                    mismatch_reason = None
            if case.expected_status == "mrca-error":
                if bijux_error is None:
                    mismatch_reason = "bijux_expected_mrca_error_missing"
                elif reference_error is None:
                    mismatch_reason = "reference_expected_mrca_error_missing"
                elif not bijux_error.get("message") or not reference_error.get(
                    "message"
                ):
                    mismatch_reason = "mrca_error_message_missing"
                else:
                    status = "passed"
                    mismatch_reason = None
            if case.expected_status == "monophyly-error":
                if bijux_error is None:
                    mismatch_reason = "bijux_expected_monophyly_error_missing"
                elif reference_error is None:
                    mismatch_reason = "reference_expected_monophyly_error_missing"
                elif not bijux_error.get("message") or not reference_error.get(
                    "message"
                ):
                    mismatch_reason = "monophyly_error_message_missing"
                else:
                    status = "passed"
                    mismatch_reason = None
            if case.expected_status == "consensus-error":
                if bijux_error is None:
                    mismatch_reason = "bijux_expected_consensus_error_missing"
                elif reference_error is None:
                    mismatch_reason = "reference_expected_consensus_error_missing"
                elif not bijux_error.get("message") or not reference_error.get(
                    "message"
                ):
                    mismatch_reason = "consensus_error_message_missing"
                else:
                    status = "passed"
                    mismatch_reason = None
            if case.expected_status == "prop-clades-error":
                if bijux_error is None:
                    mismatch_reason = "bijux_expected_prop_clades_error_missing"
                elif reference_error is None:
                    mismatch_reason = "reference_expected_prop_clades_error_missing"
                elif not bijux_error.get("message") or not reference_error.get(
                    "message"
                ):
                    mismatch_reason = "prop_clades_error_message_missing"
                else:
                    status = "passed"
                    mismatch_reason = None
            if case.expected_status == "dna-distance-error":
                if bijux_error is None:
                    mismatch_reason = "bijux_expected_dna_distance_error_missing"
                elif reference_error is None:
                    mismatch_reason = "reference_expected_dna_distance_error_missing"
                elif not bijux_error.get("message") or not reference_error.get(
                    "message"
                ):
                    mismatch_reason = "dna_distance_error_message_missing"
                else:
                    status = "passed"
                    mismatch_reason = None
            if status != "passed":
                artifact_root = _persist_failure_bundle(
                    failure_root=active_failure_root,
                    case=case,
                    case_file=case_file,
                    execution_root=execution_root,
                    execution_payload=execution_payload,
                    reference_summary=reference_summary,
                    bijux_summary=bijux_summary,
                    reference_error=reference_error,
                    bijux_error=bijux_error,
                    reference_rows=reference_rows,
                    bijux_rows=bijux_rows,
                    bijux_normalized_text=bijux_normalized_text,
                    mismatch_reason=mismatch_reason or "reference_execution_failed",
                )
                if process_stdout:
                    (artifact_root / "reference-stdout.txt").write_text(
                        process_stdout, encoding="utf-8"
                    )
                if process_stderr:
                    (artifact_root / "reference-stderr.txt").write_text(
                        process_stderr, encoding="utf-8"
                    )
            observations.append(
                ApeParityObservation(
                    case_id=case.case_id,
                    fixture_kind=case.fixture_kind,
                    fixture_id=case.fixture_id,
                    function_name=case.function_name,
                    python_function_name=case.python_function_name,
                    input_fixture=case.input_fixture,
                    tolerance=case.tolerance,
                    r_version=r_version,
                    ape_version=ape_version,
                    bijux_version=bijux_version,
                    bijux_commit=bijux_commit,
                    status=status,
                    passed=status == "passed",
                    mismatch_reason=mismatch_reason,
                    reproducible_artifact_root=artifact_root,
                    reference_summary=reference_summary,
                    bijux_summary=bijux_summary,
                    reference_error=reference_error,
                    bijux_error=bijux_error,
                )
            )
    return build_ape_parity_report(observations)
