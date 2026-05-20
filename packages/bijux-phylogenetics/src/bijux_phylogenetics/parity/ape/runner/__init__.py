from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import subprocess  # nosec B404 - parity helpers invoke repository-owned reference commands
import tempfile

from .normalization import (
    _compare_json,
    _load_json,
    _load_rows_table,
    _normalize_reference_summary,
    _optional_payload_string,
)
from .mismatch_policy import (
    _consensus_tree_mismatch_reason,
    _drop_tip_tree_mismatch_reason,
    _extract_clade_tree_mismatch_reason,
    _keep_tip_tree_mismatch_reason,
    _neighbor_joining_tree_mismatch_reason,
    _root_tree_outgroup_mismatch_reason,
    _transition_rate_rows_match,
    _tree_set_structure_mismatch_reason,
    _tree_structure_mismatch_reason,
    _unroot_tree_mismatch_reason,
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
from .reference_payloads import _load_reference_case_payload
from .ancestral_payloads import (
    _build_bijux_continuous_ancestral_rows,
    _build_bijux_discrete_ancestral_rows,
)
from .failure_artifacts import _persist_failure_bundle
from .comparative_payloads import (
    _build_bijux_brownian_covariance_rows,
    _build_bijux_branching_time_rows,
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
)
from .runtime import (
    ape_runner_path as _ape_runner_path,
    bijux_commit as _bijux_commit,
    bijux_version as _bijux_version,
    failure_root as _failure_root,
    reference_environment as _reference_environment,
    repository_root as _repository_root,
)
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
