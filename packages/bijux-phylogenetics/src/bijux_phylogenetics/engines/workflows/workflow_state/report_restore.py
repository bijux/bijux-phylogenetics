from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.phylo.alignment import (
    CodingSequenceExclusion,
)
from bijux_phylogenetics.runtime.errors import EngineWorkflowError

from ...artifacts.iqtree import IqtreeModelCandidate, IqtreeModelSelectionSummary
from ...artifacts.support import (
    BootstrapSupportNode,
    BootstrapSupportSummaryReport,
    FastTreeSupportNode,
    FastTreeSupportSummaryReport,
    ShAlrtSupportNode,
    ShAlrtSupportSummaryReport,
    WeakBackboneReport,
)
from ...common import EngineRunReport, EngineVersionInfo
from ..models import (
    AlignmentTrimmingSummary,
    CodonAwareAlignmentWorkflowReport,
    EngineWorkflowReport,
    IqtreeSupportValue,
    IqtreeWorkflowSummary,
)


def _restore_workflow_report(payload: dict[str, object]) -> EngineWorkflowReport:
    run_payload = payload["run"]
    if not isinstance(run_payload, dict):
        raise EngineWorkflowError("engine manifest contains an invalid run payload")
    version_payload = run_payload["version"]
    if not isinstance(version_payload, dict):
        raise EngineWorkflowError("engine manifest contains an invalid version payload")
    version = EngineVersionInfo(
        engine_name=str(version_payload["engine_name"]),
        executable=str(version_payload["executable"]),
        command=[str(item) for item in version_payload["command"]],
        text=str(version_payload["text"]),
    )
    run = EngineRunReport(
        engine_name=str(run_payload["engine_name"]),
        workflow=str(run_payload["workflow"]),
        executable=str(run_payload["executable"]),
        working_directory=Path(run_payload["working_directory"]),
        version=version,
        command=[str(item) for item in run_payload["command"]],
        exit_code=int(run_payload["exit_code"]),
        stdout_path=Path(run_payload["stdout_path"]),
        stderr_path=Path(run_payload["stderr_path"]),
        output_paths={
            str(key): Path(value)
            for key, value in dict(run_payload["output_paths"]).items()
        },
        warning_lines=[str(item) for item in run_payload["warning_lines"]],
        started_at_utc=str(run_payload.get("started_at_utc", "")),
        ended_at_utc=str(run_payload.get("ended_at_utc", "")),
        runtime_seconds=float(run_payload.get("runtime_seconds", 0.0)),
        timeout_seconds=(
            None
            if run_payload.get("timeout_seconds") is None
            else float(run_payload["timeout_seconds"])
        ),
        timed_out=bool(run_payload.get("timed_out", False)),
    )
    return EngineWorkflowReport(
        workflow=str(payload["workflow"]),
        engine_name=str(payload["engine_name"]),
        input_paths=[Path(item) for item in payload["input_paths"]],
        output_paths={
            str(key): Path(value)
            for key, value in dict(payload["output_paths"]).items()
        },
        run=run,
        manifest_path=Path(payload["manifest_path"]),
        input_checksums={
            str(key): str(value)
            for key, value in dict(payload.get("input_checksums", {})).items()
        },
        output_checksums={
            str(key): str(value)
            for key, value in dict(payload.get("output_checksums", {})).items()
        },
        config={
            str(key): value for key, value in dict(payload.get("config", {})).items()
        },
        selected_model=None
        if payload.get("selected_model") is None
        else str(payload["selected_model"]),
        log_likelihood=(
            None
            if payload.get("log_likelihood") is None
            else float(payload["log_likelihood"])
        ),
        iqtree_summary=(
            None
            if payload.get("iqtree_summary") is None
            else IqtreeWorkflowSummary(
                selected_model=(
                    None
                    if dict(payload["iqtree_summary"]).get("selected_model") is None
                    else str(dict(payload["iqtree_summary"])["selected_model"])
                ),
                log_likelihood=(
                    None
                    if dict(payload["iqtree_summary"]).get("log_likelihood") is None
                    else float(dict(payload["iqtree_summary"])["log_likelihood"])
                ),
                support_value_count=int(
                    dict(payload["iqtree_summary"])["support_value_count"]
                ),
                minimum_support=(
                    None
                    if dict(payload["iqtree_summary"]).get("minimum_support") is None
                    else float(dict(payload["iqtree_summary"])["minimum_support"])
                ),
                maximum_support=(
                    None
                    if dict(payload["iqtree_summary"]).get("maximum_support") is None
                    else float(dict(payload["iqtree_summary"])["maximum_support"])
                ),
                support_values=[
                    IqtreeSupportValue(
                        node=str(dict(item)["node"]),
                        descendant_taxa=[
                            str(taxon) for taxon in list(dict(item)["descendant_taxa"])
                        ],
                        support=float(dict(item)["support"]),
                        support_fraction=float(dict(item)["support_fraction"]),
                        is_backbone=bool(dict(item)["is_backbone"]),
                    )
                    for item in list(dict(payload["iqtree_summary"])["support_values"])
                ],
            )
        ),
        model_selection_summary=(
            None
            if payload.get("model_selection_summary") is None
            else IqtreeModelSelectionSummary(
                selected_model=(
                    None
                    if dict(payload["model_selection_summary"]).get("selected_model")
                    is None
                    else str(dict(payload["model_selection_summary"])["selected_model"])
                ),
                selected_criterion=(
                    None
                    if dict(payload["model_selection_summary"]).get(
                        "selected_criterion"
                    )
                    is None
                    else str(
                        dict(payload["model_selection_summary"])["selected_criterion"]
                    )
                ),
                best_model_aic=(
                    None
                    if dict(payload["model_selection_summary"]).get("best_model_aic")
                    is None
                    else str(dict(payload["model_selection_summary"])["best_model_aic"])
                ),
                best_model_aicc=(
                    None
                    if dict(payload["model_selection_summary"]).get("best_model_aicc")
                    is None
                    else str(
                        dict(payload["model_selection_summary"])["best_model_aicc"]
                    )
                ),
                best_model_bic=(
                    None
                    if dict(payload["model_selection_summary"]).get("best_model_bic")
                    is None
                    else str(dict(payload["model_selection_summary"])["best_model_bic"])
                ),
                best_score_aic=(
                    None
                    if dict(payload["model_selection_summary"]).get("best_score_aic")
                    is None
                    else float(
                        dict(payload["model_selection_summary"])["best_score_aic"]
                    )
                ),
                best_score_aicc=(
                    None
                    if dict(payload["model_selection_summary"]).get("best_score_aicc")
                    is None
                    else float(
                        dict(payload["model_selection_summary"])["best_score_aicc"]
                    )
                ),
                best_score_bic=(
                    None
                    if dict(payload["model_selection_summary"]).get("best_score_bic")
                    is None
                    else float(
                        dict(payload["model_selection_summary"])["best_score_bic"]
                    )
                ),
                candidate_count=int(
                    dict(payload["model_selection_summary"])["candidate_count"]
                ),
                candidates=[
                    IqtreeModelCandidate(
                        rank=int(dict(item)["rank"]),
                        model=str(dict(item)["model"]),
                        log_likelihood=float(dict(item)["log_likelihood"]),
                        parameter_count=(
                            None
                            if dict(item).get("parameter_count") is None
                            else int(dict(item)["parameter_count"])
                        ),
                        aic=float(dict(item)["aic"]),
                        aicc=float(dict(item)["aicc"]),
                        bic=float(dict(item)["bic"]),
                    )
                    for item in list(
                        dict(payload["model_selection_summary"])["candidates"]
                    )
                ],
                bic_near_best_models=[
                    str(item)
                    for item in list(
                        dict(payload["model_selection_summary"])["bic_near_best_models"]
                    )
                ],
            )
        ),
        bootstrap_support_summary=(
            None
            if payload.get("bootstrap_support_summary") is None
            else BootstrapSupportSummaryReport(
                tree_path=Path(dict(payload["bootstrap_support_summary"])["tree_path"]),
                internal_node_count=int(
                    dict(payload["bootstrap_support_summary"])["internal_node_count"]
                ),
                supported_node_count=int(
                    dict(payload["bootstrap_support_summary"])["supported_node_count"]
                ),
                minimum_support=(
                    None
                    if dict(payload["bootstrap_support_summary"]).get("minimum_support")
                    is None
                    else float(
                        dict(payload["bootstrap_support_summary"])["minimum_support"]
                    )
                ),
                maximum_support=(
                    None
                    if dict(payload["bootstrap_support_summary"]).get("maximum_support")
                    is None
                    else float(
                        dict(payload["bootstrap_support_summary"])["maximum_support"]
                    )
                ),
                median_support=(
                    None
                    if dict(payload["bootstrap_support_summary"]).get("median_support")
                    is None
                    else float(
                        dict(payload["bootstrap_support_summary"])["median_support"]
                    )
                ),
                weakly_supported_clade_count=int(
                    dict(payload["bootstrap_support_summary"])[
                        "weakly_supported_clade_count"
                    ]
                ),
                support_histogram={
                    str(key): int(value)
                    for key, value in dict(
                        dict(payload["bootstrap_support_summary"])["support_histogram"]
                    ).items()
                },
                nodes=[
                    BootstrapSupportNode(
                        node=str(dict(item)["node"]),
                        descendant_taxa=[
                            str(taxon) for taxon in list(dict(item)["descendant_taxa"])
                        ],
                        support=float(dict(item)["support"]),
                        support_fraction=float(dict(item)["support_fraction"]),
                        is_backbone=bool(dict(item)["is_backbone"]),
                    )
                    for item in list(
                        dict(payload["bootstrap_support_summary"])["nodes"]
                    )
                ],
                warnings=[
                    str(item)
                    for item in list(
                        dict(payload["bootstrap_support_summary"])["warnings"]
                    )
                ],
            )
        ),
        fasttree_support_summary=(
            None
            if payload.get("fasttree_support_summary") is None
            else FastTreeSupportSummaryReport(
                tree_path=Path(dict(payload["fasttree_support_summary"])["tree_path"]),
                internal_node_count=int(
                    dict(payload["fasttree_support_summary"])["internal_node_count"]
                ),
                annotated_node_count=int(
                    dict(payload["fasttree_support_summary"])["annotated_node_count"]
                ),
                minimum_local_support=(
                    None
                    if dict(payload["fasttree_support_summary"]).get(
                        "minimum_local_support"
                    )
                    is None
                    else float(
                        dict(payload["fasttree_support_summary"])[
                            "minimum_local_support"
                        ]
                    )
                ),
                maximum_local_support=(
                    None
                    if dict(payload["fasttree_support_summary"]).get(
                        "maximum_local_support"
                    )
                    is None
                    else float(
                        dict(payload["fasttree_support_summary"])[
                            "maximum_local_support"
                        ]
                    )
                ),
                median_local_support=(
                    None
                    if dict(payload["fasttree_support_summary"]).get(
                        "median_local_support"
                    )
                    is None
                    else float(
                        dict(payload["fasttree_support_summary"])[
                            "median_local_support"
                        ]
                    )
                ),
                weakly_supported_clade_count=int(
                    dict(payload["fasttree_support_summary"])[
                        "weakly_supported_clade_count"
                    ]
                ),
                support_histogram={
                    str(key): int(value)
                    for key, value in dict(
                        dict(payload["fasttree_support_summary"])["support_histogram"]
                    ).items()
                },
                approximate_method=bool(
                    dict(payload["fasttree_support_summary"])["approximate_method"]
                ),
                support_label_kind=str(
                    dict(payload["fasttree_support_summary"])["support_label_kind"]
                ),
                support_scale=str(
                    dict(payload["fasttree_support_summary"])["support_scale"]
                ),
                nodes=[
                    FastTreeSupportNode(
                        node=str(dict(item)["node"]),
                        descendant_taxa=[
                            str(taxon) for taxon in list(dict(item)["descendant_taxa"])
                        ],
                        local_support=float(dict(item)["local_support"]),
                        support_fraction=float(dict(item)["support_fraction"]),
                        is_backbone=bool(dict(item)["is_backbone"]),
                    )
                    for item in list(dict(payload["fasttree_support_summary"])["nodes"])
                ],
                warnings=[
                    str(item)
                    for item in list(
                        dict(payload["fasttree_support_summary"])["warnings"]
                    )
                ],
            )
        ),
        sh_alrt_support_summary=(
            None
            if payload.get("sh_alrt_support_summary") is None
            else ShAlrtSupportSummaryReport(
                tree_path=Path(dict(payload["sh_alrt_support_summary"])["tree_path"]),
                internal_node_count=int(
                    dict(payload["sh_alrt_support_summary"])["internal_node_count"]
                ),
                annotated_node_count=int(
                    dict(payload["sh_alrt_support_summary"])["annotated_node_count"]
                ),
                fully_scored_node_count=int(
                    dict(payload["sh_alrt_support_summary"])["fully_scored_node_count"]
                ),
                minimum_sh_alrt_support=(
                    None
                    if dict(payload["sh_alrt_support_summary"]).get(
                        "minimum_sh_alrt_support"
                    )
                    is None
                    else float(
                        dict(payload["sh_alrt_support_summary"])[
                            "minimum_sh_alrt_support"
                        ]
                    )
                ),
                maximum_sh_alrt_support=(
                    None
                    if dict(payload["sh_alrt_support_summary"]).get(
                        "maximum_sh_alrt_support"
                    )
                    is None
                    else float(
                        dict(payload["sh_alrt_support_summary"])[
                            "maximum_sh_alrt_support"
                        ]
                    )
                ),
                minimum_ufboot_support=(
                    None
                    if dict(payload["sh_alrt_support_summary"]).get(
                        "minimum_ufboot_support"
                    )
                    is None
                    else float(
                        dict(payload["sh_alrt_support_summary"])[
                            "minimum_ufboot_support"
                        ]
                    )
                ),
                maximum_ufboot_support=(
                    None
                    if dict(payload["sh_alrt_support_summary"]).get(
                        "maximum_ufboot_support"
                    )
                    is None
                    else float(
                        dict(payload["sh_alrt_support_summary"])[
                            "maximum_ufboot_support"
                        ]
                    )
                ),
                weak_sh_alrt_clade_count=int(
                    dict(payload["sh_alrt_support_summary"])["weak_sh_alrt_clade_count"]
                ),
                weak_ufboot_clade_count=int(
                    dict(payload["sh_alrt_support_summary"])["weak_ufboot_clade_count"]
                ),
                conflicting_support_signal_count=int(
                    dict(payload["sh_alrt_support_summary"])[
                        "conflicting_support_signal_count"
                    ]
                ),
                nodes=[
                    ShAlrtSupportNode(
                        node=str(dict(item)["node"]),
                        descendant_taxa=[
                            str(taxon) for taxon in list(dict(item)["descendant_taxa"])
                        ],
                        sh_alrt_support=(
                            None
                            if dict(item).get("sh_alrt_support") is None
                            else float(dict(item)["sh_alrt_support"])
                        ),
                        sh_alrt_support_fraction=(
                            None
                            if dict(item).get("sh_alrt_support_fraction") is None
                            else float(dict(item)["sh_alrt_support_fraction"])
                        ),
                        ufboot_support=(
                            None
                            if dict(item).get("ufboot_support") is None
                            else float(dict(item)["ufboot_support"])
                        ),
                        ufboot_support_fraction=(
                            None
                            if dict(item).get("ufboot_support_fraction") is None
                            else float(dict(item)["ufboot_support_fraction"])
                        ),
                        is_backbone=bool(dict(item)["is_backbone"]),
                        sh_alrt_strong=bool(dict(item)["sh_alrt_strong"]),
                        ufboot_strong=bool(dict(item)["ufboot_strong"]),
                        conflicting_support_signal=bool(
                            dict(item)["conflicting_support_signal"]
                        ),
                        support_agreement=str(dict(item)["support_agreement"]),
                    )
                    for item in list(dict(payload["sh_alrt_support_summary"])["nodes"])
                ],
                warnings=[
                    str(item)
                    for item in list(
                        dict(payload["sh_alrt_support_summary"])["warnings"]
                    )
                ],
            )
        ),
        weak_backbone_report=(
            None
            if payload.get("weak_backbone_report") is None
            else WeakBackboneReport(
                tree_path=Path(dict(payload["weak_backbone_report"])["tree_path"]),
                threshold=float(dict(payload["weak_backbone_report"])["threshold"]),
                evaluated_backbone_node_count=int(
                    dict(payload["weak_backbone_report"])[
                        "evaluated_backbone_node_count"
                    ]
                ),
                weak_backbone_node_count=int(
                    dict(payload["weak_backbone_report"])["weak_backbone_node_count"]
                ),
                weak_nodes=[
                    BootstrapSupportNode(
                        node=str(dict(item)["node"]),
                        descendant_taxa=[
                            str(taxon) for taxon in list(dict(item)["descendant_taxa"])
                        ],
                        support=float(dict(item)["support"]),
                        support_fraction=float(dict(item)["support_fraction"]),
                        is_backbone=bool(dict(item)["is_backbone"]),
                    )
                    for item in list(
                        dict(payload["weak_backbone_report"])["weak_nodes"]
                    )
                ],
                warnings=[
                    str(item)
                    for item in list(dict(payload["weak_backbone_report"])["warnings"])
                ],
            )
        ),
        trimming_summary=None
        if payload.get("trimming_summary") is None
        else AlignmentTrimmingSummary(
            mode=str(dict(payload["trimming_summary"])["mode"]),
            gap_threshold=(
                None
                if dict(payload["trimming_summary"]).get("gap_threshold") is None
                else float(dict(payload["trimming_summary"])["gap_threshold"])
            ),
            input_alignment_length=int(
                dict(payload["trimming_summary"])["input_alignment_length"]
            ),
            trimmed_alignment_length=int(
                dict(payload["trimming_summary"])["trimmed_alignment_length"]
            ),
            retained_site_count=int(
                dict(payload["trimming_summary"])["retained_site_count"]
            ),
            removed_site_count=int(
                dict(payload["trimming_summary"])["removed_site_count"]
            ),
            retained_site_fraction=float(
                dict(payload["trimming_summary"])["retained_site_fraction"]
            ),
            removed_site_fraction=float(
                dict(payload["trimming_summary"])["removed_site_fraction"]
            ),
            input_gap_fraction=float(
                dict(payload["trimming_summary"])["input_gap_fraction"]
            ),
            trimmed_gap_fraction=float(
                dict(payload["trimming_summary"])["trimmed_gap_fraction"]
            ),
            input_gap_percentage=float(
                dict(payload["trimming_summary"])["input_gap_percentage"]
            ),
            trimmed_gap_percentage=float(
                dict(payload["trimming_summary"])["trimmed_gap_percentage"]
            ),
        ),
        resumed=bool(payload.get("resumed", False)),
        notes=[str(item) for item in payload.get("notes", [])],
    )


def _restore_codon_aware_alignment_report(
    payload: dict[str, object],
) -> CodonAwareAlignmentWorkflowReport:
    run_payload = payload["run"]
    if not isinstance(run_payload, dict):
        raise EngineWorkflowError("engine manifest contains an invalid run payload")
    version_payload = run_payload["version"]
    if not isinstance(version_payload, dict):
        raise EngineWorkflowError("engine manifest contains an invalid version payload")
    version = EngineVersionInfo(
        engine_name=str(version_payload["engine_name"]),
        executable=str(version_payload["executable"]),
        command=[str(item) for item in version_payload["command"]],
        text=str(version_payload["text"]),
    )
    run = EngineRunReport(
        engine_name=str(run_payload["engine_name"]),
        workflow=str(run_payload["workflow"]),
        executable=str(run_payload["executable"]),
        working_directory=Path(run_payload["working_directory"]),
        version=version,
        command=[str(item) for item in run_payload["command"]],
        exit_code=int(run_payload["exit_code"]),
        stdout_path=Path(run_payload["stdout_path"]),
        stderr_path=Path(run_payload["stderr_path"]),
        output_paths={
            str(key): Path(value)
            for key, value in dict(run_payload["output_paths"]).items()
        },
        warning_lines=[str(item) for item in run_payload["warning_lines"]],
        started_at_utc=str(run_payload.get("started_at_utc", "")),
        ended_at_utc=str(run_payload.get("ended_at_utc", "")),
        runtime_seconds=float(run_payload.get("runtime_seconds", 0.0)),
        timeout_seconds=(
            None
            if run_payload.get("timeout_seconds") is None
            else float(run_payload["timeout_seconds"])
        ),
        timed_out=bool(run_payload.get("timed_out", False)),
    )
    return CodonAwareAlignmentWorkflowReport(
        workflow=str(payload["workflow"]),
        engine_name=str(payload["engine_name"]),
        input_path=Path(payload["input_path"]),
        output_paths={
            str(key): Path(value)
            for key, value in dict(payload["output_paths"]).items()
        },
        run=run,
        manifest_path=Path(payload["manifest_path"]),
        input_checksums={
            str(key): str(value)
            for key, value in dict(payload.get("input_checksums", {})).items()
        },
        output_checksums={
            str(key): str(value)
            for key, value in dict(payload.get("output_checksums", {})).items()
        },
        config={
            str(key): value for key, value in dict(payload.get("config", {})).items()
        },
        sequence_type=str(payload["sequence_type"]),
        genetic_code_id=int(payload.get("genetic_code_id", 1)),
        genetic_code_name=str(payload.get("genetic_code_name", "Standard")),
        input_sequence_count=int(payload.get("input_sequence_count", 0)),
        accepted_sequence_count=int(payload["accepted_sequence_count"]),
        invalid_codon_sequence_count=int(
            payload.get("invalid_codon_sequence_count", 0)
        ),
        excluded_sequences=[
            CodingSequenceExclusion(
                identifier=str(item["identifier"]),
                comparable_length=int(item["comparable_length"]),
                reason=str(item["reason"]),
                invalid_codon_count=int(item.get("invalid_codon_count", 0)),
                premature_stop_count=int(item["premature_stop_count"]),
                terminal_stop_count=int(item["terminal_stop_count"]),
                trailing_bases=int(item["trailing_bases"]),
                note=str(item["note"]),
            )
            for item in payload.get("excluded_sequences", [])
        ],
        terminal_stop_sequence_count=int(payload["terminal_stop_sequence_count"]),
        notes=[str(item) for item in payload.get("notes", [])],
        warnings=[str(item) for item in payload.get("warnings", [])],
        resumed=bool(payload.get("resumed", False)),
    )
