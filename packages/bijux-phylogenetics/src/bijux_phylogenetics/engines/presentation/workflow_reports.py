from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path

from bijux_phylogenetics.render.html import write_html_report

from ..artifacts.bootstrap import build_low_support_bootstrap_rows
from ..artifacts.fasttree import build_fasttree_low_support_rows
from ..artifacts.sh_alrt import build_conflicting_sh_alrt_support_rows
from ..common import load_engine_manifest
from ..validation import (
    classify_inference_workflow_failure,
    detect_weakly_supported_backbone,
    summarize_bootstrap_support_distribution,
    summarize_fasttree_support_distribution,
    summarize_sh_alrt_support_distribution,
    validate_bootstrap_tree_set,
    validate_inference_engine_outputs,
    validate_ml_tree_contains_expected_taxa,
)
from .model_limitations import build_model_selection_limitations_report


@dataclass(slots=True)
class InferenceWorkflowReportBuildResult:
    output_path: Path
    report_kind: str
    title: str
    manifest_path: Path
    engine_name: str
    workflow: str
    warning_count: int
    supplement_sections: list[str]
    machine_manifest: dict[str, object]


def render_inference_workflow_report(
    *, manifest_path: Path, out_path: Path
) -> InferenceWorkflowReportBuildResult:
    """Render a deterministic HTML report for one engine workflow manifest."""
    manifest = load_engine_manifest(manifest_path)
    output_paths = {
        key: Path(value) for key, value in dict(manifest["output_paths"]).items()
    }
    input_paths = [Path(path) for path in manifest["input_paths"]]
    failure = classify_inference_workflow_failure(
        workflow=str(manifest["workflow"]),
        input_paths=input_paths,
        output_paths=output_paths,
        run_exit_code=int(manifest["run"].get("exit_code", 0)),
    )
    consistency = validate_inference_engine_outputs(manifest_path)
    title = f"Bijux Inference Workflow Report: {manifest['workflow']}"
    sections = [
        (
            "workflow-summary",
            json.dumps(
                {
                    "engine_name": manifest["engine_name"],
                    "workflow": manifest["workflow"],
                    "resumed": manifest.get("resumed", False),
                    "selected_model": manifest.get("selected_model"),
                    "notes": manifest.get("notes", []),
                },
                indent=2,
                sort_keys=True,
            ),
        ),
        (
            "workflow-failure-taxonomy",
            json.dumps(asdict(failure), default=str, indent=2, sort_keys=True),
        ),
        (
            "workflow-consistency",
            json.dumps(asdict(consistency), default=str, indent=2, sort_keys=True),
        ),
        (
            "inputs",
            json.dumps(
                {
                    "input_paths": manifest["input_paths"],
                    "input_checksums": manifest.get("input_checksums", {}),
                },
                indent=2,
                sort_keys=True,
            ),
        ),
        (
            "outputs",
            json.dumps(
                {
                    "output_paths": manifest["output_paths"],
                    "output_checksums": manifest.get("output_checksums", {}),
                },
                indent=2,
                sort_keys=True,
            ),
        ),
        ("engine-run", json.dumps(manifest["run"], indent=2, sort_keys=True)),
    ]
    supplement_sections = [
        "workflow-summary",
        "workflow-failure-taxonomy",
        "workflow-consistency",
        "inputs",
        "outputs",
        "engine-run",
    ]
    if manifest["workflow"] == "model-selection":
        sections.append(
            (
                "model-selection-limitations",
                json.dumps(
                    asdict(build_model_selection_limitations_report(manifest_path)),
                    default=str,
                    indent=2,
                    sort_keys=True,
                ),
            )
        )
        supplement_sections.append("model-selection-limitations")
    if (
        manifest["workflow"] == "alignment-trimming"
        and manifest.get("trimming_summary") is not None
    ):
        sections.append(
            (
                "alignment-trimming-summary",
                json.dumps(
                    manifest["trimming_summary"],
                    indent=2,
                    sort_keys=True,
                ),
            )
        )
        supplement_sections.append("alignment-trimming-summary")
    if manifest["workflow"] == "maximum-likelihood-tree":
        sections.append(
            (
                "ml-tree-validation",
                json.dumps(
                    asdict(validate_ml_tree_contains_expected_taxa(manifest_path)),
                    default=str,
                    indent=2,
                    sort_keys=True,
                ),
            )
        )
        supplement_sections.append("ml-tree-validation")
    if manifest["workflow"] == "bootstrap-support":
        bootstrap_path = output_paths.get("bootstrap_trees")
        if bootstrap_path is not None:
            sections.append(
                (
                    "bootstrap-tree-set-validation",
                    json.dumps(
                        asdict(validate_bootstrap_tree_set(bootstrap_path)),
                        default=str,
                        indent=2,
                        sort_keys=True,
                    ),
                )
            )
            supplement_sections.append("bootstrap-tree-set-validation")
        tree_path = output_paths.get("tree") or output_paths.get("support_tree")
        if tree_path is not None:
            bootstrap_support_summary = summarize_bootstrap_support_distribution(
                tree_path
            )
            sections.append(
                (
                    "bootstrap-support-summary",
                    json.dumps(
                        asdict(bootstrap_support_summary),
                        default=str,
                        indent=2,
                        sort_keys=True,
                    ),
                )
            )
            sections.append(
                (
                    "bootstrap-support-histogram",
                    json.dumps(
                        bootstrap_support_summary.support_histogram,
                        indent=2,
                        sort_keys=True,
                    ),
                )
            )
            sections.append(
                (
                    "low-support-branches",
                    json.dumps(
                        [
                            asdict(row)
                            for row in build_low_support_bootstrap_rows(
                                bootstrap_support_summary
                            )
                        ],
                        default=str,
                        indent=2,
                        sort_keys=True,
                    ),
                )
            )
            sections.append(
                (
                    "weak-backbone",
                    json.dumps(
                        asdict(detect_weakly_supported_backbone(tree_path)),
                        default=str,
                        indent=2,
                        sort_keys=True,
                    ),
                )
            )
            supplement_sections.extend(
                [
                    "bootstrap-support-summary",
                    "bootstrap-support-histogram",
                    "low-support-branches",
                    "weak-backbone",
                ]
            )
    if manifest["workflow"] == "fast-approximate-tree":
        tree_path = output_paths.get("tree") or output_paths.get("support_tree")
        if tree_path is not None:
            fasttree_support_summary = summarize_fasttree_support_distribution(
                tree_path
            )
            sections.append(
                (
                    "fasttree-approximation-limits",
                    json.dumps(
                        {
                            "approximate_method": fasttree_support_summary.approximate_method,
                            "support_label_kind": fasttree_support_summary.support_label_kind,
                            "support_scale": fasttree_support_summary.support_scale,
                            "limitations": [
                                "FastTree uses an approximately maximum-likelihood search and should be reviewed as exploratory or large-alignment evidence rather than as a direct substitute for a fully optimized ML workflow",
                                "FastTree local support values are SH-like proportions on a 0-to-1 scale and are not bootstrap percentages",
                            ],
                        },
                        indent=2,
                        sort_keys=True,
                    ),
                )
            )
            sections.append(
                (
                    "fasttree-support-summary",
                    json.dumps(
                        asdict(fasttree_support_summary),
                        default=str,
                        indent=2,
                        sort_keys=True,
                    ),
                )
            )
            sections.append(
                (
                    "fasttree-support-histogram",
                    json.dumps(
                        fasttree_support_summary.support_histogram,
                        indent=2,
                        sort_keys=True,
                    ),
                )
            )
            sections.append(
                (
                    "fasttree-low-support-branches",
                    json.dumps(
                        [
                            asdict(row)
                            for row in build_fasttree_low_support_rows(
                                fasttree_support_summary
                            )
                        ],
                        default=str,
                        indent=2,
                        sort_keys=True,
                    ),
                )
            )
            supplement_sections.extend(
                [
                    "fasttree-approximation-limits",
                    "fasttree-support-summary",
                    "fasttree-support-histogram",
                    "fasttree-low-support-branches",
                ]
            )
    if manifest["workflow"] == "sh-alrt-support":
        bootstrap_path = output_paths.get("bootstrap_trees")
        if bootstrap_path is not None:
            sections.append(
                (
                    "bootstrap-tree-set-validation",
                    json.dumps(
                        asdict(validate_bootstrap_tree_set(bootstrap_path)),
                        default=str,
                        indent=2,
                        sort_keys=True,
                    ),
                )
            )
            supplement_sections.append("bootstrap-tree-set-validation")
        tree_path = output_paths.get("tree") or output_paths.get("support_tree")
        if tree_path is not None:
            sh_alrt_support_summary = summarize_sh_alrt_support_distribution(tree_path)
            sections.append(
                (
                    "sh-alrt-support-summary",
                    json.dumps(
                        asdict(sh_alrt_support_summary),
                        default=str,
                        indent=2,
                        sort_keys=True,
                    ),
                )
            )
            sections.append(
                (
                    "conflicting-support-branches",
                    json.dumps(
                        [
                            asdict(row)
                            for row in build_conflicting_sh_alrt_support_rows(
                                sh_alrt_support_summary
                            )
                        ],
                        default=str,
                        indent=2,
                        sort_keys=True,
                    ),
                )
            )
            sections.append(
                (
                    "weak-backbone",
                    json.dumps(
                        asdict(detect_weakly_supported_backbone(tree_path)),
                        default=str,
                        indent=2,
                        sort_keys=True,
                    ),
                )
            )
            supplement_sections.extend(
                [
                    "sh-alrt-support-summary",
                    "conflicting-support-branches",
                    "weak-backbone",
                ]
            )
    machine_manifest = {
        "report_kind": "inference-workflow",
        "title": title,
        "manifest_path": str(manifest_path),
        "engine_name": manifest["engine_name"],
        "workflow": manifest["workflow"],
        "warning_count": len(manifest["run"].get("warning_lines", [])),
        "sections": [name for name, _ in sections],
        "supplement_sections": supplement_sections,
    }
    write_html_report(
        title=title,
        sections=sections,
        out_path=out_path,
        embedded_json=machine_manifest,
    )
    return InferenceWorkflowReportBuildResult(
        output_path=out_path,
        report_kind="inference-workflow",
        title=title,
        manifest_path=manifest_path,
        engine_name=str(manifest["engine_name"]),
        workflow=str(manifest["workflow"]),
        warning_count=len(manifest["run"].get("warning_lines", [])),
        supplement_sections=supplement_sections,
        machine_manifest=machine_manifest,
    )
