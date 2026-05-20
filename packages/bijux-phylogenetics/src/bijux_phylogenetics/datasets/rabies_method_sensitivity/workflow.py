from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import os
from pathlib import Path

from bijux_phylogenetics.compare.topology import compare_tree_paths, write_tree_comparison_table
from bijux_phylogenetics.phylo.topology import root_tree_on_outgroup
from bijux_phylogenetics.engines.common import (
    EngineActiveRunRecord,
    acquire_active_engine_run,
    clear_active_engine_run,
    utc_now_text,
)
from bijux_phylogenetics.engines.inference import (
    run_tree_inference_comparison,
)
from bijux_phylogenetics.engines.workflows.alignment import (
    run_alignment_trimming,
    run_multiple_sequence_alignment,
)
from bijux_phylogenetics.io.fasta import load_fasta_alignment
from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.runtime.errors import EngineWorkflowError, PhylogeneticsError

from .config import (
    _resolve_selected_variant_dataset,
    load_rabies_method_sensitivity_panel_dataset,
)
from .models import (
    RabiesMethodSensitivityCladeRow,
    RabiesMethodSensitivityConclusionRow,
    RabiesMethodSensitivityPanelDataset,
    RabiesMethodSensitivityPanelWorkflowReport,
    RabiesMethodSensitivityPreprocessingComparisonRow,
    RabiesMethodSensitivityTaskRecord,
    RabiesMethodSensitivityVariant,
    RabiesMethodSensitivityVariantRun,
)
from .shared import _WORKFLOW_PREFIX

__all__ = ["run_rabies_method_sensitivity_panel_workflow"]


def run_rabies_method_sensitivity_panel_workflow(
    out_dir: Path,
    *,
    mafft_executable: str | Path = "mafft",
    trimal_executable: str | Path = "trimal",
    iqtree_executable: str | Path = "iqtree2",
    fasttree_executable: str | Path = "FastTree",
    iqtree_seed: int | None = None,
    iqtree_threads: int | None = None,
    bootstrap_replicates: int | None = None,
    parallel_workers: int | None = None,
    variant_ids: tuple[str, ...] | None = None,
) -> RabiesMethodSensitivityPanelWorkflowReport:
    """Run the owned method-sensitivity workflow over the packaged rabies panel."""
    dataset = _resolve_selected_variant_dataset(
        load_rabies_method_sensitivity_panel_dataset(),
        variant_ids=variant_ids,
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    started_at_utc = utc_now_text()
    resolved_seed = dataset.iqtree_seed if iqtree_seed is None else iqtree_seed
    resolved_threads = (
        dataset.iqtree_threads if iqtree_threads is None else iqtree_threads
    )
    resolved_bootstrap_replicates = (
        dataset.bootstrap_replicates
        if bootstrap_replicates is None
        else bootstrap_replicates
    )
    resolved_parallel_workers = (
        dataset.parallel_workers if parallel_workers is None else parallel_workers
    )
    if resolved_parallel_workers < 1:
        raise ValueError(
            f"parallel_workers must be at least 1, got {resolved_parallel_workers}"
        )
    execution_mode = (
        "parallel"
        if resolved_parallel_workers > 1 and len(dataset.variants) > 1
        else "serial"
    )
    execution_record_path = _workflow_execution_record_path(out_dir)
    active_record = EngineActiveRunRecord(
        engine_name="bijux",
        workflow="rabies_method_sensitivity_panel",
        executable="bijux-phylogenetics",
        working_directory=out_dir,
        manifest_path=execution_record_path,
        command=["run_rabies_method_sensitivity_panel_workflow"],
        process_id=os.getpid(),
        started_at_utc=started_at_utc,
    )
    try:
        acquire_active_engine_run(active_record)
    except EngineWorkflowError as error:
        if error.code == "engine_workflow_already_running":
            raise EngineWorkflowError(
                "rabies method-sensitivity workflow is already running for the requested output root",
                code="rabies_method_sensitivity_workflow_already_running",
                details={
                    "output_root": str(out_dir),
                    "execution_record_path": str(execution_record_path),
                    **(error.details or {}),
                },
            ) from error
        if error.code == "engine_workflow_running_marker_invalid":
            raise EngineWorkflowError(
                "rabies method-sensitivity workflow could not verify the active-run marker for the requested output root",
                code="rabies_method_sensitivity_workflow_running_marker_invalid",
                details={
                    "output_root": str(out_dir),
                    "execution_record_path": str(execution_record_path),
                    **(error.details or {}),
                },
            ) from error
        raise

    try:
        task_log_root = out_dir / "parallel-logs"
        task_log_root.mkdir(parents=True, exist_ok=True)
        task_records_by_variant: dict[str, RabiesMethodSensitivityTaskRecord] = {}
        variant_runs_by_variant: dict[str, RabiesMethodSensitivityVariantRun] = {}

        def run_variant(
            variant: RabiesMethodSensitivityVariant,
        ) -> tuple[
            RabiesMethodSensitivityTaskRecord,
            RabiesMethodSensitivityVariantRun | None,
        ]:
            variant_root = out_dir / "variants" / variant.variant_id
            variant_log_path = task_log_root / f"{variant.variant_id}.log"
            try:
                variant_run = _run_variant_workflow(
                    dataset=dataset,
                    variant=variant,
                    variant_root=variant_root,
                    mafft_executable=mafft_executable,
                    trimal_executable=trimal_executable,
                    iqtree_executable=iqtree_executable,
                    fasttree_executable=fasttree_executable,
                    iqtree_seed=resolved_seed,
                    iqtree_threads=resolved_threads,
                    bootstrap_replicates=resolved_bootstrap_replicates,
                )
            except Exception as error:
                error_code = (
                    error.code
                    if isinstance(error, PhylogeneticsError) and error.code is not None
                    else "parallel_variant_failed"
                )
                error_message = str(error)
                _write_task_log(
                    variant_log_path,
                    variant=variant,
                    execution_mode=execution_mode,
                    status="failed",
                    output_root=Path("variants") / variant.variant_id,
                    error_code=error_code,
                    error_message=error_message,
                )
                return (
                    RabiesMethodSensitivityTaskRecord(
                        variant_id=variant.variant_id,
                        label=variant.label,
                        execution_mode=execution_mode,
                        status="failed",
                        log_path=variant_log_path,
                        output_root=variant_root,
                        error_code=error_code,
                        error_message=error_message,
                    ),
                    None,
                )
            _write_task_log(
                variant_log_path,
                variant=variant,
                execution_mode=execution_mode,
                status="succeeded",
                output_root=Path("variants") / variant.variant_id,
                error_code=None,
                error_message=None,
            )
            return (
                RabiesMethodSensitivityTaskRecord(
                    variant_id=variant.variant_id,
                    label=variant.label,
                    execution_mode=execution_mode,
                    status="succeeded",
                    log_path=variant_log_path,
                    output_root=variant_root,
                    error_code=None,
                    error_message=None,
                ),
                variant_run,
            )

        if execution_mode == "serial":
            for variant in dataset.variants:
                task_record, variant_run = run_variant(variant)
                task_records_by_variant[variant.variant_id] = task_record
                if variant_run is not None:
                    variant_runs_by_variant[variant.variant_id] = variant_run
        else:
            with ThreadPoolExecutor(max_workers=resolved_parallel_workers) as executor:
                futures = {
                    executor.submit(run_variant, variant): variant.variant_id
                    for variant in dataset.variants
                }
                for future in as_completed(futures):
                    task_record, variant_run = future.result()
                    task_records_by_variant[task_record.variant_id] = task_record
                    if variant_run is not None:
                        variant_runs_by_variant[task_record.variant_id] = variant_run

        task_records = tuple(
            task_records_by_variant[variant.variant_id] for variant in dataset.variants
        )
        failed_task_records = [
            record for record in task_records if record.status != "succeeded"
        ]
        ended_at_utc = utc_now_text()
        if failed_task_records:
            _write_workflow_execution_record(
                execution_record_path,
                dataset=dataset,
                execution_mode=execution_mode,
                parallel_workers=resolved_parallel_workers,
                task_records=task_records,
                started_at_utc=started_at_utc,
                ended_at_utc=ended_at_utc,
                status="failed",
            )
            raise EngineWorkflowError(
                "rabies method-sensitivity workflow left one or more parallel tasks failed while preserving successful isolated outputs",
                code="workflow_parallel_task_failed",
                details={
                    "failed_variants": [
                        record.variant_id for record in failed_task_records
                    ],
                    "successful_variants": [
                        record.variant_id
                        for record in task_records
                        if record.status == "succeeded"
                    ],
                    "parallel_workers": resolved_parallel_workers,
                    "execution_mode": execution_mode,
                    "task_logs": {
                        record.variant_id: str(record.log_path)
                        for record in task_records
                    },
                    "execution_record_path": str(execution_record_path),
                },
            )
        variant_runs = [
            variant_runs_by_variant[variant.variant_id] for variant in dataset.variants
        ]

        preprocessing_comparison_rows = tuple(
            _build_preprocessing_comparison_rows(variant_runs)
        )
        stable_clade_rows = tuple(_aggregate_clades(variant_runs, stable_only=True))
        changed_clade_rows = tuple(_aggregate_clades(variant_runs, stable_only=False))
        conclusion_rows = tuple(
            _build_conclusion_rows(
                dataset=dataset,
                variant_runs=variant_runs,
                preprocessing_comparison_rows=preprocessing_comparison_rows,
            )
        )
        report = RabiesMethodSensitivityPanelWorkflowReport(
            dataset=dataset,
            execution_record_path=execution_record_path,
            iqtree_seed=resolved_seed,
            iqtree_threads=resolved_threads,
            bootstrap_replicates=resolved_bootstrap_replicates,
            parallel_workers=resolved_parallel_workers,
            execution_mode=execution_mode,
            task_records=task_records,
            variant_runs=tuple(variant_runs),
            preprocessing_comparison_rows=preprocessing_comparison_rows,
            stable_clade_rows=stable_clade_rows,
            changed_clade_rows=changed_clade_rows,
            conclusion_rows=conclusion_rows,
        )
        _write_workflow_execution_record(
            execution_record_path,
            dataset=dataset,
            execution_mode=execution_mode,
            parallel_workers=resolved_parallel_workers,
            task_records=task_records,
            started_at_utc=started_at_utc,
            ended_at_utc=ended_at_utc,
            status="succeeded",
        )
        return report
    finally:
        clear_active_engine_run(execution_record_path)


def _build_preprocessing_comparison_rows(
    variant_runs: list[RabiesMethodSensitivityVariantRun],
) -> list[RabiesMethodSensitivityPreprocessingComparisonRow]:
    rows: list[RabiesMethodSensitivityPreprocessingComparisonRow] = []
    for index, left in enumerate(variant_runs):
        for right in variant_runs[index + 1 :]:
            comparison = compare_tree_paths(
                left.rooted_iqtree_path, right.rooted_iqtree_path
            )
            rows.append(
                RabiesMethodSensitivityPreprocessingComparisonRow(
                    left_variant_id=left.config.variant_id,
                    right_variant_id=right.config.variant_id,
                    comparison_axis=_comparison_axis(left.config, right.config),
                    robinson_foulds_distance=comparison.robinson_foulds_distance,
                    normalized_robinson_foulds=comparison.normalized_robinson_foulds,
                    same_taxa_different_rooting=comparison.same_taxa_different_rooting,
                )
            )
    return rows


def _run_variant_workflow(
    *,
    dataset: RabiesMethodSensitivityPanelDataset,
    variant: RabiesMethodSensitivityVariant,
    variant_root: Path,
    mafft_executable: str | Path,
    trimal_executable: str | Path,
    iqtree_executable: str | Path,
    fasttree_executable: str | Path,
    iqtree_seed: int,
    iqtree_threads: int,
    bootstrap_replicates: int,
) -> RabiesMethodSensitivityVariantRun:
    alignment_path = variant_root / f"{variant.variant_id}.aln"
    trimmed_alignment_path = variant_root / f"{variant.variant_id}.trimmed.aln"
    alignment_workflow = run_multiple_sequence_alignment(
        dataset.sequences_path,
        alignment_path,
        executable=mafft_executable,
        mode=variant.alignment_mode,
    )
    trimming_workflow = run_alignment_trimming(
        alignment_path,
        trimmed_alignment_path,
        executable=trimal_executable,
        mode=variant.trimming_mode,
        gap_threshold=variant.trim_gap_threshold,
    )
    inference_comparison = run_tree_inference_comparison(
        trimmed_alignment_path,
        out_dir=variant_root / "engine-comparison",
        prefix=variant.variant_id,
        sequence_type=dataset.sequence_type,
        iqtree_executable=iqtree_executable,
        fasttree_executable=fasttree_executable,
        iqtree_seed=iqtree_seed,
        iqtree_threads=iqtree_threads,
        bootstrap_replicates=bootstrap_replicates,
    )
    fasttree_rooted, fasttree_rooting = root_tree_on_outgroup(
        inference_comparison.output_paths["fasttree_tree"],
        outgroup_taxa=list(dataset.outgroup_taxa),
    )
    iqtree_rooted, iqtree_rooting = root_tree_on_outgroup(
        inference_comparison.output_paths["iqtree_support_tree"],
        outgroup_taxa=list(dataset.outgroup_taxa),
    )
    rooted_fasttree_path = variant_root / "rooted-fasttree.nwk"
    rooted_iqtree_path = variant_root / "rooted-iqtree-support.nwk"
    write_newick(rooted_fasttree_path, fasttree_rooted)
    write_newick(rooted_iqtree_path, iqtree_rooted)
    rooted_engine_comparison = compare_tree_paths(
        rooted_fasttree_path, rooted_iqtree_path
    )
    rooted_engine_comparison_table_path = write_tree_comparison_table(
        variant_root / "rooted-engine-comparison.tsv",
        rooted_fasttree_path,
        rooted_iqtree_path,
    )
    aligned_records = load_fasta_alignment(alignment_workflow.output_paths["alignment"])
    trimmed_records = load_fasta_alignment(
        trimming_workflow.output_paths["trimmed_alignment"]
    )
    return RabiesMethodSensitivityVariantRun(
        config=variant,
        alignment_workflow=alignment_workflow,
        trimming_workflow=trimming_workflow,
        inference_comparison=inference_comparison,
        rooted_fasttree_path=rooted_fasttree_path,
        rooted_iqtree_path=rooted_iqtree_path,
        fasttree_rooting=fasttree_rooting,
        iqtree_rooting=iqtree_rooting,
        rooted_engine_comparison=rooted_engine_comparison,
        rooted_engine_comparison_table_path=rooted_engine_comparison_table_path,
        alignment_length=len(aligned_records[0].sequence),
        trimmed_alignment_length=len(trimmed_records[0].sequence),
    )


def _comparison_axis(
    left: RabiesMethodSensitivityVariant, right: RabiesMethodSensitivityVariant
) -> str:
    alignment_changed = left.alignment_mode != right.alignment_mode
    trimming_changed = left.trimming_mode != right.trimming_mode
    if alignment_changed and not trimming_changed:
        return "alignment_mode"
    if trimming_changed and not alignment_changed:
        return "trimming_mode"
    return "combined_preprocessing"


def _aggregate_clades(
    variant_runs: list[RabiesMethodSensitivityVariantRun], *, stable_only: bool
) -> list[RabiesMethodSensitivityCladeRow]:
    counts: dict[tuple[str, str], tuple[str, str, int]] = {}
    for variant in variant_runs:
        for row in variant.inference_comparison.conclusion_rows:
            is_stable = row.conclusion_class == "stable_clade"
            if is_stable != stable_only:
                continue
            key = (row.split_id, row.conclusion_class)
            evidence_class, detail, count = counts.get(
                key, (row.evidence_class, row.detail, 0)
            )
            counts[key] = (evidence_class, detail, count + 1)
    return [
        RabiesMethodSensitivityCladeRow(
            split_id=split_id,
            conclusion_class=conclusion_class,
            evidence_class=evidence_class,
            occurrence_count=count,
            variant_count=len(variant_runs),
            detail=detail,
        )
        for (split_id, conclusion_class), (evidence_class, detail, count) in sorted(
            counts.items()
        )
    ]


def _build_conclusion_rows(
    *,
    dataset: RabiesMethodSensitivityPanelDataset,
    variant_runs: list[RabiesMethodSensitivityVariantRun],
    preprocessing_comparison_rows: tuple[
        RabiesMethodSensitivityPreprocessingComparisonRow, ...
    ],
) -> list[RabiesMethodSensitivityConclusionRow]:
    selected_models = sorted(
        {
            variant.inference_comparison.selected_model
            for variant in variant_runs
            if variant.inference_comparison.selected_model
        }
    )
    max_serious_conflicts = max(
        variant.inference_comparison.conclusion_summary.serious_conflict_count
        for variant in variant_runs
    )
    stable_clade_counts = {
        variant.config.variant_id: variant.inference_comparison.conclusion_summary.stable_clade_count
        for variant in variant_runs
    }
    return [
        RabiesMethodSensitivityConclusionRow(
            conclusion_id="preprocessing_rooted_iqtree_topology",
            method_axis="alignment_and_trimming",
            stability_status=(
                "stable"
                if all(
                    row.robinson_foulds_distance == 0
                    and not row.same_taxa_different_rooting
                    for row in preprocessing_comparison_rows
                )
                else "changed"
            ),
            claim=(
                "The rooted IQ-TREE topology stayed unchanged across every declared "
                "alignment and trimming variant."
            ),
            evidence=(
                f"{len(preprocessing_comparison_rows)} rooted pairwise preprocessing "
                "comparisons returned RF distance 0 and no rooting-only disagreements."
            ),
            caution=(
                "This stability statement is limited to the compact nine-taxon rabies "
                "panel and the four declared preprocessing settings."
            ),
        ),
        RabiesMethodSensitivityConclusionRow(
            conclusion_id="preprocessing_selected_model",
            method_axis="alignment_and_trimming",
            stability_status="stable" if len(selected_models) == 1 else "changed",
            claim=(
                "The selected substitution model remained constant across the "
                "declared preprocessing matrix."
                if len(selected_models) == 1
                else "The selected substitution model changed across the declared preprocessing matrix."
            ),
            evidence=(
                f"selected models: {', '.join(selected_models)}"
                if selected_models
                else "no selected model was recorded"
            ),
            caution=(
                "Model-selection stability here reflects one short rabies nucleoprotein "
                "panel rather than a general claim about all pathogen alignments."
            ),
        ),
        RabiesMethodSensitivityConclusionRow(
            conclusion_id="rooted_engine_agreement",
            method_axis="inference_engine",
            stability_status=(
                "stable"
                if all(
                    variant.rooted_engine_comparison.robinson_foulds_distance == 0
                    and not variant.rooted_engine_comparison.same_taxa_different_rooting
                    for variant in variant_runs
                )
                else "changed"
            ),
            claim=(
                "After explicit outgroup rooting, FastTree and IQ-TREE preserved the "
                "same rooted topology in every declared preprocessing variant."
            ),
            evidence=(
                f"rooted engine comparisons over outgroup {', '.join(dataset.outgroup_taxa)} "
                "returned RF distance 0 in every variant."
            ),
            caution=(
                "Rooted agreement does not imply that every internal unrooted split or "
                "support value is interchangeable across engines."
            ),
        ),
        RabiesMethodSensitivityConclusionRow(
            conclusion_id="unrooted_engine_sensitivity",
            method_axis="inference_engine",
            stability_status="changed" if max_serious_conflicts > 0 else "stable",
            claim=(
                "Before rooting, the FastTree versus IQ-TREE comparison changed several "
                "internal clade conclusions on this rabies panel."
            ),
            evidence=(
                f"serious unrooted engine conflicts ranged up to {max_serious_conflicts} "
                f"per variant, while stable shared clades per variant ranged from "
                f"{min(stable_clade_counts.values())} to {max(stable_clade_counts.values())}."
            ),
            caution=(
                "The engine-sensitive clades are a warning against over-reading fine "
                "internal structure from one compact panel, especially when approximate "
                "and likelihood engines disagree before rooting."
            ),
        ),
    ]


def _workflow_execution_record_path(output_root: Path) -> Path:
    """Return the workflow-root execution record for one packaged sensitivity run."""
    return output_root / f"{_WORKFLOW_PREFIX}.run.json"


def _write_workflow_execution_record(
    path: Path,
    *,
    dataset: RabiesMethodSensitivityPanelDataset,
    execution_mode: str,
    parallel_workers: int,
    task_records: tuple[RabiesMethodSensitivityTaskRecord, ...],
    started_at_utc: str,
    ended_at_utc: str,
    status: str,
) -> Path:
    """Persist one execution record for one packaged sensitivity workflow run."""
    payload = {
        "dataset_id": dataset.dataset_id,
        "workflow": "rabies_method_sensitivity_panel",
        "workflow_prefix": dataset.workflow_prefix,
        "status": status,
        "started_at_utc": started_at_utc,
        "ended_at_utc": ended_at_utc,
        "parallel_workers": parallel_workers,
        "execution_mode": execution_mode,
        "variant_count": len(dataset.variants),
        "selected_variant_ids": [variant.variant_id for variant in dataset.variants],
        "successful_variants": [
            record.variant_id for record in task_records if record.status == "succeeded"
        ],
        "failed_variants": [
            record.variant_id for record in task_records if record.status != "succeeded"
        ],
        "task_records": [
            {
                "variant_id": record.variant_id,
                "label": record.label,
                "status": record.status,
                "execution_mode": record.execution_mode,
                "log_path": Path("parallel-logs", record.log_path.name).as_posix(),
                "output_root": Path("variants", record.variant_id).as_posix(),
                "error_code": record.error_code,
                "error_message": record.error_message,
            }
            for record in task_records
        ],
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def _write_task_log(
    path: Path,
    *,
    variant: RabiesMethodSensitivityVariant,
    execution_mode: str,
    status: str,
    output_root: Path,
    error_code: str | None,
    error_message: str | None,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        f"variant_id: {variant.variant_id}",
        f"label: {variant.label}",
        f"execution_mode: {execution_mode}",
        f"alignment_mode: {variant.alignment_mode}",
        f"trimming_mode: {variant.trimming_mode}",
        f"trim_gap_threshold: {_format_float(variant.trim_gap_threshold)}",
        f"status: {status}",
        f"output_root: {output_root.as_posix()}",
    ]
    if error_code is not None:
        lines.append(f"error_code: {error_code}")
    if error_message is not None:
        lines.append(f"error_message: {error_message}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _format_float(value: float) -> str:
    return format(value, ".12g")
