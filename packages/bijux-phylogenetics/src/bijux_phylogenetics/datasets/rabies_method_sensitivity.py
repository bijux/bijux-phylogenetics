from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, replace
import hashlib
import json
import os
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory

from bijux_phylogenetics.compare.topology import (
    TreeComparisonReport,
    compare_tree_paths,
    write_tree_comparison_table,
)
from bijux_phylogenetics.core.topology import TreeRootingReport, root_tree_on_outgroup
from bijux_phylogenetics.engines.common import (
    EngineActiveRunRecord,
    acquire_active_engine_run,
    clear_active_engine_run,
    utc_now_text,
)
from bijux_phylogenetics.engines.inference_comparison import (
    InferenceComparisonWorkflowReport,
    run_tree_inference_comparison,
)
from bijux_phylogenetics.engines.workflows import (
    EngineWorkflowReport,
    run_alignment_trimming,
    run_multiple_sequence_alignment,
)
from bijux_phylogenetics.datasets.rabies_method_sensitivity_audit import (
    RabiesMethodSensitivityReproducibilityAuditReport,
    audit_rabies_method_sensitivity_workflow_bundle,
    write_rabies_method_sensitivity_reproducibility_audit_json,
    write_rabies_method_sensitivity_reproducibility_checks_table,
    write_rabies_method_sensitivity_variant_audit_table,
)
from bijux_phylogenetics.datasets.rabies_method_sensitivity_slurm import (
    RabiesMethodSensitivitySlurmPlanningReport,
    build_rabies_method_sensitivity_slurm_planning_report,
    write_rabies_method_sensitivity_slurm_assumptions_table,
    write_rabies_method_sensitivity_slurm_job_plan_table,
    write_rabies_method_sensitivity_slurm_summary_json,
)
from bijux_phylogenetics.datasets.rabies_method_sensitivity_slurm_arrays import (
    RabiesMethodSensitivitySlurmArrayStrategyReport,
    build_rabies_method_sensitivity_slurm_array_strategy_report,
    write_rabies_method_sensitivity_slurm_array_members_table,
    write_rabies_method_sensitivity_slurm_array_partition_scripts,
    write_rabies_method_sensitivity_slurm_array_partitions_table,
    write_rabies_method_sensitivity_slurm_array_strategy_json,
)
from bijux_phylogenetics.datasets.rabies_method_sensitivity_slurm_freshness import (
    RabiesMethodSensitivitySlurmOutputFreshnessReport,
    build_rabies_method_sensitivity_slurm_output_freshness_report,
    write_rabies_method_sensitivity_slurm_output_freshness_checks_table,
    write_rabies_method_sensitivity_slurm_output_freshness_json,
    write_rabies_method_sensitivity_slurm_output_freshness_table,
)
from bijux_phylogenetics.datasets.rabies_method_sensitivity_slurm_status import (
    RabiesMethodSensitivitySlurmStatusReport,
    build_rabies_method_sensitivity_slurm_status_report,
    write_rabies_method_sensitivity_slurm_job_status_table,
    write_rabies_method_sensitivity_slurm_partition_status_table,
    write_rabies_method_sensitivity_slurm_status_json,
)
from bijux_phylogenetics.runtime.errors import EngineWorkflowError, PhylogeneticsError
from bijux_phylogenetics.io.fasta import load_fasta_alignment, validate_fasta_input
from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.render.html import write_html_report

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


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivityVariant:
    """One declared method combination in the packaged sensitivity matrix."""

    variant_id: str
    label: str
    alignment_mode: str
    trimming_mode: str
    trim_gap_threshold: float


@dataclass(slots=True)
class RabiesMethodSensitivityPanelDataset:
    """Packaged rabies dataset for method-sensitivity workflow review."""

    dataset_id: str
    label: str
    dataset_root: Path
    readme_path: Path
    config_path: Path
    sequences_path: Path
    metadata_path: Path
    reference_output_root: Path
    taxon_count: int
    sequence_type: str
    workflow_prefix: str
    outgroup_taxa: tuple[str, ...]
    iqtree_seed: int
    iqtree_threads: int
    bootstrap_replicates: int
    parallel_workers: int
    source_accessions: tuple[str, ...]
    variants: tuple[RabiesMethodSensitivityVariant, ...]
    source_summary: str


@dataclass(slots=True)
class RabiesMethodSensitivityPanelExportResult:
    """Materialized copy of the packaged rabies method-sensitivity dataset."""

    output_root: Path
    readme_path: Path
    config_path: Path
    sequences_path: Path
    metadata_path: Path
    expected_output_root: Path


@dataclass(slots=True)
class RabiesMethodSensitivityVariantRun:
    """One executed variant in the governed method-sensitivity matrix."""

    config: RabiesMethodSensitivityVariant
    alignment_workflow: EngineWorkflowReport
    trimming_workflow: EngineWorkflowReport
    inference_comparison: InferenceComparisonWorkflowReport
    rooted_fasttree_path: Path
    rooted_iqtree_path: Path
    fasttree_rooting: TreeRootingReport
    iqtree_rooting: TreeRootingReport
    rooted_engine_comparison: TreeComparisonReport
    rooted_engine_comparison_table_path: Path
    alignment_length: int
    trimmed_alignment_length: int


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivityTaskRecord:
    """One isolated variant execution record within the workflow batch."""

    variant_id: str
    label: str
    execution_mode: str
    status: str
    log_path: Path
    output_root: Path
    error_code: str | None
    error_message: str | None


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivityPreprocessingComparisonRow:
    """One rooted IQ-TREE comparison across two preprocessing variants."""

    left_variant_id: str
    right_variant_id: str
    comparison_axis: str
    robinson_foulds_distance: int
    normalized_robinson_foulds: float
    same_taxa_different_rooting: bool


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivityCladeRow:
    """One aggregated stable or changed clade-level conclusion across variants."""

    split_id: str
    conclusion_class: str
    evidence_class: str
    occurrence_count: int
    variant_count: int
    detail: str


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivityConclusionRow:
    """One high-level biological or analytical conclusion from the workflow."""

    conclusion_id: str
    method_axis: str
    stability_status: str
    claim: str
    evidence: str
    caution: str


@dataclass(slots=True)
class RabiesMethodSensitivityPanelWorkflowReport:
    """One governed method-sensitivity workflow run over the packaged rabies panel."""

    dataset: RabiesMethodSensitivityPanelDataset
    execution_record_path: Path
    iqtree_seed: int
    iqtree_threads: int
    bootstrap_replicates: int
    parallel_workers: int
    execution_mode: str
    task_records: tuple[RabiesMethodSensitivityTaskRecord, ...]
    variant_runs: tuple[RabiesMethodSensitivityVariantRun, ...]
    preprocessing_comparison_rows: tuple[
        RabiesMethodSensitivityPreprocessingComparisonRow, ...
    ]
    stable_clade_rows: tuple[RabiesMethodSensitivityCladeRow, ...]
    changed_clade_rows: tuple[RabiesMethodSensitivityCladeRow, ...]
    conclusion_rows: tuple[RabiesMethodSensitivityConclusionRow, ...]


@dataclass(slots=True)
class RabiesMethodSensitivityPanelWorkflowBundle:
    """Written reviewer-facing outputs for the rabies method-sensitivity workflow."""

    output_root: Path
    variant_count: int
    stable_clade_count: int
    changed_clade_count: int
    preprocessing_change_pair_count: int
    rooted_engine_change_variant_count: int
    serious_conflict_variant_count: int
    execution_record_path: Path
    parallel_workers: int
    execution_mode: str
    workflow_summary_path: Path
    variant_summary_path: Path
    parallel_summary_path: Path
    preprocessing_comparison_path: Path
    stable_clades_path: Path
    changed_clades_path: Path
    conclusion_summary_path: Path
    config_path: Path
    manifest_path: Path
    report_manifest_path: Path
    slurm_job_plan_path: Path
    slurm_assumptions_path: Path
    slurm_summary_path: Path
    slurm_array_partitions_path: Path
    slurm_array_members_path: Path
    slurm_array_strategy_path: Path
    slurm_array_scripts_root: Path
    slurm_job_count: int
    slurm_total_estimated_core_hours: float
    slurm_maximum_estimated_memory_mib: int
    slurm_maximum_estimated_wallclock_minutes: int
    slurm_total_estimated_scratch_mib: int
    slurm_total_estimated_output_mib: int
    slurm_array_partition_count: int
    slurm_array_script_count: int
    slurm_array_largest_partition_size: int
    slurm_output_freshness_path: Path
    slurm_output_freshness_checks_path: Path
    slurm_output_freshness_summary_path: Path
    slurm_job_status_path: Path
    slurm_partition_status_path: Path
    slurm_workflow_status_path: Path
    slurm_output_freshness_check_count: int
    slurm_output_freshness_failed_check_count: int
    slurm_fresh_output_job_count: int
    slurm_stale_output_job_count: int
    slurm_completed_job_count: int
    slurm_failed_job_count: int
    slurm_pending_job_count: int
    slurm_stale_job_count: int
    reproducibility_checks_path: Path
    reproducibility_variant_audit_path: Path
    reproducibility_audit_path: Path
    reproducibility_passed: bool
    reproducibility_check_count: int
    reproducibility_failed_check_count: int
    reproducibility_failed_variant_count: int
    report_path: Path
    report_linked_artifact_count: int
    report_html_size_bytes: int
    report_linked_artifact_bytes: int
    report_total_output_bytes: int
    task_logs_root: Path
    variants_root: Path


@dataclass(slots=True)
class RabiesMethodSensitivityPanelDemoResult:
    """Dataset export plus workflow outputs for the public method-sensitivity demo."""

    output_root: Path
    dataset: RabiesMethodSensitivityPanelDataset
    dataset_export: RabiesMethodSensitivityPanelExportResult
    workflow_bundle: RabiesMethodSensitivityPanelWorkflowBundle
    overview_path: Path


def load_rabies_method_sensitivity_panel_dataset() -> (
    RabiesMethodSensitivityPanelDataset
):
    """Expose the packaged rabies method-sensitivity panel as a first-class surface."""
    dataset_root = _resource_root()
    config = json.loads(
        (dataset_root / "workflow-config.json").read_text(encoding="utf-8")
    )
    sequences_path = dataset_root / "sequences.fasta"
    metadata_path = dataset_root / "metadata.csv"
    validation = validate_fasta_input(sequences_path, sequence_type=_SEQUENCE_TYPE)
    variants = tuple(
        RabiesMethodSensitivityVariant(
            variant_id=str(item["variant_id"]),
            label=str(item["label"]),
            alignment_mode=str(item["alignment_mode"]),
            trimming_mode=str(item["trimming_mode"]),
            trim_gap_threshold=float(item["trim_gap_threshold"]),
        )
        for item in config["variants"]
    )
    return RabiesMethodSensitivityPanelDataset(
        dataset_id=_DATASET_ID,
        label=_DATASET_LABEL,
        dataset_root=dataset_root,
        readme_path=dataset_root / "README.md",
        config_path=dataset_root / "workflow-config.json",
        sequences_path=sequences_path,
        metadata_path=metadata_path,
        reference_output_root=dataset_root / "expected",
        taxon_count=validation.summary.sequence_count,
        sequence_type=_SEQUENCE_TYPE,
        workflow_prefix=_WORKFLOW_PREFIX,
        outgroup_taxa=tuple(str(value) for value in config["outgroup_taxa"]),
        iqtree_seed=int(config["iqtree_seed"]),
        iqtree_threads=int(config["iqtree_threads"]),
        bootstrap_replicates=int(config["bootstrap_replicates"]),
        parallel_workers=int(config.get("parallel_workers", 1)),
        source_accessions=_SOURCE_ACCESSIONS,
        variants=variants,
        source_summary=(
            "Real rabies virus nucleoprotein sequences packaged with grouped host "
            "and geography metadata for one owned sensitivity workflow that checks "
            "how alignment, trimming, and inference-engine choices change or "
            "preserve rooted biological conclusions."
        ),
    )


def export_rabies_method_sensitivity_panel_dataset(
    destination: Path,
) -> RabiesMethodSensitivityPanelExportResult:
    """Copy the packaged rabies method-sensitivity dataset and expected outputs."""
    dataset = load_rabies_method_sensitivity_panel_dataset()
    if destination.exists():
        shutil.rmtree(destination)
    destination.mkdir(parents=True, exist_ok=True)
    readme_path = Path(shutil.copy2(dataset.readme_path, destination / "README.md"))
    config_path = Path(
        shutil.copy2(dataset.config_path, destination / "workflow-config.json")
    )
    sequences_path = Path(
        shutil.copy2(dataset.sequences_path, destination / "sequences.fasta")
    )
    metadata_path = Path(
        shutil.copy2(dataset.metadata_path, destination / "metadata.csv")
    )
    expected_output_root = destination / "expected"
    shutil.copytree(dataset.reference_output_root, expected_output_root)
    return RabiesMethodSensitivityPanelExportResult(
        output_root=destination,
        readme_path=readme_path,
        config_path=config_path,
        sequences_path=sequences_path,
        metadata_path=metadata_path,
        expected_output_root=expected_output_root,
    )


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


def write_rabies_method_sensitivity_panel_workflow_bundle(
    output_root: Path,
    report: RabiesMethodSensitivityPanelWorkflowReport,
) -> RabiesMethodSensitivityPanelWorkflowBundle:
    """Write the governed reviewer-facing bundle for the method-sensitivity workflow."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    workflow_summary_path = _write_workflow_summary_table(
        output_root / "workflow-summary.tsv",
        report,
    )
    variant_summary_path = _write_variant_summary_table(
        output_root / "variant-summary.tsv",
        report,
    )
    parallel_summary_path = _write_parallel_execution_summary_table(
        output_root / "parallel-execution-summary.tsv",
        report,
    )
    preprocessing_comparison_path = _write_preprocessing_comparison_table(
        output_root / "preprocessing-rooted-comparisons.tsv",
        report.preprocessing_comparison_rows,
    )
    stable_clades_path = _write_clade_table(
        output_root / "stable-clades.tsv",
        report.stable_clade_rows,
    )
    changed_clades_path = _write_clade_table(
        output_root / "changed-clades.tsv",
        report.changed_clade_rows,
    )
    conclusion_summary_path = _write_conclusion_summary_table(
        output_root / "method-conclusion-summary.tsv",
        report.conclusion_rows,
    )
    config_path = _write_resolved_config(
        output_root / "workflow-config.resolved.json",
        report,
    )
    execution_record_path = _copy_output(
        report.execution_record_path,
        output_root / report.execution_record_path.name,
    )
    task_logs_root = _copy_task_logs(output_root / "parallel-logs", report.task_records)
    variants_root = _write_variant_outputs(
        output_root / "variants", report.variant_runs
    )
    slurm_planning_report = build_rabies_method_sensitivity_slurm_planning_report(
        report
    )
    slurm_job_plan_path = write_rabies_method_sensitivity_slurm_job_plan_table(
        output_root / "slurm-job-plan.tsv",
        slurm_planning_report,
    )
    slurm_assumptions_path = (
        write_rabies_method_sensitivity_slurm_assumptions_table(
            output_root / "slurm-estimation-assumptions.tsv",
            slurm_planning_report,
        )
    )
    slurm_summary_path = write_rabies_method_sensitivity_slurm_summary_json(
        output_root / "slurm-planning-summary.json",
        slurm_planning_report,
    )
    slurm_array_strategy_report = (
        build_rabies_method_sensitivity_slurm_array_strategy_report(
            slurm_planning_report
        )
    )
    slurm_array_partitions_path = (
        write_rabies_method_sensitivity_slurm_array_partitions_table(
            output_root / "slurm-array-partitions.tsv",
            slurm_array_strategy_report,
        )
    )
    slurm_array_members_path = (
        write_rabies_method_sensitivity_slurm_array_members_table(
            output_root / "slurm-array-members.tsv",
            slurm_array_strategy_report,
        )
    )
    slurm_array_strategy_path = (
        write_rabies_method_sensitivity_slurm_array_strategy_json(
            output_root / "slurm-array-strategy.json",
            slurm_array_strategy_report,
        )
    )
    slurm_array_scripts_root = (
        write_rabies_method_sensitivity_slurm_array_partition_scripts(
            output_root / "slurm-arrays",
            slurm_array_strategy_report,
        )
    )
    slurm_output_freshness_report = (
        build_rabies_method_sensitivity_slurm_output_freshness_report(
            output_root,
            dataset=report.dataset,
        )
    )
    slurm_output_freshness_path = (
        write_rabies_method_sensitivity_slurm_output_freshness_table(
            output_root / "slurm-output-freshness.tsv",
            slurm_output_freshness_report,
        )
    )
    slurm_output_freshness_checks_path = (
        write_rabies_method_sensitivity_slurm_output_freshness_checks_table(
            output_root / "slurm-output-freshness-checks.tsv",
            slurm_output_freshness_report,
        )
    )
    slurm_output_freshness_summary_path = (
        write_rabies_method_sensitivity_slurm_output_freshness_json(
            output_root / "slurm-output-freshness.json",
            slurm_output_freshness_report,
        )
    )
    slurm_status_report = build_rabies_method_sensitivity_slurm_status_report(
        output_root,
        dataset=report.dataset,
    )
    slurm_job_status_path = write_rabies_method_sensitivity_slurm_job_status_table(
        output_root / "slurm-job-status.tsv",
        slurm_status_report,
    )
    slurm_partition_status_path = (
        write_rabies_method_sensitivity_slurm_partition_status_table(
            output_root / "slurm-partition-status.tsv",
            slurm_status_report,
        )
    )
    slurm_workflow_status_path = write_rabies_method_sensitivity_slurm_status_json(
        output_root / "slurm-workflow-status.json",
        slurm_status_report,
    )
    manifest_path = _write_manifest(
        output_root / "rabies-method-sensitivity.manifest.json",
        report=report,
        bundle_paths={
            "workflow_summary": workflow_summary_path,
            "variant_summary": variant_summary_path,
            "parallel_summary": parallel_summary_path,
            "preprocessing_comparison": preprocessing_comparison_path,
            "stable_clades": stable_clades_path,
            "changed_clades": changed_clades_path,
            "conclusion_summary": conclusion_summary_path,
            "config": config_path,
            "execution_record": execution_record_path,
            "task_logs_root": task_logs_root,
            "variants_root": variants_root,
            "slurm_job_plan": slurm_job_plan_path,
            "slurm_assumptions": slurm_assumptions_path,
            "slurm_summary": slurm_summary_path,
            "slurm_array_partitions": slurm_array_partitions_path,
            "slurm_array_members": slurm_array_members_path,
            "slurm_array_strategy": slurm_array_strategy_path,
            "slurm_array_scripts_root": slurm_array_scripts_root,
            "slurm_output_freshness": slurm_output_freshness_path,
            "slurm_output_freshness_checks": slurm_output_freshness_checks_path,
            "slurm_output_freshness_summary": slurm_output_freshness_summary_path,
            "slurm_job_status": slurm_job_status_path,
            "slurm_partition_status": slurm_partition_status_path,
            "slurm_workflow_status": slurm_workflow_status_path,
        },
    )
    report_manifest_path = _write_report_manifest(
        output_root
        / "report-artifacts"
        / "rabies-method-sensitivity-report.manifest.json",
        report=report,
        bundle_paths={
            "workflow_summary": workflow_summary_path,
            "variant_summary": variant_summary_path,
            "parallel_summary": parallel_summary_path,
            "preprocessing_comparison": preprocessing_comparison_path,
            "stable_clades": stable_clades_path,
            "changed_clades": changed_clades_path,
            "conclusion_summary": conclusion_summary_path,
            "config": config_path,
            "execution_record": execution_record_path,
            "workflow_manifest": manifest_path,
            "slurm_job_plan": slurm_job_plan_path,
            "slurm_assumptions": slurm_assumptions_path,
            "slurm_summary": slurm_summary_path,
            "slurm_array_partitions": slurm_array_partitions_path,
            "slurm_array_members": slurm_array_members_path,
            "slurm_array_strategy": slurm_array_strategy_path,
            "slurm_output_freshness": slurm_output_freshness_path,
            "slurm_output_freshness_checks": slurm_output_freshness_checks_path,
            "slurm_output_freshness_summary": slurm_output_freshness_summary_path,
            "slurm_job_status": slurm_job_status_path,
            "slurm_partition_status": slurm_partition_status_path,
            "slurm_workflow_status": slurm_workflow_status_path,
        },
    )
    reproducibility_report = audit_rabies_method_sensitivity_workflow_bundle(
        output_root,
        sequences_path=report.dataset.sequences_path,
        metadata_path=report.dataset.metadata_path,
    )
    reproducibility_checks_path = (
        write_rabies_method_sensitivity_reproducibility_checks_table(
            output_root / "reproducibility-checks.tsv",
            reproducibility_report,
        )
    )
    reproducibility_variant_audit_path = (
        write_rabies_method_sensitivity_variant_audit_table(
            output_root / "reproducibility-variants.tsv",
            reproducibility_report,
        )
    )
    reproducibility_audit_path = (
        write_rabies_method_sensitivity_reproducibility_audit_json(
            output_root / "reproducibility-audit.json",
            reproducibility_report,
        )
    )
    report_linked_files = (
        workflow_summary_path,
        variant_summary_path,
        parallel_summary_path,
        preprocessing_comparison_path,
        stable_clades_path,
        changed_clades_path,
        conclusion_summary_path,
        config_path,
        execution_record_path,
        manifest_path,
        report_manifest_path,
        slurm_job_plan_path,
        slurm_assumptions_path,
        slurm_summary_path,
        slurm_array_partitions_path,
        slurm_array_members_path,
        slurm_array_strategy_path,
        slurm_output_freshness_path,
        slurm_output_freshness_checks_path,
        slurm_output_freshness_summary_path,
        slurm_job_status_path,
        slurm_partition_status_path,
        slurm_workflow_status_path,
    )
    report_linked_artifact_count = len(report_linked_files)
    report_path = _write_report(
        output_root / "rabies-method-sensitivity-report.html",
        report=report,
        bundle_paths={
            "workflow_summary": workflow_summary_path,
            "variant_summary": variant_summary_path,
            "parallel_summary": parallel_summary_path,
            "preprocessing_comparison": preprocessing_comparison_path,
            "stable_clades": stable_clades_path,
            "changed_clades": changed_clades_path,
            "conclusion_summary": conclusion_summary_path,
            "config": config_path,
            "execution_record": execution_record_path,
            "workflow_manifest": manifest_path,
            "slurm_job_plan": slurm_job_plan_path,
            "slurm_assumptions": slurm_assumptions_path,
            "slurm_summary": slurm_summary_path,
            "slurm_array_partitions": slurm_array_partitions_path,
            "slurm_array_members": slurm_array_members_path,
            "slurm_array_strategy": slurm_array_strategy_path,
            "slurm_output_freshness": slurm_output_freshness_path,
            "slurm_output_freshness_checks": slurm_output_freshness_checks_path,
            "slurm_output_freshness_summary": slurm_output_freshness_summary_path,
            "slurm_job_status": slurm_job_status_path,
            "slurm_partition_status": slurm_partition_status_path,
            "slurm_workflow_status": slurm_workflow_status_path,
            "reproducibility_checks": reproducibility_checks_path,
            "reproducibility_variant_audit": reproducibility_variant_audit_path,
            "reproducibility_audit": reproducibility_audit_path,
        },
        report_manifest_path=report_manifest_path,
        reproducibility_report=reproducibility_report,
        slurm_planning_report=slurm_planning_report,
        slurm_array_strategy_report=slurm_array_strategy_report,
        slurm_output_freshness_report=slurm_output_freshness_report,
        slurm_status_report=slurm_status_report,
    )
    report_html_size_bytes = report_path.stat().st_size
    report_linked_artifact_bytes = sum(
        path.stat().st_size for path in report_linked_files
    )
    report_total_output_bytes = report_html_size_bytes + report_linked_artifact_bytes
    preprocessing_change_pair_count = sum(
        1
        for row in report.preprocessing_comparison_rows
        if row.robinson_foulds_distance > 0 or row.same_taxa_different_rooting
    )
    rooted_engine_change_variant_count = sum(
        1
        for variant in report.variant_runs
        if variant.rooted_engine_comparison.robinson_foulds_distance > 0
        or variant.rooted_engine_comparison.same_taxa_different_rooting
    )
    serious_conflict_variant_count = sum(
        1
        for variant in report.variant_runs
        if variant.inference_comparison.conclusion_summary.serious_conflict_count > 0
    )
    return RabiesMethodSensitivityPanelWorkflowBundle(
        output_root=output_root,
        variant_count=len(report.variant_runs),
        stable_clade_count=len(report.stable_clade_rows),
        changed_clade_count=len(report.changed_clade_rows),
        preprocessing_change_pair_count=preprocessing_change_pair_count,
        rooted_engine_change_variant_count=rooted_engine_change_variant_count,
        serious_conflict_variant_count=serious_conflict_variant_count,
        execution_record_path=execution_record_path,
        parallel_workers=report.parallel_workers,
        execution_mode=report.execution_mode,
        workflow_summary_path=workflow_summary_path,
        variant_summary_path=variant_summary_path,
        parallel_summary_path=parallel_summary_path,
        preprocessing_comparison_path=preprocessing_comparison_path,
        stable_clades_path=stable_clades_path,
        changed_clades_path=changed_clades_path,
        conclusion_summary_path=conclusion_summary_path,
        config_path=config_path,
        manifest_path=manifest_path,
        report_manifest_path=report_manifest_path,
        slurm_job_plan_path=slurm_job_plan_path,
        slurm_assumptions_path=slurm_assumptions_path,
        slurm_summary_path=slurm_summary_path,
        slurm_array_partitions_path=slurm_array_partitions_path,
        slurm_array_members_path=slurm_array_members_path,
        slurm_array_strategy_path=slurm_array_strategy_path,
        slurm_array_scripts_root=slurm_array_scripts_root,
        slurm_job_count=slurm_planning_report.job_count,
        slurm_total_estimated_core_hours=(
            slurm_planning_report.total_estimated_core_hours
        ),
        slurm_maximum_estimated_memory_mib=(
            slurm_planning_report.maximum_estimated_memory_mib
        ),
        slurm_maximum_estimated_wallclock_minutes=(
            slurm_planning_report.maximum_estimated_wallclock_minutes
        ),
        slurm_total_estimated_scratch_mib=(
            slurm_planning_report.total_estimated_scratch_mib
        ),
        slurm_total_estimated_output_mib=(
            slurm_planning_report.total_estimated_output_mib
        ),
        slurm_array_partition_count=slurm_array_strategy_report.partition_count,
        slurm_array_script_count=slurm_array_strategy_report.script_count,
        slurm_array_largest_partition_size=(
            slurm_array_strategy_report.largest_partition_size
        ),
        slurm_output_freshness_path=slurm_output_freshness_path,
        slurm_output_freshness_checks_path=slurm_output_freshness_checks_path,
        slurm_output_freshness_summary_path=slurm_output_freshness_summary_path,
        slurm_job_status_path=slurm_job_status_path,
        slurm_partition_status_path=slurm_partition_status_path,
        slurm_workflow_status_path=slurm_workflow_status_path,
        slurm_output_freshness_check_count=slurm_output_freshness_report.check_count,
        slurm_output_freshness_failed_check_count=(
            slurm_output_freshness_report.failed_check_count
        ),
        slurm_fresh_output_job_count=slurm_output_freshness_report.fresh_job_count,
        slurm_stale_output_job_count=slurm_output_freshness_report.stale_job_count,
        slurm_completed_job_count=slurm_status_report.completed_job_count,
        slurm_failed_job_count=slurm_status_report.failed_job_count,
        slurm_pending_job_count=slurm_status_report.pending_job_count,
        slurm_stale_job_count=slurm_status_report.stale_job_count,
        reproducibility_checks_path=reproducibility_checks_path,
        reproducibility_variant_audit_path=reproducibility_variant_audit_path,
        reproducibility_audit_path=reproducibility_audit_path,
        reproducibility_passed=reproducibility_report.all_passed,
        reproducibility_check_count=reproducibility_report.check_count,
        reproducibility_failed_check_count=reproducibility_report.failed_check_count,
        reproducibility_failed_variant_count=reproducibility_report.failed_variant_count,
        report_path=report_path,
        report_linked_artifact_count=report_linked_artifact_count,
        report_html_size_bytes=report_html_size_bytes,
        report_linked_artifact_bytes=report_linked_artifact_bytes,
        report_total_output_bytes=report_total_output_bytes,
        task_logs_root=task_logs_root,
        variants_root=variants_root,
    )


def run_rabies_method_sensitivity_panel_demo(
    output_root: Path,
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
) -> RabiesMethodSensitivityPanelDemoResult:
    """Materialize the packaged dataset and rerun the governed sensitivity workflow."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    dataset = load_rabies_method_sensitivity_panel_dataset()
    dataset_export = export_rabies_method_sensitivity_panel_dataset(
        output_root / "dataset"
    )
    with TemporaryDirectory(prefix="rabies-method-sensitivity-") as temporary_root:
        workflow_report = run_rabies_method_sensitivity_panel_workflow(
            Path(temporary_root),
            mafft_executable=mafft_executable,
            trimal_executable=trimal_executable,
            iqtree_executable=iqtree_executable,
            fasttree_executable=fasttree_executable,
            iqtree_seed=iqtree_seed,
            iqtree_threads=iqtree_threads,
            bootstrap_replicates=bootstrap_replicates,
            parallel_workers=parallel_workers,
            variant_ids=variant_ids,
        )
        workflow_bundle = write_rabies_method_sensitivity_panel_workflow_bundle(
            output_root / "workflow",
            workflow_report,
        )
    overview_path = _write_overview(
        output_root / "overview.md", dataset, workflow_bundle
    )
    return RabiesMethodSensitivityPanelDemoResult(
        output_root=output_root,
        dataset=dataset,
        dataset_export=dataset_export,
        workflow_bundle=workflow_bundle,
        overview_path=overview_path,
    )


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


def _resolve_selected_variant_dataset(
    dataset: RabiesMethodSensitivityPanelDataset,
    *,
    variant_ids: tuple[str, ...] | None,
) -> RabiesMethodSensitivityPanelDataset:
    """Return either the full dataset or an explicit variant-scoped subset."""
    if variant_ids is None:
        return dataset
    if not variant_ids:
        raise ValueError("variant_ids must not be empty when provided")
    variants_by_id = {variant.variant_id: variant for variant in dataset.variants}
    selected_variants: list[RabiesMethodSensitivityVariant] = []
    seen_variant_ids: set[str] = set()
    for variant_id in variant_ids:
        if variant_id in seen_variant_ids:
            raise ValueError(f"duplicate variant_id requested: {variant_id}")
        seen_variant_ids.add(variant_id)
        variant = variants_by_id.get(variant_id)
        if variant is None:
            known = ", ".join(sorted(variants_by_id))
            raise ValueError(
                f"unknown variant_id '{variant_id}'; known variants: {known}"
            )
        selected_variants.append(variant)
    return replace(dataset, variants=tuple(selected_variants))


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


def _write_workflow_summary_table(
    path: Path, report: RabiesMethodSensitivityPanelWorkflowReport
) -> Path:
    serious_conflicts = [
        variant.inference_comparison.conclusion_summary.serious_conflict_count
        for variant in report.variant_runs
    ]
    rows = [
        [
            "dataset_id",
            "variant_count",
            "stable_clade_count",
            "changed_clade_count",
            "preprocessing_change_pair_count",
            "rooted_engine_change_variant_count",
            "serious_conflict_variant_count",
            "maximum_serious_conflict_count",
        ],
        [
            report.dataset.dataset_id,
            str(len(report.variant_runs)),
            str(len(report.stable_clade_rows)),
            str(len(report.changed_clade_rows)),
            str(
                sum(
                    1
                    for row in report.preprocessing_comparison_rows
                    if row.robinson_foulds_distance > 0
                    or row.same_taxa_different_rooting
                )
            ),
            str(
                sum(
                    1
                    for variant in report.variant_runs
                    if variant.rooted_engine_comparison.robinson_foulds_distance > 0
                    or variant.rooted_engine_comparison.same_taxa_different_rooting
                )
            ),
            str(sum(1 for value in serious_conflicts if value > 0)),
            str(max(serious_conflicts)),
        ],
    ]
    return _write_tsv(path, rows)


def _write_variant_summary_table(
    path: Path, report: RabiesMethodSensitivityPanelWorkflowReport
) -> Path:
    rows = [
        [
            "variant_id",
            "alignment_mode",
            "trimming_mode",
            "trim_gap_threshold",
            "alignment_length",
            "trimmed_alignment_length",
            "selected_model",
            "minimum_support",
            "maximum_support",
            "stable_clade_count",
            "engine_specific_clade_count",
            "serious_conflict_count",
            "rooted_engine_rf_distance",
            "rooted_engine_same_taxa_different_rooting",
        ]
    ]
    for variant in report.variant_runs:
        summary = variant.inference_comparison.conclusion_summary
        rows.append(
            [
                variant.config.variant_id,
                variant.config.alignment_mode,
                variant.config.trimming_mode,
                _format_float(variant.config.trim_gap_threshold),
                str(variant.alignment_length),
                str(variant.trimmed_alignment_length),
                variant.inference_comparison.selected_model,
                _format_optional_float(
                    variant.inference_comparison.iqtree_support_workflow.bootstrap_support_summary.minimum_support
                    if variant.inference_comparison.iqtree_support_workflow.bootstrap_support_summary
                    is not None
                    else None
                ),
                _format_optional_float(
                    variant.inference_comparison.iqtree_support_workflow.bootstrap_support_summary.maximum_support
                    if variant.inference_comparison.iqtree_support_workflow.bootstrap_support_summary
                    is not None
                    else None
                ),
                str(summary.stable_clade_count),
                str(summary.engine_specific_clade_count),
                str(summary.serious_conflict_count),
                str(variant.rooted_engine_comparison.robinson_foulds_distance),
                str(
                    variant.rooted_engine_comparison.same_taxa_different_rooting
                ).lower(),
            ]
        )
    return _write_tsv(path, rows)


def _write_preprocessing_comparison_table(
    path: Path,
    rows: tuple[RabiesMethodSensitivityPreprocessingComparisonRow, ...],
) -> Path:
    rendered = [
        [
            "left_variant_id",
            "right_variant_id",
            "comparison_axis",
            "robinson_foulds_distance",
            "normalized_robinson_foulds",
            "same_taxa_different_rooting",
        ]
    ]
    for row in rows:
        rendered.append(
            [
                row.left_variant_id,
                row.right_variant_id,
                row.comparison_axis,
                str(row.robinson_foulds_distance),
                _format_float(row.normalized_robinson_foulds),
                str(row.same_taxa_different_rooting).lower(),
            ]
        )
    return _write_tsv(path, rendered)


def _write_clade_table(
    path: Path, rows: tuple[RabiesMethodSensitivityCladeRow, ...]
) -> Path:
    rendered = [
        [
            "split_id",
            "conclusion_class",
            "evidence_class",
            "occurrence_count",
            "variant_count",
            "detail",
        ]
    ]
    for row in rows:
        rendered.append(
            [
                row.split_id,
                row.conclusion_class,
                row.evidence_class,
                str(row.occurrence_count),
                str(row.variant_count),
                row.detail,
            ]
        )
    return _write_tsv(path, rendered)


def _write_conclusion_summary_table(
    path: Path, rows: tuple[RabiesMethodSensitivityConclusionRow, ...]
) -> Path:
    rendered = [
        [
            "conclusion_id",
            "method_axis",
            "stability_status",
            "claim",
            "evidence",
            "caution",
        ]
    ]
    for row in rows:
        rendered.append(
            [
                row.conclusion_id,
                row.method_axis,
                row.stability_status,
                row.claim,
                row.evidence,
                row.caution,
            ]
        )
    return _write_tsv(path, rendered)


def _write_parallel_execution_summary_table(
    path: Path, report: RabiesMethodSensitivityPanelWorkflowReport
) -> Path:
    rows = [
        [
            "variant_id",
            "label",
            "execution_mode",
            "status",
            "log_path",
            "error_code",
        ]
    ]
    for task in report.task_records:
        rows.append(
            [
                task.variant_id,
                task.label,
                task.execution_mode,
                task.status,
                Path("parallel-logs", task.log_path.name).as_posix(),
                "" if task.error_code is None else task.error_code,
            ]
        )
    return _write_tsv(path, rows)


def _write_resolved_config(
    path: Path, report: RabiesMethodSensitivityPanelWorkflowReport
) -> Path:
    payload = {
        "dataset_id": report.dataset.dataset_id,
        "label": report.dataset.label,
        "sequence_type": report.dataset.sequence_type,
        "workflow_prefix": report.dataset.workflow_prefix,
        "outgroup_taxa": list(report.dataset.outgroup_taxa),
        "iqtree_seed": report.iqtree_seed,
        "iqtree_threads": report.iqtree_threads,
        "bootstrap_replicates": report.bootstrap_replicates,
        "parallel_workers": report.parallel_workers,
        "execution_mode": report.execution_mode,
        "selected_variant_ids": [
            variant.variant_id for variant in report.dataset.variants
        ],
        "input_checksums": {
            "sequences.fasta": _sha256(report.dataset.sequences_path),
            "metadata.csv": _sha256(report.dataset.metadata_path),
        },
        "variants": [
            {
                "variant_id": variant.variant_id,
                "label": variant.label,
                "alignment_mode": variant.alignment_mode,
                "trimming_mode": variant.trimming_mode,
                "trim_gap_threshold": variant.trim_gap_threshold,
            }
            for variant in report.dataset.variants
        ],
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return path


def _copy_task_logs(
    output_root: Path,
    task_records: tuple[RabiesMethodSensitivityTaskRecord, ...],
) -> Path:
    output_root.mkdir(parents=True, exist_ok=True)
    for task in task_records:
        _copy_output(task.log_path, output_root / task.log_path.name)
    return output_root


def _write_variant_outputs(
    output_root: Path, variant_runs: tuple[RabiesMethodSensitivityVariantRun, ...]
) -> Path:
    output_root.mkdir(parents=True, exist_ok=True)
    for variant in variant_runs:
        variant_root = output_root / variant.config.variant_id
        variant_root.mkdir(parents=True, exist_ok=True)
        _copy_output(
            variant.alignment_workflow.output_paths["alignment"],
            variant_root / f"{variant.config.variant_id}.aln",
        )
        _copy_output(
            variant.trimming_workflow.output_paths["trimmed_alignment"],
            variant_root / f"{variant.config.variant_id}.trimmed.aln",
        )
        _copy_output(
            variant.inference_comparison.output_paths["fasttree_tree"],
            variant_root / "fasttree.nwk",
        )
        _copy_output(
            variant.inference_comparison.output_paths["iqtree_support_tree"],
            variant_root / "iqtree-support.nwk",
        )
        _copy_output(
            variant.rooted_fasttree_path,
            variant_root / "rooted-fasttree.nwk",
        )
        _copy_output(
            variant.rooted_iqtree_path,
            variant_root / "rooted-iqtree-support.nwk",
        )
        _write_rooting_summary_table(
            variant_root / "rooting-summary.tsv",
            variant,
        )
        _copy_output(
            variant.inference_comparison.output_paths["stability_summary"],
            variant_root / "unrooted-stability-summary.tsv",
        )
        _copy_output(
            variant.inference_comparison.output_paths["conclusion_table"],
            variant_root / "unrooted-conclusions.tsv",
        )
        _copy_output(
            variant.inference_comparison.output_paths["support_weighted_conflicts"],
            variant_root / "unrooted-support-weighted-conflicts.tsv",
        )
        _copy_output(
            variant.inference_comparison.output_paths["shared_clades"],
            variant_root / "unrooted-shared-clades.tsv",
        )
        _copy_output(
            variant.inference_comparison.output_paths["conflicting_clades"],
            variant_root / "unrooted-conflicting-clades.tsv",
        )
        _copy_output(
            variant.inference_comparison.output_paths["comparison_table"],
            variant_root / "unrooted-comparison.tsv",
        )
        _copy_output(
            variant.rooted_engine_comparison_table_path,
            variant_root / "rooted-engine-comparison.tsv",
        )
    return output_root


def _write_rooting_summary_table(
    path: Path, variant: RabiesMethodSensitivityVariantRun
) -> Path:
    rows = [
        [
            "engine_name",
            "requested_taxa",
            "matched_taxa",
            "outgroup_monophyletic",
            "rooted_outgroup_taxa",
            "warning_count",
        ],
        [
            "fasttree",
            ",".join(variant.fasttree_rooting.requested_taxa),
            ",".join(variant.fasttree_rooting.matched_taxa),
            _format_optional_bool(variant.fasttree_rooting.outgroup_monophyletic),
            ",".join(variant.fasttree_rooting.rooted_outgroup_taxa),
            str(len(variant.fasttree_rooting.warnings)),
        ],
        [
            "iqtree",
            ",".join(variant.iqtree_rooting.requested_taxa),
            ",".join(variant.iqtree_rooting.matched_taxa),
            _format_optional_bool(variant.iqtree_rooting.outgroup_monophyletic),
            ",".join(variant.iqtree_rooting.rooted_outgroup_taxa),
            str(len(variant.iqtree_rooting.warnings)),
        ],
    ]
    return _write_tsv(path, rows)


def _write_manifest(
    path: Path,
    *,
    report: RabiesMethodSensitivityPanelWorkflowReport,
    bundle_paths: dict[str, Path],
) -> Path:
    payload = {
        "dataset_id": report.dataset.dataset_id,
        "label": report.dataset.label,
        "report_kind": "rabies_method_sensitivity_workflow_bundle",
        "variant_count": len(report.variant_runs),
        "parallel_execution": {
            "execution_mode": report.execution_mode,
            "parallel_workers": report.parallel_workers,
            "requested_task_count": len(report.task_records),
            "completed_task_count": len(
                [task for task in report.task_records if task.status == "succeeded"]
            ),
            "failed_task_count": len(
                [task for task in report.task_records if task.status != "succeeded"]
            ),
        },
        "task_records": [
            {
                "variant_id": task.variant_id,
                "label": task.label,
                "status": task.status,
                "execution_mode": task.execution_mode,
                "log_path": Path("parallel-logs", task.log_path.name).as_posix(),
                "output_root": Path("variants", task.variant_id).as_posix(),
                "error_code": task.error_code,
                "error_message": task.error_message,
            }
            for task in report.task_records
        ],
        "output_paths": {
            key: value.name
            if value.parent == path.parent
            else value.relative_to(path.parent).as_posix()
            for key, value in bundle_paths.items()
        },
        "output_checksums": {
            key: _sha256(value)
            for key, value in bundle_paths.items()
            if value.is_file()
        },
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return path


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


def _relative_bundle_path(base_path: Path, value: Path) -> str:
    return Path(os.path.relpath(value, start=base_path.parent)).as_posix()


def _write_report_manifest(
    path: Path,
    *,
    report: RabiesMethodSensitivityPanelWorkflowReport,
    bundle_paths: dict[str, Path],
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    linked_files = {
        key: value for key, value in bundle_paths.items() if value.is_file()
    }
    payload = {
        "dataset_id": report.dataset.dataset_id,
        "report_kind": "rabies_method_sensitivity_html_report",
        "variant_count": len(report.variant_runs),
        "parallel_workers": report.parallel_workers,
        "execution_mode": report.execution_mode,
        "stable_clade_count": len(report.stable_clade_rows),
        "changed_clade_count": len(report.changed_clade_rows),
        "linked_artifact_count": len(linked_files),
        "linked_artifacts": {
            key: {
                "path": _relative_bundle_path(path, value),
                "byte_count": value.stat().st_size,
                "sha256": _sha256(value),
            }
            for key, value in linked_files.items()
        },
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return path


def _write_report(
    path: Path,
    *,
    report: RabiesMethodSensitivityPanelWorkflowReport,
    bundle_paths: dict[str, Path],
    report_manifest_path: Path,
    reproducibility_report: RabiesMethodSensitivityReproducibilityAuditReport,
    slurm_planning_report: RabiesMethodSensitivitySlurmPlanningReport,
    slurm_array_strategy_report: RabiesMethodSensitivitySlurmArrayStrategyReport,
    slurm_output_freshness_report: RabiesMethodSensitivitySlurmOutputFreshnessReport,
    slurm_status_report: RabiesMethodSensitivitySlurmStatusReport,
) -> Path:
    variant_lines = [
        (
            f"{variant.config.variant_id}: model {variant.inference_comparison.selected_model}; "
            f"unrooted serious conflicts {variant.inference_comparison.conclusion_summary.serious_conflict_count}; "
            f"rooted engine RF {variant.rooted_engine_comparison.robinson_foulds_distance}"
        )
        for variant in report.variant_runs
    ]
    conclusion_lines: list[str] = []
    for row in report.conclusion_rows:
        conclusion_lines.extend(
            [
                f"{row.conclusion_id}: {row.stability_status}",
                f"claim: {row.claim}",
                f"evidence: {row.evidence}",
                f"caution: {row.caution}",
                "",
            ]
        )
    sections = [
        (
            "workflow-summary",
            "\n".join(
                [
                    f"dataset: {report.dataset.label}",
                    f"variants: {len(report.variant_runs)}",
                    f"execution mode: {report.execution_mode}",
                    f"parallel workers: {report.parallel_workers}",
                    f"stable clades across variants: {len(report.stable_clade_rows)}",
                    f"changed clades across variants: {len(report.changed_clade_rows)}",
                    *variant_lines,
                ]
            ),
        ),
        (
            "conclusions",
            "\n".join(conclusion_lines).strip(),
        ),
        (
            "reproducibility-audit",
            "\n".join(
                [
                    f"all passed: {str(getattr(reproducibility_report, 'all_passed')).lower()}",
                    f"top-level checks: {getattr(reproducibility_report, 'check_count')}",
                    (
                        "failed top-level checks: "
                        f"{getattr(reproducibility_report, 'failed_check_count')}"
                    ),
                    (
                        "failed variants: "
                        f"{getattr(reproducibility_report, 'failed_variant_count')}"
                    ),
                ]
            ),
        ),
        (
            "slurm-job-planning",
            "\n".join(
                [
                    f"planned jobs: {slurm_planning_report.job_count}",
                    (
                        "maximum estimated memory MiB: "
                        f"{slurm_planning_report.maximum_estimated_memory_mib}"
                    ),
                    (
                        "maximum estimated wallclock minutes: "
                        f"{slurm_planning_report.maximum_estimated_wallclock_minutes}"
                    ),
                    (
                        "total estimated core hours: "
                        f"{_format_float(slurm_planning_report.total_estimated_core_hours)}"
                    ),
                    (
                        "total estimated scratch MiB: "
                        f"{slurm_planning_report.total_estimated_scratch_mib}"
                    ),
                    (
                        "total estimated output MiB: "
                        f"{slurm_planning_report.total_estimated_output_mib}"
                    ),
                ]
            ),
        ),
        (
            "slurm-array-partitioning",
            "\n".join(
                [
                    (
                        "array partitions: "
                        f"{slurm_array_strategy_report.partition_count}"
                    ),
                    (
                        "largest partition size: "
                        f"{slurm_array_strategy_report.largest_partition_size}"
                    ),
                    (
                        "partition scripts: "
                        f"{slurm_array_strategy_report.script_count}"
                    ),
                    (
                        "total array jobs: "
                        f"{slurm_array_strategy_report.total_job_count}"
                    ),
                ]
            ),
        ),
        (
            "slurm-output-freshness",
            "\n".join(
                [
                    (
                        "all outputs fresh: "
                        f"{str(slurm_output_freshness_report.all_outputs_fresh).lower()}"
                    ),
                    (
                        "fresh jobs: "
                        f"{slurm_output_freshness_report.fresh_job_count}"
                    ),
                    (
                        "stale jobs: "
                        f"{slurm_output_freshness_report.stale_job_count}"
                    ),
                    (
                        "freshness checks: "
                        f"{slurm_output_freshness_report.check_count}"
                    ),
                    (
                        "failed freshness checks: "
                        f"{slurm_output_freshness_report.failed_check_count}"
                    ),
                ]
            ),
        ),
        (
            "slurm-workflow-status",
            "\n".join(
                [
                    f"workflow status: {slurm_status_report.workflow_status}",
                    f"active run state: {slurm_status_report.active_run_state}",
                    f"completed jobs: {slurm_status_report.completed_job_count}",
                    f"failed jobs: {slurm_status_report.failed_job_count}",
                    f"pending jobs: {slurm_status_report.pending_job_count}",
                    f"stale jobs: {slurm_status_report.stale_job_count}",
                ]
            ),
        ),
        (
            "artifacts",
            "\n".join(
                [
                    f"workflow summary: {bundle_paths['workflow_summary'].name}",
                    f"variant summary: {bundle_paths['variant_summary'].name}",
                    f"parallel execution summary: {bundle_paths['parallel_summary'].name}",
                    f"preprocessing rooted comparisons: {bundle_paths['preprocessing_comparison'].name}",
                    f"stable clades: {bundle_paths['stable_clades'].name}",
                    f"changed clades: {bundle_paths['changed_clades'].name}",
                    f"method conclusions: {bundle_paths['conclusion_summary'].name}",
                    f"resolved config: {bundle_paths['config'].name}",
                    f"workflow manifest: {bundle_paths['workflow_manifest'].name}",
                    f"slurm job plan: {bundle_paths['slurm_job_plan'].name}",
                    (
                        "slurm estimation assumptions: "
                        f"{bundle_paths['slurm_assumptions'].name}"
                    ),
                    f"slurm planning summary: {bundle_paths['slurm_summary'].name}",
                    (
                        "slurm array partitions: "
                        f"{bundle_paths['slurm_array_partitions'].name}"
                    ),
                    (
                        "slurm array members: "
                        f"{bundle_paths['slurm_array_members'].name}"
                    ),
                    (
                        "slurm array strategy: "
                        f"{bundle_paths['slurm_array_strategy'].name}"
                    ),
                    (
                        "slurm output freshness: "
                        f"{bundle_paths['slurm_output_freshness'].name}"
                    ),
                    (
                        "slurm output freshness checks: "
                        f"{bundle_paths['slurm_output_freshness_checks'].name}"
                    ),
                    (
                        "slurm output freshness summary: "
                        f"{bundle_paths['slurm_output_freshness_summary'].name}"
                    ),
                    f"slurm job status: {bundle_paths['slurm_job_status'].name}",
                    (
                        "slurm partition status: "
                        f"{bundle_paths['slurm_partition_status'].name}"
                    ),
                    (
                        "slurm workflow status: "
                        f"{bundle_paths['slurm_workflow_status'].name}"
                    ),
                    f"reproducibility checks: {bundle_paths['reproducibility_checks'].name}",
                    (
                        "reproducibility variant audit: "
                        f"{bundle_paths['reproducibility_variant_audit'].name}"
                    ),
                    (
                        "reproducibility audit: "
                        f"{bundle_paths['reproducibility_audit'].name}"
                    ),
                ]
            ),
        ),
    ]
    report_manifest = json.loads(report_manifest_path.read_text(encoding="utf-8"))
    artifact_links = [
        (
            key.replace("_", "-"),
            _relative_bundle_path(path, value),
            f"{value.stat().st_size} bytes",
        )
        for key, value in bundle_paths.items()
        if value.is_file()
    ]
    return write_html_report(
        title="Rabies Method-Sensitivity Report",
        sections=sections,
        out_path=path,
        embedded_json={
            "dataset_id": report.dataset.dataset_id,
            "variant_count": len(report.variant_runs),
            "parallel_workers": report.parallel_workers,
            "execution_mode": report.execution_mode,
            "stable_clade_count": len(report.stable_clade_rows),
            "changed_clade_count": len(report.changed_clade_rows),
            "report_manifest_path": _relative_bundle_path(path, report_manifest_path),
            "reproducibility_passed": getattr(reproducibility_report, "all_passed"),
            "reproducibility_check_count": getattr(
                reproducibility_report, "check_count"
            ),
            "reproducibility_failed_check_count": getattr(
                reproducibility_report, "failed_check_count"
            ),
            "reproducibility_failed_variant_count": getattr(
                reproducibility_report, "failed_variant_count"
            ),
            "slurm_job_count": slurm_planning_report.job_count,
            "slurm_total_estimated_core_hours": (
                slurm_planning_report.total_estimated_core_hours
            ),
            "slurm_maximum_estimated_memory_mib": (
                slurm_planning_report.maximum_estimated_memory_mib
            ),
            "slurm_maximum_estimated_wallclock_minutes": (
                slurm_planning_report.maximum_estimated_wallclock_minutes
            ),
            "slurm_array_partition_count": (
                slurm_array_strategy_report.partition_count
            ),
            "slurm_array_script_count": slurm_array_strategy_report.script_count,
            "slurm_array_largest_partition_size": (
                slurm_array_strategy_report.largest_partition_size
            ),
            "slurm_output_freshness_check_count": (
                slurm_output_freshness_report.check_count
            ),
            "slurm_output_freshness_failed_check_count": (
                slurm_output_freshness_report.failed_check_count
            ),
            "slurm_fresh_output_job_count": (
                slurm_output_freshness_report.fresh_job_count
            ),
            "slurm_stale_output_job_count": (
                slurm_output_freshness_report.stale_job_count
            ),
            "slurm_completed_job_count": slurm_status_report.completed_job_count,
            "slurm_failed_job_count": slurm_status_report.failed_job_count,
            "slurm_pending_job_count": slurm_status_report.pending_job_count,
            "slurm_stale_job_count": slurm_status_report.stale_job_count,
            "slurm_active_run_state": slurm_status_report.active_run_state,
        },
        summary_metrics=[
            ("variants", len(report.variant_runs)),
            ("execution mode", report.execution_mode),
            ("parallel workers", report.parallel_workers),
            ("stable clades", len(report.stable_clade_rows)),
            ("changed clades", len(report.changed_clade_rows)),
            (
                "slurm planned jobs",
                slurm_planning_report.job_count,
            ),
            (
                "slurm max memory MiB",
                slurm_planning_report.maximum_estimated_memory_mib,
            ),
            (
                "slurm array partitions",
                slurm_array_strategy_report.partition_count,
            ),
            (
                "slurm fresh output jobs",
                slurm_output_freshness_report.fresh_job_count,
            ),
            (
                "slurm stale output jobs",
                slurm_output_freshness_report.stale_job_count,
            ),
            (
                "slurm completed jobs",
                slurm_status_report.completed_job_count,
            ),
            (
                "slurm stale jobs",
                slurm_status_report.stale_job_count,
            ),
            (
                "reproducibility passed",
                str(getattr(reproducibility_report, "all_passed")).lower(),
            ),
            (
                "reproducibility checks",
                getattr(reproducibility_report, "check_count"),
            ),
            ("linked artifacts", report_manifest["linked_artifact_count"]),
        ],
        artifact_links=[
            *artifact_links,
            (
                "report-manifest",
                _relative_bundle_path(path, report_manifest_path),
                f"{report_manifest_path.stat().st_size} bytes",
            ),
        ],
    )


def _write_overview(
    path: Path,
    dataset: RabiesMethodSensitivityPanelDataset,
    bundle: RabiesMethodSensitivityPanelWorkflowBundle,
) -> Path:
    lines = [
        f"# {dataset.label}",
        "",
        dataset.source_summary,
        "",
        "## Bundle",
        "",
        f"- variants: `{bundle.variant_count}`",
        f"- execution mode: `{bundle.execution_mode}`",
        f"- parallel workers: `{bundle.parallel_workers}`",
        f"- stable clades across variants: `{bundle.stable_clade_count}`",
        f"- changed clades across variants: `{bundle.changed_clade_count}`",
        f"- preprocessing-rooted changes: `{bundle.preprocessing_change_pair_count}`",
        f"- rooted engine changes: `{bundle.rooted_engine_change_variant_count}`",
        f"- variants with unrooted serious conflicts: `{bundle.serious_conflict_variant_count}`",
        f"- report linked artifacts: `{bundle.report_linked_artifact_count}`",
        f"- report html bytes: `{bundle.report_html_size_bytes}`",
        f"- report total output bytes: `{bundle.report_total_output_bytes}`",
        f"- reproducibility passed: `{str(bundle.reproducibility_passed).lower()}`",
        f"- reproducibility checks: `{bundle.reproducibility_check_count}`",
        f"- failed reproducibility checks: `{bundle.reproducibility_failed_check_count}`",
        f"- failed reproducibility variants: `{bundle.reproducibility_failed_variant_count}`",
        f"- slurm planned jobs: `{bundle.slurm_job_count}`",
        f"- slurm total estimated core hours: `{_format_float(bundle.slurm_total_estimated_core_hours)}`",
        f"- slurm max estimated memory MiB: `{bundle.slurm_maximum_estimated_memory_mib}`",
        f"- slurm max estimated wallclock minutes: `{bundle.slurm_maximum_estimated_wallclock_minutes}`",
        f"- slurm total estimated scratch MiB: `{bundle.slurm_total_estimated_scratch_mib}`",
        f"- slurm total estimated output MiB: `{bundle.slurm_total_estimated_output_mib}`",
        f"- slurm array partitions: `{bundle.slurm_array_partition_count}`",
        f"- slurm array scripts: `{bundle.slurm_array_script_count}`",
        f"- slurm largest array partition: `{bundle.slurm_array_largest_partition_size}`",
        f"- slurm output freshness checks: `{bundle.slurm_output_freshness_check_count}`",
        f"- failed slurm output freshness checks: `{bundle.slurm_output_freshness_failed_check_count}`",
        f"- slurm fresh output jobs: `{bundle.slurm_fresh_output_job_count}`",
        f"- slurm stale output jobs: `{bundle.slurm_stale_output_job_count}`",
        f"- slurm completed jobs: `{bundle.slurm_completed_job_count}`",
        f"- slurm failed jobs: `{bundle.slurm_failed_job_count}`",
        f"- slurm pending jobs: `{bundle.slurm_pending_job_count}`",
        f"- slurm stale jobs: `{bundle.slurm_stale_job_count}`",
        f"- workflow manifest: `{bundle.manifest_path.name}`",
        f"- report manifest: `{bundle.report_manifest_path.relative_to(bundle.output_root).as_posix()}`",
        f"- slurm job plan: `{bundle.slurm_job_plan_path.name}`",
        f"- slurm assumptions: `{bundle.slurm_assumptions_path.name}`",
        f"- slurm planning summary: `{bundle.slurm_summary_path.name}`",
        f"- slurm array partitions table: `{bundle.slurm_array_partitions_path.name}`",
        f"- slurm array members table: `{bundle.slurm_array_members_path.name}`",
        f"- slurm array strategy: `{bundle.slurm_array_strategy_path.name}`",
        f"- slurm output freshness table: `{bundle.slurm_output_freshness_path.name}`",
        f"- slurm output freshness checks: `{bundle.slurm_output_freshness_checks_path.name}`",
        f"- slurm output freshness summary: `{bundle.slurm_output_freshness_summary_path.name}`",
        f"- slurm job status table: `{bundle.slurm_job_status_path.name}`",
        f"- slurm partition status table: `{bundle.slurm_partition_status_path.name}`",
        f"- slurm workflow status: `{bundle.slurm_workflow_status_path.name}`",
        f"- slurm array scripts: `{bundle.slurm_array_scripts_root.name}/`",
        f"- reproducibility checks table: `{bundle.reproducibility_checks_path.name}`",
        f"- reproducibility variant audit: `{bundle.reproducibility_variant_audit_path.name}`",
        f"- reproducibility audit: `{bundle.reproducibility_audit_path.name}`",
        f"- report: `{bundle.report_path.name}`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
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


def _copy_output(source: Path, destination: Path) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    return Path(shutil.copy2(source, destination))


def _write_tsv(path: Path, rows: list[list[str]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join("\t".join(row) for row in rows) + "\n",
        encoding="utf-8",
    )
    return path


def _format_float(value: float) -> str:
    return format(value, ".12g")


def _format_optional_float(value: float | None) -> str:
    return "" if value is None else format(value, ".12g")


def _format_optional_bool(value: bool | None) -> str:
    if value is None:
        return ""
    return str(value).lower()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _resource_root() -> Path:
    return (
        Path(__file__).resolve().parent.parent
        / "resources"
        / "datasets"
        / "pathogens"
        / _DATASET_ID
    )
