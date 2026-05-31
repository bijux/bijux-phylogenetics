from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path

from bijux_phylogenetics.core.manifest import build_run_manifest, write_run_manifest
from bijux_phylogenetics.evidence.provenance.method_tiers import (
    fasta_to_tree_method_tier,
)
from bijux_phylogenetics.io.fasta import write_fasta_alignment
from bijux_phylogenetics.io.fasta.records import (
    repair_fasta_input,
    validate_fasta_input,
)
from bijux_phylogenetics.phylo.alignment import (
    AlignmentAlphabet,
    FastaInputValidationReport,
    FastaRepairReport,
)
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError

from ...common import build_file_checksums, write_engine_manifest
from ...validation import (
    summarize_bootstrap_support_distribution,
    validate_model_selection_against_engine_outputs,
)
from ...validation.preflight import require_external_engine_surface
from ...workflows.alignment import (
    run_alignment_trimming,
    run_multiple_sequence_alignment,
)
from ...workflows.iqtree import (
    run_bootstrap_support_estimation,
    run_maximum_likelihood_tree_inference,
    run_model_selection,
)
from .artifact_outputs import (
    write_fasta_to_tree_log,
    write_fasta_to_tree_model_table,
    write_fasta_to_tree_support_table,
)
from .contracts import FastaToTreeWorkflowReport
from .row_builders import (
    build_fasta_to_tree_model_rows,
    build_fasta_to_tree_support_rows,
)
from .stage_fingerprints import build_stage_fingerprint
from .workflow_layout import _artifact_prefix, _copy_output, _final_output_paths


def run_fasta_to_tree_workflow(
    input_path: Path,
    *,
    out_dir: Path,
    prefix: str | None = None,
    sequence_type: AlignmentAlphabet | None = None,
    mafft_executable: str | Path = "mafft",
    alignment_mode: str = "auto",
    trimal_executable: str | Path = "trimal",
    trimming_mode: str = "gap-threshold",
    iqtree_executable: str | Path = "iqtree2",
    iqtree_seed: int = 1,
    iqtree_threads: int = 1,
    trim_gap_threshold: float = 0.1,
    bootstrap_replicates: int = 1000,
    normalize_identifiers: bool = False,
    remove_invalid_records: bool = False,
    resume: bool = False,
    timeout_seconds: float | None = None,
    incomplete_run_policy: str = "reject",
) -> FastaToTreeWorkflowReport:
    """Run alignment, trimming, model selection, ML inference, and bootstrap support in one workflow."""
    started_at = datetime.now(UTC)
    require_external_engine_surface(
        workflow_id="fasta-to-tree",
        summary="Canonical FASTA-to-tree workflow over MAFFT, trimAl, and IQ-TREE.",
        required_engines=("mafft", "trimal", "iqtree"),
        executables={
            "mafft": mafft_executable,
            "trimal": trimal_executable,
            "iqtree": iqtree_executable,
        },
    )
    workflow_prefix = input_path.stem if prefix is None else prefix
    engine_artifact_dir = out_dir / "engine-artifacts" / workflow_prefix
    input_validation = validate_fasta_input(input_path, sequence_type=sequence_type)
    has_input_blockers = bool(
        input_validation.duplicate_identifiers
        or input_validation.illegal_characters
        or input_validation.empty_sequences
    )
    prepared_input_path = input_path
    repaired_input_validation: FastaInputValidationReport | None = None
    input_repair: FastaRepairReport | None = None
    if has_input_blockers and not (normalize_identifiers or remove_invalid_records):
        raise InvalidAlignmentError(
            "fasta-to-tree input contains duplicate identifiers, empty sequences, "
            "or illegal characters; use normalize_identifiers or "
            "remove_invalid_records to repair the input explicitly"
        )
    if normalize_identifiers or remove_invalid_records:
        repaired_records, input_repair = repair_fasta_input(
            input_path,
            sequence_type=sequence_type,
            normalize_identifiers=normalize_identifiers,
            remove_invalid_records=remove_invalid_records,
        )
        prepared_input_path = _artifact_prefix(
            out_dir, workflow_prefix, "input-curation"
        ).with_suffix(".fasta")
        write_fasta_alignment(prepared_input_path, repaired_records)
        input_repair.output_path = prepared_input_path
        repaired_input_validation = validate_fasta_input(
            prepared_input_path,
            sequence_type=sequence_type,
        )
        if (
            repaired_input_validation.duplicate_identifiers
            or repaired_input_validation.illegal_characters
            or repaired_input_validation.empty_sequences
        ):
            raise InvalidAlignmentError(
                "fasta-to-tree repair policy left unresolved duplicate identifiers, "
                "illegal characters, or empty sequences in the prepared input"
            )
    effective_input_validation = (
        input_validation
        if repaired_input_validation is None
        else repaired_input_validation
    )
    inferred_sequence_type = (
        effective_input_validation.sequence_type_report.selected_type
        if sequence_type is None
        else sequence_type
    )
    if inferred_sequence_type in {None, "unknown"}:
        raise InvalidAlignmentError(
            "fasta-to-tree workflow could not resolve one supported raw sequence type: "
            f"{effective_input_validation.sequence_type_report.note}"
        )
    final_outputs = _final_output_paths(out_dir, workflow_prefix)

    alignment_workflow = run_multiple_sequence_alignment(
        prepared_input_path,
        _artifact_prefix(out_dir, workflow_prefix, "alignment").with_suffix(".aln"),
        executable=mafft_executable,
        mode=alignment_mode,
        resume=resume,
        timeout_seconds=timeout_seconds,
        incomplete_run_policy=incomplete_run_policy,
    )
    trimming_workflow = run_alignment_trimming(
        alignment_workflow.output_paths["alignment"],
        _artifact_prefix(out_dir, workflow_prefix, "trimming").with_suffix(
            ".trimmed.aln"
        ),
        executable=trimal_executable,
        mode=trimming_mode,
        gap_threshold=trim_gap_threshold,
        resume=resume,
        timeout_seconds=timeout_seconds,
        incomplete_run_policy=incomplete_run_policy,
    )
    model_selection_workflow = run_model_selection(
        trimming_workflow.output_paths["trimmed_alignment"],
        out_dir=_artifact_prefix(out_dir, workflow_prefix, "model-selection").parent,
        prefix="model-selection",
        executable=iqtree_executable,
        sequence_type=inferred_sequence_type,
        resume=resume,
        seed=iqtree_seed,
        threads=iqtree_threads,
        timeout_seconds=timeout_seconds,
        incomplete_run_policy=incomplete_run_policy,
    )
    if model_selection_workflow.selected_model is None:
        raise ValueError("model-selection workflow did not expose a selected model")
    maximum_likelihood_workflow = run_maximum_likelihood_tree_inference(
        trimming_workflow.output_paths["trimmed_alignment"],
        out_dir=_artifact_prefix(out_dir, workflow_prefix, "maximum-likelihood").parent,
        model=model_selection_workflow.selected_model,
        prefix="maximum-likelihood",
        executable=iqtree_executable,
        sequence_type=inferred_sequence_type,
        resume=resume,
        seed=iqtree_seed,
        threads=iqtree_threads,
        timeout_seconds=timeout_seconds,
        incomplete_run_policy=incomplete_run_policy,
    )
    bootstrap_workflow = run_bootstrap_support_estimation(
        trimming_workflow.output_paths["trimmed_alignment"],
        out_dir=_artifact_prefix(out_dir, workflow_prefix, "bootstrap-support").parent,
        model=model_selection_workflow.selected_model,
        replicates=bootstrap_replicates,
        prefix="bootstrap-support",
        executable=iqtree_executable,
        sequence_type=inferred_sequence_type,
        resume=resume,
        seed=iqtree_seed,
        threads=iqtree_threads,
        timeout_seconds=timeout_seconds,
        incomplete_run_policy=incomplete_run_policy,
    )

    model_validation = validate_model_selection_against_engine_outputs(
        model_selection_workflow.manifest_path
    )
    _copy_output(
        alignment_workflow.output_paths["alignment"], final_outputs["alignment"]
    )
    _copy_output(
        trimming_workflow.output_paths["trimmed_alignment"],
        final_outputs["trimmed_alignment"],
    )
    _copy_output(bootstrap_workflow.output_paths["support_tree"], final_outputs["tree"])
    model_rows = build_fasta_to_tree_model_rows(
        model_selection_workflow,
        validation=model_validation,
        sequence_type=inferred_sequence_type,
        alignment_path=final_outputs["alignment"],
        trimmed_alignment_path=final_outputs["trimmed_alignment"],
    )
    write_fasta_to_tree_model_table(
        final_outputs["model_table"],
        model_rows,
        root_dir=out_dir,
    )
    support_summary = summarize_bootstrap_support_distribution(final_outputs["tree"])
    support_rows = build_fasta_to_tree_support_rows(support_summary)
    write_fasta_to_tree_support_table(final_outputs["support_table"], support_rows)

    warnings = list(
        dict.fromkeys(
            alignment_workflow.run.warning_lines
            + trimming_workflow.run.warning_lines
            + model_selection_workflow.run.warning_lines
            + maximum_likelihood_workflow.run.warning_lines
            + bootstrap_workflow.run.warning_lines
            + support_summary.warnings
            + input_validation.warnings
            + ([] if input_repair is None else input_repair.warnings)
            + (
                []
                if repaired_input_validation is None
                else repaired_input_validation.warnings
            )
        )
    )
    notes = [
        "final tree path contains the bootstrap-supported inference tree",
        "engine-specific intermediate artifacts remain under engine-artifacts/",
        f"mafft alignment mode: {alignment_mode}",
        f"trimal trimming mode: {trimming_mode}",
        f"iqtree random seed: {iqtree_seed}",
        f"iqtree threads: {iqtree_threads}",
        "raw sequence type detection: "
        f"{effective_input_validation.sequence_type_report.detected_type} "
        f"({effective_input_validation.sequence_type_report.confidence})",
        effective_input_validation.sequence_type_report.note,
    ]
    if input_repair is not None:
        notes.append(
            f"prepared input FASTA written to {prepared_input_path} before alignment"
        )
    run_manifest_arguments = [
        str(input_path),
        "--out-dir",
        str(out_dir),
        "--prefix",
        workflow_prefix,
        "--sequence-type",
        inferred_sequence_type,
        "--alignment-mode",
        alignment_mode,
        "--trimming-mode",
        trimming_mode,
        "--trim-gap-threshold",
        format(trim_gap_threshold, ".12g"),
        "--iqtree-seed",
        str(iqtree_seed),
        "--iqtree-threads",
        str(iqtree_threads),
        "--bootstrap-replicates",
        str(bootstrap_replicates),
        "--resume",
        "true" if resume else "false",
        "--incomplete-run-policy",
        incomplete_run_policy,
    ]
    if timeout_seconds is not None:
        run_manifest_arguments.extend(
            ["--timeout-seconds", format(timeout_seconds, ".12g")]
        )
    report = FastaToTreeWorkflowReport(
        workflow="fasta-to-tree",
        input_path=input_path,
        prepared_input_path=prepared_input_path,
        out_dir=out_dir,
        prefix=workflow_prefix,
        sequence_type=inferred_sequence_type,
        selected_model=model_selection_workflow.selected_model,
        alignment_mode=alignment_mode,
        trimming_mode=trimming_mode,
        trim_gap_threshold=trim_gap_threshold,
        iqtree_seed=iqtree_seed,
        iqtree_threads=iqtree_threads,
        bootstrap_replicates=bootstrap_replicates,
        started_at_utc=started_at.replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z"),
        ended_at_utc="",
        runtime_seconds=0.0,
        engine_artifact_dir=engine_artifact_dir,
        manifest_path=final_outputs["manifest"],
        run_manifest_path=final_outputs["run_manifest"],
        output_paths=final_outputs,
        config={
            "sequence_type": inferred_sequence_type,
            "alignment_mode": alignment_mode,
            "trimming_mode": trimming_mode,
            "trim_gap_threshold": trim_gap_threshold,
            "iqtree_seed": iqtree_seed,
            "iqtree_threads": iqtree_threads,
            "bootstrap_replicates": bootstrap_replicates,
            "timeout_seconds": timeout_seconds,
            "resume": resume,
            "incomplete_run_policy": incomplete_run_policy,
        },
        commands={
            "alignment": alignment_workflow.run.command,
            "trimming": trimming_workflow.run.command,
            "model_selection": model_selection_workflow.run.command,
            "maximum_likelihood": maximum_likelihood_workflow.run.command,
            "bootstrap_support": bootstrap_workflow.run.command,
        },
        engine_versions={
            "mafft": alignment_workflow.run.version.text,
            "trimal": trimming_workflow.run.version.text,
            "iqtree_model_selection": model_selection_workflow.run.version.text,
            "iqtree_maximum_likelihood": maximum_likelihood_workflow.run.version.text,
            "iqtree_bootstrap_support": bootstrap_workflow.run.version.text,
        },
        step_manifests={
            "alignment": alignment_workflow.manifest_path,
            "trimming": trimming_workflow.manifest_path,
            "model_selection": model_selection_workflow.manifest_path,
            "maximum_likelihood": maximum_likelihood_workflow.manifest_path,
            "bootstrap_support": bootstrap_workflow.manifest_path,
        },
        stage_fingerprints={},
        input_checksums=build_file_checksums(
            [input_path]
            if prepared_input_path == input_path
            else [input_path, prepared_input_path]
        ),
        output_checksums={},
        input_validation=input_validation,
        repaired_input_validation=repaired_input_validation,
        input_repair=input_repair,
        alignment_workflow=alignment_workflow,
        trimming_workflow=trimming_workflow,
        model_selection_workflow=model_selection_workflow,
        maximum_likelihood_workflow=maximum_likelihood_workflow,
        bootstrap_workflow=bootstrap_workflow,
        model_validation=model_validation,
        support_summary=support_summary,
        model_rows=model_rows,
        support_rows=support_rows,
        method_tier=fasta_to_tree_method_tier(),
        warnings=warnings,
        notes=notes,
    )
    ended_at = datetime.now(UTC)
    report.ended_at_utc = (
        ended_at.replace(microsecond=0).isoformat().replace("+00:00", "Z")
    )
    report.runtime_seconds = max(
        0.0,
        round((ended_at - started_at).total_seconds(), 6),
    )
    write_fasta_to_tree_log(final_outputs["log"], report, root_dir=out_dir)
    from bijux_phylogenetics.reports.methods import (
        write_tree_inference_methods_summary_text,
    )

    write_tree_inference_methods_summary_text(
        final_outputs["methods_summary"],
        workflow_report=report,
    )
    validation_output_checksums = (
        {}
        if prepared_input_path == input_path
        else build_file_checksums([prepared_input_path])
    )
    validation_stage = build_stage_fingerprint(
        stage="fasta_validation",
        input_checksums=build_file_checksums([input_path]),
        output_checksums=validation_output_checksums,
        config={
            "declared_sequence_type": sequence_type,
            "normalize_identifiers": normalize_identifiers,
            "remove_invalid_records": remove_invalid_records,
        },
        engine_versions={},
        upstream_fingerprints={},
        resumed=False,
    )
    alignment_stage = build_stage_fingerprint(
        stage="alignment",
        input_checksums=alignment_workflow.input_checksums,
        output_checksums=alignment_workflow.output_checksums,
        config={
            "alignment_mode": alignment_mode,
            "timeout_seconds": timeout_seconds,
        },
        engine_versions={"mafft": alignment_workflow.run.version.text},
        upstream_fingerprints={"fasta_validation": validation_stage.fingerprint},
        resumed=alignment_workflow.resumed,
    )
    trimming_stage = build_stage_fingerprint(
        stage="trimming",
        input_checksums=trimming_workflow.input_checksums,
        output_checksums=trimming_workflow.output_checksums,
        config={
            "trimming_mode": trimming_mode,
            "trim_gap_threshold": trim_gap_threshold,
            "timeout_seconds": timeout_seconds,
        },
        engine_versions={"trimal": trimming_workflow.run.version.text},
        upstream_fingerprints={"alignment": alignment_stage.fingerprint},
        resumed=trimming_workflow.resumed,
    )
    model_selection_stage = build_stage_fingerprint(
        stage="model_selection",
        input_checksums=model_selection_workflow.input_checksums,
        output_checksums=model_selection_workflow.output_checksums,
        config={
            "iqtree_seed": iqtree_seed,
            "iqtree_threads": iqtree_threads,
            "timeout_seconds": timeout_seconds,
            "sequence_type": inferred_sequence_type,
        },
        engine_versions={
            "iqtree_model_selection": model_selection_workflow.run.version.text
        },
        upstream_fingerprints={"trimming": trimming_stage.fingerprint},
        resumed=model_selection_workflow.resumed,
    )
    inference_stage = build_stage_fingerprint(
        stage="inference",
        input_checksums=maximum_likelihood_workflow.input_checksums,
        output_checksums=maximum_likelihood_workflow.output_checksums,
        config={
            "selected_model": model_selection_workflow.selected_model,
            "iqtree_seed": iqtree_seed,
            "iqtree_threads": iqtree_threads,
            "timeout_seconds": timeout_seconds,
            "sequence_type": inferred_sequence_type,
        },
        engine_versions={
            "iqtree_maximum_likelihood": (maximum_likelihood_workflow.run.version.text)
        },
        upstream_fingerprints={
            "trimming": trimming_stage.fingerprint,
            "model_selection": model_selection_stage.fingerprint,
        },
        resumed=maximum_likelihood_workflow.resumed,
    )
    support_stage = build_stage_fingerprint(
        stage="support",
        input_checksums=bootstrap_workflow.input_checksums,
        output_checksums=bootstrap_workflow.output_checksums,
        config={
            "selected_model": model_selection_workflow.selected_model,
            "bootstrap_replicates": bootstrap_replicates,
            "iqtree_seed": iqtree_seed,
            "iqtree_threads": iqtree_threads,
            "timeout_seconds": timeout_seconds,
            "sequence_type": inferred_sequence_type,
        },
        engine_versions={
            "iqtree_bootstrap_support": bootstrap_workflow.run.version.text
        },
        upstream_fingerprints={
            "trimming": trimming_stage.fingerprint,
            "model_selection": model_selection_stage.fingerprint,
        },
        resumed=bootstrap_workflow.resumed,
    )
    report_output_checksums = build_file_checksums(
        [
            final_outputs["log"],
            final_outputs["methods_summary"],
            final_outputs["model_table"],
            final_outputs["support_table"],
        ]
    )
    report_stage = build_stage_fingerprint(
        stage="report",
        input_checksums=build_file_checksums(
            [
                final_outputs["alignment"],
                final_outputs["trimmed_alignment"],
                final_outputs["tree"],
            ]
        ),
        output_checksums=report_output_checksums,
        config={
            "selected_model": model_selection_workflow.selected_model,
            "warning_count": len(warnings),
            "note_count": len(notes),
        },
        engine_versions={},
        upstream_fingerprints={
            "fasta_validation": validation_stage.fingerprint,
            "alignment": alignment_stage.fingerprint,
            "trimming": trimming_stage.fingerprint,
            "model_selection": model_selection_stage.fingerprint,
            "inference": inference_stage.fingerprint,
            "support": support_stage.fingerprint,
        },
        resumed=False,
    )
    report.stage_fingerprints = {
        "fasta_validation": validation_stage,
        "alignment": alignment_stage,
        "trimming": trimming_stage,
        "model_selection": model_selection_stage,
        "inference": inference_stage,
        "support": support_stage,
        "report": report_stage,
    }
    report.output_checksums = build_file_checksums(
        [
            final_outputs["alignment"],
            final_outputs["trimmed_alignment"],
            final_outputs["tree"],
            final_outputs["log"],
            final_outputs["methods_summary"],
            final_outputs["model_table"],
            final_outputs["support_table"],
        ]
    )
    write_engine_manifest(report.manifest_path, report)
    write_run_manifest(
        report.run_manifest_path,
        build_run_manifest(
            command="run_fasta_to_tree_workflow",
            arguments=run_manifest_arguments,
            input_paths=(
                [input_path]
                if prepared_input_path == input_path
                else [input_path, prepared_input_path]
            ),
            output_paths=[
                final_outputs["alignment"],
                final_outputs["trimmed_alignment"],
                final_outputs["tree"],
                final_outputs["log"],
                final_outputs["methods_summary"],
                final_outputs["model_table"],
                final_outputs["support_table"],
                report.manifest_path,
            ],
        ),
    )
    return report
