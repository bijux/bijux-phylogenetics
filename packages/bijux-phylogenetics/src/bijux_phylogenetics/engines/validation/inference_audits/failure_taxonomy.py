from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.io.fasta import load_fasta_alignment
from bijux_phylogenetics.io.fasta.records import validate_fasta_input
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.runtime.error_explanations import (
    explain_inference_workflow_failure,
)
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError

from .contracts import InferenceFailureTaxonomyReport


def classify_inference_workflow_failure(
    *,
    workflow: str,
    input_paths: list[Path],
    output_paths: dict[str, Path],
    run_exit_code: int | None = None,
) -> InferenceFailureTaxonomyReport:
    """Classify one workflow state into a stable inference-failure taxonomy."""
    issues: list[str] = []
    failure_category = "no_failure"
    missing_inputs = [path for path in input_paths if not path.exists()]
    invalid_fasta_paths: list[Path] = []
    if missing_inputs:
        failure_category = "input_failure"
        issues.append("one or more declared input files are missing")
    else:
        for input_path in input_paths:
            if input_path.suffix.lower() in {".fasta", ".fa", ".fas", ".faa", ".fna"}:
                try:
                    validate_fasta_input(input_path)
                    load_fasta_alignment(input_path)
                except InvalidAlignmentError:
                    failure_category = "input_failure"
                    issues.append("one or more alignment inputs are invalid")
                    invalid_fasta_paths.append(input_path)
                    break
    if run_exit_code is not None and run_exit_code != 0:
        failure_category = "timeout" if run_exit_code == 124 else "engine_failure"
        issues.append(f"engine exited with code {run_exit_code}")
    missing_output_map = {
        output_name: path
        for output_name, path in output_paths.items()
        if not path.exists()
    }
    missing_outputs = sorted(str(path) for path in missing_output_map.values())
    if missing_outputs:
        failure_category = "missing_output"
        issues.append("one or more expected outputs are missing")
    parse_failures: dict[str, Path] = {}
    invalid_outputs: dict[str, Path] = {}
    for key, path in output_paths.items():
        if not path.exists():
            continue
        if path.is_file() and not path.read_text(encoding="utf-8").strip():
            invalid_outputs[key] = path
            continue
        if path.suffix.lower() in {
            ".treefile",
            ".contree",
            ".nwk",
            ".tre",
            ".tree",
            ".ufboot",
        }:
            try:
                if key == "bootstrap_trees":
                    from bijux_phylogenetics.trees import load_tree_set

                    load_tree_set(path)
                else:
                    load_tree(path)
            except Exception:
                parse_failures[key] = path
    if invalid_outputs:
        failure_category = "invalid_output"
        issues.append("one or more outputs exist but are blank or unusable")
    if parse_failures:
        failure_category = "parse_failure"
        issues.append("one or more tree-like outputs could not be parsed")
    explanation = explain_inference_workflow_failure(
        workflow=workflow,
        input_paths=input_paths,
        output_paths=output_paths,
        run_exit_code=run_exit_code,
        missing_inputs=missing_inputs,
        invalid_fasta_paths=invalid_fasta_paths,
        missing_outputs=missing_output_map,
        invalid_outputs=invalid_outputs,
        parse_failures=parse_failures,
    )
    return InferenceFailureTaxonomyReport(
        workflow=workflow,
        failure_category=failure_category,
        failure_reason=explanation.failure_reason,
        scientific_explanation=explanation.scientific_explanation,
        likely_causes=explanation.likely_causes,
        actionable_fixes=explanation.actionable_fixes,
        evidence=explanation.evidence,
        valid=failure_category == "no_failure",
        issues=issues,
    )
