from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path

from bijux_phylogenetics.io.fasta.quality import summarize_alignment_readiness
from bijux_phylogenetics.io.fasta.records import validate_fasta_input
from bijux_phylogenetics.runtime.errors import PhylogeneticsError


@dataclass(slots=True)
class ScientificFailureExplanation:
    """User-facing scientific explanation for one structured failure."""

    failure_reason: str
    scientific_explanation: str
    likely_causes: list[str]
    actionable_fixes: list[str]
    evidence: dict[str, object]


def explanation_payload(
    explanation: ScientificFailureExplanation | None,
) -> dict[str, object]:
    """Serialize one optional failure explanation into JSON-ready data."""
    return {} if explanation is None else asdict(explanation)


def explain_phylogenetics_error(
    error: PhylogeneticsError,
    *,
    inputs: list[Path],
) -> ScientificFailureExplanation | None:
    """Attach scientific context to one structured runtime error."""
    if _has_structured_explanation(error.details):
        return ScientificFailureExplanation(
            failure_reason=str(error.details["failure_reason"]),
            scientific_explanation=str(error.details["scientific_explanation"]),
            likely_causes=[
                str(item) for item in list(error.details.get("likely_causes", []))
            ],
            actionable_fixes=[
                str(item) for item in list(error.details.get("actionable_fixes", []))
            ],
            evidence=dict(error.details.get("evidence", {})),
        )
    if error.code == "invalid_alignment_error":
        return _explain_invalid_alignment_error(inputs)
    if error.code.startswith("beast_"):
        return _explain_beast_parser_error(error)
    if error.code.startswith("mrbayes_"):
        return _explain_mrbayes_parser_error(error)
    if error.code in {
        "engine_required_output_missing",
        "engine_output_empty",
        "engine_model_result_missing",
        "engine_support_values_missing",
    }:
        return _explain_engine_output_error(error)
    if error.code == "metadata_join_error":
        return _explain_metadata_join_error(error)
    return None


def explain_inference_workflow_failure(
    *,
    workflow: str,
    input_paths: list[Path],
    output_paths: dict[str, Path],
    run_exit_code: int | None,
    missing_inputs: list[Path],
    invalid_fasta_paths: list[Path],
    missing_outputs: dict[str, Path],
    invalid_outputs: dict[str, Path],
    parse_failures: dict[str, Path],
) -> ScientificFailureExplanation:
    """Explain one engine-workflow failure in biological workflow terms."""
    if missing_inputs:
        return ScientificFailureExplanation(
            failure_reason="declared_input_missing",
            scientific_explanation=(
                "The workflow could not start from a complete biological input set "
                "because one or more declared input files were absent."
            ),
            likely_causes=[
                "the input path is wrong or points to a file that was moved",
                "a previous preparation step did not write its expected artifact",
            ],
            actionable_fixes=[
                "check the declared input paths and rerun the upstream preparation step",
                "confirm that the expected FASTA, alignment, tree, or metadata file exists before resuming",
            ],
            evidence={"missing_inputs": [str(path) for path in missing_inputs]},
        )
    if invalid_fasta_paths:
        return _explain_invalid_fasta_path(invalid_fasta_paths[0])
    if run_exit_code is not None and run_exit_code != 0:
        if workflow == "multiple-sequence-alignment":
            return _explain_alignment_engine_failure(input_paths)
        if workflow == "alignment-trimming":
            return ScientificFailureExplanation(
                failure_reason="alignment_trimming_engine_failure",
                scientific_explanation=(
                    "The trimming step exited before producing a reviewed alignment, "
                    "so the input alignment could not be reduced into a defensible retained-site matrix."
                ),
                likely_causes=[
                    "the input alignment is malformed or not equal-length FASTA",
                    "the trimming mode was too aggressive for the retained signal",
                    "the trimAl executable failed before writing outputs",
                ],
                actionable_fixes=[
                    "inspect the input alignment for equal sequence length and valid FASTA structure",
                    "rerun with a less aggressive trimming mode or a higher gap threshold",
                    "review the trimAl stderr log for parser or engine-specific complaints",
                ],
                evidence={"run_exit_code": run_exit_code},
            )
        return ScientificFailureExplanation(
            failure_reason="external_engine_exit",
            scientific_explanation=(
                "The external phylogenetics engine exited before completing a scientifically usable artifact."
            ),
            likely_causes=[
                "the engine rejected the supplied biological input",
                "the engine stopped before finishing tree inference or support estimation",
                "engine runtime settings exceeded what the current input could complete safely",
            ],
            actionable_fixes=[
                "inspect the engine stderr log for the first concrete complaint",
                "confirm that the input alignment or tree is valid for the requested workflow",
                "rerun with narrower settings or smaller support budgets if the run timed out",
            ],
            evidence={"run_exit_code": run_exit_code},
        )
    if missing_outputs:
        output_names = sorted(missing_outputs)
        if workflow == "alignment-trimming" and "trimmed_alignment" in missing_outputs:
            return ScientificFailureExplanation(
                failure_reason="trimmed_alignment_missing",
                scientific_explanation=(
                    "Trimming finished without a retained alignment artifact, which usually means "
                    "the trimmer rejected the alignment or every informative site was removed."
                ),
                likely_causes=[
                    "the alignment is malformed or contains inconsistent sequence lengths",
                    "the trimming policy removed every column",
                    "the engine stopped before writing the trimmed FASTA",
                ],
                actionable_fixes=[
                    "inspect the input alignment for equal-length rows and valid FASTA syntax",
                    "use a less aggressive trimming mode or a higher gap threshold",
                    "check the trimming stderr log for an engine-side refusal or crash",
                ],
                evidence={"missing_outputs": _serialize_output_map(missing_outputs)},
            )
        if any(name in missing_outputs for name in {"tree", "support_tree"}):
            return ScientificFailureExplanation(
                failure_reason="tree_output_missing",
                scientific_explanation=(
                    "Inference ended without writing a tree artifact, so there is no biological hypothesis to evaluate."
                ),
                likely_causes=[
                    "the engine failed before finishing tree search",
                    "the alignment or model request was invalid for the inference engine",
                    "downstream support or tree-writing steps never completed",
                ],
                actionable_fixes=[
                    "inspect the engine stderr log for the first tree-inference error",
                    "confirm that the alignment is nonempty and suitable for the requested inference workflow",
                    "rerun after checking that the selected model and support settings are accepted by the engine",
                ],
                evidence={"missing_outputs": _serialize_output_map(missing_outputs)},
            )
        if any(
            name in missing_outputs for name in {"posterior_trees", "consensus_tree"}
        ):
            return ScientificFailureExplanation(
                failure_reason="posterior_tree_artifact_missing",
                scientific_explanation=(
                    "The Bayesian run did not produce the posterior tree artifact needed to summarize phylogenetic uncertainty."
                ),
                likely_causes=[
                    "the engine ended before sampling or summarizing trees",
                    "the posterior tree file path was never written by the engine",
                    "the run failed after partial diagnostics but before tree export",
                ],
                actionable_fixes=[
                    "check the engine stderr log and any partial log files for the first failure",
                    "confirm that the configured sampling frequency and chain settings are valid",
                    "rerun only after removing incomplete outputs or using the governed clean policy",
                ],
                evidence={"missing_outputs": _serialize_output_map(missing_outputs)},
            )
        return ScientificFailureExplanation(
            failure_reason="required_output_missing",
            scientific_explanation=(
                "The workflow stopped without writing one or more required scientific artifacts."
            ),
            likely_causes=[
                "the engine finished early or failed after partial output",
                "an expected artifact path was never produced",
            ],
            actionable_fixes=[
                "inspect the workflow stderr log and missing output paths together",
                "rerun only after checking the upstream input and workflow settings",
            ],
            evidence={
                "workflow": workflow,
                "missing_outputs": _serialize_output_map(missing_outputs),
                "missing_output_names": output_names,
            },
        )
    if invalid_outputs:
        if workflow == "alignment-trimming" and "trimmed_alignment" in invalid_outputs:
            return ScientificFailureExplanation(
                failure_reason="trimmed_alignment_empty",
                scientific_explanation=(
                    "The trimming workflow produced an empty or unusable retained alignment, "
                    "so the trimming policy removed all usable alignment signal."
                ),
                likely_causes=[
                    "all sites were removed by the trimming policy",
                    "the trimmed FASTA contains blank records or no retained columns",
                    "the starting alignment was too short or too gap-heavy for the requested filter",
                ],
                actionable_fixes=[
                    "rerun with a less aggressive trimming mode or higher gap threshold",
                    "inspect the pre-trim alignment for sparse or low-information columns",
                    "confirm that the workflow starts from a real aligned matrix rather than raw sequences",
                ],
                evidence={"invalid_outputs": _serialize_output_map(invalid_outputs)},
            )
        return ScientificFailureExplanation(
            failure_reason="required_output_unusable",
            scientific_explanation=(
                "A required workflow artifact exists on disk but is blank or otherwise unusable for scientific review."
            ),
            likely_causes=[
                "the engine wrote a placeholder file without real content",
                "the workflow produced an empty alignment or report artifact",
            ],
            actionable_fixes=[
                "inspect the offending output file directly",
                "rerun after checking the upstream alignment or tree input quality",
            ],
            evidence={"invalid_outputs": _serialize_output_map(invalid_outputs)},
        )
    if parse_failures:
        return ScientificFailureExplanation(
            failure_reason="tree_output_unparsable",
            scientific_explanation=(
                "The workflow wrote a tree-like artifact, but its contents are not a valid parseable tree, "
                "so the result cannot be interpreted biologically."
            ),
            likely_causes=[
                "the engine wrote a truncated tree file",
                "support labels or Newick syntax are malformed",
                "a post-processing step interrupted tree export",
            ],
            actionable_fixes=[
                "inspect the tree artifact for truncation or malformed Newick syntax",
                "rerun the workflow after removing incomplete outputs",
                "compare the engine stderr log with the unparsable artifact path",
            ],
            evidence={"parse_failures": _serialize_output_map(parse_failures)},
        )
    return ScientificFailureExplanation(
        failure_reason="no_failure",
        scientific_explanation="The workflow state does not expose a scientific failure.",
        likely_causes=[],
        actionable_fixes=[],
        evidence={"workflow": workflow},
    )


def _has_structured_explanation(details: dict[str, object]) -> bool:
    return (
        "failure_reason" in details
        and "scientific_explanation" in details
        and "likely_causes" in details
        and "actionable_fixes" in details
        and "evidence" in details
    )


def _explain_invalid_alignment_error(
    inputs: list[Path],
) -> ScientificFailureExplanation | None:
    fasta_path = next(
        (
            path
            for path in inputs
            if path.suffix.lower() in {".fasta", ".fa", ".fas", ".faa", ".fna"}
            and path.exists()
        ),
        None,
    )
    return None if fasta_path is None else _explain_invalid_fasta_path(fasta_path)


def _explain_invalid_fasta_path(path: Path) -> ScientificFailureExplanation:
    report = validate_fasta_input(path)
    duplicate_ids = [row.identifier for row in report.duplicate_identifiers]
    illegal_rows = [
        {
            "identifier": row.identifier,
            "position": row.position,
            "character": row.character,
        }
        for row in report.illegal_characters[:5]
    ]
    empty_ids = [row.identifier for row in report.empty_sequences]
    outlier_ids = [row.identifier for row in report.length_outliers]
    likely_causes: list[str] = []
    actionable_fixes: list[str] = []
    if duplicate_ids:
        likely_causes.append(
            "one or more FASTA identifiers are duplicated and can no longer anchor unique biological samples"
        )
        actionable_fixes.append(
            "rename duplicate sequence identifiers so each biological sample appears once"
        )
    if illegal_rows:
        likely_causes.append(
            "unsupported sequence characters are present for the inferred sequence type"
        )
        actionable_fixes.append(
            "remove or repair residues that are not valid for the intended DNA, RNA, or protein alphabet"
        )
    if empty_ids:
        likely_causes.append(
            "one or more FASTA records are empty and contribute no biological signal"
        )
        actionable_fixes.append(
            "drop empty records or refill them from the source sequence data before alignment"
        )
    if outlier_ids:
        likely_causes.append(
            "sequence length outliers suggest mixed loci, truncation, or malformed records"
        )
        actionable_fixes.append(
            "review extreme sequence-length outliers for truncation, concatenation mistakes, or mixed-marker input"
        )
    if not likely_causes:
        likely_causes.append(
            "the FASTA input failed structural validation before alignment or inference"
        )
        actionable_fixes.append(
            "inspect the FASTA headers and record structure before rerunning the workflow"
        )
    return ScientificFailureExplanation(
        failure_reason="invalid_fasta_input",
        scientific_explanation=(
            "The FASTA input contains record-level biological data problems that block a defensible alignment or inference run."
        ),
        likely_causes=likely_causes,
        actionable_fixes=actionable_fixes,
        evidence={
            "path": str(path),
            "record_count": report.summary.sequence_count,
            "inferred_alphabet": report.summary.inferred_alphabet,
            "duplicate_identifier_count": len(report.duplicate_identifiers),
            "duplicate_identifiers": duplicate_ids,
            "illegal_character_count": len(report.illegal_characters),
            "illegal_characters": illegal_rows,
            "empty_sequence_count": len(report.empty_sequences),
            "empty_sequences": empty_ids,
            "length_outlier_count": len(report.length_outliers),
            "length_outlier_identifiers": outlier_ids,
        },
    )


def _explain_alignment_engine_failure(
    input_paths: list[Path],
) -> ScientificFailureExplanation:
    fasta_path = next((path for path in input_paths if path.exists()), None)
    readiness = (
        summarize_alignment_readiness(fasta_path)
        if fasta_path is not None
        and fasta_path.suffix.lower() in {".fasta", ".fa", ".fas", ".faa", ".fna"}
        else None
    )
    blockers = (
        []
        if readiness is None
        else [
            blocker
            for decision in readiness.decisions
            if decision.workflow == "maximum_likelihood"
            for blocker in decision.blockers
        ]
    )
    warnings = [] if readiness is None else list(readiness.warnings)
    likely_causes = [
        "the raw sequence set could not be aligned into one defensible equal-length matrix"
    ]
    actionable_fixes = [
        "inspect the MAFFT stderr log for the first engine complaint",
        "confirm that the FASTA contains real biological sequences rather than malformed or mixed-type records",
    ]
    if blockers:
        likely_causes.extend(blockers[:3])
    if any("not yet aligned" in blocker for blocker in blockers):
        actionable_fixes.append(
            "remove gross length outliers or mixed alphabets before retrying alignment"
        )
    return ScientificFailureExplanation(
        failure_reason="alignment_engine_failure",
        scientific_explanation=(
            "The multiple-sequence alignment step failed before producing an equal-length homologous matrix."
        ),
        likely_causes=list(dict.fromkeys(likely_causes)),
        actionable_fixes=list(dict.fromkeys(actionable_fixes)),
        evidence={
            "path": None if fasta_path is None else str(fasta_path),
            "blockers": blockers,
            "warnings": warnings,
        },
    )


def _explain_engine_output_error(
    error: PhylogeneticsError,
) -> ScientificFailureExplanation:
    details = error.details
    workflow = str(details.get("workflow", ""))
    output_name = (
        None if details.get("output_name") is None else str(details["output_name"])
    )
    missing_outputs = list(details.get("missing_outputs", []))
    if error.code == "engine_model_result_missing":
        return ScientificFailureExplanation(
            failure_reason="model_result_missing",
            scientific_explanation=(
                "The workflow produced inference sidecars but not a defensible substitution-model result, "
                "so downstream tree interpretation would rest on an undocumented model choice."
            ),
            likely_causes=[
                "the IQ-TREE report or model sidecar is missing or unparsable",
                "the engine terminated before recording its best-fit model",
            ],
            actionable_fixes=[
                "inspect the IQ-TREE report and model sidecar together",
                "rerun after confirming that the requested model-selection mode is supported by the input data",
            ],
            evidence={
                "workflow": workflow,
                "path": details.get("path"),
                "model_sidecar_path": details.get("model_sidecar_path"),
            },
        )
    if error.code == "engine_support_values_missing":
        return ScientificFailureExplanation(
            failure_reason="support_values_missing",
            scientific_explanation=(
                "The support workflow wrote a tree but did not expose parsable branch-support values, "
                "so clade certainty cannot be reviewed."
            ),
            likely_causes=[
                "the tree lacks the expected bootstrap or SH-aLRT/UFBoot labels",
                "support labels are malformed or were stripped during export",
            ],
            actionable_fixes=[
                "inspect the support tree directly for numeric branch labels",
                "rerun the support workflow and compare the tree artifact with the engine report",
            ],
            evidence={
                "workflow": workflow,
                "path": details.get("path"),
                "support_kind": details.get("support_kind"),
            },
        )
    if (
        error.code == "engine_output_empty"
        and workflow == "alignment-trimming"
        and output_name == "trimmed_alignment"
    ):
        return ScientificFailureExplanation(
            failure_reason="trimmed_alignment_empty",
            scientific_explanation=(
                "The trimming step returned an empty retained alignment, meaning the filter removed all usable alignment signal."
            ),
            likely_causes=[
                "the trimming settings were too aggressive for the input matrix",
                "the alignment is so sparse or short that no columns survived",
            ],
            actionable_fixes=[
                "use a less aggressive trimming mode or higher gap threshold",
                "inspect the input alignment for sparse columns and extreme missingness",
            ],
            evidence={"workflow": workflow, "path": details.get("path")},
        )
    if error.code == "engine_required_output_missing":
        missing_names = [
            str(row.get("output_name"))
            for row in missing_outputs
            if isinstance(row, dict) and row.get("output_name") is not None
        ]
        if workflow == "posterior-tree-inference":
            return ScientificFailureExplanation(
                failure_reason="posterior_artifact_missing",
                scientific_explanation=(
                    "The Bayesian run finished without the posterior artifact required to summarize phylogenetic uncertainty."
                ),
                likely_causes=[
                    "the run ended before tree sampling or consensus export completed",
                    "the engine wrote partial diagnostics but not the posterior tree artifact",
                ],
                actionable_fixes=[
                    "inspect the posterior log, stderr log, and missing artifact path together",
                    "clean incomplete outputs before rerunning the Bayesian workflow",
                ],
                evidence={
                    "missing_outputs": missing_outputs,
                    "missing_output_names": missing_names,
                },
            )
    return ScientificFailureExplanation(
        failure_reason="engine_artifact_failure",
        scientific_explanation=(
            "The engine workflow did not yield a complete scientifically usable artifact set."
        ),
        likely_causes=[
            "one or more required outputs are missing, blank, or malformed",
        ],
        actionable_fixes=[
            "inspect the structured error details and the engine stderr log together",
        ],
        evidence={
            "workflow": workflow,
            "output_name": output_name,
            "missing_outputs": missing_outputs,
            "path": details.get("path"),
        },
    )


def _explain_metadata_join_error(
    error: PhylogeneticsError,
) -> ScientificFailureExplanation:
    details = error.details
    return ScientificFailureExplanation(
        failure_reason="tree_trait_taxon_mismatch",
        scientific_explanation=(
            "The trait table and tree do not describe the same biological taxon set, "
            "so comparative interpretation would silently drop taxa or compare the wrong rows."
        ),
        likely_causes=[
            "one or more tree tips are missing from the trait table",
            "the trait table contains taxa not present in the tree",
        ],
        actionable_fixes=[
            "add the missing tree taxa to the trait table or prune the tree intentionally",
            "remove extra trait-table taxa that are not represented in the tree",
        ],
        evidence=dict(details),
    )


def _explain_beast_parser_error(
    error: PhylogeneticsError,
) -> ScientificFailureExplanation:
    details = error.details
    code = error.code
    artifact_kind = str(details.get("artifact_kind", "beast-artifact"))
    path = details.get("path")
    if code.endswith("_missing_file"):
        explanation = "The expected BEAST artifact file is absent, so posterior diagnostics cannot be reconstructed."
        causes = ["the BEAST run did not write the expected output file"]
        fixes = ["confirm the file path and rerun the BEAST step that should create it"]
        reason = "beast_artifact_missing_file"
    elif code.endswith("_missing_header"):
        explanation = "The BEAST tabular artifact is missing its header section, so sampled parameters cannot be identified safely."
        causes = ["the file is truncated or not a real BEAST tabular export"]
        fixes = [
            "inspect the first lines of the file and rerun the BEAST export if the header is missing"
        ]
        reason = "beast_tabular_header_missing"
    elif code.endswith("_missing_state_column"):
        explanation = "The BEAST posterior log lacks the state column that anchors sampled rows to MCMC iteration order."
        causes = [
            "the log file is malformed or not the expected posterior parameter log"
        ]
        fixes = [
            "use the BEAST posterior log with its full Sample or state column intact"
        ]
        reason = "beast_state_column_missing"
    elif code.endswith("_missing_rows"):
        explanation = "The BEAST artifact contains no sampled rows or tree entries, so there is no posterior evidence to summarize."
        causes = [
            "the file is truncated after the header or contains no posterior samples"
        ]
        fixes = [
            "rerun the BEAST analysis and confirm that sampling progressed beyond burn-in"
        ]
        reason = "beast_sample_rows_missing"
    elif code.endswith("_missing_entries"):
        explanation = "The BEAST posterior tree file contains no tree entries, so topology uncertainty cannot be summarized."
        causes = ["the trees block is empty or missing from the posterior tree file"]
        fixes = [
            "inspect the BEAST trees file for a populated trees block and rerun if it is empty"
        ]
        reason = "beast_tree_entries_missing"
    else:
        explanation = "The BEAST artifact contains a malformed sampled value or malformed tree content, so posterior summaries would be unreliable."
        causes = ["one sampled value or tree entry is malformed or truncated"]
        fixes = [
            "inspect the reported row, column, or tree entry and regenerate the artifact from BEAST"
        ]
        reason = "beast_artifact_malformed"
    return ScientificFailureExplanation(
        failure_reason=reason,
        scientific_explanation=explanation,
        likely_causes=causes,
        actionable_fixes=fixes,
        evidence={
            "artifact_kind": artifact_kind,
            "path": path,
            **{
                key: value
                for key, value in details.items()
                if key not in {"artifact_kind", "path"}
            },
        },
    )


def _explain_mrbayes_parser_error(
    error: PhylogeneticsError,
) -> ScientificFailureExplanation:
    details = error.details
    code = error.code
    artifact_kind = str(details.get("artifact_kind", "mrbayes-artifact"))
    path = details.get("path")
    if code.endswith("_missing_file"):
        explanation = "The expected MrBayes artifact file is absent, so posterior diagnostics cannot be reconstructed."
        causes = ["the MrBayes run did not write the expected output file"]
        fixes = [
            "confirm the file path and rerun the MrBayes step that should create it"
        ]
        reason = "mrbayes_artifact_missing_file"
    elif code.endswith("_missing_header"):
        explanation = "The MrBayes tabular artifact is missing its header section, so sampled parameters cannot be identified safely."
        causes = ["the file is truncated or not a real MrBayes table export"]
        fixes = [
            "inspect the first lines of the file and rerun the MrBayes export if the header is missing"
        ]
        reason = "mrbayes_tabular_header_missing"
    elif code.endswith("_missing_generation_column"):
        explanation = "The MrBayes table lacks the Gen column that anchors sampled rows to MCMC generation order."
        causes = [
            "the file is malformed or not the expected MrBayes trace or MCMC table"
        ]
        fixes = ["use the original MrBayes table with its Gen column intact"]
        reason = "mrbayes_generation_column_missing"
    elif code.endswith("_missing_rows"):
        explanation = "The MrBayes artifact contains no sampled rows, so convergence or parameter summaries cannot be computed."
        causes = [
            "the file is truncated after the header or contains no posterior samples"
        ]
        fixes = [
            "rerun MrBayes and confirm that sampling progressed beyond its initial setup"
        ]
        reason = "mrbayes_sample_rows_missing"
    else:
        explanation = "The MrBayes artifact contains a malformed sampled value or malformed summary content, so posterior summaries would be unreliable."
        causes = ["one sampled value is missing, non-numeric, or truncated"]
        fixes = [
            "inspect the reported row and column and regenerate the artifact from MrBayes"
        ]
        reason = "mrbayes_artifact_malformed"
    return ScientificFailureExplanation(
        failure_reason=reason,
        scientific_explanation=explanation,
        likely_causes=causes,
        actionable_fixes=fixes,
        evidence={
            "artifact_kind": artifact_kind,
            "path": path,
            **{
                key: value
                for key, value in details.items()
                if key not in {"artifact_kind", "path"}
            },
        },
    )


def _serialize_output_map(paths: dict[str, Path]) -> list[dict[str, str]]:
    return [
        {"output_name": output_name, "path": str(path)}
        for output_name, path in sorted(paths.items())
    ]
