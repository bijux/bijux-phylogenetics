from __future__ import annotations

from dataclasses import dataclass, field, replace
from datetime import datetime
import json
import math
from pathlib import Path
import shutil
import tempfile

from bijux_phylogenetics.datasets.shared_fixtures import (
    SharedBeastPosteriorFixture,
)
from bijux_phylogenetics.runtime.errors import EngineWorkflowError

from ..common import build_file_checksums, utc_now_text
from ..workflows.models import EngineWorkflowReport

_CANONICAL_ENGINE_NAMES = {
    "mafft": "MAFFT",
    "trimal": "trimAl",
    "iqtree": "IQ-TREE",
    "fasttree": "FastTree",
    "mrbayes": "MrBayes",
    "beast": "BEAST",
}


def _parse_utc_timestamp(text: str) -> datetime:
    return datetime.fromisoformat(text.replace("Z", "+00:00"))


def _runtime_seconds(*, started_at_utc: str, ended_at_utc: str) -> float:
    return max(
        (
            _parse_utc_timestamp(ended_at_utc) - _parse_utc_timestamp(started_at_utc)
        ).total_seconds(),
        0.0,
    )


def _canonical_engine_name(engine_name: str) -> str:
    return _CANONICAL_ENGINE_NAMES.get(engine_name.lower(), engine_name)


def _labeled_output_checksums(
    output_paths: dict[str, Path],
    *,
    stored_checksums: dict[str, str] | None = None,
) -> dict[str, str]:
    checksum_payload = (
        build_file_checksums(list(output_paths.values()))
        if stored_checksums is None or not stored_checksums
        else dict(stored_checksums)
    )
    checksums_by_path = {str(key): value for key, value in checksum_payload.items()}
    return {label: checksums_by_path[str(path)] for label, path in output_paths.items()}


@dataclass(slots=True)
class ExternalEngineValidationCase:
    engine_name: str
    validation_name: str
    validation_mode: str
    manifest_path: Path | None
    executable: str | None
    version_text: str | None
    command: list[str]
    exit_code: int | None
    runtime_seconds: float | None
    output_paths: dict[str, Path]
    output_checksums: dict[str, str]
    notes: list[str] = field(default_factory=list)


@dataclass(slots=True)
class ExternalEngineValidationMatrixReport:
    generated_at_utc: str
    cases: list[ExternalEngineValidationCase]


def build_external_engine_validation_case(
    validation_name: str,
    report: EngineWorkflowReport,
    *,
    validation_mode: str = "workflow-run",
    notes: list[str] | None = None,
) -> ExternalEngineValidationCase:
    """Normalize one workflow report into a reviewer-facing validation row."""
    output_checksums = _labeled_output_checksums(
        dict(report.output_paths),
        stored_checksums=(
            None if not report.output_checksums else dict(report.output_checksums)
        ),
    )
    return ExternalEngineValidationCase(
        engine_name=_canonical_engine_name(report.engine_name),
        validation_name=validation_name,
        validation_mode=validation_mode,
        manifest_path=report.manifest_path,
        executable=report.run.executable,
        version_text=report.run.version.text,
        command=list(report.run.command),
        exit_code=report.run.exit_code,
        runtime_seconds=_runtime_seconds(
            started_at_utc=report.run.started_at_utc,
            ended_at_utc=report.run.ended_at_utc,
        ),
        output_paths=dict(report.output_paths),
        output_checksums=output_checksums,
        notes=[*report.notes, *(notes or [])],
    )


def build_beast_artifact_validation_case(
    validation_name: str,
    *,
    xml_path: Path,
    log_path: Path,
    tree_path: Path,
    burnin_fraction: float = 0.1,
) -> ExternalEngineValidationCase:
    """Build one validation row from governed real BEAST analysis artifacts."""
    from bijux_phylogenetics.bayesian.beast.logs import (
        summarize_beast_log,
    )
    from bijux_phylogenetics.bayesian.beast.posterior_trees import (
        parse_beast_posterior_tree_samples,
    )
    from bijux_phylogenetics.bayesian.beast.xml_analysis import (
        summarize_beast_analysis_xml,
    )

    xml_report = summarize_beast_analysis_xml(xml_path)
    if not xml_report.valid:
        raise EngineWorkflowError(
            f"BEAST analysis XML failed validation for matrix case '{validation_name}': {xml_path}"
        )
    log_summary = summarize_beast_log(log_path, burnin_fraction=burnin_fraction)
    tree_report = parse_beast_posterior_tree_samples(
        tree_path,
        burnin_fraction=burnin_fraction,
    )
    output_paths = {
        "analysis_xml": xml_path,
        "posterior_log": log_path,
        "posterior_trees": tree_path,
    }
    return ExternalEngineValidationCase(
        engine_name="BEAST",
        validation_name=validation_name,
        validation_mode="fixture-parse",
        manifest_path=None,
        executable=None,
        version_text=xml_report.beast_version,
        command=[],
        exit_code=None,
        runtime_seconds=None,
        output_paths=output_paths,
        output_checksums=_labeled_output_checksums(output_paths),
        notes=[
            f"beast xml taxon count: {xml_report.taxon_count}",
            f"beast log kept rows after burn-in: {log_summary.kept_row_count}",
            f"beast posterior trees kept after burn-in: {tree_report.kept_tree_count}",
        ],
    )


def build_governed_beast_fixture_validation_case(
    validation_name: str,
    fixture: SharedBeastPosteriorFixture,
    *,
    burnin_fraction: float | None = None,
) -> ExternalEngineValidationCase:
    """Build one reviewer-facing validation row from the governed BEAST posterior corpus."""
    from bijux_phylogenetics.bayesian.beast.logs import (
        summarize_beast_log,
    )
    from bijux_phylogenetics.bayesian.beast.posterior_trees import (
        parse_beast_posterior_tree_samples,
        summarize_beast_posterior_trees,
    )
    from bijux_phylogenetics.bayesian.posterior_sets.tree_sets import (
        summarize_maximum_clade_credibility_tree,
    )

    selected_burnin = (
        fixture.recommended_burnin_fraction
        if burnin_fraction is None
        else burnin_fraction
    )
    reference = fixture.load_reference()
    if selected_burnin not in reference.burnin_reference:
        supported = ", ".join(
            format(value, ".15g") for value in sorted(reference.burnin_reference)
        )
        raise EngineWorkflowError(
            "governed BEAST posterior fixture reference is unavailable for "
            f"burn-in fraction {selected_burnin}; expected one of: {supported}"
        )
    case = build_beast_artifact_validation_case(
        validation_name,
        xml_path=fixture.analysis_xml_path,
        log_path=fixture.posterior_log_path,
        tree_path=fixture.posterior_trees_path,
        burnin_fraction=selected_burnin,
    )
    expected_counts = reference.burnin_reference[selected_burnin]
    log_summary = summarize_beast_log(
        fixture.posterior_log_path,
        burnin_fraction=selected_burnin,
    )
    tree_report = parse_beast_posterior_tree_samples(
        fixture.posterior_trees_path,
        burnin_fraction=selected_burnin,
    )
    if log_summary.burnin_row_count != expected_counts.burnin_row_count:
        raise EngineWorkflowError(
            "governed BEAST posterior fixture burn-in rows do not match the "
            f"checked reference for {selected_burnin}: "
            f"observed {log_summary.burnin_row_count}, "
            f"expected {expected_counts.burnin_row_count}"
        )
    if log_summary.kept_row_count != expected_counts.kept_row_count:
        raise EngineWorkflowError(
            "governed BEAST posterior fixture retained log rows do not match the "
            f"checked reference for {selected_burnin}: "
            f"observed {log_summary.kept_row_count}, "
            f"expected {expected_counts.kept_row_count}"
        )
    if tree_report.burnin_tree_count != expected_counts.burnin_tree_count:
        raise EngineWorkflowError(
            "governed BEAST posterior fixture burn-in trees do not match the "
            f"checked reference for {selected_burnin}: "
            f"observed {tree_report.burnin_tree_count}, "
            f"expected {expected_counts.burnin_tree_count}"
        )
    if tree_report.kept_tree_count != expected_counts.kept_tree_count:
        raise EngineWorkflowError(
            "governed BEAST posterior fixture retained trees do not match the "
            f"checked reference for {selected_burnin}: "
            f"observed {tree_report.kept_tree_count}, "
            f"expected {expected_counts.kept_tree_count}"
        )
    with tempfile.TemporaryDirectory(prefix="bijux-beast-reference-") as temp_dir:
        temp_tree_path = Path(temp_dir) / fixture.posterior_trees_path.name
        shutil.copyfile(fixture.posterior_trees_path, temp_tree_path)
        _consensus_tree, consensus = summarize_beast_posterior_trees(
            temp_tree_path,
            burnin_fraction=selected_burnin,
        )
        _mcc_tree, mcc = summarize_maximum_clade_credibility_tree(
            temp_tree_path,
            burnin_fraction=selected_burnin,
        )
    if selected_burnin == reference.consensus_reference.burnin_fraction:
        _assert_beast_reference_close(
            observed=consensus.minimum_posterior_probability,
            expected=reference.consensus_reference.minimum_posterior_probability,
            label="minimum posterior probability",
        )
        _assert_beast_reference_close(
            observed=consensus.maximum_posterior_probability,
            expected=reference.consensus_reference.maximum_posterior_probability,
            label="maximum posterior probability",
        )
        if (
            consensus.annotated_node_count
            != reference.consensus_reference.annotated_node_count
        ):
            raise EngineWorkflowError(
                "governed BEAST consensus annotation counts do not match the "
                f"checked reference for {selected_burnin}"
            )
        if consensus.consensus_newick != reference.consensus_reference.newick:
            raise EngineWorkflowError(
                "governed BEAST consensus topology does not match the checked "
                f"reference for {selected_burnin}"
            )
        if mcc.selected_tree_index != reference.mcc_reference.selected_tree_index:
            raise EngineWorkflowError(
                "governed BEAST maximum clade credibility tree selection does not "
                f"match the checked reference for {selected_burnin}"
            )
        _assert_beast_reference_close(
            observed=mcc.clade_credibility_score,
            expected=reference.mcc_reference.clade_credibility_score,
            label="clade credibility score",
        )
        if mcc.mcc_newick != reference.mcc_reference.newick:
            raise EngineWorkflowError(
                "governed BEAST maximum clade credibility topology does not match "
                f"the checked reference for {selected_burnin}"
            )
        parameter_summaries = {
            row.parameter: row for row in log_summary.parameter_summaries
        }
        for parameter, expected in reference.parameter_reference.items():
            observed = parameter_summaries[parameter]
            _assert_beast_reference_close(
                observed=observed.effective_sample_size,
                expected=expected.effective_sample_size,
                label=f"{parameter} ESS",
            )
            _assert_beast_reference_close(
                observed=observed.mean,
                expected=expected.mean,
                label=f"{parameter} mean",
            )
            _assert_beast_reference_close(
                observed=observed.median,
                expected=expected.median,
                label=f"{parameter} median",
            )
            _assert_beast_reference_close(
                observed=observed.hpd_95_lower,
                expected=expected.hpd_95_lower,
                label=f"{parameter} HPD lower",
            )
            _assert_beast_reference_close(
                observed=observed.hpd_95_upper,
                expected=expected.hpd_95_upper,
                label=f"{parameter} HPD upper",
            )
    return replace(
        case,
        notes=[
            *case.notes,
            f"governed beast fixture id: {fixture.fixture_id}",
            f"governed beast shared taxa: {', '.join(fixture.shared_taxa)}",
            (
                "governed beast consensus and maximum clade credibility summaries "
                "match the checked posterior reference bundle"
            ),
            (
                "governed beast posterior parameter ESS and interval summaries "
                "match the checked reference bundle"
            ),
        ],
    )


def _assert_beast_reference_close(
    *,
    observed: float | None,
    expected: float,
    label: str,
) -> None:
    if observed is None or not math.isclose(
        observed, expected, rel_tol=1e-9, abs_tol=1e-6
    ):
        raise EngineWorkflowError(
            f"governed BEAST reference mismatch for {label}: observed {observed}, expected {expected}"
        )


def build_external_engine_validation_matrix(
    cases: list[ExternalEngineValidationCase],
) -> ExternalEngineValidationMatrixReport:
    """Assemble one ordered external-engine validation matrix."""
    return ExternalEngineValidationMatrixReport(
        generated_at_utc=utc_now_text(),
        cases=list(cases),
    )


def merge_external_engine_validation_matrices(
    reports: list[ExternalEngineValidationMatrixReport],
) -> ExternalEngineValidationMatrixReport:
    """Merge multiple ordered validation matrices into one reviewer-facing report."""
    merged_cases: list[ExternalEngineValidationCase] = []
    for report in reports:
        merged_cases.extend(report.cases)
    return ExternalEngineValidationMatrixReport(
        generated_at_utc=utc_now_text(),
        cases=merged_cases,
    )


def write_external_engine_validation_matrix(
    path: Path,
    report: ExternalEngineValidationMatrixReport,
) -> Path:
    """Persist one portable JSON validation matrix for reviewer inspection."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "generated_at_utc": report.generated_at_utc,
        "case_count": len(report.cases),
        "engine_names": sorted({case.engine_name for case in report.cases}),
        "cases": [
            {
                "engine_name": case.engine_name,
                "validation_name": case.validation_name,
                "validation_mode": case.validation_mode,
                "manifest_path": (
                    None if case.manifest_path is None else str(case.manifest_path)
                ),
                "executable": case.executable,
                "version_text": case.version_text,
                "command": case.command,
                "exit_code": case.exit_code,
                "runtime_seconds": case.runtime_seconds,
                "output_paths": {
                    key: str(value) for key, value in case.output_paths.items()
                },
                "output_checksums": dict(case.output_checksums),
                "notes": list(case.notes),
            }
            for case in report.cases
        ],
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return path
