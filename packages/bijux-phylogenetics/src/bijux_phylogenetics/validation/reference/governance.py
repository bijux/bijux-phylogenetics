from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from bijux_phylogenetics.core.dataset import audit_dataset_inputs
from bijux_phylogenetics.diagnostics.validation import validate_tree_path
from bijux_phylogenetics.io.fasta.records import summarise_fasta
from bijux_phylogenetics.phylo.taxa import build_taxon_workflow_loss_report
from bijux_phylogenetics.render.tree_figure_package import build_tree_figure_package

from .core_suites import (
    validate_alignment_quality_reference_fixtures,
    validate_dataset_audit_reference_fixtures,
    validate_taxon_naming_reference_fixtures,
    validate_tree_reference_fixtures,
)
from .models import (
    CoreWorkflowFailureCase,
    CoreWorkflowValidationReport,
    CoreWorkflowValidationRow,
    LevelOneReleaseGateDecision,
    LevelOneReleaseGateReport,
    ReferenceValidationSuiteReport,
    WorkflowMaturityClassification,
)
from .publication_suites import validate_figure_reference_fixtures
from .report_regression import validate_report_regression_fixtures
from .shared import default_fixtures_root, fixture, temp_reference_dir


def _core_workflow_validation_rows(
    suites: list[ReferenceValidationSuiteReport],
) -> list[CoreWorkflowValidationRow]:
    suite_map = {suite.suite: suite for suite in suites}
    definitions = [
        (
            "tree-review",
            ["tree-validation-reference", "report-regression-reference"],
            ["tree validation report", "tree inspection report", "tree html golden"],
            ["single-tree diagnostics do not replace multi-tree uncertainty analysis"],
        ),
        (
            "taxon-audit",
            ["taxon-naming-reference", "report-regression-reference"],
            [
                "taxon audit report",
                "mapping conflict audit",
                "taxonomy html section contract",
            ],
            [
                "naming audits still rely on the supplied synonym table and cannot infer authority completeness"
            ],
        ),
        (
            "alignment-review",
            ["alignment-quality-reference", "report-regression-reference"],
            [
                "alignment quality report",
                "duplicate policy report",
                "alignment html section contract",
            ],
            [
                "alignment quality signals are deterministic heuristics, not biological truth"
            ],
        ),
        (
            "dataset-audit",
            ["dataset-audit-reference", "report-regression-reference"],
            [
                "dataset audit report",
                "exclusion table",
                "dataset html section contract",
            ],
            [
                "dataset fixture coverage is strongest for mismatch traceability rather than for large cohorts"
            ],
        ),
        (
            "figure-package",
            ["figure-correctness-reference"],
            ["support-label render audit", "topology-preserving figure metadata"],
            [
                "figure correctness is validated through render metadata rather than screenshot diffs"
            ],
        ),
        (
            "phylo-inputs-review",
            [
                "tree-validation-reference",
                "alignment-quality-reference",
                "report-regression-reference",
            ],
            [
                "phylo-inputs report",
                "alignment linkage report",
                "phylo html section contract",
            ],
            [
                "tree and alignment fixtures are validated separately before they are linked together"
            ],
        ),
    ]
    rows: list[CoreWorkflowValidationRow] = []
    for workflow, suite_names, outputs, limitations in definitions:
        fixture_count = sum(suite_map[name].fixture_count for name in suite_names)
        passed = all(suite_map[name].passed for name in suite_names)
        rows.append(
            CoreWorkflowValidationRow(
                workflow=workflow,
                fixture_suite_names=suite_names,
                fixture_count=fixture_count,
                expected_outputs=outputs,
                limitations=limitations,
                passed=passed,
                notes=[
                    f"draws fixture coverage from {', '.join(suite_names)}",
                    "expected outputs and limitations are explicitly listed for reviewer-facing traceability",
                ],
            )
        )
    return rows


def build_core_workflow_failure_gallery(
    *, fixtures_root: Path | None = None
) -> list[CoreWorkflowFailureCase]:
    """Document known failure and warning cases for core workflows."""
    root = default_fixtures_root() if fixtures_root is None else fixtures_root
    cases: list[CoreWorkflowFailureCase] = []

    duplicate_tree = fixture(root, "trees", "example_tree_duplicate.nwk")
    try:
        validate_tree_path(duplicate_tree)
    except Exception as error:
        cases.append(
            CoreWorkflowFailureCase(
                workflow="tree-review",
                fixture_name=duplicate_tree.name,
                outcome_kind="error",
                observed_code=getattr(error, "code", type(error).__name__),
                observed_summary=str(error),
                passed=getattr(error, "code", "") == "duplicate_taxon_error",
            )
        )

    invalid_alignment = fixture(
        root, "alignments", "example_alignment_invalid_lengths.fasta"
    )
    try:
        summarise_fasta(invalid_alignment)
    except Exception as error:
        cases.append(
            CoreWorkflowFailureCase(
                workflow="alignment-review",
                fixture_name=invalid_alignment.name,
                outcome_kind="error",
                observed_code=getattr(error, "code", type(error).__name__),
                observed_summary=str(error),
                passed=getattr(error, "code", "") == "invalid_alignment_error",
            )
        )

    dataset_report = audit_dataset_inputs(
        fixture(root, "trees", "example_taxon_workflow_tree.nwk"),
        fixture(root, "metadata", "example_taxon_workflow_metadata.csv"),
        fixture(root, "metadata", "example_taxon_workflow_traits.csv"),
        alignment_path=fixture(
            root, "alignments", "example_taxon_workflow_alignment.fasta"
        ),
    )
    cases.append(
        CoreWorkflowFailureCase(
            workflow="dataset-audit",
            fixture_name="example_taxon_workflow_*",
            outcome_kind="blocked",
            observed_code=dataset_report.readiness_decision,
            observed_summary="; ".join(dataset_report.blockers),
            passed=dataset_report.readiness_decision == "blocked"
            and len(dataset_report.blockers) >= 2,
        )
    )

    figure_report = build_tree_figure_package(
        fixture(root, "trees", "example_tree_support_invalid.nwk"),
        out_dir=temp_reference_dir("bijux-failure-gallery-figure"),
        show_support_values=True,
    )
    cases.append(
        CoreWorkflowFailureCase(
            workflow="figure-package",
            fixture_name="example_tree_support_invalid.nwk",
            outcome_kind="warning",
            observed_code="support_withheld",
            observed_summary="; ".join(figure_report.audit.support_audit.warnings),
            passed=figure_report.audit.support_audit.validated is False
            and figure_report.render.rendered_support_count == 0,
        )
    )
    return cases


def classify_core_workflow_maturity(
    *, fixtures_root: Path | None = None
) -> list[WorkflowMaturityClassification]:
    """Assign a maturity label to each Level 1 workflow."""
    _ = fixtures_root
    return [
        WorkflowMaturityClassification(
            workflow="tree-review",
            maturity="production-capable",
            rationale=[
                "tree validation and regression fixtures cover stable success and failure cases",
                "reviewer-facing HTML and machine manifests are deterministic",
            ],
            outstanding_risks=[
                "multi-tree or posterior uncertainty remains a separate review surface"
            ],
        ),
        WorkflowMaturityClassification(
            workflow="taxon-audit",
            maturity="usable",
            rationale=[
                "taxon naming, rank, namespace, and synonym conflicts are fixture-backed",
                "taxonomy reports expose explicit reviewer warnings and conflict rows",
            ],
            outstanding_risks=[
                "taxonomy completeness still depends on the supplied synonym authority tables"
            ],
        ),
        WorkflowMaturityClassification(
            workflow="alignment-review",
            maturity="production-capable",
            rationale=[
                "alignment diagnostics are validated across clean, duplicate-heavy, ambiguity-heavy, and missingness-heavy fixtures",
                "reviewer-facing diagnostics are bundled in stable report contracts",
            ],
            outstanding_risks=[
                "quality scores remain heuristic and should be interpreted with domain context"
            ],
        ),
        WorkflowMaturityClassification(
            workflow="dataset-audit",
            maturity="usable",
            rationale=[
                "dataset mismatch, exclusion, and blocked-analysis logic are pinned to a checked-in workflow example",
                "dataset reports expose ledgers, findings, exclusions, and reviewer checklists",
            ],
            outstanding_risks=[
                "fixture coverage currently emphasizes traceability over broad real-world dataset diversity"
            ],
        ),
        WorkflowMaturityClassification(
            workflow="figure-package",
            maturity="usable",
            rationale=[
                "figure packages preserve render metadata and support-label audits across safe and unsafe fixtures",
                "invalid support labels are explicitly withheld instead of silently rendered",
            ],
            outstanding_risks=[
                "visual regression is validated through render metadata rather than screenshot goldens"
            ],
        ),
        WorkflowMaturityClassification(
            workflow="phylo-inputs-review",
            maturity="production-capable",
            rationale=[
                "tree and alignment fixtures are independently validated and then surfaced together in a reviewer report",
                "stable section contracts keep combined-input reports reviewable across releases",
            ],
            outstanding_risks=[
                "combined tree-alignment trust still inherits the limitations of the underlying tree and alignment heuristics"
            ],
        ),
    ]


def build_core_workflow_validation_report(
    *, fixtures_root: Path | None = None
) -> CoreWorkflowValidationReport:
    """Build the Level 1 workflow validation report across goals 91 to 99."""
    root = default_fixtures_root() if fixtures_root is None else fixtures_root
    suites = [
        validate_tree_reference_fixtures(fixtures_root=root),
        validate_taxon_naming_reference_fixtures(fixtures_root=root),
        validate_alignment_quality_reference_fixtures(fixtures_root=root),
        validate_dataset_audit_reference_fixtures(fixtures_root=root),
        validate_figure_reference_fixtures(fixtures_root=root),
        validate_report_regression_fixtures(fixtures_root=root),
    ]
    workflows = _core_workflow_validation_rows(suites)
    failure_gallery = build_core_workflow_failure_gallery(fixtures_root=root)
    maturity_classifications = classify_core_workflow_maturity(fixtures_root=root)
    total_fixture_count = sum(suite.fixture_count for suite in suites)
    passed_fixture_count = sum(suite.passed_fixture_count for suite in suites)
    failed_fixture_count = total_fixture_count - passed_fixture_count
    limitations = sorted(
        {limitation for suite in suites for limitation in suite.limitations}
    )
    return CoreWorkflowValidationReport(
        suites=suites,
        workflows=workflows,
        failure_gallery=failure_gallery,
        maturity_classifications=maturity_classifications,
        total_fixture_count=total_fixture_count,
        passed_fixture_count=passed_fixture_count,
        failed_fixture_count=failed_fixture_count,
        limitations=limitations,
    )


def build_level_one_release_gate_report(
    *, fixtures_root: Path | None = None
) -> LevelOneReleaseGateReport:
    """Audit the checked-in Level 1 workflow end to end for reviewer traceability."""
    root = default_fixtures_root() if fixtures_root is None else fixtures_root
    validation = build_core_workflow_validation_report(fixtures_root=root)
    tree_path = fixture(root, "trees", "example_taxon_workflow_tree.nwk")
    metadata_path = fixture(root, "metadata", "example_taxon_workflow_metadata.csv")
    traits_path = fixture(root, "metadata", "example_taxon_workflow_traits.csv")
    alignment_path = fixture(
        root, "alignments", "example_taxon_workflow_alignment.fasta"
    )
    filtered_alignment_path = fixture(
        root, "alignments", "example_taxon_workflow_filtered_alignment.fasta"
    )
    inference_tree_path = fixture(root, "trees", "example_taxon_workflow_inference.nwk")
    reported_taxa_path = fixture(
        root, "metadata", "example_taxon_workflow_reported.csv"
    )

    dataset = audit_dataset_inputs(
        tree_path,
        metadata_path,
        traits_path,
        alignment_path=alignment_path,
    )
    taxon_loss = build_taxon_workflow_loss_report(
        tree_path,
        metadata_path,
        traits_path,
        alignment_path=alignment_path,
        filtered_alignment_path=filtered_alignment_path,
        inference_tree_path=inference_tree_path,
        reported_taxa_path=reported_taxa_path,
    )

    retained_taxa = sorted(
        row.taxon for row in taxon_loss.rows if row.first_loss_stage is None
    )
    excluded_taxa = sorted(
        row.taxon for row in taxon_loss.rows if row.first_loss_stage is not None
    )
    exclusion_causes = {
        row.taxon: [event.stage for event in row.loss_events]
        for row in taxon_loss.rows
        if row.loss_events
    }
    taxon_first_loss_stage = {
        row.taxon: row.first_loss_stage for row in taxon_loss.rows
    }
    decision = "blocked"
    rationale = [
        f"dataset readiness is {dataset.readiness_decision}",
        f"fixture validation passed {validation.passed_fixture_count} of {validation.total_fixture_count} checks",
        "every excluded taxon has an explicit first-loss stage and loss-event trail",
    ]
    if validation.failed_fixture_count > 0:
        rationale.append("one or more validation fixtures are currently failing")
    if dataset.readiness_decision != "blocked" and validation.failed_fixture_count == 0:
        decision = "pass"
    gate = LevelOneReleaseGateDecision(
        decision=decision,
        rationale=rationale,
        retained_taxa=retained_taxa,
        excluded_taxa=excluded_taxa,
        blocked_analyses=dataset.blocked_analyses,
        allowed_analyses=dataset.allowed_analyses,
        reviewer_visible_warnings=dataset.warnings,
    )
    return LevelOneReleaseGateReport(
        fixtures_root=root,
        validation=validation,
        dataset_readiness_decision=dataset.readiness_decision,
        dataset_blockers=dataset.blockers,
        dataset_warnings=dataset.warnings,
        exclusion_causes=exclusion_causes,
        taxon_first_loss_stage=taxon_first_loss_stage,
        gate=gate,
    )


def write_core_workflow_validation_json(
    path: Path,
    *,
    fixtures_root: Path | None = None,
) -> Path:
    """Write the aggregate workflow validation report as deterministic JSON."""
    report = build_core_workflow_validation_report(fixtures_root=fixtures_root)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(asdict(report), default=str, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def write_level_one_release_gate_json(
    path: Path,
    *,
    fixtures_root: Path | None = None,
) -> Path:
    """Write the Level 1 release gate report as deterministic JSON."""
    report = build_level_one_release_gate_report(fixtures_root=fixtures_root)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(asdict(report), default=str, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path
