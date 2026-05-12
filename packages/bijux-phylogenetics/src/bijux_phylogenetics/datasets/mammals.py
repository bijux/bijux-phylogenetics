from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import shutil

from bijux_phylogenetics.ancestral import (
    ContinuousAncestralReport,
    DiscreteAncestralReport,
    reconstruct_continuous_ancestral_states,
    reconstruct_discrete_ancestral_states,
    summarize_continuous_ancestral_report,
    summarize_discrete_ancestral_report,
    write_continuous_ancestral_exclusion_table,
    write_continuous_ancestral_summary_table,
    write_continuous_ancestral_uncertainty_table,
    write_discrete_ancestral_exclusion_table,
    write_discrete_ancestral_probability_table,
    write_discrete_ancestral_summary_table,
)
from bijux_phylogenetics.comparative import (
    BrownianTraitEvolutionSummaryReport,
    OUTraitEvolutionSummaryReport,
    PGLSResult,
    PhylogeneticSignalSummaryReport,
    run_pgls,
    summarize_brownian_trait_evolution,
    summarize_ou_trait_evolution,
    summarize_phylogenetic_signal,
    write_brownian_trait_evolution_exclusion_table,
    write_brownian_trait_evolution_summary_table,
    write_ou_trait_evolution_exclusion_table,
    write_ou_trait_evolution_summary_table,
    write_phylogenetic_signal_permutation_table,
    write_phylogenetic_signal_summary_table,
    write_pgls_lambda_profile_table,
)
from bijux_phylogenetics.core.metadata import load_taxon_table, write_taxon_rows

_DATASET_ID = "primate_comparative"
_DATASET_LABEL = "Primate comparative mammal dataset"
_TAXON_COLUMN = "species"
_CONTINUOUS_TRAIT = "longevity"
_PGLS_PREDICTOR = "social_group_size"
_DISCRETE_TRAIT = "mating_system"
_SIGNAL_PERMUTATIONS = 11
_SIGNAL_SEED = 7


@dataclass(slots=True)
class PrimateComparativeDataset:
    """Packaged real mammal dataset for comparative workflow review."""

    dataset_id: str
    label: str
    dataset_root: Path
    tree_path: Path
    traits_path: Path
    reference_output_root: Path
    taxon_column: str
    taxon_count: int
    continuous_traits: tuple[str, ...]
    categorical_traits: tuple[str, ...]
    workflow_continuous_trait: str
    workflow_pgls_predictor: str
    workflow_discrete_trait: str
    source_locator: str


@dataclass(slots=True)
class PrimateComparativeDatasetExportResult:
    """Materialized copy of the packaged primate dataset."""

    output_root: Path
    readme_path: Path
    tree_path: Path
    traits_path: Path
    expected_output_root: Path


@dataclass(slots=True)
class PrimateComparativeWorkflowReport:
    """One governed comparative workflow run over the packaged primate dataset."""

    dataset: PrimateComparativeDataset
    pgls: PGLSResult
    brownian: BrownianTraitEvolutionSummaryReport
    ou: OUTraitEvolutionSummaryReport
    signal: PhylogeneticSignalSummaryReport
    continuous_ancestral: ContinuousAncestralReport
    discrete_ancestral: DiscreteAncestralReport


@dataclass(slots=True)
class PrimateComparativeWorkflowBundle:
    """Written comparative workflow outputs for the packaged primate dataset."""

    output_root: Path
    summary_path: Path
    pgls_lambda_profile_path: Path
    brownian_summary_path: Path
    brownian_exclusion_path: Path
    ou_summary_path: Path
    ou_exclusion_path: Path
    signal_summary_path: Path
    signal_permutations_path: Path
    continuous_ancestral_summary_path: Path
    continuous_ancestral_uncertainty_path: Path
    continuous_ancestral_exclusion_path: Path
    discrete_ancestral_summary_path: Path
    discrete_ancestral_probability_path: Path
    discrete_ancestral_exclusion_path: Path


@dataclass(slots=True)
class PrimateComparativeDemoResult:
    """Dataset export plus workflow outputs for the public primate demo."""

    output_root: Path
    dataset: PrimateComparativeDataset
    dataset_export: PrimateComparativeDatasetExportResult
    workflow_bundle: PrimateComparativeWorkflowBundle
    overview_path: Path


def load_primate_comparative_dataset() -> PrimateComparativeDataset:
    """Expose the packaged real primate dataset as a first-class runtime surface."""
    dataset_root = _resource_root()
    traits_path = dataset_root / "traits.csv"
    table = load_taxon_table(traits_path, taxon_column=_TAXON_COLUMN)
    return PrimateComparativeDataset(
        dataset_id=_DATASET_ID,
        label=_DATASET_LABEL,
        dataset_root=dataset_root,
        tree_path=dataset_root / "tree.nwk",
        traits_path=traits_path,
        reference_output_root=dataset_root / "expected",
        taxon_column=_TAXON_COLUMN,
        taxon_count=table.row_count,
        continuous_traits=(
            "body_mass",
            "gestation",
            "home_range",
            "longevity",
            "social_group_size",
        ),
        categorical_traits=("family", "sex_dimorphism", "mating_system"),
        workflow_continuous_trait=_CONTINUOUS_TRAIT,
        workflow_pgls_predictor=_PGLS_PREDICTOR,
        workflow_discrete_trait=_DISCRETE_TRAIT,
        source_locator=(
            "evidence-book/studies/primate-longevity-signal/datasets/"
            "reference_primate.csv"
        ),
    )


def export_primate_comparative_dataset(
    destination: Path,
) -> PrimateComparativeDatasetExportResult:
    """Copy the packaged primate dataset and reference outputs to one directory."""
    dataset = load_primate_comparative_dataset()
    if destination.exists():
        shutil.rmtree(destination)
    destination.mkdir(parents=True, exist_ok=True)
    readme_path = shutil.copy2(dataset.dataset_root / "README.md", destination / "README.md")
    tree_path = shutil.copy2(dataset.tree_path, destination / "tree.nwk")
    traits_path = shutil.copy2(dataset.traits_path, destination / "traits.csv")
    expected_output_root = destination / "expected"
    shutil.copytree(dataset.reference_output_root, expected_output_root)
    return PrimateComparativeDatasetExportResult(
        output_root=destination,
        readme_path=Path(readme_path),
        tree_path=Path(tree_path),
        traits_path=Path(traits_path),
        expected_output_root=expected_output_root,
    )


def run_primate_comparative_workflow() -> PrimateComparativeWorkflowReport:
    """Run the owned comparative workflow over the packaged real primate dataset."""
    dataset = load_primate_comparative_dataset()
    pgls = run_pgls(
        dataset.tree_path,
        dataset.traits_path,
        response=dataset.workflow_continuous_trait,
        predictors=[dataset.workflow_pgls_predictor],
        taxon_column=dataset.taxon_column,
        lambda_value="estimate",
    )
    brownian = summarize_brownian_trait_evolution(
        dataset.tree_path,
        dataset.traits_path,
        trait=dataset.workflow_continuous_trait,
        taxon_column=dataset.taxon_column,
    )
    ou = summarize_ou_trait_evolution(
        dataset.tree_path,
        dataset.traits_path,
        trait=dataset.workflow_continuous_trait,
        taxon_column=dataset.taxon_column,
    )
    signal = summarize_phylogenetic_signal(
        dataset.tree_path,
        dataset.traits_path,
        trait=dataset.workflow_continuous_trait,
        taxon_column=dataset.taxon_column,
        permutations=_SIGNAL_PERMUTATIONS,
        seed=_SIGNAL_SEED,
    )
    continuous_ancestral = reconstruct_continuous_ancestral_states(
        dataset.tree_path,
        dataset.traits_path,
        trait=dataset.workflow_continuous_trait,
        taxon_column=dataset.taxon_column,
        model="brownian",
    )
    discrete_ancestral = reconstruct_discrete_ancestral_states(
        dataset.tree_path,
        dataset.traits_path,
        trait=dataset.workflow_discrete_trait,
        taxon_column=dataset.taxon_column,
        model="equal-rates",
    )
    return PrimateComparativeWorkflowReport(
        dataset=dataset,
        pgls=pgls,
        brownian=brownian,
        ou=ou,
        signal=signal,
        continuous_ancestral=continuous_ancestral,
        discrete_ancestral=discrete_ancestral,
    )


def write_primate_comparative_workflow_bundle(
    output_root: Path, report: PrimateComparativeWorkflowReport
) -> PrimateComparativeWorkflowBundle:
    """Write the governed comparative workflow outputs for the packaged dataset."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    summary_path = _write_workflow_summary_table(output_root / "workflow-summary.tsv", report)
    pgls_lambda_profile_path = write_pgls_lambda_profile_table(
        output_root / "pgls-lambda-profile.tsv",
        report.pgls.lambda_fit,
    )
    brownian_summary_path = write_brownian_trait_evolution_summary_table(
        output_root / "brownian-summary.tsv",
        report.brownian,
    )
    brownian_exclusion_path = write_brownian_trait_evolution_exclusion_table(
        output_root / "brownian-excluded.tsv",
        report.brownian,
    )
    ou_summary_path = write_ou_trait_evolution_summary_table(
        output_root / "ou-summary.tsv",
        report.ou,
    )
    ou_exclusion_path = write_ou_trait_evolution_exclusion_table(
        output_root / "ou-excluded.tsv",
        report.ou,
    )
    signal_summary_path = write_phylogenetic_signal_summary_table(
        output_root / "signal-summary.tsv",
        report.signal,
    )
    signal_permutations_path = write_phylogenetic_signal_permutation_table(
        output_root / "signal-permutations.tsv",
        report.signal,
    )
    continuous_ancestral_summary_path = write_continuous_ancestral_summary_table(
        output_root / "continuous-ancestral-summary.tsv",
        report.continuous_ancestral,
    )
    continuous_ancestral_uncertainty_path = write_continuous_ancestral_uncertainty_table(
        output_root / "continuous-ancestral-uncertainty.tsv",
        report.continuous_ancestral,
    )
    continuous_ancestral_exclusion_path = write_continuous_ancestral_exclusion_table(
        output_root / "continuous-ancestral-excluded.tsv",
        report.continuous_ancestral,
    )
    discrete_ancestral_summary_path = write_discrete_ancestral_summary_table(
        output_root / "discrete-ancestral-summary.tsv",
        report.discrete_ancestral,
    )
    discrete_ancestral_probability_path = write_discrete_ancestral_probability_table(
        output_root / "discrete-ancestral-probabilities.tsv",
        report.discrete_ancestral,
    )
    discrete_ancestral_exclusion_path = write_discrete_ancestral_exclusion_table(
        output_root / "discrete-ancestral-excluded.tsv",
        report.discrete_ancestral,
    )
    return PrimateComparativeWorkflowBundle(
        output_root=output_root,
        summary_path=summary_path,
        pgls_lambda_profile_path=pgls_lambda_profile_path,
        brownian_summary_path=brownian_summary_path,
        brownian_exclusion_path=brownian_exclusion_path,
        ou_summary_path=ou_summary_path,
        ou_exclusion_path=ou_exclusion_path,
        signal_summary_path=signal_summary_path,
        signal_permutations_path=signal_permutations_path,
        continuous_ancestral_summary_path=continuous_ancestral_summary_path,
        continuous_ancestral_uncertainty_path=continuous_ancestral_uncertainty_path,
        continuous_ancestral_exclusion_path=continuous_ancestral_exclusion_path,
        discrete_ancestral_summary_path=discrete_ancestral_summary_path,
        discrete_ancestral_probability_path=discrete_ancestral_probability_path,
        discrete_ancestral_exclusion_path=discrete_ancestral_exclusion_path,
    )


def run_primate_comparative_demo(output_root: Path) -> PrimateComparativeDemoResult:
    """Materialize the packaged dataset and its governed workflow outputs."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    report = run_primate_comparative_workflow()
    dataset_export = export_primate_comparative_dataset(output_root / "dataset")
    workflow_bundle = write_primate_comparative_workflow_bundle(
        output_root / "workflow",
        report,
    )
    overview_path = _write_demo_overview(output_root / "overview.md", workflow_bundle)
    return PrimateComparativeDemoResult(
        output_root=output_root,
        dataset=report.dataset,
        dataset_export=dataset_export,
        workflow_bundle=workflow_bundle,
        overview_path=overview_path,
    )


def _resource_root() -> Path:
    return (
        Path(__file__).resolve().parents[1]
        / "resources"
        / "datasets"
        / "mammals"
        / "primate_comparative"
    )


def _write_workflow_summary_table(
    path: Path, report: PrimateComparativeWorkflowReport
) -> Path:
    slope = next(
        coefficient
        for coefficient in report.pgls.coefficients
        if coefficient.name == report.dataset.workflow_pgls_predictor
    )
    continuous_summary = summarize_continuous_ancestral_report(
        report.continuous_ancestral
    )
    discrete_summary = summarize_discrete_ancestral_report(report.discrete_ancestral)
    continuous_root = next(
        estimate
        for estimate in report.continuous_ancestral.estimates
        if not estimate.is_tip and estimate.node == continuous_summary.root_node
    )
    discrete_root = next(
        estimate
        for estimate in report.discrete_ancestral.estimates
        if not estimate.is_tip and estimate.node == discrete_summary.root_node
    )
    return write_taxon_rows(
        path,
        columns=[
            "dataset_id",
            "taxon_count",
            "taxon_column",
            "continuous_trait",
            "pgls_predictor",
            "pgls_lambda",
            "pgls_predictor_estimate",
            "pgls_predictor_p_value",
            "brownian_sigma_squared",
            "brownian_aicc",
            "ou_alpha",
            "ou_theta",
            "ou_sigma_squared",
            "ou_aicc",
            "signal_blombergs_k",
            "signal_pagels_lambda",
            "signal_permutation_p_value",
            "continuous_root_node",
            "continuous_root_estimate",
            "discrete_trait",
            "discrete_root_node",
            "discrete_root_state",
            "discrete_root_confidence",
        ],
        rows=[
            {
                "dataset_id": report.dataset.dataset_id,
                "taxon_count": str(report.dataset.taxon_count),
                "taxon_column": report.dataset.taxon_column,
                "continuous_trait": report.dataset.workflow_continuous_trait,
                "pgls_predictor": report.dataset.workflow_pgls_predictor,
                "pgls_lambda": format(report.pgls.lambda_value, ".15g"),
                "pgls_predictor_estimate": format(slope.estimate, ".15g"),
                "pgls_predictor_p_value": format(slope.p_value, ".15g"),
                "brownian_sigma_squared": format(report.brownian.sigma_squared, ".15g"),
                "brownian_aicc": format(report.brownian.aicc, ".15g"),
                "ou_alpha": format(report.ou.alpha, ".15g"),
                "ou_theta": format(report.ou.theta, ".15g"),
                "ou_sigma_squared": format(report.ou.sigma_squared, ".15g"),
                "ou_aicc": format(report.ou.aicc, ".15g"),
                "signal_blombergs_k": format(report.signal.blombergs_k.k, ".15g"),
                "signal_pagels_lambda": format(
                    report.signal.pagels_lambda.lambda_value, ".15g"
                ),
                "signal_permutation_p_value": format(
                    report.signal.signal_test.p_value, ".15g"
                ),
                "continuous_root_node": continuous_root.node,
                "continuous_root_estimate": format(
                    continuous_root.estimate,
                    ".15g",
                ),
                "discrete_trait": report.dataset.workflow_discrete_trait,
                "discrete_root_node": discrete_root.node,
                "discrete_root_state": discrete_root.most_likely_state,
                "discrete_root_confidence": format(discrete_root.confidence, ".15g"),
            }
        ],
    )


def _write_demo_overview(path: Path, bundle: PrimateComparativeWorkflowBundle) -> Path:
    lines = [
        "# Primate Comparative Demo",
        "",
        "This demo materializes the packaged mammal dataset and regenerates the",
        "governed comparative workflow outputs that ship with the package.",
        "",
        "## Workflow Outputs",
        "",
        f"- workflow summary: `{bundle.summary_path}`",
        f"- PGLS lambda profile: `{bundle.pgls_lambda_profile_path}`",
        f"- Brownian summary: `{bundle.brownian_summary_path}`",
        f"- OU summary: `{bundle.ou_summary_path}`",
        f"- signal summary: `{bundle.signal_summary_path}`",
        f"- continuous ancestral summary: `{bundle.continuous_ancestral_summary_path}`",
        f"- discrete ancestral summary: `{bundle.discrete_ancestral_summary_path}`",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")
    return path
