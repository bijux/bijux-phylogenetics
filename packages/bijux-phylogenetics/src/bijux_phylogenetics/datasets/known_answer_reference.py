from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory

from bijux_phylogenetics.ancestral.continuous import (
    ContinuousAncestralReport,
    reconstruct_continuous_ancestral_states,
    write_continuous_ancestral_summary_table,
    write_continuous_ancestral_uncertainty_table,
)
from bijux_phylogenetics.ancestral.discrete import (
    DiscreteAncestralReport,
    reconstruct_discrete_ancestral_states,
    write_discrete_ancestral_probability_table,
    write_discrete_ancestral_summary_table,
)
from bijux_phylogenetics.comparative.brownian_trait_evolution import (
    BrownianTraitEvolutionSummaryReport,
    summarize_brownian_trait_evolution,
    write_brownian_trait_evolution_summary_table,
)
from bijux_phylogenetics.compare.topology import TreeComparisonReport, compare_tree_paths
from bijux_phylogenetics.core.metadata import write_taxon_rows
from bijux_phylogenetics.distance import DistanceTreeBuildReport, build_distance_tree
from bijux_phylogenetics.io.fasta import load_fasta_alignment, validate_fasta_input
from bijux_phylogenetics.io.newick import dumps_newick, write_newick
from bijux_phylogenetics.io.trees import load_tree

_DATASET_ID = "known_answer_reference_panel"
_DATASET_LABEL = "Known-answer simulation reference panel"
_SEQUENCE_TYPE = "dna"
_DISTANCE_METHOD = "neighbor-joining"
_DISTANCE_MODEL = "p-distance"
_CONTINUOUS_TRAIT = "value"
_DISCRETE_TRAIT = "state"


@dataclass(slots=True)
class KnownAnswerContinuousNodeTruth:
    """One true continuous node value stored with the packaged simulation panel."""

    node: str
    node_name: str | None
    is_tip: bool
    descendant_taxa: list[str]
    true_value: float


@dataclass(slots=True)
class KnownAnswerDiscreteNodeTruth:
    """One true discrete node state stored with the packaged simulation panel."""

    node: str
    node_name: str | None
    is_tip: bool
    descendant_taxa: list[str]
    true_state: str


@dataclass(slots=True)
class KnownAnswerParameterRecoveryRow:
    """One recovered continuous parameter compared directly against truth."""

    parameter: str
    true_value: float
    estimated_value: float
    absolute_error: float
    relative_error: float
    interpretation: str


@dataclass(slots=True)
class KnownAnswerContinuousNodeRecoveryRow:
    """One internal-node continuous ancestral estimate compared against truth."""

    node: str
    descendant_taxa: list[str]
    true_value: float
    estimated_value: float
    absolute_error: float
    standard_error: float
    lower_95_interval: float
    upper_95_interval: float
    confidence: float


@dataclass(slots=True)
class KnownAnswerDiscreteNodeRecoveryRow:
    """One internal-node discrete ancestral estimate compared against truth."""

    node: str
    descendant_taxa: list[str]
    true_state: str
    estimated_state: str
    true_state_probability: float
    confidence: float
    correct: bool
    ambiguous: bool


@dataclass(slots=True)
class KnownAnswerReferenceDataset:
    """Packaged deterministic simulation panel with stored truth artifacts."""

    dataset_id: str
    label: str
    dataset_root: Path
    true_tree_path: Path
    alignment_path: Path
    continuous_traits_path: Path
    discrete_traits_path: Path
    true_parameters_path: Path
    true_continuous_nodes_path: Path
    true_discrete_nodes_path: Path
    reference_output_root: Path
    taxon_count: int
    sequence_length: int
    sequence_type: str
    distance_method: str
    distance_model: str
    source_summary: str


@dataclass(slots=True)
class KnownAnswerReferenceExportResult:
    """Materialized copy of the packaged simulation dataset."""

    output_root: Path
    readme_path: Path
    true_tree_path: Path
    alignment_path: Path
    continuous_traits_path: Path
    discrete_traits_path: Path
    true_parameters_path: Path
    true_continuous_nodes_path: Path
    true_discrete_nodes_path: Path
    expected_output_root: Path


@dataclass(slots=True)
class KnownAnswerReferenceWorkflowReport:
    """Recovery workflow run over the packaged known-answer simulation panel."""

    dataset: KnownAnswerReferenceDataset
    distance_tree_build: DistanceTreeBuildReport
    distance_tree_newick: str
    tree_recovery: TreeComparisonReport
    brownian_fit: BrownianTraitEvolutionSummaryReport
    continuous_ancestral: ContinuousAncestralReport
    discrete_ancestral: DiscreteAncestralReport
    parameter_recovery_rows: list[KnownAnswerParameterRecoveryRow]
    continuous_node_recovery_rows: list[KnownAnswerContinuousNodeRecoveryRow]
    discrete_node_recovery_rows: list[KnownAnswerDiscreteNodeRecoveryRow]


@dataclass(slots=True)
class KnownAnswerReferenceWorkflowBundle:
    """Written recovery outputs for the packaged known-answer simulation panel."""

    output_root: Path
    rooted_topology_equal: bool
    same_unrooted_topology: bool
    same_taxa_different_rooting: bool
    robinson_foulds_distance: int
    continuous_internal_node_mean_absolute_error: float
    discrete_internal_node_accuracy: float
    discrete_mean_true_state_probability: float
    workflow_summary_path: Path
    distance_tree_path: Path
    tree_recovery_path: Path
    parameter_recovery_path: Path
    brownian_fit_summary_path: Path
    continuous_ancestral_summary_path: Path
    continuous_ancestral_uncertainty_path: Path
    continuous_node_recovery_path: Path
    discrete_ancestral_summary_path: Path
    discrete_ancestral_probability_path: Path
    discrete_node_recovery_path: Path


@dataclass(slots=True)
class KnownAnswerReferenceDemoResult:
    """Dataset export plus recovery workflow outputs for the public simulation demo."""

    output_root: Path
    dataset: KnownAnswerReferenceDataset
    dataset_export: KnownAnswerReferenceExportResult
    workflow_bundle: KnownAnswerReferenceWorkflowBundle
    overview_path: Path


def load_known_answer_reference_dataset() -> KnownAnswerReferenceDataset:
    """Expose the packaged deterministic simulation panel as a first-class surface."""
    dataset_root = _resource_root()
    true_tree_path = dataset_root / "true-tree.nwk"
    alignment_path = dataset_root / "simulated-alignment.fasta"
    continuous_traits_path = dataset_root / "continuous-traits.tsv"
    discrete_traits_path = dataset_root / "discrete-traits.tsv"
    true_parameters_path = dataset_root / "true-parameters.tsv"
    true_continuous_nodes_path = dataset_root / "true-continuous-nodes.tsv"
    true_discrete_nodes_path = dataset_root / "true-discrete-nodes.tsv"
    validate_fasta_input(alignment_path, sequence_type=_SEQUENCE_TYPE)
    records = load_fasta_alignment(alignment_path)
    tree = load_tree(true_tree_path)
    return KnownAnswerReferenceDataset(
        dataset_id=_DATASET_ID,
        label=_DATASET_LABEL,
        dataset_root=dataset_root,
        true_tree_path=true_tree_path,
        alignment_path=alignment_path,
        continuous_traits_path=continuous_traits_path,
        discrete_traits_path=discrete_traits_path,
        true_parameters_path=true_parameters_path,
        true_continuous_nodes_path=true_continuous_nodes_path,
        true_discrete_nodes_path=true_discrete_nodes_path,
        reference_output_root=dataset_root / "expected",
        taxon_count=tree.tip_count,
        sequence_length=len(records[0].sequence),
        sequence_type=_SEQUENCE_TYPE,
        distance_method=_DISTANCE_METHOD,
        distance_model=_DISTANCE_MODEL,
        source_summary=(
            "Deterministic owned simulation panel with one birth-death tree, one "
            "JC-like DNA alignment, one Brownian continuous trait, and one "
            "symmetric discrete trait, packaged with full node-level truth ledgers "
            "for recovery review."
        ),
    )


def export_known_answer_reference_dataset(
    destination: Path,
) -> KnownAnswerReferenceExportResult:
    """Copy the packaged known-answer simulation dataset and reference outputs."""
    dataset = load_known_answer_reference_dataset()
    if destination.exists():
        shutil.rmtree(destination)
    destination.mkdir(parents=True, exist_ok=True)
    readme_path = shutil.copy2(dataset.dataset_root / "README.md", destination / "README.md")
    true_tree_path = shutil.copy2(dataset.true_tree_path, destination / "true-tree.nwk")
    alignment_path = shutil.copy2(
        dataset.alignment_path,
        destination / "simulated-alignment.fasta",
    )
    continuous_traits_path = shutil.copy2(
        dataset.continuous_traits_path,
        destination / "continuous-traits.tsv",
    )
    discrete_traits_path = shutil.copy2(
        dataset.discrete_traits_path,
        destination / "discrete-traits.tsv",
    )
    true_parameters_path = shutil.copy2(
        dataset.true_parameters_path,
        destination / "true-parameters.tsv",
    )
    true_continuous_nodes_path = shutil.copy2(
        dataset.true_continuous_nodes_path,
        destination / "true-continuous-nodes.tsv",
    )
    true_discrete_nodes_path = shutil.copy2(
        dataset.true_discrete_nodes_path,
        destination / "true-discrete-nodes.tsv",
    )
    expected_output_root = destination / "expected"
    shutil.copytree(dataset.reference_output_root, expected_output_root)
    return KnownAnswerReferenceExportResult(
        output_root=destination,
        readme_path=Path(readme_path),
        true_tree_path=Path(true_tree_path),
        alignment_path=Path(alignment_path),
        continuous_traits_path=Path(continuous_traits_path),
        discrete_traits_path=Path(discrete_traits_path),
        true_parameters_path=Path(true_parameters_path),
        true_continuous_nodes_path=Path(true_continuous_nodes_path),
        true_discrete_nodes_path=Path(true_discrete_nodes_path),
        expected_output_root=expected_output_root,
    )


def run_known_answer_reference_workflow() -> KnownAnswerReferenceWorkflowReport:
    """Run the governed recovery workflow over the packaged known-answer panel."""
    dataset = load_known_answer_reference_dataset()
    true_parameters = _load_true_parameter_map(dataset.true_parameters_path)
    continuous_truth = _load_true_continuous_nodes(dataset.true_continuous_nodes_path)
    discrete_truth = _load_true_discrete_nodes(dataset.true_discrete_nodes_path)

    distance_tree, distance_tree_build = build_distance_tree(
        dataset.alignment_path,
        method=dataset.distance_method,
        model=dataset.distance_model,
    )
    distance_tree_newick = dumps_newick(distance_tree)
    with TemporaryDirectory(prefix="known-answer-reference-") as temporary_root:
        built_tree_path = write_newick(
            Path(temporary_root) / "recovered-distance-tree.nwk",
            distance_tree,
        )
        tree_recovery = compare_tree_paths(built_tree_path, dataset.true_tree_path)

    brownian_fit = summarize_brownian_trait_evolution(
        dataset.true_tree_path,
        dataset.continuous_traits_path,
        trait=_CONTINUOUS_TRAIT,
    )
    continuous_ancestral = reconstruct_continuous_ancestral_states(
        dataset.true_tree_path,
        dataset.continuous_traits_path,
        trait=_CONTINUOUS_TRAIT,
        model="brownian",
    )
    discrete_ancestral = reconstruct_discrete_ancestral_states(
        dataset.true_tree_path,
        dataset.discrete_traits_path,
        trait=_DISCRETE_TRAIT,
        model="equal-rates",
    )
    parameter_recovery_rows = _build_parameter_recovery_rows(
        true_parameters=true_parameters,
        brownian_fit=brownian_fit,
    )
    continuous_node_recovery_rows = _build_continuous_node_recovery_rows(
        true_nodes=continuous_truth,
        report=continuous_ancestral,
    )
    discrete_node_recovery_rows = _build_discrete_node_recovery_rows(
        true_nodes=discrete_truth,
        report=discrete_ancestral,
    )
    return KnownAnswerReferenceWorkflowReport(
        dataset=dataset,
        distance_tree_build=distance_tree_build,
        distance_tree_newick=distance_tree_newick,
        tree_recovery=tree_recovery,
        brownian_fit=brownian_fit,
        continuous_ancestral=continuous_ancestral,
        discrete_ancestral=discrete_ancestral,
        parameter_recovery_rows=parameter_recovery_rows,
        continuous_node_recovery_rows=continuous_node_recovery_rows,
        discrete_node_recovery_rows=discrete_node_recovery_rows,
    )


def write_known_answer_reference_workflow_bundle(
    output_root: Path,
    report: KnownAnswerReferenceWorkflowReport,
) -> KnownAnswerReferenceWorkflowBundle:
    """Write the governed recovery outputs for the packaged simulation panel."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    continuous_mae = _mean(
        row.absolute_error for row in report.continuous_node_recovery_rows
    )
    discrete_accuracy = _mean(
        1.0 if row.correct else 0.0 for row in report.discrete_node_recovery_rows
    )
    discrete_mean_true_probability = _mean(
        row.true_state_probability for row in report.discrete_node_recovery_rows
    )
    workflow_summary_path = _write_workflow_summary_table(
        output_root / "workflow-summary.tsv",
        report=report,
        continuous_mae=continuous_mae,
        discrete_accuracy=discrete_accuracy,
        discrete_mean_true_probability=discrete_mean_true_probability,
    )
    distance_tree_path = _write_distance_tree(
        output_root / "recovered-distance-tree.nwk",
        report.distance_tree_newick,
    )
    tree_recovery_path = _write_tree_recovery_table(
        output_root / "tree-recovery.tsv",
        report,
    )
    parameter_recovery_path = _write_parameter_recovery_table(
        output_root / "parameter-recovery.tsv",
        report.parameter_recovery_rows,
    )
    brownian_fit_summary_path = write_brownian_trait_evolution_summary_table(
        output_root / "brownian-fit-summary.tsv",
        report.brownian_fit,
    )
    continuous_ancestral_summary_path = write_continuous_ancestral_summary_table(
        output_root / "continuous-ancestral-summary.tsv",
        report.continuous_ancestral,
    )
    continuous_ancestral_uncertainty_path = write_continuous_ancestral_uncertainty_table(
        output_root / "continuous-ancestral-uncertainty.tsv",
        report.continuous_ancestral,
    )
    continuous_node_recovery_path = _write_continuous_node_recovery_table(
        output_root / "continuous-node-recovery.tsv",
        report.continuous_node_recovery_rows,
    )
    discrete_ancestral_summary_path = write_discrete_ancestral_summary_table(
        output_root / "discrete-ancestral-summary.tsv",
        report.discrete_ancestral,
    )
    discrete_ancestral_probability_path = write_discrete_ancestral_probability_table(
        output_root / "discrete-ancestral-probabilities.tsv",
        report.discrete_ancestral,
    )
    discrete_node_recovery_path = _write_discrete_node_recovery_table(
        output_root / "discrete-node-recovery.tsv",
        report.discrete_node_recovery_rows,
    )
    return KnownAnswerReferenceWorkflowBundle(
        output_root=output_root,
        rooted_topology_equal=report.tree_recovery.topology_equal,
        same_unrooted_topology=report.tree_recovery.same_unrooted_topology,
        same_taxa_different_rooting=report.tree_recovery.same_taxa_different_rooting,
        robinson_foulds_distance=report.tree_recovery.robinson_foulds_distance,
        continuous_internal_node_mean_absolute_error=continuous_mae,
        discrete_internal_node_accuracy=discrete_accuracy,
        discrete_mean_true_state_probability=discrete_mean_true_probability,
        workflow_summary_path=workflow_summary_path,
        distance_tree_path=distance_tree_path,
        tree_recovery_path=tree_recovery_path,
        parameter_recovery_path=parameter_recovery_path,
        brownian_fit_summary_path=brownian_fit_summary_path,
        continuous_ancestral_summary_path=continuous_ancestral_summary_path,
        continuous_ancestral_uncertainty_path=continuous_ancestral_uncertainty_path,
        continuous_node_recovery_path=continuous_node_recovery_path,
        discrete_ancestral_summary_path=discrete_ancestral_summary_path,
        discrete_ancestral_probability_path=discrete_ancestral_probability_path,
        discrete_node_recovery_path=discrete_node_recovery_path,
    )


def run_known_answer_reference_demo(output_root: Path) -> KnownAnswerReferenceDemoResult:
    """Materialize the packaged simulation dataset and rerun the recovery outputs."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    report = run_known_answer_reference_workflow()
    dataset_export = export_known_answer_reference_dataset(output_root / "dataset")
    workflow_bundle = write_known_answer_reference_workflow_bundle(
        output_root / "workflow",
        report,
    )
    overview_path = _write_overview(output_root / "overview.md", report, workflow_bundle)
    return KnownAnswerReferenceDemoResult(
        output_root=output_root,
        dataset=report.dataset,
        dataset_export=dataset_export,
        workflow_bundle=workflow_bundle,
        overview_path=overview_path,
    )


def _build_parameter_recovery_rows(
    *,
    true_parameters: dict[str, str],
    brownian_fit: BrownianTraitEvolutionSummaryReport,
) -> list[KnownAnswerParameterRecoveryRow]:
    true_root_state = float(true_parameters["continuous_root_state"])
    true_sigma_squared = float(true_parameters["continuous_sigma_squared"])
    return [
        _parameter_recovery_row(
            parameter="continuous_root_state",
            true_value=true_root_state,
            estimated_value=brownian_fit.root_state,
            interpretation="Brownian root-state estimate recovered from observed tip values on the true tree.",
        ),
        _parameter_recovery_row(
            parameter="continuous_sigma_squared",
            true_value=true_sigma_squared,
            estimated_value=brownian_fit.sigma_squared,
            interpretation="Brownian evolutionary rate recovered from observed tip values on the true tree.",
        ),
    ]


def _parameter_recovery_row(
    *,
    parameter: str,
    true_value: float,
    estimated_value: float,
    interpretation: str,
) -> KnownAnswerParameterRecoveryRow:
    absolute_error = abs(estimated_value - true_value)
    denominator = abs(true_value)
    relative_error = 0.0 if denominator == 0.0 else absolute_error / denominator
    return KnownAnswerParameterRecoveryRow(
        parameter=parameter,
        true_value=true_value,
        estimated_value=estimated_value,
        absolute_error=absolute_error,
        relative_error=relative_error,
        interpretation=interpretation,
    )


def _build_continuous_node_recovery_rows(
    *,
    true_nodes: list[KnownAnswerContinuousNodeTruth],
    report: ContinuousAncestralReport,
) -> list[KnownAnswerContinuousNodeRecoveryRow]:
    truth_by_node = {
        row.node: row
        for row in true_nodes
        if not row.is_tip
    }
    return [
        KnownAnswerContinuousNodeRecoveryRow(
            node=estimate.node,
            descendant_taxa=estimate.descendant_taxa,
            true_value=truth_by_node[estimate.node].true_value,
            estimated_value=estimate.estimate,
            absolute_error=abs(
                estimate.estimate - truth_by_node[estimate.node].true_value
            ),
            standard_error=estimate.standard_error,
            lower_95_interval=estimate.lower_95_interval,
            upper_95_interval=estimate.upper_95_interval,
            confidence=estimate.confidence,
        )
        for estimate in report.estimates
        if not estimate.is_tip
    ]


def _build_discrete_node_recovery_rows(
    *,
    true_nodes: list[KnownAnswerDiscreteNodeTruth],
    report: DiscreteAncestralReport,
) -> list[KnownAnswerDiscreteNodeRecoveryRow]:
    truth_by_node = {
        row.node: row
        for row in true_nodes
        if not row.is_tip
    }
    return [
        KnownAnswerDiscreteNodeRecoveryRow(
            node=estimate.node,
            descendant_taxa=estimate.descendant_taxa,
            true_state=truth_by_node[estimate.node].true_state,
            estimated_state=estimate.most_likely_state,
            true_state_probability=estimate.state_probabilities.get(
                truth_by_node[estimate.node].true_state,
                0.0,
            ),
            confidence=estimate.confidence,
            correct=estimate.most_likely_state == truth_by_node[estimate.node].true_state,
            ambiguous=estimate.ambiguous,
        )
        for estimate in report.estimates
        if not estimate.is_tip
    ]


def _write_workflow_summary_table(
    path: Path,
    *,
    report: KnownAnswerReferenceWorkflowReport,
    continuous_mae: float,
    discrete_accuracy: float,
    discrete_mean_true_probability: float,
) -> Path:
    rows = [
        "\t".join(
            [
                "dataset_id",
                "taxon_count",
                "sequence_length",
                "distance_method",
                "distance_model",
                "rooted_topology_equal",
                "same_unrooted_topology",
                "same_taxa_different_rooting",
                "robinson_foulds_distance",
                "continuous_root_absolute_error",
                "continuous_sigma_squared_absolute_error",
                "continuous_internal_node_mean_absolute_error",
                "discrete_internal_node_accuracy",
                "discrete_mean_true_state_probability",
            ]
        ),
        "\t".join(
            [
                report.dataset.dataset_id,
                str(report.dataset.taxon_count),
                str(report.dataset.sequence_length),
                report.dataset.distance_method,
                report.dataset.distance_model,
                str(report.tree_recovery.topology_equal).lower(),
                str(report.tree_recovery.same_unrooted_topology).lower(),
                str(report.tree_recovery.same_taxa_different_rooting).lower(),
                str(report.tree_recovery.robinson_foulds_distance),
                _format_number(report.parameter_recovery_rows[0].absolute_error),
                _format_number(report.parameter_recovery_rows[1].absolute_error),
                _format_number(continuous_mae),
                _format_number(discrete_accuracy),
                _format_number(discrete_mean_true_probability),
            ]
        ),
    ]
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def _write_distance_tree(path: Path, newick: str) -> Path:
    path.write_text(f"{newick}\n", encoding="utf-8")
    return path


def _write_tree_recovery_table(
    path: Path,
    report: KnownAnswerReferenceWorkflowReport,
) -> Path:
    topology = report.tree_recovery
    build = report.distance_tree_build
    return write_taxon_rows(
        path,
        columns=[
            "method",
            "model",
            "taxon_count",
            "pair_count",
            "rooted_topology_equal",
            "same_unrooted_topology",
            "same_taxa_different_rooting",
            "robinson_foulds_distance",
            "rooted_robinson_foulds_distance",
            "unrooted_robinson_foulds_distance",
            "normalized_robinson_foulds",
            "rooted_normalized_robinson_foulds",
            "unrooted_normalized_robinson_foulds",
        ],
        rows=[
            {
                "method": build.method,
                "model": build.model,
                "taxon_count": str(build.taxon_count),
                "pair_count": str(build.pair_count),
                "rooted_topology_equal": str(topology.topology_equal).lower(),
                "same_unrooted_topology": str(
                    topology.same_unrooted_topology
                ).lower(),
                "same_taxa_different_rooting": str(
                    topology.same_taxa_different_rooting
                ).lower(),
                "robinson_foulds_distance": str(topology.robinson_foulds_distance),
                "rooted_robinson_foulds_distance": str(
                    topology.rooted_robinson_foulds_distance
                ),
                "unrooted_robinson_foulds_distance": str(
                    topology.unrooted_robinson_foulds_distance
                ),
                "normalized_robinson_foulds": _format_number(
                    topology.normalized_robinson_foulds
                ),
                "rooted_normalized_robinson_foulds": _format_number(
                    topology.rooted_normalized_robinson_foulds
                ),
                "unrooted_normalized_robinson_foulds": _format_number(
                    topology.unrooted_normalized_robinson_foulds
                ),
            }
        ],
    )


def _write_parameter_recovery_table(
    path: Path,
    rows: list[KnownAnswerParameterRecoveryRow],
) -> Path:
    return write_taxon_rows(
        path,
        columns=[
            "parameter",
            "true_value",
            "estimated_value",
            "absolute_error",
            "relative_error",
            "interpretation",
        ],
        rows=[
            {
                "parameter": row.parameter,
                "true_value": _format_number(row.true_value),
                "estimated_value": _format_number(row.estimated_value),
                "absolute_error": _format_number(row.absolute_error),
                "relative_error": _format_number(row.relative_error),
                "interpretation": row.interpretation,
            }
            for row in rows
        ],
    )


def _write_continuous_node_recovery_table(
    path: Path,
    rows: list[KnownAnswerContinuousNodeRecoveryRow],
) -> Path:
    return write_taxon_rows(
        path,
        columns=[
            "node",
            "descendant_taxa",
            "true_value",
            "estimated_value",
            "absolute_error",
            "standard_error",
            "lower_95_interval",
            "upper_95_interval",
            "confidence",
        ],
        rows=[
            {
                "node": row.node,
                "descendant_taxa": ",".join(row.descendant_taxa),
                "true_value": _format_number(row.true_value),
                "estimated_value": _format_number(row.estimated_value),
                "absolute_error": _format_number(row.absolute_error),
                "standard_error": _format_number(row.standard_error),
                "lower_95_interval": _format_number(row.lower_95_interval),
                "upper_95_interval": _format_number(row.upper_95_interval),
                "confidence": _format_number(row.confidence),
            }
            for row in rows
        ],
    )


def _write_discrete_node_recovery_table(
    path: Path,
    rows: list[KnownAnswerDiscreteNodeRecoveryRow],
) -> Path:
    return write_taxon_rows(
        path,
        columns=[
            "node",
            "descendant_taxa",
            "true_state",
            "estimated_state",
            "true_state_probability",
            "confidence",
            "correct",
            "ambiguous",
        ],
        rows=[
            {
                "node": row.node,
                "descendant_taxa": ",".join(row.descendant_taxa),
                "true_state": row.true_state,
                "estimated_state": row.estimated_state,
                "true_state_probability": _format_number(
                    row.true_state_probability
                ),
                "confidence": _format_number(row.confidence),
                "correct": str(row.correct).lower(),
                "ambiguous": str(row.ambiguous).lower(),
            }
            for row in rows
        ],
    )


def _write_overview(
    path: Path,
    report: KnownAnswerReferenceWorkflowReport,
    bundle: KnownAnswerReferenceWorkflowBundle,
) -> Path:
    lines = [
        "# Known-Answer Simulation Demo",
        "",
        f"- dataset id: `{report.dataset.dataset_id}`",
        f"- taxon count: `{report.dataset.taxon_count}`",
        f"- alignment length: `{report.dataset.sequence_length}`",
        f"- distance recovery preserves rooted topology: `{str(bundle.rooted_topology_equal).lower()}`",
        f"- distance recovery preserves unrooted topology: `{str(bundle.same_unrooted_topology).lower()}`",
        f"- continuous internal-node mean absolute error: `{_format_number(bundle.continuous_internal_node_mean_absolute_error)}`",
        f"- discrete internal-node accuracy: `{_format_number(bundle.discrete_internal_node_accuracy)}`",
        "",
        "Generated outputs:",
        "",
        f"- workflow summary: `{bundle.workflow_summary_path.name}`",
        f"- recovered distance tree: `{bundle.distance_tree_path.name}`",
        f"- tree recovery ledger: `{bundle.tree_recovery_path.name}`",
        f"- parameter recovery ledger: `{bundle.parameter_recovery_path.name}`",
        f"- continuous node recovery ledger: `{bundle.continuous_node_recovery_path.name}`",
        f"- discrete node recovery ledger: `{bundle.discrete_node_recovery_path.name}`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _load_true_parameter_map(path: Path) -> dict[str, str]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {
            row["parameter"]: row["value"]
            for row in reader
        }


def _load_true_continuous_nodes(path: Path) -> list[KnownAnswerContinuousNodeTruth]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [
            KnownAnswerContinuousNodeTruth(
                node=row["node"],
                node_name=row["node_name"] or None,
                is_tip=row["is_tip"].strip().lower() == "true",
                descendant_taxa=_split_descendant_taxa(row["descendant_taxa"]),
                true_value=float(row["true_value"]),
            )
            for row in reader
        ]


def _load_true_discrete_nodes(path: Path) -> list[KnownAnswerDiscreteNodeTruth]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [
            KnownAnswerDiscreteNodeTruth(
                node=row["node"],
                node_name=row["node_name"] or None,
                is_tip=row["is_tip"].strip().lower() == "true",
                descendant_taxa=_split_descendant_taxa(row["descendant_taxa"]),
                true_state=row["true_state"],
            )
            for row in reader
        ]


def _split_descendant_taxa(value: str) -> list[str]:
    if not value:
        return []
    return value.split(",")


def _mean(values) -> float:
    materialized = list(values)
    if not materialized:
        return 0.0
    return sum(materialized) / len(materialized)


def _format_number(value: float) -> str:
    return format(value, ".15g")


def _resource_root() -> Path:
    return (
        Path(__file__).resolve().parent.parent
        / "resources"
        / "datasets"
        / "simulation"
        / _DATASET_ID
    )
