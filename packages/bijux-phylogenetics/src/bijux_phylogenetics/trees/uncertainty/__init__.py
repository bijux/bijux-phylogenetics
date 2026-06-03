from __future__ import annotations

from importlib import import_module

from .bootstrap import (
    summarize_bootstrap_tree_set as summarize_bootstrap_tree_set,
)
from .bootstrap import (
    write_bootstrap_tree_set_artifacts as write_bootstrap_tree_set_artifacts,
)
from .bootstrap import (
    write_bootstrap_tree_set_summary_table as write_bootstrap_tree_set_summary_table,
)
from .bootstrap import (
    write_bootstrap_unstable_branch_table as write_bootstrap_unstable_branch_table,
)
from .comparisons import (
    compare_bootstrap_and_posterior_uncertainty as compare_bootstrap_and_posterior_uncertainty,
)
from .comparisons import (
    compare_posterior_tree_sets as compare_posterior_tree_sets,
)
from .gene_tree_conflicts import (
    summarize_gene_tree_conflicts as summarize_gene_tree_conflicts,
)
from .gene_tree_conflicts import (
    write_gene_tree_conflict_artifacts as write_gene_tree_conflict_artifacts,
)
from .gene_tree_conflicts import (
    write_gene_tree_conflict_quartet_table as write_gene_tree_conflict_quartet_table,
)
from .gene_tree_conflicts import (
    write_gene_tree_conflict_summary_table as write_gene_tree_conflict_summary_table,
)
from .instability import (
    detect_unstable_clades as detect_unstable_clades,
)
from .instability import (
    detect_unstable_taxa as detect_unstable_taxa,
)
from .instability import (
    summarize_clade_credibility_conflicts as summarize_clade_credibility_conflicts,
)
from .instability import (
    summarize_uncertainty_aware_conclusions as summarize_uncertainty_aware_conclusions,
)
from .instability import (
    write_clade_credibility_conflict_table as write_clade_credibility_conflict_table,
)
from .instability import (
    write_uncertainty_conclusion_table as write_uncertainty_conclusion_table,
)
from .instability import (
    write_unstable_clade_table as write_unstable_clade_table,
)
from .methods_text import (
    TreeSetUncertaintyMethodReport as TreeSetUncertaintyMethodReport,
)
from .methods_text import (
    TreeSetUncertaintyMethodsSummaryTextResult as TreeSetUncertaintyMethodsSummaryTextResult,
)
from .methods_text import (
    build_tree_set_uncertainty_method_report as build_tree_set_uncertainty_method_report,
)
from .methods_text import (
    build_tree_set_uncertainty_methods_summary_text as build_tree_set_uncertainty_methods_summary_text,
)
from .methods_text import (
    write_tree_set_uncertainty_methods_summary_text as write_tree_set_uncertainty_methods_summary_text,
)
from .models import (
    BootstrapPosteriorCladeComparison as BootstrapPosteriorCladeComparison,
)
from .models import (
    BootstrapPosteriorSupportComparisonReport as BootstrapPosteriorSupportComparisonReport,
)
from .models import (
    BootstrapTreeSetArtifactReport as BootstrapTreeSetArtifactReport,
)
from .models import (
    BootstrapTreeSetSummaryReport as BootstrapTreeSetSummaryReport,
)
from .models import (
    BootstrapUnstableBranch as BootstrapUnstableBranch,
)
from .models import (
    CladeCredibilityConflict as CladeCredibilityConflict,
)
from .models import (
    CladeCredibilityConflictReport as CladeCredibilityConflictReport,
)
from .models import (
    CladeFrequencyDelta as CladeFrequencyDelta,
)
from .models import (
    ConsensusThresholdSensitivityReport as ConsensusThresholdSensitivityReport,
)
from .models import (
    ConsensusThresholdSensitivityRow as ConsensusThresholdSensitivityRow,
)
from .models import (
    GeneTreeConflictArtifactReport as GeneTreeConflictArtifactReport,
)
from .models import (
    GeneTreeConflictQuartetSummary as GeneTreeConflictQuartetSummary,
)
from .models import (
    GeneTreeConflictReferenceTree as GeneTreeConflictReferenceTree,
)
from .models import (
    GeneTreeConflictSummaryReport as GeneTreeConflictSummaryReport,
)
from .models import (
    PosteriorTopologicalDiversityComparisonReport as PosteriorTopologicalDiversityComparisonReport,
)
from .models import (
    PosteriorTopologicalDiversitySummary as PosteriorTopologicalDiversitySummary,
)
from .models import (
    PosteriorTopologyDiversityReport as PosteriorTopologyDiversityReport,
)
from .models import (
    PosteriorTopologyMode as PosteriorTopologyMode,
)
from .models import (
    PosteriorTopologyMultimodalityReport as PosteriorTopologyMultimodalityReport,
)
from .models import (
    PosteriorTreeSetComparisonReport as PosteriorTreeSetComparisonReport,
)
from .models import (
    RogueTaxonDetectionReport as RogueTaxonDetectionReport,
)
from .models import (
    RogueTaxonScoreRow as RogueTaxonScoreRow,
)
from .models import (
    TaxonPlacementSignature as TaxonPlacementSignature,
)
from .models import (
    TreeDistanceDistributionRow as TreeDistanceDistributionRow,
)
from .models import (
    TreeSetBenchmarkRow as TreeSetBenchmarkRow,
)
from .models import (
    TreeSetMaturityGateCheck as TreeSetMaturityGateCheck,
)
from .models import (
    TreeSetMaturityGateReport as TreeSetMaturityGateReport,
)
from .models import (
    TreeSetScalingBenchmarkReport as TreeSetScalingBenchmarkReport,
)
from .models import (
    TreeSetStorageRiskReport as TreeSetStorageRiskReport,
)
from .models import (
    TreeSetThinningSensitivityReport as TreeSetThinningSensitivityReport,
)
from .models import (
    TreeSetThinningSensitivityRow as TreeSetThinningSensitivityRow,
)
from .models import (
    TreeTopologyCluster as TreeTopologyCluster,
)
from .models import (
    TreeTopologyClusterReport as TreeTopologyClusterReport,
)
from .models import (
    UncertaintyAwareCladeConclusion as UncertaintyAwareCladeConclusion,
)
from .models import (
    UncertaintyAwareConclusionSummaryReport as UncertaintyAwareConclusionSummaryReport,
)
from .models import (
    UnstableClade as UnstableClade,
)
from .models import (
    UnstableCladeReport as UnstableCladeReport,
)
from .models import (
    UnstableTaxaReport as UnstableTaxaReport,
)
from .models import (
    UnstableTaxon as UnstableTaxon,
)
from .rogue_taxa import (
    detect_rogue_taxa as detect_rogue_taxa,
)
from .rogue_taxa import (
    write_rogue_taxon_table as write_rogue_taxon_table,
)
from .sensitivity import (
    assess_tree_set_maturity as assess_tree_set_maturity,
)
from .sensitivity import (
    assess_tree_set_storage_risk as assess_tree_set_storage_risk,
)
from .sensitivity import (
    assess_tree_set_thinning_sensitivity as assess_tree_set_thinning_sensitivity,
)
from .sensitivity import (
    benchmark_tree_set_uncertainty as benchmark_tree_set_uncertainty,
)
from .sensitivity import (
    compare_consensus_thresholds as compare_consensus_thresholds,
)
from .topology_diversity import (
    cluster_trees_by_topology as cluster_trees_by_topology,
)
from .topology_diversity import (
    compare_posterior_topological_diversity as compare_posterior_topological_diversity,
)
from .topology_diversity import (
    detect_posterior_topology_multimodality as detect_posterior_topology_multimodality,
)
from .topology_diversity import (
    summarize_posterior_topology_diversity as summarize_posterior_topology_diversity,
)
from .topology_diversity import (
    write_topology_cluster_table as write_topology_cluster_table,
)
from .topology_diversity import (
    write_tree_distance_distribution_table as write_tree_distance_distribution_table,
)

__all__ = [
    "BootstrapPosteriorCladeComparison",
    "BootstrapPosteriorSupportComparisonReport",
    "BootstrapTreeSetArtifactReport",
    "BootstrapTreeSetSummaryReport",
    "BootstrapUnstableBranch",
    "CladeCredibilityConflict",
    "CladeCredibilityConflictReport",
    "CladeFrequencyDelta",
    "ConsensusThresholdSensitivityReport",
    "ConsensusThresholdSensitivityRow",
    "GeneTreeConflictArtifactReport",
    "GeneTreeConflictQuartetSummary",
    "GeneTreeConflictReferenceTree",
    "GeneTreeConflictSummaryReport",
    "PosteriorTopologicalDiversityComparisonReport",
    "PosteriorTopologicalDiversitySummary",
    "PosteriorTopologyDiversityReport",
    "PosteriorTopologyMode",
    "PosteriorTopologyMultimodalityReport",
    "PosteriorTreeSetComparisonReport",
    "RogueTaxonDetectionReport",
    "RogueTaxonScoreRow",
    "TaxonPlacementSignature",
    "TreeDistanceDistributionRow",
    "TreeSetBenchmarkRow",
    "TreeSetMaturityGateCheck",
    "TreeSetMaturityGateReport",
    "TreeSetScalingBenchmarkReport",
    "TreeSetStorageRiskReport",
    "TreeSetThinningSensitivityReport",
    "TreeSetThinningSensitivityRow",
    "TreeSetUncertaintyCaptionDraft",
    "TreeSetUncertaintyFigurePackageResult",
    "TreeSetUncertaintyLegendEntry",
    "TreeSetUncertaintyMethodReport",
    "TreeSetUncertaintyMethodsSummaryTextResult",
    "TreeSetUncertaintyPublicationAudit",
    "TreeTopologyCluster",
    "TreeTopologyClusterReport",
    "UncertaintyAwareCladeConclusion",
    "UncertaintyAwareConclusionSummaryReport",
    "UnstableClade",
    "UnstableCladeReport",
    "UnstableTaxaReport",
    "UnstableTaxon",
    "assess_tree_set_maturity",
    "assess_tree_set_storage_risk",
    "assess_tree_set_thinning_sensitivity",
    "benchmark_tree_set_uncertainty",
    "build_tree_set_uncertainty_method_report",
    "build_tree_set_uncertainty_methods_summary_text",
    "build_tree_set_uncertainty_figure_package",
    "cluster_trees_by_topology",
    "compare_bootstrap_and_posterior_uncertainty",
    "compare_consensus_thresholds",
    "compare_posterior_topological_diversity",
    "compare_posterior_tree_sets",
    "detect_rogue_taxa",
    "detect_posterior_topology_multimodality",
    "detect_unstable_clades",
    "detect_unstable_taxa",
    "summarize_gene_tree_conflicts",
    "summarize_bootstrap_tree_set",
    "summarize_clade_credibility_conflicts",
    "summarize_posterior_topology_diversity",
    "summarize_uncertainty_aware_conclusions",
    "write_bootstrap_tree_set_artifacts",
    "write_bootstrap_tree_set_summary_table",
    "write_bootstrap_unstable_branch_table",
    "write_clade_credibility_conflict_table",
    "write_gene_tree_conflict_artifacts",
    "write_gene_tree_conflict_quartet_table",
    "write_gene_tree_conflict_summary_table",
    "write_rogue_taxon_table",
    "write_topology_cluster_table",
    "write_tree_distance_distribution_table",
    "write_tree_set_uncertainty_methods_summary_text",
    "write_uncertainty_conclusion_table",
    "write_unstable_clade_table",
]

_LAZY_FIGURE_PACKAGE_EXPORTS = {
    "TreeSetUncertaintyCaptionDraft",
    "TreeSetUncertaintyFigurePackageResult",
    "TreeSetUncertaintyLegendEntry",
    "TreeSetUncertaintyPublicationAudit",
    "build_tree_set_uncertainty_figure_package",
}


def __getattr__(name: str):
    if name in _LAZY_FIGURE_PACKAGE_EXPORTS:
        module = import_module(".figure_package", __name__)
        return getattr(module, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
