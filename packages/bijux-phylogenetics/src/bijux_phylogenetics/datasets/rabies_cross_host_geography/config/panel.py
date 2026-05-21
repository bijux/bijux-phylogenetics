# ruff: noqa: F401, F403, F405
from __future__ import annotations

from dataclasses import dataclass, replace
from hashlib import sha256
from html import escape
from pathlib import Path
from tempfile import TemporaryDirectory

from Bio import Phylo

from bijux_phylogenetics.ancestral.discrete import (
    DiscreteAncestralReport,
    reconstruct_discrete_ancestral_states,
)
from bijux_phylogenetics.ancestral.tree_set import (
    DiscreteAncestralTreeSetReport,
    summarize_discrete_ancestral_tree_set,
)
from bijux_phylogenetics.biogeography import (
    BiogeographyReportPackageResult,
    build_biogeography_report_package,
)
from bijux_phylogenetics.trees import (
    CladeTableReport,
    CladeTableRow,
    extract_tree_clades,
    write_clade_table,
)
from bijux_phylogenetics.comparative.pgls import (
    PGLSResult,
    run_pgls,
    write_pgls_model_matrix_table,
)
from bijux_phylogenetics.comparative.pgls.categorical_contrasts import (
    PGLSCategoricalContrastReport,
    summarize_pgls_categorical_contrasts,
    write_pgls_categorical_contrast_table,
)
from bijux_phylogenetics.comparative.pgls.lambda_fit import (
    write_pgls_lambda_profile_table,
)
from bijux_phylogenetics.comparative.pgls.posterior_tree import (
    PosteriorTreePGLSReport,
    run_posterior_tree_pgls,
)
from bijux_phylogenetics.comparative.reporting.analysis_package import (
    ComparativeAnalysisSummaryRow,
    ComparativeCoefficientTableRow,
    ComparativeInterpretationRow,
    ComparativeResidualTableRow,
    ComparativeSignalTableRow,
    summarize_comparative_analysis,
    summarize_comparative_audit,
    summarize_comparative_coefficients,
    summarize_comparative_interpretation,
    summarize_comparative_residuals,
    summarize_comparative_signal,
    write_comparative_audit_table,
    write_comparative_coefficient_table,
    write_comparative_contrast_table,
    write_comparative_interpretation_table,
    write_comparative_model_comparison_table,
    write_comparative_residual_table,
    write_comparative_signal_table,
    write_comparative_summary_table,
)
from bijux_phylogenetics.comparative.reporting import (
    ComparativeMethodReport,
    build_comparative_method_report,
)
from bijux_phylogenetics.compare.presentation import (
    ComparisonReportBuildResult,
    build_tree_comparison_report,
)
from bijux_phylogenetics.compare.topology import write_tree_comparison_table
from bijux_phylogenetics.phylo.alignment import (
    AlignmentQualityReport,
    SequenceQualityRankingReport,
)
from bijux_phylogenetics.phylo.topology import (
    TreeRootingReport,
    root_tree_on_outgroup,
    write_tree_rooting_report,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.datasets.rabies_method_sensitivity import (
    RabiesMethodSensitivityPanelWorkflowReport,
    run_rabies_method_sensitivity_panel_workflow,
)
from bijux_phylogenetics.diagnostics.conclusion_stability import (
    ConclusionStabilityReport,
    build_ancestral_state_stability_rows,
    build_comparative_coefficient_stability_rows,
    build_conclusion_stability_report,
    build_key_clade_stability_rows,
    build_support_value_stability_rows,
    write_ancestral_state_stability_table,
    write_comparative_coefficient_stability_table,
    write_conclusion_stability_report_html,
    write_conclusion_stability_summary_table,
    write_key_clade_stability_table,
    write_support_value_stability_table,
)
from bijux_phylogenetics.engines.inference import (
    FastaToTreeWorkflowReport,
    run_fasta_to_tree_workflow,
)
from bijux_phylogenetics.ecology import (
    HostSwitchingReport,
    summarize_host_switching,
    write_host_state_node_table,
    write_host_switch_branch_table,
    write_host_switch_count_table,
    write_host_switch_exclusion_table,
    write_host_switch_fit_table,
    write_host_switch_summary_table,
    write_unsupported_host_switch_claim_table,
)
from bijux_phylogenetics.io.biopython import tree_from_biophylo, tree_to_biophylo
from bijux_phylogenetics.io.fasta.quality import (
    build_alignment_quality_report,
    build_sequence_quality_ranking,
)
from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.simulation import write_tree_set
from bijux_phylogenetics.trees import (
    BootstrapTreeSetArtifactReport,
    BootstrapTreeSetSummaryReport,
    compute_clade_frequency_table,
    write_bootstrap_tree_set_artifacts,
)

