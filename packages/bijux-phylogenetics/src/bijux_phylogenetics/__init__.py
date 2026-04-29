"""Phylogenetic analysis, reporting, and evidence tools for Bijux."""

from importlib import metadata

from .compare.topology import TreeComparisonReport, compare_tree_paths
from .core.topology import reroot_tree_by_midpoint, root_tree_on_outgroup, unroot_tree
from .diagnostics.assumptions import (
    BranchLengthUnitReport,
    StandardizedSupportLabel,
    TreeAssumptionReport,
    assess_tree_assumptions,
    inspect_branch_length_units,
    standardize_support_labels,
)
from .distance import (
    DistanceMatrixAsymmetry,
    ImportedDistanceMatrixReport,
    ImportedDistanceTreeBuildReport,
    DistanceTreeBuildReport,
    DistanceTreeTopologyComparison,
    GeneticDistanceMatrix,
    build_distance_tree,
    build_tree_from_imported_distance_matrix,
    compare_distance_tree_topologies,
    compute_pairwise_genetic_distance_matrix,
    load_imported_distance_matrix,
    validate_imported_distance_matrix,
    write_genetic_distance_matrix,
)
from .diagnostics.validation import (
    TreeDiagnosticReport,
    TreeInspectionReport,
    TreeValidationReport,
    diagnose_tree_path,
    inspect_tree_path,
    validate_tree_path,
)
from .evidence.bundles import EvidenceBundleReport, bundle_directory
from .identity import CLI_NAME, IDENTITY, IMPORT_NAME, PACKAGE_NAME, PRODUCT_NAME, UMBRELLA_COMMAND
from .io.fasta import (
    build_alignment_quality_report,
    compute_pairwise_sequence_identity_matrix,
    inspect_coding_alignment,
    summarise_fasta,
    trim_columns_above_missingness_threshold,
    translate_coding_alignment,
    trim_alignment,
)
from .render import AnnotationStrip, TreeFigurePackageResult, TreeRenderResult, build_tree_figure_package, render_tree_svg
from .reports.service import ReportBuildResult, render_phylo_inputs_report, render_phylogenetics_report
from .reports.service import DistanceReportBuildResult, distance_method_limitations, render_distance_report

try:
    __version__ = metadata.version("bijux-phylogenetics")
except metadata.PackageNotFoundError:
    __version__ = "0.1.0"

__all__ = [
    "EvidenceBundleReport",
    "CLI_NAME",
    "BranchLengthUnitReport",
    "DistanceMatrixAsymmetry",
    "DistanceReportBuildResult",
    "DistanceTreeBuildReport",
    "DistanceTreeTopologyComparison",
    "IDENTITY",
    "ImportedDistanceMatrixReport",
    "ImportedDistanceTreeBuildReport",
    "IMPORT_NAME",
    "GeneticDistanceMatrix",
    "PACKAGE_NAME",
    "PRODUCT_NAME",
    "ReportBuildResult",
    "AnnotationStrip",
    "StandardizedSupportLabel",
    "TreeAssumptionReport",
    "TreeDiagnosticReport",
    "TreeFigurePackageResult",
    "TreeComparisonReport",
    "TreeInspectionReport",
    "TreeRenderResult",
    "TreeValidationReport",
    "UMBRELLA_COMMAND",
    "__version__",
    "build_alignment_quality_report",
    "build_tree_figure_package",
    "build_tree_from_imported_distance_matrix",
    "bundle_directory",
    "compare_tree_paths",
    "compare_distance_tree_topologies",
    "compute_pairwise_sequence_identity_matrix",
    "compute_pairwise_genetic_distance_matrix",
    "distance_method_limitations",
    "diagnose_tree_path",
    "inspect_branch_length_units",
    "inspect_coding_alignment",
    "inspect_tree_path",
    "load_imported_distance_matrix",
    "render_distance_report",
    "render_phylo_inputs_report",
    "render_phylogenetics_report",
    "render_tree_svg",
    "reroot_tree_by_midpoint",
    "root_tree_on_outgroup",
    "assess_tree_assumptions",
    "standardize_support_labels",
    "summarise_fasta",
    "build_distance_tree",
    "trim_columns_above_missingness_threshold",
    "translate_coding_alignment",
    "trim_alignment",
    "unroot_tree",
    "validate_imported_distance_matrix",
    "validate_tree_path",
    "write_genetic_distance_matrix",
]
