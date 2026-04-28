"""Phylogenetic analysis, reporting, and evidence tools for Bijux."""

from importlib import metadata

from .compare.topology import TreeComparisonReport, compare_tree_paths
from .core.topology import reroot_tree_by_midpoint, root_tree_on_outgroup
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
    translate_coding_alignment,
    trim_alignment,
)
from .reports.service import ReportBuildResult, render_phylo_inputs_report, render_phylogenetics_report

try:
    __version__ = metadata.version("bijux-phylogenetics")
except metadata.PackageNotFoundError:
    __version__ = "0.1.0"

__all__ = [
    "EvidenceBundleReport",
    "CLI_NAME",
    "IDENTITY",
    "IMPORT_NAME",
    "PACKAGE_NAME",
    "PRODUCT_NAME",
    "ReportBuildResult",
    "TreeDiagnosticReport",
    "TreeComparisonReport",
    "TreeInspectionReport",
    "TreeValidationReport",
    "UMBRELLA_COMMAND",
    "__version__",
    "build_alignment_quality_report",
    "bundle_directory",
    "compare_tree_paths",
    "compute_pairwise_sequence_identity_matrix",
    "diagnose_tree_path",
    "inspect_coding_alignment",
    "inspect_tree_path",
    "render_phylo_inputs_report",
    "render_phylogenetics_report",
    "reroot_tree_by_midpoint",
    "root_tree_on_outgroup",
    "summarise_fasta",
    "translate_coding_alignment",
    "trim_alignment",
    "validate_tree_path",
]
