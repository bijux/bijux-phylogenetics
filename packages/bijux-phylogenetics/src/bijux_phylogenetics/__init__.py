"""Phylogenetic analysis, reporting, and evidence tools for Bijux."""

from importlib import metadata

from .compare.topology import TreeComparisonReport, compare_tree_paths
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
from .reports.service import ReportBuildResult, render_phylogenetics_report

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
    "bundle_directory",
    "compare_tree_paths",
    "diagnose_tree_path",
    "inspect_tree_path",
    "render_phylogenetics_report",
    "validate_tree_path",
]
