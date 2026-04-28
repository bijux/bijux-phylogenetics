"""Phylogenetic analysis, reporting, and evidence tools for Bijux."""

from importlib import metadata

from .compare.topology import TreeComparisonReport, compare_tree_paths
from .diagnostics.validation import TreeInspectionReport, TreeValidationReport, inspect_tree_path, validate_tree_path
from .evidence.bundles import EvidenceBundleReport, bundle_directory
from .reports.service import ReportBuildResult, render_phylogenetics_report

try:
    __version__ = metadata.version("bijux-phylogenetics")
except metadata.PackageNotFoundError:
    __version__ = "0.1.0"

__all__ = [
    "EvidenceBundleReport",
    "ReportBuildResult",
    "TreeComparisonReport",
    "TreeInspectionReport",
    "TreeValidationReport",
    "__version__",
    "bundle_directory",
    "compare_tree_paths",
    "inspect_tree_path",
    "render_phylogenetics_report",
    "validate_tree_path",
]

