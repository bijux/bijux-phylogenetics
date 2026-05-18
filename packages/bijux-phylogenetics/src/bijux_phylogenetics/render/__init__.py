"""Rendering helpers for phylogenetics reports and tree figures."""

from .package import TreeFigurePackageResult, build_tree_figure_package
from .svg import AnnotationStrip, TreeRenderResult, render_tree_svg
from .time_tree_svg import (
    TimeTreeNodeInterval,
    TimeTreeRenderResult,
    render_time_tree_svg,
)
from .trait_tree_package import (
    AnnotatedTraitTreePackageResult,
    build_annotated_trait_tree_package,
)

__all__ = [
    "AnnotationStrip",
    "AnnotatedTraitTreePackageResult",
    "TimeTreeNodeInterval",
    "TimeTreeRenderResult",
    "TreeFigurePackageResult",
    "TreeRenderResult",
    "build_annotated_trait_tree_package",
    "build_tree_figure_package",
    "render_time_tree_svg",
    "render_tree_svg",
]
