"""Rendering helpers for phylogenetics reports and tree figures."""

from .package import TreeFigurePackageResult, build_tree_figure_package
from .svg import AnnotationStrip, TreeRenderResult, render_tree_svg
from .trait_tree_package import (
    AnnotatedTraitTreePackageResult,
    build_annotated_trait_tree_package,
)

__all__ = [
    "AnnotationStrip",
    "AnnotatedTraitTreePackageResult",
    "TreeFigurePackageResult",
    "TreeRenderResult",
    "build_annotated_trait_tree_package",
    "build_tree_figure_package",
    "render_tree_svg",
]
