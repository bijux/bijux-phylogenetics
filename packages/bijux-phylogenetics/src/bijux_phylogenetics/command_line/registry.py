from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class CommandSpec:
    """Stable metadata for a CLI command."""

    name: str
    domain: str
    summary: str
    inputs: tuple[str, ...]
    outputs: tuple[str, ...]


COMMAND_SPECS: tuple[CommandSpec, ...] = (
    CommandSpec(
        name="env",
        domain="runtime",
        summary="Inspect runtime dependency availability.",
        inputs=(),
        outputs=("environment-report",),
    ),
    CommandSpec(
        name="phylo",
        domain="engines",
        summary="Inspect external engine availability and workflow readiness.",
        inputs=("workflow-id",),
        outputs=("engine-preflight-report",),
    ),
    CommandSpec(
        name="metadata",
        domain="metadata-linkage",
        summary="Inspect metadata tables keyed by taxon.",
        inputs=("metadata-table",),
        outputs=("metadata-table-report",),
    ),
    CommandSpec(
        name="traits",
        domain="comparative-inputs",
        summary="Validate and link trait tables keyed by taxon.",
        inputs=("traits-table",),
        outputs=("traits-report",),
    ),
    CommandSpec(
        name="prune",
        domain="tree-io",
        summary="Prune a tree to taxa present in a metadata or traits table.",
        inputs=("tree", "taxon-table"),
        outputs=("pruned-tree", "pruned-taxa"),
    ),
    CommandSpec(
        name="alignment",
        domain="comparative-inputs",
        summary="Inspect and link FASTA alignments against tree taxa.",
        inputs=("alignment",),
        outputs=("alignment-report",),
    ),
    CommandSpec(
        name="comparative",
        domain="comparative-analysis",
        summary="Assess, summarize, and model numeric traits on a phylogeny.",
        inputs=("tree", "traits-table"),
        outputs=("comparative-report",),
    ),
    CommandSpec(
        name="ancestral",
        domain="ancestral-state",
        summary="Reconstruct, compare, render, and report ancestral trait states.",
        inputs=("tree", "traits-table"),
        outputs=("ancestral-state-report",),
    ),
    CommandSpec(
        name="biogeography",
        domain="biogeography",
        summary="Model ancestral geographic states and report region transition evidence.",
        inputs=("tree", "traits-table"),
        outputs=("biogeography-report",),
    ),
    CommandSpec(
        name="host-association",
        domain="host-association",
        summary="Reconstruct host-state evolution and report host-switch evidence on a phylogeny.",
        inputs=("tree", "traits-table"),
        outputs=("host-association-report",),
    ),
    CommandSpec(
        name="ecological-niche",
        domain="ecological-niche",
        summary="Model ecological niche evolution and report niche-transition evidence on a phylogeny.",
        inputs=("tree", "traits-table"),
        outputs=("ecological-niche-report",),
    ),
    CommandSpec(
        name="phylogeography",
        domain="geographic-reconstruction",
        summary="Reconstruct geographic coordinates or regional states and render map-based geographic review outputs on a phylogeny.",
        inputs=("tree", "traits-table"),
        outputs=("geographic-reconstruction-report",),
    ),
    CommandSpec(
        name="discrete-evolution",
        domain="discrete-state-evolution",
        summary="Model, compare, render, and report discrete-state evolution on a phylogeny.",
        inputs=("tree", "traits-table"),
        outputs=("discrete-state-evolution-report",),
    ),
    CommandSpec(
        name="diversification",
        domain="diversification-analysis",
        summary="Estimate diversification rates, compare models, and report macroevolution summaries on time trees.",
        inputs=("tree", "metadata-table", "traits-table"),
        outputs=("diversification-report",),
    ),
    CommandSpec(
        name="distance",
        domain="distance-analysis",
        summary="Validate, build, and report on explicit distance matrices.",
        inputs=("distance-matrix",),
        outputs=("distance-analysis-report",),
    ),
    CommandSpec(
        name="tree-set",
        domain="tree-uncertainty",
        summary="Summarize consensus, topology diversity, and instability across a tree set.",
        inputs=("tree-set",),
        outputs=("tree-set-report",),
    ),
    CommandSpec(
        name="simulate",
        domain="simulation",
        summary="Simulate trees, traits, and alignments under explicit stochastic models.",
        inputs=("simulation-input",),
        outputs=("simulation-output",),
    ),
    CommandSpec(
        name="benchmark",
        domain="scientific-validation",
        summary="Benchmark validation, comparison, rendering, reporting, diagnostics, trimming, distance, and readiness across size scales.",
        inputs=(),
        outputs=("benchmark-report",),
    ),
    CommandSpec(
        name="parity",
        domain="reference-validation",
        summary="Validate core numerical phylogenetics methods against checked-in external reference outputs.",
        inputs=(),
        outputs=("reference-parity-report",),
    ),
    CommandSpec(
        name="inspect",
        domain="tree-io",
        summary="Inspect a tree and report high-level summary metrics.",
        inputs=("tree",),
        outputs=("tree-inspection-report",),
    ),
    CommandSpec(
        name="validate",
        domain="tree-diagnostics",
        summary="Validate tree structure and branch-length hygiene.",
        inputs=("tree",),
        outputs=("tree-validation-report",),
    ),
    CommandSpec(
        name="normalize",
        domain="tree-io",
        summary="Normalize a tree into canonical output formats.",
        inputs=("tree",),
        outputs=("normalized-tree",),
    ),
    CommandSpec(
        name="normalize-taxa",
        domain="taxa",
        summary="Apply an explicit taxon normalization policy and emit a mapping file.",
        inputs=("tree",),
        outputs=("normalized-tree", "taxon-mapping"),
    ),
    CommandSpec(
        name="taxonomy",
        domain="taxa",
        summary="Audit synonym, namespace, loss, and stability behavior for phylogenetic taxa.",
        inputs=("tree", "taxon-table"),
        outputs=("taxon-audit-report",),
    ),
    CommandSpec(
        name="topology",
        domain="tree-io",
        summary="Apply explicit rooting and ordering transforms to a tree.",
        inputs=("tree",),
        outputs=("transformed-tree",),
    ),
    CommandSpec(
        name="compare",
        domain="tree-comparison",
        summary="Compare two trees over their shared taxa.",
        inputs=("left-tree", "right-tree"),
        outputs=("tree-comparison-report",),
    ),
    CommandSpec(
        name="annotate",
        domain="metadata-linkage",
        summary="Check trait or metadata linkage against tree tips.",
        inputs=("tree", "metadata-table"),
        outputs=("metadata-linkage-report",),
    ),
    CommandSpec(
        name="diagnose",
        domain="tree-diagnostics",
        summary="Run combined inspection and validation diagnostics.",
        inputs=("tree",),
        outputs=("tree-diagnostic-report",),
    ),
    CommandSpec(
        name="render",
        domain="reporting",
        summary="Render an HTML report for a tree and optional metadata.",
        inputs=("tree", "metadata-table"),
        outputs=("html-report",),
    ),
    CommandSpec(
        name="report",
        domain="reporting",
        summary="Build an evidence-first phylogenetics HTML report.",
        inputs=("tree", "alignment", "traits-table", "metadata-table"),
        outputs=("html-report", "provenance-manifest"),
    ),
    CommandSpec(
        name="demo",
        domain="examples",
        summary="Run the repository capability demo workflow.",
        inputs=(),
        outputs=("demo-artifacts",),
    ),
    CommandSpec(
        name="evidence",
        domain="evidence",
        summary="Bundle evidence packs and govern the repository evidence-book.",
        inputs=("run-root", "study-id", "evidence-id"),
        outputs=("evidence-bundle", "evidence-book-index", "reviewer-summary"),
    ),
    CommandSpec(
        name="adapter",
        domain="engines",
        summary="Inspect or invoke configured phylogenetics adapters.",
        inputs=("adapter-name",),
        outputs=("adapter-report",),
    ),
)


def get_command_spec(name: str) -> CommandSpec:
    """Return the registered specification for a command."""
    for spec in COMMAND_SPECS:
        if spec.name == name:
            return spec
    raise KeyError(name)
