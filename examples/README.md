# Examples

This directory holds runnable example inputs and workflows for
`bijux-phylogenetics`.

Use repository-owned example material here when documenting CLI behavior or
teaching a workflow. Keep large generated outputs under `artifacts/` unless the
task explicitly governs a tracked destination in this tree.

Run the public capability walkthrough with:

```bash
make demo
```

Distance iteration example:

```bash
uv run bijux-phylogenetics alignment distance-quality path/to/alignment.fasta --model kimura-2-parameter --json
uv run bijux-phylogenetics alignment distance-suitability path/to/alignment.fasta --model kimura-2-parameter --json
uv run bijux-phylogenetics alignment bootstrap-tree path/to/alignment.fasta --method neighbor-joining --replicates 200 --support-out artifacts/distance-support.tsv --tree-set-out artifacts/distance-bootstrap.trees --json
uv run bijux-phylogenetics alignment distance-support-summary path/to/alignment.fasta --method neighbor-joining --replicates 50 --json
uv run bijux-phylogenetics alignment distance-models path/to/alignment.fasta --json
uv run bijux-phylogenetics alignment distance-gap-sensitivity path/to/alignment.fasta --model p-distance --json
uv run bijux-phylogenetics alignment distance-method-report path/to/alignment.fasta --method neighbor-joining --replicates 50 --json
uv run bijux-phylogenetics alignment distance-maturity path/to/alignment.fasta --method neighbor-joining --replicates 50 --json
uv run bijux-phylogenetics alignment distance-bundle path/to/alignment.fasta --method neighbor-joining --replicates 200 --out-dir artifacts/distance-bundle --json
uv run bijux-phylogenetics distance quality path/to/exported-distances.tsv --json
```

Comparative iteration example:

```bash
uv run bijux-phylogenetics comparative validate-reference --json
uv run bijux-phylogenetics comparative maturity path/to/tree.nwk path/to/traits.tsv --formula "height_cm ~ body_mass + habitat" --lambda-value 1.0 --json
uv run bijux-phylogenetics comparative report path/to/tree.nwk path/to/traits.tsv --formula "height_cm ~ body_mass + habitat" --out artifacts/comparative-report.html --json
```
