from __future__ import annotations

MAXIMUM_LIKELIHOOD_WRAPPER_REFERENCE_PAYLOADS: dict[str, dict[str, object]] = {
    "iqtree-jc69-maximum-likelihood": {
        "engine_name": "iqtree",
        "reference_surface": "maximum-likelihood-tree",
        "alignment_path": "packages/bijux-phylogenetics/tests/fixtures/alignments/jc69_likelihood_nni_alignment_4_taxa.fasta",
        "model_name": "JC",
        "tree_newick": "(A:0.0000010000,B:0.0000010000,(C:0.0000010000,D:0.0000010000):9.9999986410);",
        "log_likelihood": -33.271,
        "runtime_seconds": 0.018,
        "peak_memory_bytes": 8110080,
        "memory_observation_kind": "macos-time-maximum-resident-set-size",
        "support_rows": [],
        "notes": [
            "cached from a local iqtree2 2.3.6 fixed-model JC run on the governed four-taxon likelihood fixture",
            "the cached wall-clock runtime uses the engine-reported total wall-clock time and the cached peak memory uses macos /usr/bin/time maximum resident set size",
        ],
    },
    "iqtree-jc69-ultrafast-bootstrap-support": {
        "engine_name": "iqtree",
        "reference_surface": "ultrafast-bootstrap-support",
        "alignment_path": "packages/bijux-phylogenetics/tests/fixtures/alignments/jc69_likelihood_nni_alignment_4_taxa.fasta",
        "model_name": "JC",
        "tree_newick": "(A:0.0000023479,B:0.0000023479,(C:0.0000023479,D:0.0000023479)100:10.0000000000);",
        "log_likelihood": -33.271,
        "runtime_seconds": 0.129,
        "peak_memory_bytes": 10502144,
        "memory_observation_kind": "macos-time-maximum-resident-set-size",
        "support_rows": [
            {
                "split_taxa": ["A", "B"],
                "support_percent": 100.0,
            }
        ],
        "notes": [
            "cached from a local iqtree2 2.3.6 JC ultrafast-bootstrap run with 1000 replicates on the governed four-taxon likelihood fixture",
            "the cached support split uses one canonical unrooted bipartition so the same internal edge can be matched against rooted native bootstrap clades",
        ],
    },
    "fasttree-jc69-approximate-maximum-likelihood": {
        "engine_name": "fasttree",
        "reference_surface": "approximate-maximum-likelihood-tree",
        "alignment_path": "packages/bijux-phylogenetics/tests/fixtures/alignments/jc69_likelihood_nni_alignment_4_taxa.fasta",
        "model_name": "Jukes-Cantor+CAT",
        "tree_newick": "((A:0.0,B:0.0):2.995048178,(C:0.0,D:0.0):2.995048178);",
        "log_likelihood": -33.275,
        "runtime_seconds": 0.006905541056767,
        "peak_memory_bytes": 2260992,
        "memory_observation_kind": "macos-time-maximum-resident-set-size",
        "support_rows": [],
        "notes": [
            "cached from a local FastTree 2.1.11 nucleotide run on the governed four-taxon likelihood fixture",
            "the cached wall-clock runtime uses one Python perf_counter measurement because FastTree rounds its own reported total time to 0.00 seconds on this small fixture",
            "FastTree emits one approximate CAT likelihood surface, so likelihood and support comparisons must remain expectation-guided rather than exact parity claims",
        ],
    },
}
