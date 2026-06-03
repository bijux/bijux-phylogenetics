from __future__ import annotations

from pathlib import Path

from tests.support.fake_external_engines import write_executable


def fake_ape_rscript(
    path: Path,
    *,
    ape_available: bool = True,
    summary_overrides: dict[str, object] | None = None,
    normalized_tree_overrides: dict[str, str] | None = None,
) -> Path:
    summary_payload = repr(summary_overrides or {})
    normalized_tree_payload = repr(normalized_tree_overrides or {})
    return write_executable(
        path,
        f"""#!/usr/bin/env python3
import csv
import json
import math
import sys
from pathlib import Path

from Bio import Phylo

TABULAR_CASES = {{
    "dna-dnabin-structure-clean": {{
        "summary": {{
            "sequence_count": 4,
            "alignment_length": 8,
            "state_count": 32,
        }},
        "rows_name": "dnabin.tsv",
        "rows": [
            {{"identifier": "A", "position": 1, "state": "a"}},
            {{"identifier": "A", "position": 2, "state": "a"}},
            {{"identifier": "A", "position": 3, "state": "a"}},
            {{"identifier": "A", "position": 4, "state": "a"}},
            {{"identifier": "A", "position": 5, "state": "a"}},
            {{"identifier": "A", "position": 6, "state": "a"}},
            {{"identifier": "A", "position": 7, "state": "c"}},
            {{"identifier": "A", "position": 8, "state": "c"}},
            {{"identifier": "B", "position": 1, "state": "a"}},
            {{"identifier": "B", "position": 2, "state": "a"}},
            {{"identifier": "B", "position": 3, "state": "a"}},
            {{"identifier": "B", "position": 4, "state": "a"}},
            {{"identifier": "B", "position": 5, "state": "a"}},
            {{"identifier": "B", "position": 6, "state": "a"}},
            {{"identifier": "B", "position": 7, "state": "c"}},
            {{"identifier": "B", "position": 8, "state": "t"}},
            {{"identifier": "C", "position": 1, "state": "t"}},
            {{"identifier": "C", "position": 2, "state": "t"}},
            {{"identifier": "C", "position": 3, "state": "t"}},
            {{"identifier": "C", "position": 4, "state": "t"}},
            {{"identifier": "C", "position": 5, "state": "a"}},
            {{"identifier": "C", "position": 6, "state": "a"}},
            {{"identifier": "C", "position": 7, "state": "c"}},
            {{"identifier": "C", "position": 8, "state": "c"}},
            {{"identifier": "D", "position": 1, "state": "t"}},
            {{"identifier": "D", "position": 2, "state": "t"}},
            {{"identifier": "D", "position": 3, "state": "t"}},
            {{"identifier": "D", "position": 4, "state": "t"}},
            {{"identifier": "D", "position": 5, "state": "a"}},
            {{"identifier": "D", "position": 6, "state": "a"}},
            {{"identifier": "D", "position": 7, "state": "c"}},
            {{"identifier": "D", "position": 8, "state": "t"}},
        ],
    }},
    "dna-dnabin-structure-lowercase": {{
        "summary": {{
            "sequence_count": 3,
            "alignment_length": 6,
            "state_count": 18,
        }},
        "rows_name": "dnabin.tsv",
        "rows": [
            {{"identifier": "lower_a", "position": 1, "state": "a"}},
            {{"identifier": "lower_a", "position": 2, "state": "c"}},
            {{"identifier": "lower_a", "position": 3, "state": "g"}},
            {{"identifier": "lower_a", "position": 4, "state": "t"}},
            {{"identifier": "lower_a", "position": 5, "state": "a"}},
            {{"identifier": "lower_a", "position": 6, "state": "a"}},
            {{"identifier": "lower_b", "position": 1, "state": "a"}},
            {{"identifier": "lower_b", "position": 2, "state": "c"}},
            {{"identifier": "lower_b", "position": 3, "state": "g"}},
            {{"identifier": "lower_b", "position": 4, "state": "t"}},
            {{"identifier": "lower_b", "position": 5, "state": "t"}},
            {{"identifier": "lower_b", "position": 6, "state": "a"}},
            {{"identifier": "lower_c", "position": 1, "state": "a"}},
            {{"identifier": "lower_c", "position": 2, "state": "t"}},
            {{"identifier": "lower_c", "position": 3, "state": "g"}},
            {{"identifier": "lower_c", "position": 4, "state": "t"}},
            {{"identifier": "lower_c", "position": 5, "state": "n"}},
            {{"identifier": "lower_c", "position": 6, "state": "n"}},
        ],
    }},
    "dna-dnabin-structure-gaps": {{
        "summary": {{
            "sequence_count": 4,
            "alignment_length": 6,
            "state_count": 24,
        }},
        "rows_name": "dnabin.tsv",
        "rows": [
            {{"identifier": "A", "position": 1, "state": "a"}},
            {{"identifier": "A", "position": 2, "state": "c"}},
            {{"identifier": "A", "position": 3, "state": "g"}},
            {{"identifier": "A", "position": 4, "state": "t"}},
            {{"identifier": "A", "position": 5, "state": "a"}},
            {{"identifier": "A", "position": 6, "state": "a"}},
            {{"identifier": "B", "position": 1, "state": "a"}},
            {{"identifier": "B", "position": 2, "state": "c"}},
            {{"identifier": "B", "position": 3, "state": "g"}},
            {{"identifier": "B", "position": 4, "state": "t"}},
            {{"identifier": "B", "position": 5, "state": "t"}},
            {{"identifier": "B", "position": 6, "state": "a"}},
            {{"identifier": "C", "position": 1, "state": "a"}},
            {{"identifier": "C", "position": 2, "state": "t"}},
            {{"identifier": "C", "position": 3, "state": "g"}},
            {{"identifier": "C", "position": 4, "state": "t"}},
            {{"identifier": "C", "position": 5, "state": "c"}},
            {{"identifier": "C", "position": 6, "state": "g"}},
            {{"identifier": "D", "position": 1, "state": "a"}},
            {{"identifier": "D", "position": 2, "state": "-"}},
            {{"identifier": "D", "position": 3, "state": "g"}},
            {{"identifier": "D", "position": 4, "state": "t"}},
            {{"identifier": "D", "position": 5, "state": "n"}},
            {{"identifier": "D", "position": 6, "state": "a"}},
        ],
    }},
    "dna-dnabin-structure-ambiguity": {{
        "summary": {{
            "sequence_count": 3,
            "alignment_length": 6,
            "state_count": 18,
        }},
        "rows_name": "dnabin.tsv",
        "rows": [
            {{"identifier": "A", "position": 1, "state": "a"}},
            {{"identifier": "A", "position": 2, "state": "c"}},
            {{"identifier": "A", "position": 3, "state": "g"}},
            {{"identifier": "A", "position": 4, "state": "t"}},
            {{"identifier": "A", "position": 5, "state": "n"}},
            {{"identifier": "A", "position": 6, "state": "?"}},
            {{"identifier": "B", "position": 1, "state": "a"}},
            {{"identifier": "B", "position": 2, "state": "c"}},
            {{"identifier": "B", "position": 3, "state": "g"}},
            {{"identifier": "B", "position": 4, "state": "t"}},
            {{"identifier": "B", "position": 5, "state": "r"}},
            {{"identifier": "B", "position": 6, "state": "?"}},
            {{"identifier": "C", "position": 1, "state": "a"}},
            {{"identifier": "C", "position": 2, "state": "c"}},
            {{"identifier": "C", "position": 3, "state": "g"}},
            {{"identifier": "C", "position": 4, "state": "t"}},
            {{"identifier": "C", "position": 5, "state": "-"}},
            {{"identifier": "C", "position": 6, "state": "?"}},
        ],
    }},
    "dna-base-frequency-lowercase": {{
        "summary": {{
            "sequence_count": 3,
            "alignment_length": 6,
            "state_count": 17,
        }},
        "rows_name": "base-frequency.tsv",
        "rows": [
            {{"state": "a", "frequency": 6 / 18}},
            {{"state": "c", "frequency": 2 / 18}},
            {{"state": "g", "frequency": 3 / 18}},
            {{"state": "t", "frequency": 5 / 18}},
            {{"state": "r", "frequency": 0.0}},
            {{"state": "m", "frequency": 0.0}},
            {{"state": "w", "frequency": 0.0}},
            {{"state": "s", "frequency": 0.0}},
            {{"state": "k", "frequency": 0.0}},
            {{"state": "y", "frequency": 0.0}},
            {{"state": "v", "frequency": 0.0}},
            {{"state": "h", "frequency": 0.0}},
            {{"state": "d", "frequency": 0.0}},
            {{"state": "b", "frequency": 0.0}},
            {{"state": "n", "frequency": 2 / 18}},
            {{"state": "-", "frequency": 0.0}},
            {{"state": "?", "frequency": 0.0}},
        ],
    }},
    "dna-base-frequency-ambiguity": {{
        "summary": {{
            "sequence_count": 3,
            "alignment_length": 6,
            "state_count": 17,
        }},
        "rows_name": "base-frequency.tsv",
        "rows": [
            {{"state": "a", "frequency": 3 / 18}},
            {{"state": "c", "frequency": 3 / 18}},
            {{"state": "g", "frequency": 3 / 18}},
            {{"state": "t", "frequency": 3 / 18}},
            {{"state": "r", "frequency": 1 / 18}},
            {{"state": "m", "frequency": 0.0}},
            {{"state": "w", "frequency": 0.0}},
            {{"state": "s", "frequency": 0.0}},
            {{"state": "k", "frequency": 0.0}},
            {{"state": "y", "frequency": 0.0}},
            {{"state": "v", "frequency": 0.0}},
            {{"state": "h", "frequency": 0.0}},
            {{"state": "d", "frequency": 0.0}},
            {{"state": "b", "frequency": 0.0}},
            {{"state": "n", "frequency": 1 / 18}},
            {{"state": "-", "frequency": 1 / 18}},
            {{"state": "?", "frequency": 3 / 18}},
        ],
    }},
    "dna-base-frequency-missing-data": {{
        "summary": {{
            "sequence_count": 3,
            "alignment_length": 6,
            "state_count": 17,
        }},
        "rows_name": "base-frequency.tsv",
        "rows": [
            {{"state": "a", "frequency": 5 / 18}},
            {{"state": "c", "frequency": 3 / 18}},
            {{"state": "g", "frequency": 3 / 18}},
            {{"state": "t", "frequency": 3 / 18}},
            {{"state": "r", "frequency": 0.0}},
            {{"state": "m", "frequency": 0.0}},
            {{"state": "w", "frequency": 0.0}},
            {{"state": "s", "frequency": 0.0}},
            {{"state": "k", "frequency": 0.0}},
            {{"state": "y", "frequency": 0.0}},
            {{"state": "v", "frequency": 0.0}},
            {{"state": "h", "frequency": 0.0}},
            {{"state": "d", "frequency": 0.0}},
            {{"state": "b", "frequency": 0.0}},
            {{"state": "n", "frequency": 3 / 18}},
            {{"state": "-", "frequency": 0.0}},
            {{"state": "?", "frequency": 1 / 18}},
        ],
    }},
    "dna-base-frequency-all-gap-missing": {{
        "summary": {{
            "sequence_count": 3,
            "alignment_length": 6,
            "state_count": 17,
        }},
        "rows_name": "base-frequency.tsv",
        "rows": [
            {{"state": "a", "frequency": 0.0}},
            {{"state": "c", "frequency": 0.0}},
            {{"state": "g", "frequency": 0.0}},
            {{"state": "t", "frequency": 0.0}},
            {{"state": "r", "frequency": 0.0}},
            {{"state": "m", "frequency": 0.0}},
            {{"state": "w", "frequency": 0.0}},
            {{"state": "s", "frequency": 0.0}},
            {{"state": "k", "frequency": 0.0}},
            {{"state": "y", "frequency": 0.0}},
            {{"state": "v", "frequency": 0.0}},
            {{"state": "h", "frequency": 0.0}},
            {{"state": "d", "frequency": 0.0}},
            {{"state": "b", "frequency": 0.0}},
            {{"state": "n", "frequency": 0.0}},
            {{"state": "-", "frequency": 8 / 18}},
            {{"state": "?", "frequency": 10 / 18}},
        ],
    }},
    "dna-segregating-sites-lowercase": {{
        "summary": {{
            "sequence_count": 3,
            "alignment_length": 6,
            "segregating_site_count": 2,
        }},
        "rows_name": "segregating-sites.tsv",
        "rows": [
            {{"position": 2}},
            {{"position": 5}},
        ],
    }},
    "dna-segregating-sites-invariant": {{
        "summary": {{
            "sequence_count": 3,
            "alignment_length": 6,
            "segregating_site_count": 0,
        }},
        "rows_name": "segregating-sites.tsv",
        "rows": [],
    }},
    "dna-segregating-sites-one-variable": {{
        "summary": {{
            "sequence_count": 3,
            "alignment_length": 6,
            "segregating_site_count": 1,
        }},
        "rows_name": "segregating-sites.tsv",
        "rows": [
            {{"position": 4}},
        ],
    }},
    "dna-segregating-sites-gaps": {{
        "summary": {{
            "sequence_count": 4,
            "alignment_length": 6,
            "segregating_site_count": 3,
        }},
        "rows_name": "segregating-sites.tsv",
        "rows": [
            {{"position": 2}},
            {{"position": 5}},
            {{"position": 6}},
        ],
    }},
    "dna-segregating-sites-ambiguity": {{
        "summary": {{
            "sequence_count": 3,
            "alignment_length": 6,
            "segregating_site_count": 1,
        }},
        "rows_name": "segregating-sites.tsv",
        "rows": [
            {{"position": 5}},
        ],
    }},
    "dna-segregating-sites-missing-data": {{
        "summary": {{
            "sequence_count": 3,
            "alignment_length": 6,
            "segregating_site_count": 0,
        }},
        "rows_name": "segregating-sites.tsv",
        "rows": [],
    }},
    "dna-segregating-sites-all-gap-missing": {{
        "summary": {{
            "sequence_count": 3,
            "alignment_length": 6,
            "segregating_site_count": 0,
        }},
        "rows_name": "segregating-sites.tsv",
        "rows": [],
    }},
    "dna-raw-distance-clean": {{
        "summary": {{
            "sequence_count": 4,
            "alignment_length": 8,
            "pairwise_deletion": False,
        }},
        "rows_name": "distance-matrix.tsv",
        "rows": [
            {{"left_identifier": "A", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "B", "distance": 0.125}},
            {{"left_identifier": "A", "right_identifier": "C", "distance": 0.5}},
            {{"left_identifier": "A", "right_identifier": "D", "distance": 0.625}},
            {{"left_identifier": "B", "right_identifier": "A", "distance": 0.125}},
            {{"left_identifier": "B", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "C", "distance": 0.625}},
            {{"left_identifier": "B", "right_identifier": "D", "distance": 0.5}},
            {{"left_identifier": "C", "right_identifier": "A", "distance": 0.5}},
            {{"left_identifier": "C", "right_identifier": "B", "distance": 0.625}},
            {{"left_identifier": "C", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "D", "distance": 0.125}},
            {{"left_identifier": "D", "right_identifier": "A", "distance": 0.625}},
            {{"left_identifier": "D", "right_identifier": "B", "distance": 0.5}},
            {{"left_identifier": "D", "right_identifier": "C", "distance": 0.125}},
            {{"left_identifier": "D", "right_identifier": "D", "distance": 0.0}},
        ],
    }},
    "dna-raw-distance-gaps": {{
        "summary": {{
            "sequence_count": 4,
            "alignment_length": 6,
            "pairwise_deletion": True,
        }},
        "rows_name": "distance-matrix.tsv",
        "rows": [
            {{"left_identifier": "A", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "B", "distance": 1 / 6}},
            {{"left_identifier": "A", "right_identifier": "C", "distance": 1 / 2}},
            {{"left_identifier": "A", "right_identifier": "D", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "A", "distance": 1 / 6}},
            {{"left_identifier": "B", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "C", "distance": 1 / 2}},
            {{"left_identifier": "B", "right_identifier": "D", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "A", "distance": 1 / 2}},
            {{"left_identifier": "C", "right_identifier": "B", "distance": 1 / 2}},
            {{"left_identifier": "C", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "D", "distance": 1 / 4}},
            {{"left_identifier": "D", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "D", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "D", "right_identifier": "C", "distance": 1 / 4}},
            {{"left_identifier": "D", "right_identifier": "D", "distance": 0.0}},
        ],
    }},
    "dna-raw-distance-gaps-complete-deletion": {{
        "summary": {{
            "sequence_count": 4,
            "alignment_length": 6,
            "pairwise_deletion": False,
        }},
        "rows_name": "distance-matrix.tsv",
        "rows": [
            {{"left_identifier": "A", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "C", "distance": 1 / 4}},
            {{"left_identifier": "A", "right_identifier": "D", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "C", "distance": 1 / 4}},
            {{"left_identifier": "B", "right_identifier": "D", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "A", "distance": 1 / 4}},
            {{"left_identifier": "C", "right_identifier": "B", "distance": 1 / 4}},
            {{"left_identifier": "C", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "D", "distance": 1 / 4}},
            {{"left_identifier": "D", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "D", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "D", "right_identifier": "C", "distance": 1 / 4}},
            {{"left_identifier": "D", "right_identifier": "D", "distance": 0.0}},
        ],
    }},
    "dna-raw-distance-ambiguity": {{
        "summary": {{
            "sequence_count": 3,
            "alignment_length": 6,
            "pairwise_deletion": True,
        }},
        "rows_name": "distance-matrix.tsv",
        "rows": [
            {{"left_identifier": "A", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "C", "distance": 0.0}},
        ],
    }},
    "dna-raw-distance-identical": {{
        "summary": {{
            "sequence_count": 4,
            "alignment_length": 8,
            "pairwise_deletion": False,
        }},
        "rows_name": "distance-matrix.tsv",
        "rows": [
            {{"left_identifier": "A", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "C", "distance": 1 / 8}},
            {{"left_identifier": "A", "right_identifier": "D", "distance": 1 / 8}},
            {{"left_identifier": "B", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "C", "distance": 1 / 8}},
            {{"left_identifier": "B", "right_identifier": "D", "distance": 1 / 8}},
            {{"left_identifier": "C", "right_identifier": "A", "distance": 1 / 8}},
            {{"left_identifier": "C", "right_identifier": "B", "distance": 1 / 8}},
            {{"left_identifier": "C", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "D", "distance": 1 / 8}},
            {{"left_identifier": "D", "right_identifier": "A", "distance": 1 / 8}},
            {{"left_identifier": "D", "right_identifier": "B", "distance": 1 / 8}},
            {{"left_identifier": "D", "right_identifier": "C", "distance": 1 / 8}},
            {{"left_identifier": "D", "right_identifier": "D", "distance": 0.0}},
        ],
    }},
    "dna-raw-distance-high-divergence": {{
        "summary": {{
            "sequence_count": 3,
            "alignment_length": 4,
            "pairwise_deletion": False,
        }},
        "rows_name": "distance-matrix.tsv",
        "rows": [
            {{"left_identifier": "A", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "B", "distance": 1.0}},
            {{"left_identifier": "A", "right_identifier": "C", "distance": 0.25}},
            {{"left_identifier": "B", "right_identifier": "A", "distance": 1.0}},
            {{"left_identifier": "B", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "C", "distance": 0.75}},
            {{"left_identifier": "C", "right_identifier": "A", "distance": 0.25}},
            {{"left_identifier": "C", "right_identifier": "B", "distance": 0.75}},
            {{"left_identifier": "C", "right_identifier": "C", "distance": 0.0}},
        ],
    }},
    "dna-raw-distance-missing-data": {{
        "summary": {{
            "sequence_count": 3,
            "alignment_length": 6,
            "pairwise_deletion": True,
        }},
        "rows_name": "distance-matrix.tsv",
        "rows": [
            {{"left_identifier": "A", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "C", "distance": 0.0}},
        ],
    }},
    "dna-jc69-distance-clean": {{
        "summary": {{
            "sequence_count": 4,
            "alignment_length": 8,
            "pairwise_deletion": False,
        }},
        "rows_name": "distance-matrix.tsv",
        "rows": [
            {{"left_identifier": "A", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "B", "distance": 0.136741167595466}},
            {{"left_identifier": "A", "right_identifier": "C", "distance": 0.823959216501082}},
            {{"left_identifier": "A", "right_identifier": "D", "distance": 1.343819601921041}},
            {{"left_identifier": "B", "right_identifier": "A", "distance": 0.136741167595466}},
            {{"left_identifier": "B", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "C", "distance": 1.343819601921041}},
            {{"left_identifier": "B", "right_identifier": "D", "distance": 0.823959216501082}},
            {{"left_identifier": "C", "right_identifier": "A", "distance": 0.823959216501082}},
            {{"left_identifier": "C", "right_identifier": "B", "distance": 1.343819601921041}},
            {{"left_identifier": "C", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "D", "distance": 0.136741167595466}},
            {{"left_identifier": "D", "right_identifier": "A", "distance": 1.343819601921041}},
            {{"left_identifier": "D", "right_identifier": "B", "distance": 0.823959216501082}},
            {{"left_identifier": "D", "right_identifier": "C", "distance": 0.136741167595466}},
            {{"left_identifier": "D", "right_identifier": "D", "distance": 0.0}},
        ],
    }},
    "dna-jc69-distance-gaps": {{
        "summary": {{
            "sequence_count": 4,
            "alignment_length": 6,
            "pairwise_deletion": True,
        }},
        "rows_name": "distance-matrix.tsv",
        "rows": [
            {{"left_identifier": "A", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "B", "distance": 0.188485821210679}},
            {{"left_identifier": "A", "right_identifier": "C", "distance": 0.823959216501082}},
            {{"left_identifier": "A", "right_identifier": "D", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "A", "distance": 0.188485821210679}},
            {{"left_identifier": "B", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "C", "distance": 0.823959216501082}},
            {{"left_identifier": "B", "right_identifier": "D", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "A", "distance": 0.823959216501082}},
            {{"left_identifier": "C", "right_identifier": "B", "distance": 0.823959216501082}},
            {{"left_identifier": "C", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "D", "distance": 0.304098831081123}},
            {{"left_identifier": "D", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "D", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "D", "right_identifier": "C", "distance": 0.304098831081123}},
            {{"left_identifier": "D", "right_identifier": "D", "distance": 0.0}},
        ],
    }},
    "dna-jc69-distance-gaps-complete-deletion": {{
        "summary": {{
            "sequence_count": 4,
            "alignment_length": 6,
            "pairwise_deletion": False,
        }},
        "rows_name": "distance-matrix.tsv",
        "rows": [
            {{"left_identifier": "A", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "C", "distance": 0.304098831081123}},
            {{"left_identifier": "A", "right_identifier": "D", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "C", "distance": 0.304098831081123}},
            {{"left_identifier": "B", "right_identifier": "D", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "A", "distance": 0.304098831081123}},
            {{"left_identifier": "C", "right_identifier": "B", "distance": 0.304098831081123}},
            {{"left_identifier": "C", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "D", "distance": 0.304098831081123}},
            {{"left_identifier": "D", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "D", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "D", "right_identifier": "C", "distance": 0.304098831081123}},
            {{"left_identifier": "D", "right_identifier": "D", "distance": 0.0}},
        ],
    }},
    "dna-jc69-distance-ambiguity": {{
        "summary": {{
            "sequence_count": 3,
            "alignment_length": 6,
            "pairwise_deletion": True,
        }},
        "rows_name": "distance-matrix.tsv",
        "rows": [
            {{"left_identifier": "A", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "C", "distance": 0.0}},
        ],
    }},
    "dna-jc69-distance-identical": {{
        "summary": {{
            "sequence_count": 4,
            "alignment_length": 8,
            "pairwise_deletion": False,
        }},
        "rows_name": "distance-matrix.tsv",
        "rows": [
            {{"left_identifier": "A", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "C", "distance": 0.136741167595466}},
            {{"left_identifier": "A", "right_identifier": "D", "distance": 0.136741167595466}},
            {{"left_identifier": "B", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "C", "distance": 0.136741167595466}},
            {{"left_identifier": "B", "right_identifier": "D", "distance": 0.136741167595466}},
            {{"left_identifier": "C", "right_identifier": "A", "distance": 0.136741167595466}},
            {{"left_identifier": "C", "right_identifier": "B", "distance": 0.136741167595466}},
            {{"left_identifier": "C", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "D", "distance": 0.136741167595466}},
            {{"left_identifier": "D", "right_identifier": "A", "distance": 0.136741167595466}},
            {{"left_identifier": "D", "right_identifier": "B", "distance": 0.136741167595466}},
            {{"left_identifier": "D", "right_identifier": "C", "distance": 0.136741167595466}},
            {{"left_identifier": "D", "right_identifier": "D", "distance": 0.0}},
        ],
    }},
    "dna-jc69-distance-high-divergence": {{
        "summary": {{
            "sequence_count": 3,
            "alignment_length": 4,
            "pairwise_deletion": False,
        }},
        "rows_name": "distance-matrix.tsv",
        "rows": [
            {{"left_identifier": "A", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "B", "distance": "", "distance_status": "undefined"}},
            {{"left_identifier": "A", "right_identifier": "C", "distance": 0.304098831081123}},
            {{"left_identifier": "B", "right_identifier": "A", "distance": "", "distance_status": "undefined"}},
            {{"left_identifier": "B", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "C", "distance": "", "distance_status": "infinite"}},
            {{"left_identifier": "C", "right_identifier": "A", "distance": 0.304098831081123}},
            {{"left_identifier": "C", "right_identifier": "B", "distance": "", "distance_status": "infinite"}},
            {{"left_identifier": "C", "right_identifier": "C", "distance": 0.0}},
        ],
    }},
    "dna-jc69-distance-missing-data": {{
        "summary": {{
            "sequence_count": 3,
            "alignment_length": 6,
            "pairwise_deletion": True,
        }},
        "rows_name": "distance-matrix.tsv",
        "rows": [
            {{"left_identifier": "A", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "C", "distance": 0.0}},
        ],
    }},
    "dna-k80-distance-clean": {{
        "summary": {{
            "sequence_count": 4,
            "alignment_length": 8,
            "pairwise_deletion": False,
        }},
        "rows_name": "distance-matrix.tsv",
        "rows": [
            {{"left_identifier": "A", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "B", "distance": 0.14384103622589}},
            {{"left_identifier": "A", "right_identifier": "C", "distance": "", "distance_status": "infinite"}},
            {{"left_identifier": "A", "right_identifier": "D", "distance": "", "distance_status": "infinite"}},
            {{"left_identifier": "B", "right_identifier": "A", "distance": 0.14384103622589}},
            {{"left_identifier": "B", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "C", "distance": "", "distance_status": "infinite"}},
            {{"left_identifier": "B", "right_identifier": "D", "distance": "", "distance_status": "infinite"}},
            {{"left_identifier": "C", "right_identifier": "A", "distance": "", "distance_status": "infinite"}},
            {{"left_identifier": "C", "right_identifier": "B", "distance": "", "distance_status": "infinite"}},
            {{"left_identifier": "C", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "D", "distance": 0.14384103622589}},
            {{"left_identifier": "D", "right_identifier": "A", "distance": "", "distance_status": "infinite"}},
            {{"left_identifier": "D", "right_identifier": "B", "distance": "", "distance_status": "infinite"}},
            {{"left_identifier": "D", "right_identifier": "C", "distance": 0.14384103622589}},
            {{"left_identifier": "D", "right_identifier": "D", "distance": 0.0}},
        ],
    }},
    "dna-k80-distance-gaps": {{
        "summary": {{
            "sequence_count": 4,
            "alignment_length": 6,
            "pairwise_deletion": True,
        }},
        "rows_name": "distance-matrix.tsv",
        "rows": [
            {{"left_identifier": "A", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "B", "distance": 0.192527055424018}},
            {{"left_identifier": "A", "right_identifier": "C", "distance": 0.997246011641068}},
            {{"left_identifier": "A", "right_identifier": "D", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "A", "distance": 0.192527055424018}},
            {{"left_identifier": "B", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "C", "distance": "", "distance_status": "infinite"}},
            {{"left_identifier": "B", "right_identifier": "D", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "A", "distance": 0.997246011641068}},
            {{"left_identifier": "C", "right_identifier": "B", "distance": "", "distance_status": "infinite"}},
            {{"left_identifier": "C", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "D", "distance": 0.346573590279973}},
            {{"left_identifier": "D", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "D", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "D", "right_identifier": "C", "distance": 0.346573590279973}},
            {{"left_identifier": "D", "right_identifier": "D", "distance": 0.0}},
        ],
    }},
    "dna-k80-distance-gaps-complete-deletion": {{
        "summary": {{
            "sequence_count": 4,
            "alignment_length": 6,
            "pairwise_deletion": False,
        }},
        "rows_name": "distance-matrix.tsv",
        "rows": [
            {{"left_identifier": "A", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "C", "distance": 0.346573590279973}},
            {{"left_identifier": "A", "right_identifier": "D", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "C", "distance": 0.346573590279973}},
            {{"left_identifier": "B", "right_identifier": "D", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "A", "distance": 0.346573590279973}},
            {{"left_identifier": "C", "right_identifier": "B", "distance": 0.346573590279973}},
            {{"left_identifier": "C", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "D", "distance": 0.346573590279973}},
            {{"left_identifier": "D", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "D", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "D", "right_identifier": "C", "distance": 0.346573590279973}},
            {{"left_identifier": "D", "right_identifier": "D", "distance": 0.0}},
        ],
    }},
    "dna-k80-distance-ambiguity": {{
        "summary": {{
            "sequence_count": 3,
            "alignment_length": 6,
            "pairwise_deletion": True,
        }},
        "rows_name": "distance-matrix.tsv",
        "rows": [
            {{"left_identifier": "A", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "C", "distance": 0.0}},
        ],
    }},
    "dna-k80-distance-identical": {{
        "summary": {{
            "sequence_count": 4,
            "alignment_length": 8,
            "pairwise_deletion": False,
        }},
        "rows_name": "distance-matrix.tsv",
        "rows": [
            {{"left_identifier": "A", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "C", "distance": 0.14384103622589}},
            {{"left_identifier": "A", "right_identifier": "D", "distance": 0.138686214425207}},
            {{"left_identifier": "B", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "C", "distance": 0.14384103622589}},
            {{"left_identifier": "B", "right_identifier": "D", "distance": 0.138686214425207}},
            {{"left_identifier": "C", "right_identifier": "A", "distance": 0.14384103622589}},
            {{"left_identifier": "C", "right_identifier": "B", "distance": 0.14384103622589}},
            {{"left_identifier": "C", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "D", "distance": 0.138686214425207}},
            {{"left_identifier": "D", "right_identifier": "A", "distance": 0.138686214425207}},
            {{"left_identifier": "D", "right_identifier": "B", "distance": 0.138686214425207}},
            {{"left_identifier": "D", "right_identifier": "C", "distance": 0.138686214425207}},
            {{"left_identifier": "D", "right_identifier": "D", "distance": 0.0}},
        ],
    }},
    "dna-k80-distance-high-divergence": {{
        "summary": {{
            "sequence_count": 3,
            "alignment_length": 4,
            "pairwise_deletion": False,
        }},
        "rows_name": "distance-matrix.tsv",
        "rows": [
            {{"left_identifier": "A", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "B", "distance": "", "distance_status": "undefined"}},
            {{"left_identifier": "A", "right_identifier": "C", "distance": 0.317127831365877}},
            {{"left_identifier": "B", "right_identifier": "A", "distance": "", "distance_status": "undefined"}},
            {{"left_identifier": "B", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "C", "distance": "", "distance_status": "undefined"}},
            {{"left_identifier": "C", "right_identifier": "A", "distance": 0.317127831365877}},
            {{"left_identifier": "C", "right_identifier": "B", "distance": "", "distance_status": "undefined"}},
            {{"left_identifier": "C", "right_identifier": "C", "distance": 0.0}},
        ],
    }},
    "dna-k80-distance-missing-data": {{
        "summary": {{
            "sequence_count": 3,
            "alignment_length": 6,
            "pairwise_deletion": True,
        }},
        "rows_name": "distance-matrix.tsv",
        "rows": [
            {{"left_identifier": "A", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "C", "distance": 0.0}},
        ],
    }},
    "dna-f81-distance-clean": {{
        "summary": {{
            "sequence_count": 4,
            "alignment_length": 8,
            "pairwise_deletion": False,
        }},
        "rows_name": "distance-matrix.tsv",
        "rows": [
            {{"left_identifier": "A", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "B", "distance": 0.139677632499716}},
            {{"left_identifier": "A", "right_identifier": "C", "distance": 1.02539386295172}},
            {{"left_identifier": "A", "right_identifier": "D", "distance": "", "distance_status": "undefined"}},
            {{"left_identifier": "B", "right_identifier": "A", "distance": 0.139677632499716}},
            {{"left_identifier": "B", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "C", "distance": "", "distance_status": "undefined"}},
            {{"left_identifier": "B", "right_identifier": "D", "distance": 1.02539386295172}},
            {{"left_identifier": "C", "right_identifier": "A", "distance": 1.02539386295172}},
            {{"left_identifier": "C", "right_identifier": "B", "distance": "", "distance_status": "undefined"}},
            {{"left_identifier": "C", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "D", "distance": 0.139677632499716}},
            {{"left_identifier": "D", "right_identifier": "A", "distance": "", "distance_status": "undefined"}},
            {{"left_identifier": "D", "right_identifier": "B", "distance": 1.02539386295172}},
            {{"left_identifier": "D", "right_identifier": "C", "distance": 0.139677632499716}},
            {{"left_identifier": "D", "right_identifier": "D", "distance": 0.0}},
        ],
    }},
    "dna-f81-distance-gaps": {{
        "summary": {{
            "sequence_count": 4,
            "alignment_length": 6,
            "pairwise_deletion": True,
        }},
        "rows_name": "distance-matrix.tsv",
        "rows": [
            {{"left_identifier": "A", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "B", "distance": 0.189450794670797}},
            {{"left_identifier": "A", "right_identifier": "C", "distance": 0.850269988792838}},
            {{"left_identifier": "A", "right_identifier": "D", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "A", "distance": 0.189450794670797}},
            {{"left_identifier": "B", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "C", "distance": 0.850269988792838}},
            {{"left_identifier": "B", "right_identifier": "D", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "A", "distance": 0.850269988792838}},
            {{"left_identifier": "C", "right_identifier": "B", "distance": 0.850269988792838}},
            {{"left_identifier": "C", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "D", "distance": 0.306764262035593}},
            {{"left_identifier": "D", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "D", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "D", "right_identifier": "C", "distance": 0.306764262035593}},
            {{"left_identifier": "D", "right_identifier": "D", "distance": 0.0}},
        ],
    }},
    "dna-f81-distance-gaps-complete-deletion": {{
        "summary": {{
            "sequence_count": 4,
            "alignment_length": 6,
            "pairwise_deletion": False,
        }},
        "rows_name": "distance-matrix.tsv",
        "rows": [
            {{"left_identifier": "A", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "C", "distance": 0.306764262035593}},
            {{"left_identifier": "A", "right_identifier": "D", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "C", "distance": 0.306764262035593}},
            {{"left_identifier": "B", "right_identifier": "D", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "A", "distance": 0.306764262035593}},
            {{"left_identifier": "C", "right_identifier": "B", "distance": 0.306764262035593}},
            {{"left_identifier": "C", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "D", "distance": 0.306764262035593}},
            {{"left_identifier": "D", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "D", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "D", "right_identifier": "C", "distance": 0.306764262035593}},
            {{"left_identifier": "D", "right_identifier": "D", "distance": 0.0}},
        ],
    }},
    "dna-f81-distance-ambiguity": {{
        "summary": {{
            "sequence_count": 3,
            "alignment_length": 6,
            "pairwise_deletion": True,
        }},
        "rows_name": "distance-matrix.tsv",
        "rows": [
            {{"left_identifier": "A", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "C", "distance": 0.0}},
        ],
    }},
    "dna-f81-distance-identical": {{
        "summary": {{
            "sequence_count": 4,
            "alignment_length": 8,
            "pairwise_deletion": False,
        }},
        "rows_name": "distance-matrix.tsv",
        "rows": [
            {{"left_identifier": "A", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "C", "distance": 0.136845675929744}},
            {{"left_identifier": "A", "right_identifier": "D", "distance": 0.136845675929744}},
            {{"left_identifier": "B", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "C", "distance": 0.136845675929744}},
            {{"left_identifier": "B", "right_identifier": "D", "distance": 0.136845675929744}},
            {{"left_identifier": "C", "right_identifier": "A", "distance": 0.136845675929744}},
            {{"left_identifier": "C", "right_identifier": "B", "distance": 0.136845675929744}},
            {{"left_identifier": "C", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "D", "distance": 0.136845675929744}},
            {{"left_identifier": "D", "right_identifier": "A", "distance": 0.136845675929744}},
            {{"left_identifier": "D", "right_identifier": "B", "distance": 0.136845675929744}},
            {{"left_identifier": "D", "right_identifier": "C", "distance": 0.136845675929744}},
            {{"left_identifier": "D", "right_identifier": "D", "distance": 0.0}},
        ],
    }},
    "dna-f81-distance-high-divergence": {{
        "summary": {{
            "sequence_count": 3,
            "alignment_length": 4,
            "pairwise_deletion": False,
        }},
        "rows_name": "distance-matrix.tsv",
        "rows": [
            {{"left_identifier": "A", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "B", "distance": "", "distance_status": "undefined"}},
            {{"left_identifier": "A", "right_identifier": "C", "distance": 0.35103770986336}},
            {{"left_identifier": "B", "right_identifier": "A", "distance": "", "distance_status": "undefined"}},
            {{"left_identifier": "B", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "C", "distance": "", "distance_status": "undefined"}},
            {{"left_identifier": "C", "right_identifier": "A", "distance": 0.35103770986336}},
            {{"left_identifier": "C", "right_identifier": "B", "distance": "", "distance_status": "undefined"}},
            {{"left_identifier": "C", "right_identifier": "C", "distance": 0.0}},
        ],
    }},
    "dna-f81-distance-missing-data": {{
        "summary": {{
            "sequence_count": 3,
            "alignment_length": 6,
            "pairwise_deletion": True,
        }},
        "rows_name": "distance-matrix.tsv",
        "rows": [
            {{"left_identifier": "A", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "C", "distance": 0.0}},
        ],
    }},
    "dna-tn93-distance-clean": {{
        "summary": {{
            "sequence_count": 4,
            "alignment_length": 8,
            "pairwise_deletion": False,
        }},
        "rows_name": "distance-matrix.tsv",
        "rows": [
            {{"left_identifier": "A", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "B", "distance": "", "distance_status": "undefined"}},
            {{"left_identifier": "A", "right_identifier": "C", "distance": "", "distance_status": "undefined"}},
            {{"left_identifier": "A", "right_identifier": "D", "distance": "", "distance_status": "undefined"}},
            {{"left_identifier": "B", "right_identifier": "A", "distance": "", "distance_status": "undefined"}},
            {{"left_identifier": "B", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "C", "distance": "", "distance_status": "undefined"}},
            {{"left_identifier": "B", "right_identifier": "D", "distance": "", "distance_status": "undefined"}},
            {{"left_identifier": "C", "right_identifier": "A", "distance": "", "distance_status": "undefined"}},
            {{"left_identifier": "C", "right_identifier": "B", "distance": "", "distance_status": "undefined"}},
            {{"left_identifier": "C", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "D", "distance": "", "distance_status": "undefined"}},
            {{"left_identifier": "D", "right_identifier": "A", "distance": "", "distance_status": "undefined"}},
            {{"left_identifier": "D", "right_identifier": "B", "distance": "", "distance_status": "undefined"}},
            {{"left_identifier": "D", "right_identifier": "C", "distance": "", "distance_status": "undefined"}},
            {{"left_identifier": "D", "right_identifier": "D", "distance": 0.0}},
        ],
    }},
    "dna-tn93-distance-gaps": {{
        "summary": {{
            "sequence_count": 4,
            "alignment_length": 6,
            "pairwise_deletion": True,
        }},
        "rows_name": "distance-matrix.tsv",
        "rows": [
            {{"left_identifier": "A", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "B", "distance": 0.194518264955326}},
            {{"left_identifier": "A", "right_identifier": "C", "distance": "", "distance_status": "undefined"}},
            {{"left_identifier": "A", "right_identifier": "D", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "A", "distance": 0.194518264955326}},
            {{"left_identifier": "B", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "C", "distance": "", "distance_status": "undefined"}},
            {{"left_identifier": "B", "right_identifier": "D", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "A", "distance": "", "distance_status": "undefined"}},
            {{"left_identifier": "C", "right_identifier": "B", "distance": "", "distance_status": "undefined"}},
            {{"left_identifier": "C", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "D", "distance": 0.627121810119612}},
            {{"left_identifier": "D", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "D", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "D", "right_identifier": "C", "distance": 0.627121810119612}},
            {{"left_identifier": "D", "right_identifier": "D", "distance": 0.0}},
        ],
    }},
    "dna-tn93-distance-gaps-complete-deletion": {{
        "summary": {{
            "sequence_count": 4,
            "alignment_length": 6,
            "pairwise_deletion": False,
        }},
        "rows_name": "distance-matrix.tsv",
        "rows": [
            {{"left_identifier": "A", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "C", "distance": 0.627121810119612}},
            {{"left_identifier": "A", "right_identifier": "D", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "C", "distance": 0.627121810119612}},
            {{"left_identifier": "B", "right_identifier": "D", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "A", "distance": 0.627121810119612}},
            {{"left_identifier": "C", "right_identifier": "B", "distance": 0.627121810119612}},
            {{"left_identifier": "C", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "D", "distance": 0.627121810119612}},
            {{"left_identifier": "D", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "D", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "D", "right_identifier": "C", "distance": 0.627121810119612}},
            {{"left_identifier": "D", "right_identifier": "D", "distance": 0.0}},
        ],
    }},
    "dna-tn93-distance-ambiguity": {{
        "summary": {{
            "sequence_count": 3,
            "alignment_length": 6,
            "pairwise_deletion": True,
        }},
        "rows_name": "distance-matrix.tsv",
        "rows": [
            {{"left_identifier": "A", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "C", "distance": 0.0}},
        ],
    }},
    "dna-tn93-distance-identical": {{
        "summary": {{
            "sequence_count": 4,
            "alignment_length": 8,
            "pairwise_deletion": False,
        }},
        "rows_name": "distance-matrix.tsv",
        "rows": [
            {{"left_identifier": "A", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "C", "distance": 0.182459298648674}},
            {{"left_identifier": "A", "right_identifier": "D", "distance": 0.138858236087742}},
            {{"left_identifier": "B", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "C", "distance": 0.182459298648674}},
            {{"left_identifier": "B", "right_identifier": "D", "distance": 0.138858236087742}},
            {{"left_identifier": "C", "right_identifier": "A", "distance": 0.182459298648674}},
            {{"left_identifier": "C", "right_identifier": "B", "distance": 0.182459298648674}},
            {{"left_identifier": "C", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "D", "distance": 0.138858236087742}},
            {{"left_identifier": "D", "right_identifier": "A", "distance": 0.138858236087742}},
            {{"left_identifier": "D", "right_identifier": "B", "distance": 0.138858236087742}},
            {{"left_identifier": "D", "right_identifier": "C", "distance": 0.138858236087742}},
            {{"left_identifier": "D", "right_identifier": "D", "distance": 0.0}},
        ],
    }},
    "dna-tn93-distance-high-divergence": {{
        "summary": {{
            "sequence_count": 3,
            "alignment_length": 4,
            "pairwise_deletion": False,
        }},
        "rows_name": "distance-matrix.tsv",
        "rows": [
            {{"left_identifier": "A", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "B", "distance": "", "distance_status": "undefined"}},
            {{"left_identifier": "A", "right_identifier": "C", "distance": "", "distance_status": "undefined"}},
            {{"left_identifier": "B", "right_identifier": "A", "distance": "", "distance_status": "undefined"}},
            {{"left_identifier": "B", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "C", "distance": "", "distance_status": "undefined"}},
            {{"left_identifier": "C", "right_identifier": "A", "distance": "", "distance_status": "undefined"}},
            {{"left_identifier": "C", "right_identifier": "B", "distance": "", "distance_status": "undefined"}},
            {{"left_identifier": "C", "right_identifier": "C", "distance": 0.0}},
        ],
    }},
    "dna-tn93-distance-missing-data": {{
        "summary": {{
            "sequence_count": 3,
            "alignment_length": 6,
            "pairwise_deletion": True,
        }},
        "rows_name": "distance-matrix.tsv",
        "rows": [
            {{"left_identifier": "A", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "A", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "B", "right_identifier": "C", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "A", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "B", "distance": 0.0}},
            {{"left_identifier": "C", "right_identifier": "C", "distance": 0.0}},
        ],
    }},
    "ace-continuous-balanced-rooted-ultrametric": {{
        "summary": {{
            "trait": "response",
            "taxon_count": 4,
            "excluded_taxon_count": 0,
            "dropped_missing_taxa": [],
            "dropped_non_numeric_taxa": [],
            "internal_node_count": 3,
            "method": "pic",
            "tree_is_ultrametric": True,
            "minimum_root_to_tip_depth": 0.3,
            "maximum_root_to_tip_depth": 0.3,
        }},
        "rows_name": "continuous-ancestral.tsv",
        "rows": [
            {{"node_id": 5, "node": "A|B|C|D", "estimate": 2.80555555555555, "standard_error": 0.670820393249937, "lower_95_interval": 1.49077174469068, "upper_95_interval": 4.12033936642043}},
            {{"node_id": 6, "node": "A|B", "estimate": 2.25, "standard_error": 0.447213595499958, "lower_95_interval": 1.37347745942342, "upper_95_interval": 3.12652254057658}},
            {{"node_id": 7, "node": "C|D", "estimate": 3.25, "standard_error": 0.632455532033676, "lower_95_interval": 2.01040993539088, "upper_95_interval": 4.48959006460912}},
        ],
    }},
    "ace-continuous-pectinate-non-ultrametric": {{
        "summary": {{
            "trait": "response",
            "taxon_count": 4,
            "excluded_taxon_count": 0,
            "dropped_missing_taxa": [],
            "dropped_non_numeric_taxa": [],
            "internal_node_count": 3,
            "method": "pic",
            "tree_is_ultrametric": False,
            "minimum_root_to_tip_depth": 0.1,
            "maximum_root_to_tip_depth": 0.3,
        }},
        "rows_name": "continuous-ancestral.tsv",
        "rows": [
            {{"node_id": 5, "node": "A|B|C|D", "estimate": 3.38461538461539, "standard_error": 0.509901951359279, "lower_95_interval": 2.38522592430451, "upper_95_interval": 4.38400484492626}},
            {{"node_id": 6, "node": "A|B|C", "estimate": 2.4, "standard_error": 0.5, "lower_95_interval": 1.42001800772997, "upper_95_interval": 3.37998199227003}},
            {{"node_id": 7, "node": "A|B", "estimate": 2.25, "standard_error": 0.447213595499958, "lower_95_interval": 1.37347745942342, "upper_95_interval": 3.12652254057658}},
        ],
    }},
    "ace-continuous-balanced-six-taxon": {{
        "summary": {{
            "trait": "response_growth",
            "taxon_count": 6,
            "excluded_taxon_count": 0,
            "dropped_missing_taxa": [],
            "dropped_non_numeric_taxa": [],
            "internal_node_count": 5,
            "method": "pic",
            "tree_is_ultrametric": True,
            "minimum_root_to_tip_depth": 3.0,
            "maximum_root_to_tip_depth": 3.0,
        }},
        "rows_name": "continuous-ancestral.tsv",
        "rows": [
            {{"node_id": 7, "node": "A|B|C|D|E|F", "estimate": 2.86764705882353, "standard_error": 2.06155281280883, "lower_95_interval": -1.17292220650902, "upper_95_interval": 6.90821632415608}},
            {{"node_id": 8, "node": "A|B|C|D", "estimate": 2.25, "standard_error": 1.73205080756888, "lower_95_interval": -1.14475720222852, "upper_95_interval": 5.64475720222852}},
            {{"node_id": 9, "node": "A|B", "estimate": 1.75, "standard_error": 1.41421356237309, "lower_95_interval": -1.02180764869936, "upper_95_interval": 4.52180764869936}},
            {{"node_id": 10, "node": "C|D", "estimate": 2.75, "standard_error": 1.41421356237309, "lower_95_interval": -0.021807648699356, "upper_95_interval": 5.52180764869936}},
            {{"node_id": 11, "node": "E|F", "estimate": 3.75, "standard_error": 1.41421356237309, "lower_95_interval": 0.978192351300644, "upper_95_interval": 6.52180764869936}},
        ],
    }},
    "ace-continuous-missing-values-pruned": {{
        "summary": {{
            "trait": "response_growth",
            "taxon_count": 4,
            "excluded_taxon_count": 2,
            "dropped_missing_taxa": ["B"],
            "dropped_non_numeric_taxa": ["C"],
            "internal_node_count": 3,
            "method": "pic",
            "tree_is_ultrametric": True,
            "minimum_root_to_tip_depth": 3.0,
            "maximum_root_to_tip_depth": 3.0,
        }},
        "rows_name": "continuous-ancestral.tsv",
        "rows": [
            {{"node_id": 5, "node": "A|D|E|F", "estimate": 2.91666666666667, "standard_error": 2.12132034355964, "lower_95_interval": -1.24104480638237, "upper_95_interval": 7.0743781397157}},
            {{"node_id": 6, "node": "A|D", "estimate": 2.25, "standard_error": 2.0, "lower_95_interval": -1.66992796908011, "upper_95_interval": 6.16992796908011}},
            {{"node_id": 7, "node": "E|F", "estimate": 3.75, "standard_error": 1.41421356237309, "lower_95_interval": 0.978192351300644, "upper_95_interval": 6.52180764869936}},
        ],
    }},
    "ace-discrete-binary-balanced-rooted-ultrametric": {{
        "summary": {{
            "trait": "presence",
            "taxon_count": 4,
            "excluded_taxon_count": 0,
            "dropped_missing_taxa": [],
            "internal_node_count": 3,
            "model": "equal-rates",
            "state_count": 2,
            "state_labels": ["0", "1"],
            "log_likelihood": -1.757639727067767,
            "parameter_count": 1,
            "aic": 5.515279454135534,
            "overparameterized": False,
            "baseline_model": None,
            "baseline_delta_aic": None,
            "preferred_model_by_aic": None,
            "transition_rate_rows": [
                {{"source_state": "0", "target_state": "1", "transition_allowed": True, "step_distance": 1, "rate": 1.9279052278659354}},
                {{"source_state": "1", "target_state": "0", "transition_allowed": True, "step_distance": 1, "rate": 1.9279052278659354}},
            ],
        }},
        "rows_name": "discrete-ancestral.tsv",
        "rows": [
            {{"node_id": 5, "node": "A|B|C|D", "state": 0, "posterior_probability": 0.443328715428089, "most_likely_state": 1, "max_posterior_probability": 0.556671284571911}},
            {{"node_id": 5, "node": "A|B|C|D", "state": 1, "posterior_probability": 0.556671284571911, "most_likely_state": 1, "max_posterior_probability": 0.556671284571911}},
            {{"node_id": 6, "node": "A|B", "state": 0, "posterior_probability": 0.944173762264708, "most_likely_state": 0, "max_posterior_probability": 0.944173762264708}},
            {{"node_id": 6, "node": "A|B", "state": 1, "posterior_probability": 0.055826237735292, "most_likely_state": 0, "max_posterior_probability": 0.944173762264708}},
            {{"node_id": 7, "node": "C|D", "state": 0, "posterior_probability": 0.197937244108908, "most_likely_state": 1, "max_posterior_probability": 0.802062755891092}},
            {{"node_id": 7, "node": "C|D", "state": 1, "posterior_probability": 0.802062755891092, "most_likely_state": 1, "max_posterior_probability": 0.802062755891092}},
        ],
    }},
    "ace-discrete-multistate-balanced-rooted-ultrametric": {{
        "summary": {{
            "trait": "region",
            "taxon_count": 4,
            "excluded_taxon_count": 0,
            "dropped_missing_taxa": [],
            "internal_node_count": 3,
            "model": "equal-rates",
            "state_count": 3,
            "state_labels": ["island", "north", "south"],
            "log_likelihood": -2.9885756563927695,
            "parameter_count": 1,
            "aic": 7.977151312785539,
            "overparameterized": False,
            "baseline_model": None,
            "baseline_delta_aic": None,
            "preferred_model_by_aic": None,
            "transition_rate_rows": [
                {{"source_state": "island", "target_state": "north", "transition_allowed": True, "step_distance": 1, "rate": 2.006603696572599}},
                {{"source_state": "island", "target_state": "south", "transition_allowed": True, "step_distance": 2, "rate": 2.006603696572599}},
                {{"source_state": "north", "target_state": "island", "transition_allowed": True, "step_distance": 1, "rate": 2.006603696572599}},
                {{"source_state": "north", "target_state": "south", "transition_allowed": True, "step_distance": 1, "rate": 2.006603696572599}},
                {{"source_state": "south", "target_state": "island", "transition_allowed": True, "step_distance": 2, "rate": 2.006603696572599}},
                {{"source_state": "south", "target_state": "north", "transition_allowed": True, "step_distance": 1, "rate": 2.006603696572599}},
            ],
        }},
        "rows_name": "discrete-ancestral.tsv",
        "rows": [
            {{"node_id": 5, "node": "A|B|C|D", "state": "island", "posterior_probability": 0.296861633587191, "most_likely_state": "north", "max_posterior_probability": 0.406276732825617}},
            {{"node_id": 5, "node": "A|B|C|D", "state": "north", "posterior_probability": 0.406276732825617, "most_likely_state": "north", "max_posterior_probability": 0.406276732825617}},
            {{"node_id": 5, "node": "A|B|C|D", "state": "south", "posterior_probability": 0.296861633587191, "most_likely_state": "north", "max_posterior_probability": 0.406276732825617}},
            {{"node_id": 6, "node": "A|B", "state": "island", "posterior_probability": 0.047363699539525, "most_likely_state": "north", "max_posterior_probability": 0.905272600920949}},
            {{"node_id": 6, "node": "A|B", "state": "north", "posterior_probability": 0.905272600920949, "most_likely_state": "north", "max_posterior_probability": 0.905272600920949}},
            {{"node_id": 6, "node": "A|B", "state": "south", "posterior_probability": 0.047363699539525, "most_likely_state": "north", "max_posterior_probability": 0.905272600920949}},
            {{"node_id": 7, "node": "C|D", "state": "island", "posterior_probability": 0.376356180907049, "most_likely_state": "south", "max_posterior_probability": 0.376356180907049}},
            {{"node_id": 7, "node": "C|D", "state": "north", "posterior_probability": 0.247287638185903, "most_likely_state": "south", "max_posterior_probability": 0.376356180907049}},
            {{"node_id": 7, "node": "C|D", "state": "south", "posterior_probability": 0.376356180907049, "most_likely_state": "south", "max_posterior_probability": 0.376356180907049}},
        ],
    }},
    "ace-discrete-multistate-pectinate-non-ultrametric": {{
        "summary": {{
            "trait": "region",
            "taxon_count": 4,
            "excluded_taxon_count": 0,
            "dropped_missing_taxa": [],
            "internal_node_count": 3,
            "model": "equal-rates",
            "state_count": 3,
            "state_labels": ["island", "north", "south"],
            "log_likelihood": -3.218012508155108,
            "parameter_count": 1,
            "aic": 8.436025016310216,
            "overparameterized": False,
            "baseline_model": None,
            "baseline_delta_aic": None,
            "preferred_model_by_aic": None,
            "transition_rate_rows": [
                {{"source_state": "island", "target_state": "north", "transition_allowed": True, "step_distance": 1, "rate": 3.6620407772788095}},
                {{"source_state": "island", "target_state": "south", "transition_allowed": True, "step_distance": 2, "rate": 3.6620407772788095}},
                {{"source_state": "north", "target_state": "island", "transition_allowed": True, "step_distance": 1, "rate": 3.6620407772788095}},
                {{"source_state": "north", "target_state": "south", "transition_allowed": True, "step_distance": 1, "rate": 3.6620407772788095}},
                {{"source_state": "south", "target_state": "island", "transition_allowed": True, "step_distance": 2, "rate": 3.6620407772788095}},
                {{"source_state": "south", "target_state": "north", "transition_allowed": True, "step_distance": 1, "rate": 3.6620407772788095}},
            ],
        }},
        "rows_name": "discrete-ancestral.tsv",
        "rows": [
            {{"node_id": 5, "node": "A|B|C|D", "state": "island", "posterior_probability": 0.499153982225924, "most_likely_state": "island", "max_posterior_probability": 0.499153982225924}},
            {{"node_id": 5, "node": "A|B|C|D", "state": "north", "posterior_probability": 0.235194582068628, "most_likely_state": "island", "max_posterior_probability": 0.499153982225924}},
            {{"node_id": 5, "node": "A|B|C|D", "state": "south", "posterior_probability": 0.265651435705448, "most_likely_state": "island", "max_posterior_probability": 0.499153982225924}},
            {{"node_id": 6, "node": "A|B|C", "state": "island", "posterior_probability": 0.241962771577934, "most_likely_state": "south", "max_posterior_probability": 0.439932318949848}},
            {{"node_id": 6, "node": "A|B|C", "state": "north", "posterior_probability": 0.318104909472218, "most_likely_state": "south", "max_posterior_probability": 0.439932318949848}},
            {{"node_id": 6, "node": "A|B|C", "state": "south", "posterior_probability": 0.439932318949848, "most_likely_state": "south", "max_posterior_probability": 0.439932318949848}},
            {{"node_id": 7, "node": "A|B", "state": "island", "posterior_probability": 0.125211497754402, "most_likely_state": "north", "max_posterior_probability": 0.719120150854377}},
            {{"node_id": 7, "node": "A|B", "state": "north", "posterior_probability": 0.719120150854377, "most_likely_state": "north", "max_posterior_probability": 0.719120150854377}},
            {{"node_id": 7, "node": "A|B", "state": "south", "posterior_probability": 0.155668351391221, "most_likely_state": "north", "max_posterior_probability": 0.719120150854377}},
        ],
    }},
    "ace-discrete-missing-values-pruned": {{
        "summary": {{
            "trait": "habitat",
            "taxon_count": 3,
            "excluded_taxon_count": 1,
            "dropped_missing_taxa": ["D"],
            "internal_node_count": 2,
            "model": "equal-rates",
            "state_count": 2,
            "state_labels": ["forest", "tundra"],
            "log_likelihood": -1.1455737733979303,
            "parameter_count": 1,
            "aic": 4.291147546795861,
            "overparameterized": False,
            "baseline_model": None,
            "baseline_delta_aic": None,
            "preferred_model_by_aic": None,
            "transition_rate_rows": [
                {{"source_state": "forest", "target_state": "tundra", "transition_allowed": True, "step_distance": 1, "rate": 2.2396993279589243}},
                {{"source_state": "tundra", "target_state": "forest", "transition_allowed": True, "step_distance": 1, "rate": 2.2396993279589243}},
            ],
        }},
        "rows_name": "discrete-ancestral.tsv",
        "rows": [
            {{"node_id": 4, "node": "A|B|C", "state": "forest", "posterior_probability": 0.560666992538983, "most_likely_state": "forest", "max_posterior_probability": 0.560666992538983}},
            {{"node_id": 4, "node": "A|B|C", "state": "tundra", "posterior_probability": 0.439333007461017, "most_likely_state": "forest", "max_posterior_probability": 0.560666992538983}},
            {{"node_id": 5, "node": "A|B", "state": "forest", "posterior_probability": 0.943307532034481, "most_likely_state": "forest", "max_posterior_probability": 0.943307532034481}},
            {{"node_id": 5, "node": "A|B", "state": "tundra", "posterior_probability": 0.056692467965519, "most_likely_state": "forest", "max_posterior_probability": 0.943307532034481}},
        ],
    }},
    "ace-discrete-sym-balanced-rooted-ultrametric": {{
        "summary": {{
            "trait": "region",
            "taxon_count": 4,
            "excluded_taxon_count": 0,
            "dropped_missing_taxa": [],
            "internal_node_count": 3,
            "model": "symmetric",
            "state_count": 3,
            "state_labels": ["island", "north", "south"],
            "log_likelihood": -2.7877249232794385,
            "parameter_count": 3,
            "aic": 11.575449846558877,
            "overparameterized": False,
            "baseline_model": "equal-rates",
            "baseline_delta_aic": 3.5982985333733435,
            "preferred_model_by_aic": "equal-rates",
            "transition_rate_rows": [
                {{"source_state": "island", "target_state": "north", "transition_allowed": True, "step_distance": 1, "rate": 1.2594253892995917}},
                {{"source_state": "island", "target_state": "south", "transition_allowed": True, "step_distance": 2, "rate": 28.400649574079157}},
                {{"source_state": "north", "target_state": "island", "transition_allowed": True, "step_distance": 1, "rate": 1.2593828279036432}},
                {{"source_state": "north", "target_state": "south", "transition_allowed": True, "step_distance": 1, "rate": 1.2594240474115468}},
                {{"source_state": "south", "target_state": "island", "transition_allowed": True, "step_distance": 2, "rate": 28.400649574079157}},
                {{"source_state": "south", "target_state": "north", "transition_allowed": True, "step_distance": 1, "rate": 1.2594240474115468}},
            ],
        }},
        "rows_name": "discrete-ancestral.tsv",
        "rows": [
            {{"node_id": 5, "node": "A|B|C|D", "state": "island", "posterior_probability": 0.300605046187356, "most_likely_state": "north", "max_posterior_probability": 0.398789910678949}},
            {{"node_id": 5, "node": "A|B|C|D", "state": "north", "posterior_probability": 0.398789910678949, "most_likely_state": "north", "max_posterior_probability": 0.398789910678949}},
            {{"node_id": 5, "node": "A|B|C|D", "state": "south", "posterior_probability": 0.300605043133695, "most_likely_state": "north", "max_posterior_probability": 0.398789910678949}},
            {{"node_id": 6, "node": "A|B", "state": "island", "posterior_probability": 0.024694666608751, "most_likely_state": "north", "max_posterior_probability": 0.950610668557755}},
            {{"node_id": 6, "node": "A|B", "state": "north", "posterior_probability": 0.950610668557755, "most_likely_state": "north", "max_posterior_probability": 0.950610668557755}},
            {{"node_id": 6, "node": "A|B", "state": "south", "posterior_probability": 0.024694664833495, "most_likely_state": "north", "max_posterior_probability": 0.950610668557755}},
            {{"node_id": 7, "node": "C|D", "state": "island", "posterior_probability": 0.411901434831837, "most_likely_state": "island", "max_posterior_probability": 0.411901434831837}},
            {{"node_id": 7, "node": "C|D", "state": "north", "posterior_probability": 0.176197130750738, "most_likely_state": "island", "max_posterior_probability": 0.411901434831837}},
            {{"node_id": 7, "node": "C|D", "state": "south", "posterior_probability": 0.411901434417425, "most_likely_state": "island", "max_posterior_probability": 0.411901434831837}},
        ],
    }},
    "ace-discrete-sym-pectinate-non-ultrametric": {{
        "summary": {{
            "trait": "region",
            "taxon_count": 4,
            "excluded_taxon_count": 0,
            "dropped_missing_taxa": [],
            "internal_node_count": 3,
            "model": "symmetric",
            "state_count": 3,
            "state_labels": ["island", "north", "south"],
            "log_likelihood": -3.102887409018013,
            "parameter_count": 3,
            "aic": 12.205774818036026,
            "overparameterized": False,
            "baseline_model": "equal-rates",
            "baseline_delta_aic": 3.7697498013258155,
            "preferred_model_by_aic": "equal-rates",
            "transition_rate_rows": [
                {{"source_state": "island", "target_state": "north", "transition_allowed": True, "step_distance": 1, "rate": 4.5399929762484854e-05}},
                {{"source_state": "island", "target_state": "south", "transition_allowed": True, "step_distance": 2, "rate": 8.96363692114651}},
                {{"source_state": "north", "target_state": "island", "transition_allowed": True, "step_distance": 1, "rate": 4.5399929762484854e-05}},
                {{"source_state": "north", "target_state": "south", "transition_allowed": True, "step_distance": 1, "rate": 5.616352384948936}},
                {{"source_state": "south", "target_state": "island", "transition_allowed": True, "step_distance": 2, "rate": 8.96363692114651}},
                {{"source_state": "south", "target_state": "north", "transition_allowed": True, "step_distance": 1, "rate": 5.616352384948936}},
            ],
        }},
        "rows_name": "discrete-ancestral.tsv",
        "rows": [
            {{"node_id": 5, "node": "A|B|C|D", "state": "island", "posterior_probability": 0.489718578919389, "most_likely_state": "island", "max_posterior_probability": 0.489718578919389}},
            {{"node_id": 5, "node": "A|B|C|D", "state": "north", "posterior_probability": 0.146079461872196, "most_likely_state": "island", "max_posterior_probability": 0.489718578919389}},
            {{"node_id": 5, "node": "A|B|C|D", "state": "south", "posterior_probability": 0.364201959208415, "most_likely_state": "island", "max_posterior_probability": 0.489718578919389}},
            {{"node_id": 6, "node": "A|B|C", "state": "island", "posterior_probability": 0.231890137639492, "most_likely_state": "south", "max_posterior_probability": 0.434011898424568}},
            {{"node_id": 6, "node": "A|B|C", "state": "north", "posterior_probability": 0.33409796393594, "most_likely_state": "south", "max_posterior_probability": 0.434011898424568}},
            {{"node_id": 6, "node": "A|B|C", "state": "south", "posterior_probability": 0.434011898424568, "most_likely_state": "south", "max_posterior_probability": 0.434011898424568}},
            {{"node_id": 7, "node": "A|B", "state": "island", "posterior_probability": 0.033878505842044, "most_likely_state": "north", "max_posterior_probability": 0.794971398041303}},
            {{"node_id": 7, "node": "A|B", "state": "north", "posterior_probability": 0.794971398041303, "most_likely_state": "north", "max_posterior_probability": 0.794971398041303}},
            {{"node_id": 7, "node": "A|B", "state": "south", "posterior_probability": 0.171150096116653, "most_likely_state": "north", "max_posterior_probability": 0.794971398041303}},
        ],
    }},
    "ace-discrete-sym-balanced-six-taxon": {{
        "summary": {{
            "trait": "region",
            "taxon_count": 6,
            "excluded_taxon_count": 0,
            "dropped_missing_taxa": [],
            "internal_node_count": 5,
            "model": "symmetric",
            "state_count": 3,
            "state_labels": ["island", "north", "south"],
            "log_likelihood": -3.4150348666333974,
            "parameter_count": 3,
            "aic": 12.830069733266795,
            "overparameterized": False,
            "baseline_model": "equal-rates",
            "baseline_delta_aic": 1.89034299315602,
            "preferred_model_by_aic": "equal-rates",
            "transition_rate_rows": [
                {{"source_state": "island", "target_state": "north", "transition_allowed": True, "step_distance": 1, "rate": 0.06263343139952941}},
                {{"source_state": "island", "target_state": "south", "transition_allowed": True, "step_distance": 2, "rate": 5.74955899176828}},
                {{"source_state": "north", "target_state": "island", "transition_allowed": True, "step_distance": 1, "rate": 0.06263155732731363}},
                {{"source_state": "north", "target_state": "south", "transition_allowed": True, "step_distance": 1, "rate": 0.06263843808777063}},
                {{"source_state": "south", "target_state": "island", "transition_allowed": True, "step_distance": 2, "rate": 5.74955899176828}},
                {{"source_state": "south", "target_state": "north", "transition_allowed": True, "step_distance": 1, "rate": 0.06263843808777063}},
            ],
        }},
        "rows_name": "discrete-ancestral.tsv",
        "rows": [
            {{"node_id": 7, "node": "A|B|C|D|E|F", "state": "island", "posterior_probability": 0.180590722637636, "most_likely_state": "north", "max_posterior_probability": 0.63881851370543}},
            {{"node_id": 7, "node": "A|B|C|D|E|F", "state": "north", "posterior_probability": 0.63881851370543, "most_likely_state": "north", "max_posterior_probability": 0.63881851370543}},
            {{"node_id": 7, "node": "A|B|C|D|E|F", "state": "south", "posterior_probability": 0.180590763656934, "most_likely_state": "north", "max_posterior_probability": 0.63881851370543}},
            {{"node_id": 8, "node": "A|B|C|D", "state": "island", "posterior_probability": 0.013305152752442, "most_likely_state": "north", "max_posterior_probability": 0.973389688184691}},
            {{"node_id": 8, "node": "A|B|C|D", "state": "north", "posterior_probability": 0.973389688184691, "most_likely_state": "north", "max_posterior_probability": 0.973389688184691}},
            {{"node_id": 8, "node": "A|B|C|D", "state": "south", "posterior_probability": 0.013305159062868, "most_likely_state": "north", "max_posterior_probability": 0.973389688184691}},
            {{"node_id": 9, "node": "A|B", "state": "island", "posterior_probability": 0.001115019136557, "most_likely_state": "north", "max_posterior_probability": 0.997769961090887}},
            {{"node_id": 9, "node": "A|B", "state": "north", "posterior_probability": 0.997769961090887, "most_likely_state": "north", "max_posterior_probability": 0.997769961090887}},
            {{"node_id": 9, "node": "A|B", "state": "south", "posterior_probability": 0.001115019772557, "most_likely_state": "north", "max_posterior_probability": 0.997769961090887}},
            {{"node_id": 10, "node": "C|D", "state": "island", "posterior_probability": 0.001115019136557, "most_likely_state": "north", "max_posterior_probability": 0.997769961090887}},
            {{"node_id": 10, "node": "C|D", "state": "north", "posterior_probability": 0.997769961090887, "most_likely_state": "north", "max_posterior_probability": 0.997769961090887}},
            {{"node_id": 10, "node": "C|D", "state": "south", "posterior_probability": 0.001115019772557, "most_likely_state": "north", "max_posterior_probability": 0.997769961090887}},
            {{"node_id": 11, "node": "E|F", "state": "island", "posterior_probability": 0.482874993955513, "most_likely_state": "south", "max_posterior_probability": 0.482875011977419}},
            {{"node_id": 11, "node": "E|F", "state": "north", "posterior_probability": 0.034249994067067, "most_likely_state": "south", "max_posterior_probability": 0.482875011977419}},
            {{"node_id": 11, "node": "E|F", "state": "south", "posterior_probability": 0.482875011977419, "most_likely_state": "south", "max_posterior_probability": 0.482875011977419}},
        ],
    }},
    "ace-discrete-sym-missing-values-pruned": {{
        "summary": {{
            "trait": "region",
            "taxon_count": 5,
            "excluded_taxon_count": 1,
            "dropped_missing_taxa": ["F"],
            "internal_node_count": 4,
            "model": "symmetric",
            "state_count": 3,
            "state_labels": ["island", "north", "south"],
            "log_likelihood": -3.704217157959887,
            "parameter_count": 3,
            "aic": 13.408434315919774,
            "overparameterized": False,
            "baseline_model": "equal-rates",
            "baseline_delta_aic": 3.8949246834485027,
            "preferred_model_by_aic": "equal-rates",
            "transition_rate_rows": [
                {{"source_state": "island", "target_state": "north", "transition_allowed": True, "step_distance": 1, "rate": 0.1219349630011045}},
                {{"source_state": "island", "target_state": "south", "transition_allowed": True, "step_distance": 2, "rate": 0.12193622593762252}},
                {{"source_state": "north", "target_state": "island", "transition_allowed": True, "step_distance": 1, "rate": 0.1219349630011045}},
                {{"source_state": "north", "target_state": "south", "transition_allowed": True, "step_distance": 1, "rate": 0.24764718612573017}},
                {{"source_state": "south", "target_state": "island", "transition_allowed": True, "step_distance": 2, "rate": 0.12193622593762252}},
                {{"source_state": "south", "target_state": "north", "transition_allowed": True, "step_distance": 1, "rate": 0.24764718612573017}},
            ],
        }},
        "rows_name": "discrete-ancestral.tsv",
        "rows": [
            {{"node_id": 6, "node": "A|B|C|D|E", "state": "island", "posterior_probability": 0.277744577273476, "most_likely_state": "south", "max_posterior_probability": 0.361127715527266}},
            {{"node_id": 6, "node": "A|B|C|D|E", "state": "north", "posterior_probability": 0.361127707199258, "most_likely_state": "south", "max_posterior_probability": 0.361127715527266}},
            {{"node_id": 6, "node": "A|B|C|D|E", "state": "south", "posterior_probability": 0.361127715527266, "most_likely_state": "south", "max_posterior_probability": 0.361127715527266}},
            {{"node_id": 7, "node": "A|B|C|D", "state": "island", "posterior_probability": 0.08195148273393, "most_likely_state": "south", "max_posterior_probability": 0.459024261477649}},
            {{"node_id": 7, "node": "A|B|C|D", "state": "north", "posterior_probability": 0.459024255788421, "most_likely_state": "south", "max_posterior_probability": 0.459024261477649}},
            {{"node_id": 7, "node": "A|B|C|D", "state": "south", "posterior_probability": 0.459024261477649, "most_likely_state": "south", "max_posterior_probability": 0.459024261477649}},
            {{"node_id": 8, "node": "A|B", "state": "island", "posterior_probability": 0.016055890641513, "most_likely_state": "north", "max_posterior_probability": 0.881986522646983}},
            {{"node_id": 8, "node": "A|B", "state": "north", "posterior_probability": 0.881986522646983, "most_likely_state": "north", "max_posterior_probability": 0.881986522646983}},
            {{"node_id": 8, "node": "A|B", "state": "south", "posterior_probability": 0.101957586711504, "most_likely_state": "north", "max_posterior_probability": 0.881986522646983}},
            {{"node_id": 9, "node": "C|D", "state": "island", "posterior_probability": 0.016055892347348, "most_likely_state": "south", "max_posterior_probability": 0.881986519985111}},
            {{"node_id": 9, "node": "C|D", "state": "north", "posterior_probability": 0.101957587667541, "most_likely_state": "south", "max_posterior_probability": 0.881986519985111}},
            {{"node_id": 9, "node": "C|D", "state": "south", "posterior_probability": 0.881986519985111, "most_likely_state": "south", "max_posterior_probability": 0.881986519985111}},
        ],
    }},
    "ace-discrete-ard-binary-balanced-rooted-ultrametric": {{
        "summary": {{
            "trait": "habitat",
            "taxon_count": 4,
            "excluded_taxon_count": 0,
            "dropped_missing_taxa": [],
            "internal_node_count": 3,
            "model": "all-rates-different",
            "state_count": 2,
            "state_labels": ["forest", "tundra"],
            "log_likelihood": -1.7423030956462013,
            "parameter_count": 2,
            "aic": 7.484606191292403,
            "overparameterized": False,
            "baseline_model": "equal-rates",
            "baseline_delta_aic": 1.9693267371568686,
            "preferred_model_by_aic": "equal-rates",
            "transition_rate_rows": [
                {{"source_state": "forest", "target_state": "tundra", "transition_allowed": True, "step_distance": 1, "rate": 2.2937761295508943}},
                {{"source_state": "tundra", "target_state": "forest", "transition_allowed": True, "step_distance": 1, "rate": 1.6512202819475246}},
            ],
        }},
        "rows_name": "discrete-ancestral.tsv",
        "rows": [
            {{"node_id": 5, "node": "A|B|C|D", "state": "forest", "posterior_probability": 0.504463640200633, "most_likely_state": "forest", "max_posterior_probability": 0.504463640200633}},
            {{"node_id": 5, "node": "A|B|C|D", "state": "tundra", "posterior_probability": 0.495536359799367, "most_likely_state": "forest", "max_posterior_probability": 0.504463640200633}},
            {{"node_id": 6, "node": "A|B", "state": "forest", "posterior_probability": 0.955807851977392, "most_likely_state": "forest", "max_posterior_probability": 0.955807851977392}},
            {{"node_id": 6, "node": "A|B", "state": "tundra", "posterior_probability": 0.044192148022608, "most_likely_state": "forest", "max_posterior_probability": 0.955807851977392}},
            {{"node_id": 7, "node": "C|D", "state": "forest", "posterior_probability": 0.242967749034647, "most_likely_state": "tundra", "max_posterior_probability": 0.757032250965353}},
            {{"node_id": 7, "node": "C|D", "state": "tundra", "posterior_probability": 0.757032250965353, "most_likely_state": "tundra", "max_posterior_probability": 0.757032250965353}},
        ],
    }},
    "ace-discrete-ard-multistate-pectinate-non-ultrametric": {{
        "summary": {{
            "trait": "region",
            "taxon_count": 4,
            "excluded_taxon_count": 0,
            "dropped_missing_taxa": [],
            "internal_node_count": 3,
            "model": "all-rates-different",
            "state_count": 3,
            "state_labels": ["island", "north", "south"],
            "log_likelihood": -2.5693089407614753,
            "parameter_count": 6,
            "aic": 17.13861788152295,
            "overparameterized": True,
            "baseline_model": "equal-rates",
            "baseline_delta_aic": 8.702592865212734,
            "preferred_model_by_aic": "equal-rates",
            "transition_rate_rows": [
                {{"source_state": "island", "target_state": "north", "transition_allowed": True, "step_distance": 1, "rate": 4.5399929762484854e-05}},
                {{"source_state": "island", "target_state": "south", "transition_allowed": True, "step_distance": 2, "rate": 7.0560196734205665}},
                {{"source_state": "north", "target_state": "island", "transition_allowed": True, "step_distance": 1, "rate": 4.5399929762484854e-05}},
                {{"source_state": "north", "target_state": "south", "transition_allowed": True, "step_distance": 1, "rate": 4.5399929762484854e-05}},
                {{"source_state": "south", "target_state": "island", "transition_allowed": True, "step_distance": 2, "rate": 0.6517960181064918}},
                {{"source_state": "south", "target_state": "north", "transition_allowed": True, "step_distance": 1, "rate": 6.141245382921997}},
            ],
        }},
        "rows_name": "discrete-ancestral.tsv",
        "rows": [
            {{"node_id": 5, "node": "A|B|C|D", "state": "island", "posterior_probability": 0.935019548662187, "most_likely_state": "island", "max_posterior_probability": 0.935019548662187}},
            {{"node_id": 5, "node": "A|B|C|D", "state": "north", "posterior_probability": 2.5778e-10, "most_likely_state": "island", "max_posterior_probability": 0.935019548662187}},
            {{"node_id": 5, "node": "A|B|C|D", "state": "south", "posterior_probability": 0.064980451080033, "most_likely_state": "island", "max_posterior_probability": 0.935019548662187}},
            {{"node_id": 6, "node": "A|B|C", "state": "island", "posterior_probability": 0.274103686147979, "most_likely_state": "south", "max_posterior_probability": 0.725891882827243}},
            {{"node_id": 6, "node": "A|B|C", "state": "north", "posterior_probability": 4.431024779e-06, "most_likely_state": "south", "max_posterior_probability": 0.725891882827243}},
            {{"node_id": 6, "node": "A|B|C", "state": "south", "posterior_probability": 0.725891882827243, "most_likely_state": "south", "max_posterior_probability": 0.725891882827243}},
            {{"node_id": 7, "node": "A|B", "state": "island", "posterior_probability": 0.022934190650633, "most_likely_state": "north", "max_posterior_probability": 0.719609497603507}},
            {{"node_id": 7, "node": "A|B", "state": "north", "posterior_probability": 0.719609497603507, "most_likely_state": "north", "max_posterior_probability": 0.719609497603507}},
            {{"node_id": 7, "node": "A|B", "state": "south", "posterior_probability": 0.25745631174586, "most_likely_state": "north", "max_posterior_probability": 0.719609497603507}},
        ],
    }},
    "ace-discrete-ard-balanced-six-taxon": {{
        "summary": {{
            "trait": "region",
            "taxon_count": 6,
            "excluded_taxon_count": 0,
            "dropped_missing_taxa": [],
            "internal_node_count": 5,
            "model": "all-rates-different",
            "state_count": 3,
            "state_labels": ["island", "north", "south"],
            "log_likelihood": -2.9812150238347677,
            "parameter_count": 6,
            "aic": 17.962430047669535,
            "overparameterized": True,
            "baseline_model": "equal-rates",
            "baseline_delta_aic": 7.022703307958688,
            "preferred_model_by_aic": "equal-rates",
            "transition_rate_rows": [
                {{"source_state": "island", "target_state": "north", "transition_allowed": True, "step_distance": 1, "rate": 0.542389523475714}},
                {{"source_state": "island", "target_state": "south", "transition_allowed": True, "step_distance": 2, "rate": 1.4963870483597863}},
                {{"source_state": "north", "target_state": "island", "transition_allowed": True, "step_distance": 1, "rate": 4.5399929762484854e-05}},
                {{"source_state": "north", "target_state": "south", "transition_allowed": True, "step_distance": 1, "rate": 0.011915663508013425}},
                {{"source_state": "south", "target_state": "island", "transition_allowed": True, "step_distance": 2, "rate": 1.575766742871299}},
                {{"source_state": "south", "target_state": "north", "transition_allowed": True, "step_distance": 1, "rate": 4.5399929762484854e-05}},
            ],
        }},
        "rows_name": "discrete-ancestral.tsv",
        "rows": [
            {{"node_id": 7, "node": "A|B|C|D|E|F", "state": "island", "posterior_probability": 0.532202345475697, "most_likely_state": "island", "max_posterior_probability": 0.532202345475697}},
            {{"node_id": 7, "node": "A|B|C|D|E|F", "state": "north", "posterior_probability": 0.056349835943089, "most_likely_state": "island", "max_posterior_probability": 0.532202345475697}},
            {{"node_id": 7, "node": "A|B|C|D|E|F", "state": "south", "posterior_probability": 0.411447818581214, "most_likely_state": "island", "max_posterior_probability": 0.532202345475697}},
            {{"node_id": 8, "node": "A|B|C|D", "state": "island", "posterior_probability": 0.134323458951379, "most_likely_state": "north", "max_posterior_probability": 0.806889166517314}},
            {{"node_id": 8, "node": "A|B|C|D", "state": "north", "posterior_probability": 0.806889166517314, "most_likely_state": "north", "max_posterior_probability": 0.806889166517314}},
            {{"node_id": 8, "node": "A|B|C|D", "state": "south", "posterior_probability": 0.058787374531307, "most_likely_state": "north", "max_posterior_probability": 0.806889166517314}},
            {{"node_id": 9, "node": "A|B", "state": "island", "posterior_probability": 0.064094056411562, "most_likely_state": "north", "max_posterior_probability": 0.914130423464103}},
            {{"node_id": 9, "node": "A|B", "state": "north", "posterior_probability": 0.914130423464103, "most_likely_state": "north", "max_posterior_probability": 0.914130423464103}},
            {{"node_id": 9, "node": "A|B", "state": "south", "posterior_probability": 0.021775520124336, "most_likely_state": "north", "max_posterior_probability": 0.914130423464103}},
            {{"node_id": 10, "node": "C|D", "state": "island", "posterior_probability": 0.064094056411562, "most_likely_state": "north", "max_posterior_probability": 0.914130423464103}},
            {{"node_id": 10, "node": "C|D", "state": "north", "posterior_probability": 0.914130423464103, "most_likely_state": "north", "max_posterior_probability": 0.914130423464103}},
            {{"node_id": 10, "node": "C|D", "state": "south", "posterior_probability": 0.021775520124336, "most_likely_state": "north", "max_posterior_probability": 0.914130423464103}},
            {{"node_id": 11, "node": "E|F", "state": "island", "posterior_probability": 0.389284284230715, "most_likely_state": "south", "max_posterior_probability": 0.610579721722943}},
            {{"node_id": 11, "node": "E|F", "state": "north", "posterior_probability": 0.000135994046341, "most_likely_state": "south", "max_posterior_probability": 0.610579721722943}},
            {{"node_id": 11, "node": "E|F", "state": "south", "posterior_probability": 0.610579721722943, "most_likely_state": "south", "max_posterior_probability": 0.610579721722943}},
        ],
    }},
    "ace-discrete-ard-missing-values-pruned": {{
        "summary": {{
            "trait": "region",
            "taxon_count": 5,
            "excluded_taxon_count": 1,
            "dropped_missing_taxa": ["F"],
            "internal_node_count": 4,
            "model": "all-rates-different",
            "state_count": 3,
            "state_labels": ["island", "north", "south"],
            "log_likelihood": -3.6947644572775626,
            "parameter_count": 6,
            "aic": 19.389528914555125,
            "overparameterized": True,
            "baseline_model": "equal-rates",
            "baseline_delta_aic": 9.87601928248385,
            "preferred_model_by_aic": "equal-rates",
            "transition_rate_rows": [
                {{"source_state": "island", "target_state": "north", "transition_allowed": True, "step_distance": 1, "rate": 0.16373825344651594}},
                {{"source_state": "island", "target_state": "south", "transition_allowed": True, "step_distance": 2, "rate": 0.16374043672392322}},
                {{"source_state": "north", "target_state": "island", "transition_allowed": True, "step_distance": 1, "rate": 0.11881834258757133}},
                {{"source_state": "north", "target_state": "south", "transition_allowed": True, "step_distance": 1, "rate": 0.23789662659738825}},
                {{"source_state": "south", "target_state": "island", "transition_allowed": True, "step_distance": 2, "rate": 0.11881838590700566}},
                {{"source_state": "south", "target_state": "north", "transition_allowed": True, "step_distance": 1, "rate": 0.237896582835569}},
            ],
        }},
        "rows_name": "discrete-ancestral.tsv",
        "rows": [
            {{"node_id": 6, "node": "A|B|C|D|E", "state": "island", "posterior_probability": 0.327896912638208, "most_likely_state": "north", "max_posterior_probability": 0.336051604924961}},
            {{"node_id": 6, "node": "A|B|C|D|E", "state": "north", "posterior_probability": 0.336051604924961, "most_likely_state": "north", "max_posterior_probability": 0.336051604924961}},
            {{"node_id": 6, "node": "A|B|C|D|E", "state": "south", "posterior_probability": 0.336051482436831, "most_likely_state": "north", "max_posterior_probability": 0.336051604924961}},
            {{"node_id": 7, "node": "A|B|C|D", "state": "island", "posterior_probability": 0.128573641653441, "most_likely_state": "south", "max_posterior_probability": 0.435713690096739}},
            {{"node_id": 7, "node": "A|B|C|D", "state": "north", "posterior_probability": 0.435712668249821, "most_likely_state": "south", "max_posterior_probability": 0.435713690096739}},
            {{"node_id": 7, "node": "A|B|C|D", "state": "south", "posterior_probability": 0.435713690096739, "most_likely_state": "south", "max_posterior_probability": 0.435713690096739}},
            {{"node_id": 8, "node": "A|B", "state": "island", "posterior_probability": 0.031645589136504, "most_likely_state": "north", "max_posterior_probability": 0.874529916631945}},
            {{"node_id": 8, "node": "A|B", "state": "north", "posterior_probability": 0.874529916631945, "most_likely_state": "north", "max_posterior_probability": 0.874529916631945}},
            {{"node_id": 8, "node": "A|B", "state": "south", "posterior_probability": 0.093824494231551, "most_likely_state": "north", "max_posterior_probability": 0.874529916631945}},
            {{"node_id": 9, "node": "C|D", "state": "island", "posterior_probability": 0.031646086796991, "most_likely_state": "south", "max_posterior_probability": 0.874529753151947}},
            {{"node_id": 9, "node": "C|D", "state": "north", "posterior_probability": 0.093824160051062, "most_likely_state": "south", "max_posterior_probability": 0.874529753151947}},
            {{"node_id": 9, "node": "C|D", "state": "south", "posterior_probability": 0.874529753151947, "most_likely_state": "south", "max_posterior_probability": 0.874529753151947}},
        ],
    }},
    "pic-balanced-rooted-ultrametric": {{
        "summary": {{
            "trait": "response",
            "taxon_count": 4,
            "contrast_count": 3,
            "tree_is_ultrametric": True,
            "minimum_root_to_tip_depth": 0.3,
            "maximum_root_to_tip_depth": 0.3,
        }},
        "rows_name": "independent-contrasts.tsv",
        "rows": [
            {{"node_id": 5, "node": "A|B|C|D", "left_taxa": "A|B", "right_taxa": "C|D", "contrast": -1.49071198499986, "expected_variance": 0.45}},
            {{"node_id": 6, "node": "A|B", "left_taxa": "A", "right_taxa": "B", "contrast": -3.3541019662496847, "expected_variance": 0.2}},
            {{"node_id": 7, "node": "C|D", "left_taxa": "C", "right_taxa": "D", "contrast": -2.3717082451262845, "expected_variance": 0.4}},
        ],
    }},
    "pic-pectinate-non-ultrametric": {{
        "summary": {{
            "trait": "response",
            "taxon_count": 4,
            "contrast_count": 3,
            "tree_is_ultrametric": False,
            "minimum_root_to_tip_depth": 0.1,
            "maximum_root_to_tip_depth": 0.3,
        }},
        "rows_name": "independent-contrasts.tsv",
        "rows": [
            {{"node_id": 5, "node": "A|B|C|D", "left_taxa": "A|B|C", "right_taxa": "D", "contrast": -3.137858162210944, "expected_variance": 0.26}},
            {{"node_id": 6, "node": "A|B|C", "left_taxa": "A|B", "right_taxa": "C", "contrast": -0.5, "expected_variance": 0.25}},
            {{"node_id": 7, "node": "A|B", "left_taxa": "A", "right_taxa": "B", "contrast": -3.3541019662496847, "expected_variance": 0.2}},
        ],
    }},
    "pic-balanced-six-taxon": {{
        "summary": {{
            "trait": "response_growth",
            "taxon_count": 6,
            "contrast_count": 5,
            "tree_is_ultrametric": True,
            "minimum_root_to_tip_depth": 3.0,
            "maximum_root_to_tip_depth": 3.0,
        }},
        "rows_name": "independent-contrasts.tsv",
        "rows": [
            {{"node_id": 7, "node": "A|B|C|D|E|F", "left_taxa": "A|B|C|D", "right_taxa": "E|F", "contrast": -0.7276068751089989, "expected_variance": 4.25}},
            {{"node_id": 8, "node": "A|B|C|D", "left_taxa": "A|B", "right_taxa": "C|D", "contrast": -0.5773502691896258, "expected_variance": 3.0}},
            {{"node_id": 9, "node": "A|B", "left_taxa": "A", "right_taxa": "B", "contrast": -0.35355339059327373, "expected_variance": 2.0}},
            {{"node_id": 10, "node": "C|D", "left_taxa": "C", "right_taxa": "D", "contrast": -0.35355339059327373, "expected_variance": 2.0}},
            {{"node_id": 11, "node": "E|F", "left_taxa": "E", "right_taxa": "F", "contrast": -0.35355339059327373, "expected_variance": 2.0}},
        ],
    }},
    "dna-translation-valid-frame": {{
        "summary": {{
            "sequence_count": 3,
            "translated_length": 3,
            "stop_codon_count": 0,
            "dropped_trailing_nucleotide_count": 0,
            "warning_count": 0,
            "warnings": [],
        }},
        "rows_name": "translation.tsv",
        "rows": [
            {{"identifier": "valid_a", "amino_acid_sequence": "MEL"}},
            {{"identifier": "valid_b", "amino_acid_sequence": "MKM"}},
            {{"identifier": "valid_c", "amino_acid_sequence": "MTW"}},
        ],
    }},
    "dna-translation-ambiguous-codon": {{
        "summary": {{
            "sequence_count": 1,
            "translated_length": 3,
            "stop_codon_count": 0,
            "dropped_trailing_nucleotide_count": 0,
            "warning_count": 0,
            "warnings": [],
        }},
        "rows_name": "translation.tsv",
        "rows": [
            {{"identifier": "ambiguous_codon", "amino_acid_sequence": "MXG"}},
        ],
    }},
    "dna-translation-internal-stop": {{
        "summary": {{
            "sequence_count": 1,
            "translated_length": 3,
            "stop_codon_count": 1,
            "dropped_trailing_nucleotide_count": 0,
            "warning_count": 0,
            "warnings": [],
        }},
        "rows_name": "translation.tsv",
        "rows": [
            {{"identifier": "internal_stop", "amino_acid_sequence": "M*W"}},
        ],
    }},
    "dna-translation-terminal-stop": {{
        "summary": {{
            "sequence_count": 1,
            "translated_length": 3,
            "stop_codon_count": 1,
            "dropped_trailing_nucleotide_count": 0,
            "warning_count": 0,
            "warnings": [],
        }},
        "rows_name": "translation.tsv",
        "rows": [
            {{"identifier": "terminal_stop", "amino_acid_sequence": "ME*"}},
        ],
    }},
    "dna-translation-frame-error-truncation": {{
        "summary": {{
            "sequence_count": 1,
            "translated_length": 2,
            "stop_codon_count": 0,
            "dropped_trailing_nucleotide_count": 2,
            "warning_count": 1,
            "warnings": ["sequence length not a multiple of 3: 2 nucleotides dropped"],
        }},
        "rows_name": "translation.tsv",
        "rows": [
            {{"identifier": "frame_error", "amino_acid_sequence": "ME"}},
        ],
    }},
    "dna-translation-vertebrate-mitochondrial": {{
        "summary": {{
            "sequence_count": 1,
            "translated_length": 3,
            "stop_codon_count": 0,
            "dropped_trailing_nucleotide_count": 0,
            "warning_count": 0,
            "warnings": [],
        }},
        "rows_name": "translation.tsv",
        "rows": [
            {{"identifier": "mitochondrial_triplet", "amino_acid_sequence": "MWG"}},
        ],
    }},
    "rtree-rooted-six-taxon-uniform": {{
        "summary": {{
            "simulation_model": "random-tree",
            "reference_function": "ape::rtree",
            "tree_count": 64,
            "tip_count": 6,
            "seed": 17,
            "branch_length_model": "uniform",
            "population_size": None,
            "rooted": True,
            "binary": True,
            "pooled_branch_count": 640,
            "envelope_metric_count": 6,
        }},
        "rows_name": "simulation-envelope.tsv",
        "rows": [
            {{"metric": "tree_height_branch_length", "sample_scope": "tree", "observation_count": 64, "mean": 2.200019969538129, "standard_deviation": 0.568199526954071, "minimum": 1.180148204554461, "median": 2.172595851103772, "maximum": 3.926725072850143}},
            {{"metric": "total_branch_length", "sample_scope": "tree", "observation_count": 64, "mean": 5.011133250674653, "standard_deviation": 0.8717600443239, "minimum": 3.130211330959545, "median": 5.101709492928771, "maximum": 6.815146136377609}},
            {{"metric": "branch_length", "sample_scope": "edge", "observation_count": 640, "mean": 0.501113325067465, "standard_deviation": 0.29599061860717, "minimum": 0.000746097731662, "median": 0.511955271949145, "maximum": 0.999895000159003}},
            {{"metric": "cherry_count", "sample_scope": "tree", "observation_count": 64, "mean": 2.0625, "standard_deviation": 0.496078370824611, "minimum": 1.0, "median": 2.0, "maximum": 3.0}},
            {{"metric": "sackin_imbalance_index", "sample_scope": "tree", "observation_count": 64, "mean": 17.375, "standard_deviation": 1.268611445636527, "minimum": 16.0, "median": 17.0, "maximum": 20.0}},
            {{"metric": "normalized_colless_imbalance", "sample_scope": "tree", "observation_count": 64, "mean": 0.490625, "standard_deviation": 0.242202620495733, "minimum": 0.2, "median": 0.5, "maximum": 1.0}},
        ],
    }},
    "rtree-rooted-twelve-taxon-uniform": {{
        "summary": {{
            "simulation_model": "random-tree",
            "reference_function": "ape::rtree",
            "tree_count": 128,
            "tip_count": 12,
            "seed": 29,
            "branch_length_model": "uniform",
            "population_size": None,
            "rooted": True,
            "binary": True,
            "pooled_branch_count": 2816,
            "envelope_metric_count": 6,
        }},
        "rows_name": "simulation-envelope.tsv",
        "rows": [
            {{"metric": "tree_height_branch_length", "sample_scope": "tree", "observation_count": 128, "mean": 3.293935035043397, "standard_deviation": 0.680477733621439, "minimum": 1.859634394572088, "median": 3.194875305517843, "maximum": 5.170210547581812}},
            {{"metric": "total_branch_length", "sample_scope": "tree", "observation_count": 128, "mean": 10.81044224466413, "standard_deviation": 1.392459143501874, "minimum": 7.638574808773966, "median": 10.691554469854701, "maximum": 14.658633860473989}},
            {{"metric": "branch_length", "sample_scope": "edge", "observation_count": 2816, "mean": 0.491383738393823, "standard_deviation": 0.284450712368707, "minimum": 4.7054522384e-05, "median": 0.483198671141144, "maximum": 0.999667835574992}},
            {{"metric": "cherry_count", "sample_scope": "tree", "observation_count": 128, "mean": 3.953125, "standard_deviation": 0.727360113269211, "minimum": 2.0, "median": 4.0, "maximum": 6.0}},
            {{"metric": "sackin_imbalance_index", "sample_scope": "tree", "observation_count": 128, "mean": 50.0078125, "standard_deviation": 4.400809751039432, "minimum": 44.0, "median": 49.0, "maximum": 66.0}},
            {{"metric": "normalized_colless_imbalance", "sample_scope": "tree", "observation_count": 128, "mean": 0.305539772727273, "standard_deviation": 0.121501849154168, "minimum": 0.072727272727273, "median": 0.290909090909091, "maximum": 0.727272727272727}},
        ],
    }},
    "rcoal-rooted-six-taxon": {{
        "summary": {{
            "simulation_model": "coalescent",
            "reference_function": "ape::rcoal",
            "tree_count": 64,
            "tip_count": 6,
            "seed": 17,
            "branch_length_model": "coalescent-waiting-times",
            "population_size": 1.0,
            "rooted": True,
            "binary": True,
            "pooled_branch_count": 640,
            "envelope_metric_count": 6,
        }},
        "rows_name": "simulation-envelope.tsv",
        "rows": [
            {{"metric": "tree_height_branch_length", "sample_scope": "tree", "observation_count": 64, "mean": 1.774295612977311, "standard_deviation": 1.094448031695996, "minimum": 0.315287735719744, "median": 1.44689642194404, "maximum": 5.657108945358354}},
            {{"metric": "total_branch_length", "sample_scope": "tree", "observation_count": 64, "mean": 4.836969469523078, "standard_deviation": 2.502517797978043, "minimum": 1.082371620932966, "median": 4.114422515274727, "maximum": 12.314041485523623}},
            {{"metric": "branch_length", "sample_scope": "edge", "observation_count": 640, "mean": 0.483696946952307, "standard_deviation": 0.721376576530333, "minimum": 0.000823616521307, "median": 0.216835144077217, "maximum": 5.359401701040186}},
            {{"metric": "cherry_count", "sample_scope": "tree", "observation_count": 64, "mean": 1.921875, "standard_deviation": 0.620097963530763, "minimum": 1.0, "median": 2.0, "maximum": 3.0}},
            {{"metric": "sackin_imbalance_index", "sample_scope": "tree", "observation_count": 64, "mean": 17.671875, "standard_deviation": 1.541657868781203, "minimum": 16.0, "median": 17.0, "maximum": 20.0}},
            {{"metric": "normalized_colless_imbalance", "sample_scope": "tree", "observation_count": 64, "mean": 0.5484375, "standard_deviation": 0.299474637646913, "minimum": 0.2, "median": 0.5, "maximum": 1.0}},
        ],
    }},
    "rcoal-rooted-twelve-taxon": {{
        "summary": {{
            "simulation_model": "coalescent",
            "reference_function": "ape::rcoal",
            "tree_count": 128,
            "tip_count": 12,
            "seed": 29,
            "branch_length_model": "coalescent-waiting-times",
            "population_size": 1.0,
            "rooted": True,
            "binary": True,
            "pooled_branch_count": 2816,
            "envelope_metric_count": 6,
        }},
        "rows_name": "simulation-envelope.tsv",
        "rows": [
            {{"metric": "tree_height_branch_length", "sample_scope": "tree", "observation_count": 128, "mean": 1.848760485535123, "standard_deviation": 1.116849418772524, "minimum": 0.439983410386754, "median": 1.595097706740242, "maximum": 6.195802225859135}},
            {{"metric": "total_branch_length", "sample_scope": "tree", "observation_count": 128, "mean": 6.164387781665651, "standard_deviation": 2.586531949750664, "minimum": 1.826587536526831, "median": 5.652029095234272, "maximum": 15.697139126436669}},
            {{"metric": "branch_length", "sample_scope": "edge", "observation_count": 2816, "mean": 0.280199444621166, "standard_deviation": 0.519106863790402, "minimum": 2.2898700613e-05, "median": 0.109624730396113, "maximum": 6.147730825333207}},
            {{"metric": "cherry_count", "sample_scope": "tree", "observation_count": 128, "mean": 4.0859375, "standard_deviation": 0.70740882528687, "minimum": 3.0, "median": 4.0, "maximum": 6.0}},
            {{"metric": "sackin_imbalance_index", "sample_scope": "tree", "observation_count": 128, "mean": 50.640625, "standard_deviation": 4.582122827835915, "minimum": 44.0, "median": 50.0, "maximum": 66.0}},
            {{"metric": "normalized_colless_imbalance", "sample_scope": "tree", "observation_count": 128, "mean": 0.319602272727273, "standard_deviation": 0.129378069193442, "minimum": 0.072727272727273, "median": 0.318181818181818, "maximum": 0.727272727272727}},
        ],
    }},
}}

TREE_CASES = {{
    "neighbor-joining-analytical-three-taxon": {{
        "summary": {{
            "tree_count": 1,
            "tip_count": 3,
            "internal_node_count": 1,
            "edge_count": 3,
            "rooted": False,
            "tip_labels": ["A", "B", "C"],
            "branch_length_count": 3,
        }},
        "newick": "(A:0,B:0.125,C:0.5);",
    }},
    "neighbor-joining-ultrametric-four-taxon": {{
        "summary": {{
            "tree_count": 1,
            "tip_count": 4,
            "internal_node_count": 2,
            "edge_count": 5,
            "rooted": False,
            "tip_labels": ["A", "B", "C", "D"],
            "branch_length_count": 5,
        }},
        "newick": "((A:1,B:1):4,C:1,D:1);",
    }},
    "neighbor-joining-nonultrametric-four-taxon": {{
        "summary": {{
            "tree_count": 1,
            "tip_count": 4,
            "internal_node_count": 2,
            "edge_count": 5,
            "rooted": False,
            "tip_labels": ["A", "B", "C", "D"],
            "branch_length_count": 5,
        }},
        "newick": "((A:1,B:1):4.5,C:0.5,D:1.5);",
    }},
}}

ERROR_CASES = {{
    "dna-raw-distance-unequal-length-invalid": {{
        "error_type": "DnaDistanceError",
        "message": "DNA sequences in list not of the same length.",
    }},
    "dna-jc69-distance-unequal-length-invalid": {{
        "error_type": "DnaDistanceError",
        "message": "DNA sequences in list not of the same length.",
    }},
    "dna-k80-distance-unequal-length-invalid": {{
        "error_type": "DnaDistanceError",
        "message": "DNA sequences in list not of the same length.",
    }},
    "dna-f81-distance-unequal-length-invalid": {{
        "error_type": "DnaDistanceError",
        "message": "DNA sequences in list not of the same length.",
    }},
    "dna-tn93-distance-unequal-length-invalid": {{
        "error_type": "DnaDistanceError",
        "message": "DNA sequences in list not of the same length.",
    }},
}}

SUMMARY_OVERRIDES = {summary_payload}
NORMALIZED_TREE_OVERRIDES = {normalized_tree_payload}
APE_AVAILABLE = {str(ape_available)}

def write_json(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\\n", encoding="utf-8")

def write_tsv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\\t")
        writer.writeheader()
        writer.writerows(rows)

def normalize_label(value):
    return "" if value is None else str(value)

def parse_support_label(value):
    text = normalize_label(value)
    if not text:
        return ""
    if "/" in text:
        try:
            return float(text.split("/")[-1])
        except ValueError:
            return ""
    try:
        return float(text)
    except ValueError:
        return ""

def descendant_taxa(clade):
    return sorted(terminal.name for terminal in clade.get_terminals() if terminal.name)

def clade_rows(tree, tree_index_value):
    rows = []
    for clade in tree.find_clades(order="preorder"):
        taxa = descendant_taxa(clade)
        node_kind = "tip" if clade.is_terminal() else "internal"
        if clade == tree.root:
            node_kind = "root"
        if clade.is_terminal():
            node_label = normalize_label(clade.name)
        else:
            node_label = normalize_label(clade.name if clade.name is not None else getattr(clade, "confidence", None))
        rows.append(
            {{
                "tree_index": tree_index_value,
                "node_kind": node_kind,
                "clade_id": "|".join(taxa),
                "node_label": node_label,
                "taxon_count": len(taxa),
                "taxa": "|".join(taxa),
                "support": parse_support_label(node_label),
                "branch_length": "" if clade.branch_length is None else clade.branch_length,
            }}
        )
    node_order = {{"root": 0, "internal": 1, "tip": 2}}
    return sorted(
        rows,
        key=lambda row: (
            0 if row["tree_index"] == "" else int(row["tree_index"]),
            node_order.get(row["node_kind"], 9),
            row["clade_id"],
            row["node_label"],
        ),
    )

def is_rooted_tree(tree):
    return len(getattr(tree.root, "clades", [])) == 2

def iter_internal_clades_preorder(clade):
    yield clade
    for child in getattr(clade, "clades", []):
        if child.is_terminal():
            continue
        yield from iter_internal_clades_preorder(child)

def node_depth_lookup(tree):
    lookup = {{id(tree.root): 0.0}}
    def walk(clade):
        base_depth = lookup[id(clade)]
        for child in getattr(clade, "clades", []):
            lookup[id(child)] = base_depth + float(child.branch_length or 0.0)
            walk(child)
    walk(tree.root)
    return lookup

def node_depth_rows(tree):
    depths = node_depth_lookup(tree)
    rows = []
    tip_clades = list(tree.get_terminals())
    internal_clades = list(iter_internal_clades_preorder(tree.root))
    for node_id, clade in enumerate(tip_clades, start=1):
        rows.append(
            {{
                "node_id": node_id,
                "node_kind": "tip",
                "node_label": normalize_label(clade.name),
                "descendant_taxa": "|".join(descendant_taxa(clade)),
                "branch_length_depth": depths[id(clade)],
                "branch_length": ""
                if clade.branch_length is None
                else clade.branch_length,
            }}
        )
    for offset, clade in enumerate(internal_clades, start=1):
        node_id = len(tip_clades) + offset
        node_label = normalize_label(
            clade.name if clade.name is not None else getattr(clade, "confidence", None)
        )
        rows.append(
            {{
                "node_id": node_id,
                "node_kind": "root" if clade is tree.root else "internal",
                "node_label": node_label,
                "descendant_taxa": "|".join(descendant_taxa(clade)),
                "branch_length_depth": depths[id(clade)],
                "branch_length": ""
                if clade.branch_length is None
                else clade.branch_length,
            }}
        )
    return rows

def branching_time_rows(tree):
    node_depth_table = node_depth_rows(tree)
    tip_depths = [
        row["branch_length_depth"]
        for row in node_depth_table
        if row["node_kind"] == "tip"
    ]
    root_age = max(tip_depths)
    return [
        {{
            "node_id": row["node_id"],
            "node_kind": row["node_kind"],
            "node_label": row["node_label"],
            "descendant_taxa": row["descendant_taxa"],
            "node_depth": row["branch_length_depth"],
            "branching_time": root_age - row["branch_length_depth"],
        }}
        for row in node_depth_table
        if row["node_kind"] != "tip"
    ]

def ultrametric_rows(tree):
    node_depth_table = node_depth_rows(tree)
    return [
        {{
            "node_id": row["node_id"],
            "tip_label": row["node_label"],
            "root_to_tip_depth": row["branch_length_depth"],
            "deviation_from_mean_depth": None,
            "deviation_from_min_depth": None,
            "deviation_from_max_depth": None,
            "is_offending_taxon": False,
        }}
        for row in node_depth_table
        if row["node_kind"] == "tip"
    ]

def signature_id(taxa):
    return "|".join(sorted(set(taxa)))

def rooted_topology_signatures(tree):
    total_tip_count = len(list(tree.get_terminals()))
    signatures = {{}}
    for clade in tree.find_clades(order="preorder"):
        if clade is tree.root or clade.is_terminal():
            continue
        taxa = descendant_taxa(clade)
        if len(taxa) <= 1 or len(taxa) >= total_tip_count:
            continue
        signatures[signature_id(taxa)] = taxa
    return signatures

def canonical_unrooted_signature(taxa, all_taxa):
    selected = sorted(set(taxa))
    complement = sorted(taxon for taxon in all_taxa if taxon not in selected)
    if len(selected) < len(complement):
        return selected
    if len(complement) < len(selected):
        return complement
    if "|".join(selected) <= "|".join(complement):
        return selected
    return complement

def unrooted_topology_signatures(tree):
    all_taxa = descendant_taxa(tree.root)
    signatures = {{}}
    for clade in tree.find_clades(order="preorder"):
        if clade is tree.root or clade.is_terminal():
            continue
        selected = canonical_unrooted_signature(descendant_taxa(clade), all_taxa)
        if len(selected) <= 1 or len(selected) >= len(all_taxa):
            continue
        signatures[signature_id(selected)] = selected
    return signatures

def topology_distance_rows(left_tree, right_tree, rf_mode):
    if rf_mode == "rooted":
        left_signatures = rooted_topology_signatures(left_tree)
        right_signatures = rooted_topology_signatures(right_tree)
        split_kind = "clade"
    else:
        left_signatures = unrooted_topology_signatures(left_tree)
        right_signatures = unrooted_topology_signatures(right_tree)
        split_kind = "split"
    left_ids = set(left_signatures)
    right_ids = set(right_signatures)
    shared_ids = left_ids & right_ids
    left_only_ids = left_ids - right_ids
    right_only_ids = right_ids - left_ids
    all_ids = sorted(
        left_ids | right_ids,
        key=lambda split_id: (len(split_id.split("|")), split_id),
    )
    rows = [
        {{
            "split_id": split_id,
            "split_kind": split_kind,
            "comparison_status": (
                "shared"
                if split_id in shared_ids
                else "left_only"
                if split_id in left_only_ids
                else "right_only"
            ),
            "taxon_count": len((left_signatures.get(split_id) or right_signatures[split_id])),
            "descendant_taxa": "|".join(left_signatures.get(split_id) or right_signatures[split_id]),
            "left_present": split_id in left_ids,
            "right_present": split_id in right_ids,
        }}
        for split_id in all_ids
    ]
    return rows, len(left_ids), len(right_ids), len(shared_ids), len(left_only_ids), len(right_only_ids)

def informative_clades(tree, shared_taxa):
    total_tip_count = len(shared_taxa)
    clades = []
    for clade in tree.find_clades(order="preorder"):
        if clade is tree.root or clade.is_terminal():
            continue
        taxa = frozenset(descendant_taxa(clade))
        if 1 < len(taxa) < total_tip_count and taxa <= shared_taxa:
            clades.append(taxa)
    return clades

def split_counts(tree_set):
    counts = {{}}
    for tree in tree_set:
        for split_id in unrooted_topology_signatures(tree):
            counts[split_id] = counts.get(split_id, 0) + 1
    return counts

def canonical_bipartition(taxa, all_taxa):
    return "|".join(canonical_unrooted_signature(sorted(taxa), sorted(all_taxa)))

def clade_support_status(supporting_tree_count, tree_count, node_kind, unscored_reason=None):
    if node_kind == "root":
        return (
            "fixed",
            "the root spans the full compatible taxon set and is present in every comparison tree",
        )
    if supporting_tree_count is None:
        if unscored_reason == "absent-root-split":
            return (
                "not-counted",
                "ape::prop.clades leaves this root-adjacent split unscored when the comparison tree set never realizes the matching bipartition",
            )
        return (
            "not-counted",
            "ape::prop.clades leaves this root-adjacent clade unscored because its complement is a singleton tip",
        )
    if supporting_tree_count == 0:
        return (
            "absent",
            "the reference clade is absent from the comparison tree set",
        )
    if supporting_tree_count == tree_count:
        return (
            "fixed",
            "the reference clade is present in every comparison tree",
        )
    return (
        "partial-support",
        "the reference clade is present in only a subset of comparison trees",
    )

def prop_clades_rows(reference_tree, comparison_trees):
    reference_taxa = frozenset(terminal.name for terminal in reference_tree.get_terminals())
    comparison_taxa = frozenset(terminal.name for terminal in comparison_trees[0].get_terminals())
    if any(
        frozenset(terminal.name for terminal in tree.get_terminals()) != comparison_taxa
        for tree in comparison_trees[1:]
    ):
        raise ValueError(
            "reference tree support mapping requires all comparison trees to share the exact same taxon set"
        )
    if reference_taxa != comparison_taxa:
        raise ValueError(
            "reference tree and comparison tree set must share the exact same taxon set"
        )

    tree_count = len(comparison_trees)
    clade_counts = {{}}
    for tree in comparison_trees:
        for clade in informative_clades(tree, comparison_taxa):
            clade_id = "|".join(sorted(clade))
            clade_counts[clade_id] = clade_counts.get(clade_id, 0) + 1
    split_count_lookup = split_counts(comparison_trees)
    depths = node_depth_lookup(reference_tree)
    rows = []
    supported_clade_count = 0
    absent_clade_count = 0
    unscored_clade_count = 0
    root_children = list(getattr(reference_tree.root, "clades", []))
    tip_count = len(list(reference_tree.get_terminals()))
    internal_clades = list(iter_internal_clades_preorder(reference_tree.root))

    for offset, clade in enumerate(internal_clades, start=1):
        node_id = tip_count + offset
        taxa = descendant_taxa(clade)
        clade_id = "|".join(taxa)
        node_kind = "root" if clade is reference_tree.root else "internal"
        if clade is reference_tree.root:
            supporting_tree_count = tree_count
            clade_frequency = 1.0
            support_percent = 100.0
            unscored_reason = None
        elif len(taxa) == len(reference_taxa) - 1:
            supporting_tree_count = None
            clade_frequency = ""
            support_percent = ""
            unscored_reason = "singleton-complement"
            unscored_clade_count += 1
        elif clade in root_children:
            split_support = split_count_lookup.get(
                canonical_bipartition(taxa, reference_taxa),
                0,
            )
            if split_support == 0:
                supporting_tree_count = None
                clade_frequency = ""
                support_percent = ""
                unscored_reason = "absent-root-split"
                unscored_clade_count += 1
            else:
                supporting_tree_count = split_support
                clade_frequency = supporting_tree_count / tree_count
                support_percent = clade_frequency * 100.0
                unscored_reason = None
                supported_clade_count += 1
        else:
            supporting_tree_count = clade_counts.get(clade_id, 0)
            clade_frequency = supporting_tree_count / tree_count
            support_percent = clade_frequency * 100.0
            unscored_reason = None
            if supporting_tree_count == 0:
                absent_clade_count += 1
            else:
                supported_clade_count += 1

        support_status, explanation = clade_support_status(
            supporting_tree_count,
            tree_count,
            node_kind,
            unscored_reason,
        )
        node_label = normalize_label(
            clade.name if clade.name is not None else getattr(clade, "confidence", None)
        )
        rows.append(
            {{
                "node_id": node_id,
                "node_kind": node_kind,
                "node_label": node_label,
                "descendant_taxa": clade_id,
                "supporting_tree_count": ""
                if supporting_tree_count is None
                else supporting_tree_count,
                "clade_frequency": clade_frequency,
                "support_percent": support_percent,
                "support_status": support_status,
                "explanation": explanation,
                "reference_branch_length": ""
                if clade.branch_length is None
                else clade.branch_length,
                "reference_root_depth": depths[id(clade)],
            }}
        )

    summary = {{
        "tree_count": tree_count,
        "shared_taxa": sorted(reference_taxa),
        "shared_taxon_count": len(reference_taxa),
        "internal_node_count": len(internal_clades),
        "supported_clade_count": supported_clade_count,
        "absent_clade_count": absent_clade_count,
        "unscored_clade_count": unscored_clade_count,
    }}
    return summary, rows

def matrix_rank(matrix, tolerance=1e-12):
    working = [list(map(float, row)) for row in matrix]
    row_count = len(working)
    column_count = len(working[0]) if working else 0
    rank = 0
    pivot_row = 0
    for pivot_column in range(column_count):
        candidate_row = max(
            range(pivot_row, row_count),
            key=lambda index: abs(working[index][pivot_column]),
            default=None,
        )
        if candidate_row is None:
            break
        pivot_value = working[candidate_row][pivot_column]
        if abs(pivot_value) <= tolerance:
            continue
        working[pivot_row], working[candidate_row] = (
            working[candidate_row],
            working[pivot_row],
        )
        pivot = working[pivot_row][pivot_column]
        working[pivot_row] = [value / pivot for value in working[pivot_row]]
        for row_index in range(row_count):
            if row_index == pivot_row:
                continue
            factor = working[row_index][pivot_column]
            if abs(factor) <= tolerance:
                continue
            working[row_index] = [
                row_value - factor * pivot_value
                for row_value, pivot_value in zip(
                    working[row_index], working[pivot_row], strict=True
                )
            ]
        rank += 1
        pivot_row += 1
        if pivot_row == row_count:
            break
    return rank

def invert_matrix(matrix):
    size = len(matrix)
    augmented = [
        [float(value) for value in row] + [1.0 if row_index == column_index else 0.0 for column_index in range(size)]
        for row_index, row in enumerate(matrix)
    ]
    for pivot_index in range(size):
        pivot_row = max(
            range(pivot_index, size),
            key=lambda row_index: abs(augmented[row_index][pivot_index]),
        )
        pivot_value = augmented[pivot_row][pivot_index]
        if abs(pivot_value) <= 1e-12:
            raise ValueError("matrix is singular")
        if pivot_row != pivot_index:
            augmented[pivot_index], augmented[pivot_row] = augmented[pivot_row], augmented[pivot_index]
        pivot_value = augmented[pivot_index][pivot_index]
        augmented[pivot_index] = [value / pivot_value for value in augmented[pivot_index]]
        for row_index in range(size):
            if row_index == pivot_index:
                continue
            factor = augmented[row_index][pivot_index]
            if abs(factor) <= 1e-15:
                continue
            augmented[row_index] = [
                row_value - factor * pivot_value
                for row_value, pivot_value in zip(
                    augmented[row_index], augmented[pivot_index], strict=True
                )
            ]
    return [row[size:] for row in augmented]

def matrix_infinity_norm(matrix):
    return max((sum(abs(value) for value in row) for row in matrix), default=0.0)

def symmetric_matrix_eigenvalues(matrix, tolerance=1e-15, max_iterations=10000):
    size = len(matrix)
    if size == 0:
        return []
    if size == 1:
        return [float(matrix[0][0])]
    working = [list(map(float, row)) for row in matrix]
    for _ in range(max_iterations):
        pivot_row = 0
        pivot_column = 1
        pivot_value = 0.0
        for row_index in range(size):
            for column_index in range(row_index + 1, size):
                candidate = abs(working[row_index][column_index])
                if candidate > pivot_value:
                    pivot_row = row_index
                    pivot_column = column_index
                    pivot_value = candidate
        if pivot_value <= tolerance:
            return [working[index][index] for index in range(size)]
        app = working[pivot_row][pivot_row]
        aqq = working[pivot_column][pivot_column]
        apq = working[pivot_row][pivot_column]
        tau = (aqq - app) / (2.0 * apq)
        tangent = (
            math.copysign(1.0, tau) / (abs(tau) + math.sqrt(1.0 + tau * tau))
            if abs(tau) > tolerance
            else 1.0
        )
        cosine = 1.0 / math.sqrt(1.0 + tangent * tangent)
        sine = tangent * cosine
        for index in range(size):
            if index in (pivot_row, pivot_column):
                continue
            left = working[index][pivot_row]
            right = working[index][pivot_column]
            working[index][pivot_row] = working[pivot_row][index] = cosine * left - sine * right
            working[index][pivot_column] = working[pivot_column][index] = sine * left + cosine * right
        working[pivot_row][pivot_row] = (
            cosine * cosine * app
            - 2.0 * sine * cosine * apq
            + sine * sine * aqq
        )
        working[pivot_column][pivot_column] = (
            sine * sine * app
            + 2.0 * sine * cosine * apq
            + cosine * cosine * aqq
        )
        working[pivot_row][pivot_column] = 0.0
        working[pivot_column][pivot_row] = 0.0
    raise ValueError("symmetric eigenvalue iteration did not converge")

def symmetric_matrix_condition_number(matrix, tolerance=1e-12):
    singular_values = sorted(
        abs(value)
        for value in symmetric_matrix_eigenvalues(matrix, tolerance=tolerance)
    )
    if not singular_values:
        return 0.0
    if singular_values[0] <= tolerance:
        return math.inf
    return singular_values[-1] / singular_values[0]

def matrix_log_determinant(matrix):
    size = len(matrix)
    working = [list(map(float, row)) for row in matrix]
    sign = 1.0
    log_abs_det = 0.0
    for pivot_index in range(size):
        pivot_row = max(
            range(pivot_index, size),
            key=lambda row_index: abs(working[row_index][pivot_index]),
        )
        pivot_value = working[pivot_row][pivot_index]
        if abs(pivot_value) <= 1e-12:
            raise ValueError("matrix determinant is zero")
        if pivot_row != pivot_index:
            working[pivot_index], working[pivot_row] = working[pivot_row], working[pivot_index]
            sign *= -1.0
        pivot_value = working[pivot_index][pivot_index]
        if pivot_value < 0:
            sign *= -1.0
        log_abs_det += math.log(abs(pivot_value))
        for row_index in range(pivot_index + 1, size):
            factor = working[row_index][pivot_index] / pivot_value
            if abs(factor) <= 1e-15:
                continue
            for column_index in range(pivot_index, size):
                working[row_index][column_index] -= factor * working[pivot_index][column_index]
    if sign <= 0:
        raise ValueError("matrix determinant is not positive")
    return log_abs_det

def ancestor_depths(tree):
    lookup = {{}}
    def walk(clade, depth, path):
        current = dict(path)
        current[id(clade)] = depth
        if clade.is_terminal():
            lookup[clade.name] = current
            return
        for child in clade.clades:
            walk(child, depth + (child.branch_length or 0.0), current)
    walk(tree.root, 0.0, {{}})
    return lookup

runner_path = Path(sys.argv[1])
case_path = Path(sys.argv[2])
output_root = Path(sys.argv[3])
output_root.mkdir(parents=True, exist_ok=True)
case_payload = json.loads(case_path.read_text(encoding="utf-8"))
execution_path = output_root / "reference-execution.json"

if not APE_AVAILABLE:
    write_json(
        execution_path,
        {{
            "status": "unavailable",
            "mismatch_reason": "ape_package_unavailable",
            "message": "ape is not installed in the fake reference environment",
            "case_id": case_payload["case_id"],
            "function_name": case_payload["function_name"],
            "input_fixture": case_payload["input_fixture"],
            "r_version": "4.6.0",
            "ape_version": None,
        }},
    )
    raise SystemExit(0)

case_id = case_payload["case_id"]
if case_payload["operation"] in {{
    "read-tree-structure",
    "write-tree-structure",
    "root-tree-outgroup",
    "unroot-tree",
    "drop-tree-taxa",
    "keep-tree-taxa",
    "extract-tree-clade",
    "read-tree-set-structure",
    "write-tree-set-structure",
}}:
    try:
        if case_payload["operation"] in {{"read-tree-structure", "write-tree-structure", "root-tree-outgroup", "unroot-tree", "drop-tree-taxa", "keep-tree-taxa", "extract-tree-clade"}}:
            newick_path = output_root / "normalized-tree.nwk"
            if case_payload["operation"] == "root-tree-outgroup":
                outgroup_taxa = tuple(case_payload.get("outgroup_taxa", []))
                if outgroup_taxa == ("Z",):
                    raise ValueError("specified outgroup not in labels of the tree")
                if outgroup_taxa == ("B", "D"):
                    raise ValueError("the specified outgroup is not monophyletic")
                if outgroup_taxa == ("D",):
                    newick_text = "(((A:0.2,B:0.2):0.7,C:0.1):0,D:0.1);\\n"
                elif outgroup_taxa == ("C", "D"):
                    newick_text = "((A:0.2,B:0.2):0.7,(C:0.1,D:0.1):0);\\n"
                else:
                    newick_text = Path(case_payload["input_fixture"]).read_text(encoding="utf-8")
                newick_path.write_text(newick_text, encoding="utf-8")
                tree = Phylo.read(newick_path, "newick")
            elif case_payload["operation"] == "unroot-tree":
                if case_id == "unroot-tree-balanced-rooted":
                    newick_text = "(A:0.1,B:0.1,(C:0.2,D:0.2):0.3);\\n"
                elif case_id == "unroot-tree-rootable":
                    newick_text = "(A:0.2,B:0.2,(C:0.1,D:0.1):0.7);\\n"
                elif case_id == "unroot-tree-after-outgroup-rooting":
                    newick_text = "((A:0.2,B:0.2):0.7,C:0.1,D:0.1);\\n"
                else:
                    newick_text = Path(case_payload["input_fixture"]).read_text(encoding="utf-8")
                newick_path.write_text(newick_text, encoding="utf-8")
                tree = Phylo.read(newick_path, "newick")
            elif case_payload["operation"] == "drop-tree-taxa":
                if case_id == "drop-tip-rooted-single":
                    newick_text = "((A:0.1,B:0.1):0.2,C:0.3);\\n"
                    dropped_taxa = ["D"]
                    absent_requested_taxa = []
                elif case_id == "drop-tip-rooted-multiple":
                    newick_text = "(A:0.3,C:0.3);\\n"
                    dropped_taxa = ["B", "D"]
                    absent_requested_taxa = []
                elif case_id == "drop-tip-root-change-after-outgroup-rooting":
                    newick_text = "((A:0.2,B:0.2):0.7,C:0.1);\\n"
                    dropped_taxa = ["D"]
                    absent_requested_taxa = []
                elif case_id == "drop-tip-unrooted-three-tip":
                    newick_text = "(A:0.1,B:0.2,C:0.3);\\n"
                    dropped_taxa = ["D"]
                    absent_requested_taxa = []
                elif case_id == "drop-tip-unrooted-two-tip":
                    newick_text = "(A:0.1,B:0.2);\\n"
                    dropped_taxa = ["C", "D"]
                    absent_requested_taxa = []
                else:
                    newick_text = Path(case_payload["input_fixture"]).read_text(encoding="utf-8")
                    dropped_taxa = []
                    absent_requested_taxa = ["Z"]
                newick_path.write_text(newick_text, encoding="utf-8")
                tree = Phylo.read(newick_path, "newick")
            elif case_payload["operation"] == "keep-tree-taxa":
                if case_id == "keep-tip-rooted-selected-two":
                    newick_text = "(A:0.3,C:0.3);\\n"
                    requested_taxa = ["A", "C"]
                    dropped_taxa = ["B", "D"]
                elif case_id == "keep-tip-rooted-order-insensitive":
                    newick_text = "(A:0.3,C:0.3);\\n"
                    requested_taxa = ["A", "C"]
                    dropped_taxa = ["B", "D"]
                elif case_id == "keep-tip-root-change-after-outgroup-rooting":
                    newick_text = "((A:0.2,B:0.2):0.7,C:0.1);\\n"
                    requested_taxa = ["A", "B", "C"]
                    dropped_taxa = ["D"]
                elif case_id == "keep-tip-unrooted-three-tip":
                    newick_text = "(A:0.1,B:0.2,C:0.3);\\n"
                    requested_taxa = ["A", "B", "C"]
                    dropped_taxa = ["D"]
                elif case_id == "keep-tip-unrooted-two-tip":
                    newick_text = "(A:0.1,B:0.2);\\n"
                    requested_taxa = ["A", "B"]
                    dropped_taxa = ["C", "D"]
                else:
                    newick_text = Path(case_payload["input_fixture"]).read_text(encoding="utf-8")
                    requested_taxa = sorted(set(case_payload.get("requested_taxa", [])))
                    dropped_taxa = []
                newick_path.write_text(newick_text, encoding="utf-8")
                tree = Phylo.read(newick_path, "newick")
            elif case_payload["operation"] == "extract-tree-clade":
                if case_id == "extract-clade-tip-node-invalid":
                    raise ValueError("node number must be greater than the number of tips")
                elif case_id == "extract-clade-node-out-of-bounds":
                    raise IndexError("subscript out of bounds")
                elif case_id == "extract-clade-root":
                    newick_text = "((A:0.1,B:0.1)Mammals:0.2,(C:0.2,D:0.2)Birds:0.1)Root;\\n"
                    requested_node_id = 5
                    matched_node_id = 5
                    matched_node_name = "Root"
                elif case_id == "extract-clade-mammals":
                    newick_text = "(A:0.1,B:0.1)Mammals;\\n"
                    requested_node_id = 6
                    matched_node_id = 6
                    matched_node_name = "Mammals"
                else:
                    newick_text = "(C:0.2,D:0.2)Birds;\\n"
                    requested_node_id = 7
                    matched_node_id = 7
                    matched_node_name = "Birds"
                newick_path.write_text(newick_text, encoding="utf-8")
                tree = Phylo.read(newick_path, "newick")
            else:
                tree = Phylo.read(case_payload["input_fixture"], "newick")
            summary = {{
                "tree_count": 1,
                "tip_count": len(tree.get_terminals()),
                "internal_node_count": len(tree.get_nonterminals()),
                "edge_count": len(tree.get_terminals()) + len(tree.get_nonterminals()) - 1,
                "rooted": is_rooted_tree(tree),
                "tip_labels": [terminal.name for terminal in tree.get_terminals()],
                "branch_length_count": sum(
                    1 for clade in tree.find_clades(order="preorder")
                    if clade is not tree.root and clade.branch_length is not None
                ),
            }}
            if case_payload["operation"] == "drop-tree-taxa":
                summary["dropped_taxa"] = dropped_taxa
                summary["absent_requested_taxa"] = absent_requested_taxa
            if case_payload["operation"] == "keep-tree-taxa":
                summary["requested_taxa"] = requested_taxa
                summary["dropped_taxa"] = dropped_taxa
            if case_payload["operation"] == "extract-tree-clade":
                summary["requested_node_id"] = requested_node_id
                summary["matched_node_id"] = matched_node_id
                summary["matched_node_name"] = matched_node_name
            rows = clade_rows(tree, "")
            summary.update(SUMMARY_OVERRIDES)
            summary_path = output_root / "summary.json"
            clades_path = output_root / "clades.tsv"
            write_json(summary_path, summary)
            write_tsv(clades_path, rows)
            if case_payload["operation"] in {{"root-tree-outgroup", "unroot-tree", "drop-tree-taxa", "keep-tree-taxa", "extract-tree-clade"}}:
                pass
            elif case_id == "read-tree-quoted-taxon-labels":
                newick_path.write_text(
                    "('Homo_sapiens':0.1,'Mus_musculus':0.2,'A.B-1':0.3);\\n",
                    encoding="utf-8",
                )
            elif case_id in NORMALIZED_TREE_OVERRIDES:
                newick_path.write_text(
                    NORMALIZED_TREE_OVERRIDES[case_id],
                    encoding="utf-8",
                )
            else:
                newick_path.write_text(
                    Path(case_payload["input_fixture"]).read_text(encoding="utf-8"),
                    encoding="utf-8",
                )
            outputs = {{
                "summary_json": str(summary_path),
                "clades": str(clades_path),
                "normalized_tree": str(newick_path),
            }}
        else:
            trees = list(Phylo.parse(case_payload["input_fixture"], "newick"))
            if not trees:
                raise ValueError("tree set contains no trees")
            rows = []
            for index, tree in enumerate(trees, start=1):
                rows.extend(clade_rows(tree, index))
            shared_tip_labels = sorted(set(t.name for t in trees[0].get_terminals()))
            summary = {{
                "tree_count": len(trees),
                "source_format": "newick",
                "tree_indices": list(range(1, len(trees) + 1)),
                "shared_tip_labels": shared_tip_labels,
                "unique_tip_label_count": len(shared_tip_labels),
            }}
            summary.update(SUMMARY_OVERRIDES)
            summary_path = output_root / "summary.json"
            clades_path = output_root / "clades.tsv"
            tree_set_path = output_root / "normalized-tree-set.nwk"
            write_json(summary_path, summary)
            write_tsv(clades_path, rows)
            tree_set_path.write_text(
                Path(case_payload["input_fixture"]).read_text(encoding="utf-8"),
                encoding="utf-8",
            )
            outputs = {{
                "summary_json": str(summary_path),
                "clades": str(clades_path),
                "normalized_tree_set": str(tree_set_path),
            }}
    except Exception as error:
        write_json(
            execution_path,
            {{
                "status": "failed",
                "mismatch_reason": "reference_execution_failed",
                "error_type": "TreeRootingError" if case_payload["operation"] == "root-tree-outgroup" else "TreeParseError",
                "message": str(error),
                "case_id": case_payload["case_id"],
                "function_name": case_payload["function_name"],
                "input_fixture": case_payload["input_fixture"],
                "r_version": "4.6.0",
                "ape_version": "5.0.0",
            }},
        )
        raise SystemExit(0)

    write_json(
        execution_path,
        {{
            "status": "ok",
            "case_id": case_payload["case_id"],
            "function_name": case_payload["function_name"],
            "input_fixture": case_payload["input_fixture"],
            "r_version": "4.6.0",
            "ape_version": "5.0.0",
            "outputs": outputs,
        }},
    )
    raise SystemExit(0)

if case_payload["operation"] == "get-tree-mrca":
    try:
        if case_id == "get-mrca-missing-tip":
            raise ValueError("missing value where TRUE/FALSE needed")
        if case_id == "get-mrca-balanced-two-tip":
            summary = {{
                "requested_taxa": ["A", "B"],
                "unique_requested_taxa": ["A", "B"],
                "duplicate_requested_taxa": [],
                "matched_node_id": 6,
                "matched_node_name": "",
                "matched_taxa": ["A", "B"],
                "matched_extra_taxa": [],
                "matched_tip_count": 2,
                "is_root": False,
            }}
        elif case_id == "get-mrca-balanced-full-tip-set":
            summary = {{
                "requested_taxa": ["A", "B", "C", "D"],
                "unique_requested_taxa": ["A", "B", "C", "D"],
                "duplicate_requested_taxa": [],
                "matched_node_id": 5,
                "matched_node_name": "",
                "matched_taxa": ["A", "B", "C", "D"],
                "matched_extra_taxa": [],
                "matched_tip_count": 4,
                "is_root": True,
            }}
        elif case_id == "get-mrca-balanced-duplicate-request":
            summary = {{
                "requested_taxa": ["A", "A", "B"],
                "unique_requested_taxa": ["A", "B"],
                "duplicate_requested_taxa": ["A"],
                "matched_node_id": 6,
                "matched_node_name": "",
                "matched_taxa": ["A", "B"],
                "matched_extra_taxa": [],
                "matched_tip_count": 2,
                "is_root": False,
            }}
        elif case_id == "get-mrca-pectinate-many-tip":
            summary = {{
                "requested_taxa": ["A", "B", "C"],
                "unique_requested_taxa": ["A", "B", "C"],
                "duplicate_requested_taxa": [],
                "matched_node_id": 6,
                "matched_node_name": "",
                "matched_taxa": ["A", "B", "C"],
                "matched_extra_taxa": [],
                "matched_tip_count": 3,
                "is_root": False,
            }}
        elif case_id == "get-mrca-rooted-polytomy":
            summary = {{
                "requested_taxa": ["A", "B", "C"],
                "unique_requested_taxa": ["A", "B", "C"],
                "duplicate_requested_taxa": [],
                "matched_node_id": 6,
                "matched_node_name": "",
                "matched_taxa": ["A", "B", "C"],
                "matched_extra_taxa": [],
                "matched_tip_count": 3,
                "is_root": False,
            }}
        else:
            summary = {{
                "requested_taxa": ["A", "B", "C"],
                "unique_requested_taxa": ["A", "B", "C"],
                "duplicate_requested_taxa": [],
                "matched_node_id": 6,
                "matched_node_name": "",
                "matched_taxa": ["A", "B", "C"],
                "matched_extra_taxa": [],
                "matched_tip_count": 3,
                "is_root": False,
            }}
    except Exception as error:
        write_json(
            execution_path,
            {{
                "status": "failed",
                "mismatch_reason": "reference_execution_failed",
                "error_type": "TreeMrcaError",
                "message": str(error),
                "case_id": case_payload["case_id"],
                "function_name": case_payload["function_name"],
                "input_fixture": case_payload["input_fixture"],
                "r_version": "4.6.0",
                "ape_version": "5.0.0",
            }},
        )
        raise SystemExit(0)

    summary.update(SUMMARY_OVERRIDES)
    summary_path = output_root / "summary.json"
    write_json(summary_path, summary)
    write_json(
        execution_path,
        {{
            "status": "ok",
            "case_id": case_payload["case_id"],
            "function_name": case_payload["function_name"],
            "input_fixture": case_payload["input_fixture"],
            "r_version": "4.6.0",
            "ape_version": "5.0.0",
            "outputs": {{"summary_json": str(summary_path)}},
        }},
    )
    raise SystemExit(0)

if case_payload["operation"] == "assess-tree-monophyly":
    try:
        if case_id == "is-monophyletic-all-missing-rerooted":
            raise ValueError("specified outgroup not in labels of the tree")
        if case_id == "is-monophyletic-rooted-two-tip":
            summary = {{
                "requested_taxa": ["A", "B"],
                "unique_requested_taxa": ["A", "B"],
                "duplicate_requested_taxa": [],
                "missing_requested_taxa": [],
                "present_requested_taxa": ["A", "B"],
                "reroot": False,
                "rooted": True,
                "monophyletic": True,
                "complementary_clade_used": False,
                "matched_node_id": 6,
                "matched_node_name": "",
                "matched_taxa": ["A", "B"],
                "matched_extra_taxa": [],
                "matched_tip_count": 2,
                "is_root": False,
            }}
        elif case_id == "is-monophyletic-rooted-three-tip-reroot-false":
            summary = {{
                "requested_taxa": ["A", "B", "C"],
                "unique_requested_taxa": ["A", "B", "C"],
                "duplicate_requested_taxa": [],
                "missing_requested_taxa": [],
                "present_requested_taxa": ["A", "B", "C"],
                "reroot": False,
                "rooted": True,
                "monophyletic": False,
                "complementary_clade_used": False,
                "matched_node_id": 5,
                "matched_node_name": "",
                "matched_taxa": ["A", "B", "C", "D"],
                "matched_extra_taxa": ["D"],
                "matched_tip_count": 4,
                "is_root": True,
            }}
        elif case_id == "is-monophyletic-rooted-three-tip-reroot-true":
            summary = {{
                "requested_taxa": ["A", "B", "C"],
                "unique_requested_taxa": ["A", "B", "C"],
                "duplicate_requested_taxa": [],
                "missing_requested_taxa": [],
                "present_requested_taxa": ["A", "B", "C"],
                "reroot": True,
                "rooted": True,
                "monophyletic": True,
                "complementary_clade_used": True,
                "matched_node_id": 5,
                "matched_node_name": "",
                "matched_taxa": ["A", "B", "C", "D"],
                "matched_extra_taxa": ["D"],
                "matched_tip_count": 4,
                "is_root": True,
            }}
        elif case_id == "is-monophyletic-rooted-full-tip-set":
            summary = {{
                "requested_taxa": ["A", "B", "C", "D"],
                "unique_requested_taxa": ["A", "B", "C", "D"],
                "duplicate_requested_taxa": [],
                "missing_requested_taxa": [],
                "present_requested_taxa": ["A", "B", "C", "D"],
                "reroot": False,
                "rooted": True,
                "monophyletic": True,
                "complementary_clade_used": False,
                "matched_node_id": 5,
                "matched_node_name": "",
                "matched_taxa": ["A", "B", "C", "D"],
                "matched_extra_taxa": [],
                "matched_tip_count": 4,
                "is_root": True,
            }}
        elif case_id == "is-monophyletic-rooted-mixed-missing":
            summary = {{
                "requested_taxa": ["A", "Z"],
                "unique_requested_taxa": ["A", "Z"],
                "duplicate_requested_taxa": [],
                "missing_requested_taxa": ["Z"],
                "present_requested_taxa": ["A"],
                "reroot": False,
                "rooted": True,
                "monophyletic": True,
                "complementary_clade_used": False,
                "matched_node_id": 1,
                "matched_node_name": "A",
                "matched_taxa": ["A"],
                "matched_extra_taxa": [],
                "matched_tip_count": 1,
                "is_root": False,
            }}
        elif case_id == "is-monophyletic-unrooted-two-tip":
            summary = {{
                "requested_taxa": ["A", "B"],
                "unique_requested_taxa": ["A", "B"],
                "duplicate_requested_taxa": [],
                "missing_requested_taxa": [],
                "present_requested_taxa": ["A", "B"],
                "reroot": True,
                "rooted": False,
                "monophyletic": False,
                "complementary_clade_used": False,
                "matched_node_id": 5,
                "matched_node_name": "",
                "matched_taxa": ["A", "B", "C", "D"],
                "matched_extra_taxa": ["C", "D"],
                "matched_tip_count": 4,
                "is_root": True,
            }}
        elif case_id == "is-monophyletic-unrooted-three-tip":
            summary = {{
                "requested_taxa": ["A", "B", "C"],
                "unique_requested_taxa": ["A", "B", "C"],
                "duplicate_requested_taxa": [],
                "missing_requested_taxa": [],
                "present_requested_taxa": ["A", "B", "C"],
                "reroot": True,
                "rooted": False,
                "monophyletic": True,
                "complementary_clade_used": True,
                "matched_node_id": 5,
                "matched_node_name": "",
                "matched_taxa": ["A", "B", "C", "D"],
                "matched_extra_taxa": ["D"],
                "matched_tip_count": 4,
                "is_root": True,
            }}
        elif case_id == "is-monophyletic-after-outgroup-rooting":
            summary = {{
                "requested_taxa": ["A", "B", "C"],
                "unique_requested_taxa": ["A", "B", "C"],
                "duplicate_requested_taxa": [],
                "missing_requested_taxa": [],
                "present_requested_taxa": ["A", "B", "C"],
                "reroot": False,
                "rooted": True,
                "monophyletic": True,
                "complementary_clade_used": False,
                "matched_node_id": 6,
                "matched_node_name": "",
                "matched_taxa": ["A", "B", "C"],
                "matched_extra_taxa": [],
                "matched_tip_count": 3,
                "is_root": False,
            }}
        else:
            summary = {{
                "requested_taxa": ["A", "B", "C"],
                "unique_requested_taxa": ["A", "B", "C"],
                "duplicate_requested_taxa": [],
                "missing_requested_taxa": [],
                "present_requested_taxa": ["A", "B", "C"],
                "reroot": False,
                "rooted": True,
                "monophyletic": True,
                "complementary_clade_used": False,
                "matched_node_id": 6,
                "matched_node_name": "",
                "matched_taxa": ["A", "B", "C"],
                "matched_extra_taxa": [],
                "matched_tip_count": 3,
                "is_root": False,
            }}
    except Exception as error:
        write_json(
            execution_path,
            {{
                "status": "failed",
                "mismatch_reason": "reference_execution_failed",
                "error_type": "TreeMonophylyError",
                "message": str(error),
                "case_id": case_payload["case_id"],
                "function_name": case_payload["function_name"],
                "input_fixture": case_payload["input_fixture"],
                "r_version": "4.6.0",
                "ape_version": "5.0.0",
            }},
        )
        raise SystemExit(0)

    summary.update(SUMMARY_OVERRIDES)
    summary_path = output_root / "summary.json"
    write_json(summary_path, summary)
    write_json(
        execution_path,
        {{
            "status": "ok",
            "case_id": case_payload["case_id"],
            "function_name": case_payload["function_name"],
            "input_fixture": case_payload["input_fixture"],
            "r_version": "4.6.0",
            "ape_version": "5.0.0",
            "outputs": {{"summary_json": str(summary_path)}},
        }},
    )
    raise SystemExit(0)

if case_payload["operation"] == "tree-consensus":
    if case_id == "consensus-mismatched-taxon-set":
        write_json(
            execution_path,
            {{
                "status": "failed",
                "mismatch_reason": "reference_execution_failed",
                "error_type": "ConsensusTreeError",
                "message": "consensus requires all trees to share the exact same taxon set",
                "case_id": case_payload["case_id"],
                "function_name": case_payload["function_name"],
                "input_fixture": case_payload["input_fixture"],
                "r_version": "4.6.0",
                "ape_version": "5.0.0",
            }},
        )
        raise SystemExit(0)

    if case_id == "consensus-majority-conflicting-four-taxon":
        newick_text = "((A:0.1,B:0.1)66.6666666666667:0.2,(C:0.1,D:0.1)66.6666666666667:0.2);\\n"
        rows = [
            {{"clade": "A|B", "tree_count": 2, "frequency": 2 / 3}},
            {{"clade": "A|C", "tree_count": 1, "frequency": 1 / 3}},
            {{"clade": "B|D", "tree_count": 1, "frequency": 1 / 3}},
            {{"clade": "C|D", "tree_count": 2, "frequency": 2 / 3}},
        ]
        summary = {{
            "tree_count": 3,
            "shared_taxa": ["A", "B", "C", "D"],
            "shared_taxon_count": 4,
            "tip_count": 4,
            "rooted": False,
            "consensus_method": "majority-rule",
            "consensus_threshold": 0.5,
            "included_clade_count": 1,
            "clade_frequency_count": 4,
        }}
    elif case_id == "consensus-strict-conflicting-four-taxon":
        newick_text = "(A:0.1,B:0.1,C:0.1,D:0.1);\\n"
        rows = [
            {{"clade": "A|B", "tree_count": 2, "frequency": 2 / 3}},
            {{"clade": "A|C", "tree_count": 1, "frequency": 1 / 3}},
            {{"clade": "B|D", "tree_count": 1, "frequency": 1 / 3}},
            {{"clade": "C|D", "tree_count": 2, "frequency": 2 / 3}},
        ]
        summary = {{
            "tree_count": 3,
            "shared_taxa": ["A", "B", "C", "D"],
            "shared_taxon_count": 4,
            "tip_count": 4,
            "rooted": False,
            "consensus_method": "strict",
            "consensus_threshold": 1.0,
            "included_clade_count": 0,
            "clade_frequency_count": 4,
        }}
    elif case_id == "consensus-majority-posterior-six-taxon":
        newick_text = "(A:1,B:1,(C:1,D:1)60:1.66666666666667,(E:1,F:1)60:1.66666666666667);\\n"
        rows = [
            {{"clade": "A|B", "tree_count": 2, "frequency": 0.4}},
            {{"clade": "A|B|C|D", "tree_count": 2, "frequency": 0.4}},
            {{"clade": "A|B|D|E", "tree_count": 1, "frequency": 0.2}},
            {{"clade": "A|B|E|F", "tree_count": 2, "frequency": 0.4}},
            {{"clade": "A|C", "tree_count": 1, "frequency": 0.2}},
            {{"clade": "A|D", "tree_count": 1, "frequency": 0.2}},
            {{"clade": "A|E", "tree_count": 1, "frequency": 0.2}},
            {{"clade": "B|D", "tree_count": 1, "frequency": 0.2}},
            {{"clade": "B|E", "tree_count": 1, "frequency": 0.2}},
            {{"clade": "B|F", "tree_count": 1, "frequency": 0.2}},
            {{"clade": "C|D", "tree_count": 3, "frequency": 0.6}},
            {{"clade": "C|F", "tree_count": 1, "frequency": 0.2}},
            {{"clade": "E|F", "tree_count": 3, "frequency": 0.6}},
        ]
        summary = {{
            "tree_count": 5,
            "shared_taxa": ["A", "B", "C", "D", "E", "F"],
            "shared_taxon_count": 6,
            "tip_count": 6,
            "rooted": False,
            "consensus_method": "majority-rule",
            "consensus_threshold": 0.5,
            "included_clade_count": 2,
            "clade_frequency_count": 13,
        }}
    else:
        raise ValueError(f"unsupported fake consensus case: {{case_id}}")

    summary.update(SUMMARY_OVERRIDES)
    summary_path = output_root / "summary.json"
    rows_path = output_root / "clade-frequencies.tsv"
    newick_path = output_root / "normalized-tree.nwk"
    write_json(summary_path, summary)
    write_tsv(rows_path, rows)
    newick_path.write_text(newick_text, encoding="utf-8")
    write_json(
        execution_path,
        {{
            "status": "ok",
            "case_id": case_payload["case_id"],
            "function_name": case_payload["function_name"],
            "input_fixture": case_payload["input_fixture"],
            "r_version": "4.6.0",
            "ape_version": "5.0.0",
            "outputs": {{
                "summary_json": str(summary_path),
                "clade_frequencies": str(rows_path),
                "normalized_tree": str(newick_path),
            }},
        }},
    )
    raise SystemExit(0)

if case_payload["operation"] == "tree-clade-support":
    try:
        reference_tree = Phylo.read(case_payload["reference_tree_path"], "newick")
        comparison_trees = list(Phylo.parse(case_payload["input_fixture"], "newick"))
        if not comparison_trees:
            raise ValueError("tree set contains no trees")
        summary, rows = prop_clades_rows(reference_tree, comparison_trees)
    except Exception as error:
        write_json(
            execution_path,
            {{
                "status": "failed",
                "mismatch_reason": "reference_execution_failed",
                "error_type": "PropCladesError",
                "message": str(error),
                "case_id": case_payload["case_id"],
                "function_name": case_payload["function_name"],
                "input_fixture": case_payload["input_fixture"],
                "r_version": "4.6.0",
                "ape_version": "5.0.0",
            }},
        )
        raise SystemExit(0)

    summary.update(SUMMARY_OVERRIDES)
    summary_path = output_root / "summary.json"
    rows_path = output_root / "support-table.tsv"
    write_json(summary_path, summary)
    write_tsv(rows_path, rows)
    write_json(
        execution_path,
        {{
            "status": "ok",
            "case_id": case_payload["case_id"],
            "function_name": case_payload["function_name"],
            "input_fixture": case_payload["input_fixture"],
            "r_version": "4.6.0",
            "ape_version": "5.0.0",
            "outputs": {{
                "summary_json": str(summary_path),
                "support_table": str(rows_path),
            }},
        }},
    )
    raise SystemExit(0)

if case_payload["operation"] == "tree-tip-distance":
    tree = Phylo.read(case_payload["input_fixture"], "newick")
    tip_labels = [terminal.name for terminal in tree.get_terminals()]
    rows = []
    for left in tip_labels:
        for right in tip_labels:
            rows.append(
                {{
                    "left_identifier": left,
                    "right_identifier": right,
                    "distance": tree.distance(left, right),
                }}
            )
    summary = {{
        "tip_count": len(tip_labels),
        "rooted": is_rooted_tree(tree),
        "tip_labels": tip_labels,
        "pair_count": len(rows),
        "diagonal_zero": True,
        "symmetric": True,
        "complete_branch_lengths": True,
        "missing_branch_length_policy": "error",
    }}
    summary.update(SUMMARY_OVERRIDES)
    summary_path = output_root / "summary.json"
    rows_path = output_root / "tip-distance-long.tsv"
    write_json(summary_path, summary)
    write_tsv(rows_path, rows)
    write_json(
        execution_path,
        {{
            "status": "ok",
            "case_id": case_payload["case_id"],
            "function_name": case_payload["function_name"],
            "input_fixture": case_payload["input_fixture"],
            "r_version": "4.6.0",
            "ape_version": "5.0.0",
            "outputs": {{
                "summary_json": str(summary_path),
                "tip_distance_long": str(rows_path),
            }},
        }},
    )
    raise SystemExit(0)

if case_payload["operation"] == "tree-topology-distance":
    trees = list(Phylo.parse(case_payload["input_fixture"], "newick"))
    if len(trees) != 2:
        raise ValueError("ape topology-distance parity fixtures must contain exactly two trees")
    left_tree, right_tree = trees
    rf_mode = case_payload.get("rf_mode", "rooted")
    (
        rows,
        left_split_count,
        right_split_count,
        shared_split_count,
        left_only_split_count,
        right_only_split_count,
    ) = topology_distance_rows(left_tree, right_tree, rf_mode)
    summary = {{
        "tip_count": len(descendant_taxa(left_tree.root)),
        "shared_taxa": descendant_taxa(left_tree.root),
        "left_only_taxa": [],
        "right_only_taxa": [],
        "taxon_overlap_policy": "require-identical",
        "rf_mode": rf_mode,
        "rooted_left": is_rooted_tree(left_tree),
        "rooted_right": is_rooted_tree(right_tree),
        "polytomy_present_left": any(len(clade.clades) > 2 for clade in left_tree.find_clades()),
        "polytomy_present_right": any(len(clade.clades) > 2 for clade in right_tree.find_clades()),
        "left_split_count": left_split_count,
        "right_split_count": right_split_count,
        "shared_split_count": shared_split_count,
        "left_only_split_count": left_only_split_count,
        "right_only_split_count": right_only_split_count,
        "robinson_foulds_distance": left_only_split_count + right_only_split_count,
        "normalized_robinson_foulds": (
            0.0
            if (left_split_count + right_split_count) == 0
            else (left_only_split_count + right_only_split_count)
            / (left_split_count + right_split_count)
        ),
        "topology_equal": (left_only_split_count + right_only_split_count) == 0,
    }}
    summary.update(SUMMARY_OVERRIDES)
    summary_path = output_root / "summary.json"
    rows_path = output_root / "split-table.tsv"
    write_json(summary_path, summary)
    write_tsv(rows_path, rows)
    write_json(
        execution_path,
        {{
            "status": "ok",
            "case_id": case_payload["case_id"],
            "function_name": case_payload["function_name"],
            "input_fixture": case_payload["input_fixture"],
            "r_version": "4.6.0",
            "ape_version": "5.0.0",
            "outputs": {{
                "summary_json": str(summary_path),
                "split_table": str(rows_path),
            }},
        }},
    )
    raise SystemExit(0)

if case_payload["operation"] == "tree-brownian-covariance":
    tree = Phylo.read(case_payload["input_fixture"], "newick")
    tip_labels = [terminal.name for terminal in tree.get_terminals()]
    depth_lookup = ancestor_depths(tree)
    covariance = []
    for left in tip_labels:
        row = []
        left_path = depth_lookup[left]
        for right in tip_labels:
            right_path = depth_lookup[right]
            shared_depth = max(
                left_path[node_id]
                for node_id in set(left_path) & set(right_path)
            )
            row.append(shared_depth)
        covariance.append(row)
    root_depths = [covariance[index][index] for index in range(len(tip_labels))]
    branch_lengths = [
        clade.branch_length
        for clade in tree.find_clades(order="preorder")
        if clade is not tree.root and clade.branch_length is not None
    ]
    covariance_rank = matrix_rank(covariance)
    singular = covariance_rank < len(covariance)
    try:
        raw_log_determinant = matrix_log_determinant(covariance)
        positive_definite = True
    except ValueError:
        raw_log_determinant = None
        positive_definite = False
    condition_number = None if singular else symmetric_matrix_condition_number(covariance)
    near_singular = singular or (
        condition_number is not None and condition_number >= 1e12
    )
    matrix_rows = []
    long_rows = []
    for row_index, left in enumerate(tip_labels):
        matrix_row = {{"taxon": left}}
        for column_index, right in enumerate(tip_labels):
            value = covariance[row_index][column_index]
            matrix_row[right] = value
            long_rows.append(
                {{
                    "left_taxon": left,
                    "right_taxon": right,
                    "shared_ancestry_covariance": value,
                }}
            )
        matrix_rows.append(matrix_row)
    summary = {{
        "tip_count": len(tip_labels),
        "rooted": is_rooted_tree(tree),
        "tip_labels": tip_labels,
        "pair_count": len(long_rows),
        "tree_is_ultrametric": max(root_depths) - min(root_depths) <= 1e-12,
        "minimum_root_to_tip_depth": min(root_depths),
        "maximum_root_to_tip_depth": max(root_depths),
        "minimum_branch_length": min(branch_lengths),
        "maximum_branch_length": max(branch_lengths),
        "matrix_dimension": len(covariance),
        "matrix_rank": covariance_rank,
        "singular": singular,
        "near_singular": near_singular,
        "positive_definite": positive_definite,
        "condition_number": condition_number,
        "raw_log_determinant": raw_log_determinant,
    }}
    summary.update(SUMMARY_OVERRIDES)
    summary_path = output_root / "summary.json"
    matrix_path = output_root / "covariance-matrix.tsv"
    rows_path = output_root / "covariance-long.tsv"
    write_json(summary_path, summary)
    write_tsv(matrix_path, matrix_rows)
    write_tsv(rows_path, long_rows)
    write_json(
        execution_path,
        {{
            "status": "ok",
            "case_id": case_payload["case_id"],
            "function_name": case_payload["function_name"],
            "input_fixture": case_payload["input_fixture"],
            "r_version": "4.6.0",
            "ape_version": "5.0.0",
            "outputs": {{
                "summary_json": str(summary_path),
                "covariance_matrix": str(matrix_path),
                "covariance_long": str(rows_path),
            }},
        }},
    )
    raise SystemExit(0)

if case_payload["operation"] == "tree-node-depth":
    tree = Phylo.read(case_payload["input_fixture"], "newick")
    rows = node_depth_rows(tree)
    tip_depths = [
        row["branch_length_depth"] for row in rows if row["node_kind"] == "tip"
    ]
    internal_depths = [
        row["branch_length_depth"] for row in rows if row["node_kind"] != "tip"
    ]
    summary = {{
        "node_count": len(rows),
        "tip_count": len(tip_depths),
        "internal_node_count": len(internal_depths),
        "rooted": is_rooted_tree(tree),
        "tip_labels": [terminal.name for terminal in tree.get_terminals()],
        "tree_is_ultrametric": (
            abs(max(tip_depths) - min(tip_depths)) <= 1e-12 if tip_depths else True
        ),
        "zero_branch_length_count": sum(
            1
            for clade in tree.find_clades(order="preorder")
            if clade is not tree.root and clade.branch_length == 0.0
        ),
        "minimum_tip_depth": min(tip_depths),
        "maximum_tip_depth": max(tip_depths),
        "minimum_internal_depth": min(internal_depths),
        "maximum_internal_depth": max(internal_depths),
    }}
    summary.update(SUMMARY_OVERRIDES)
    summary_path = output_root / "summary.json"
    rows_path = output_root / "node-depths.tsv"
    write_json(summary_path, summary)
    write_tsv(rows_path, rows)
    write_json(
        execution_path,
        {{
            "status": "ok",
            "case_id": case_payload["case_id"],
            "function_name": case_payload["function_name"],
            "input_fixture": case_payload["input_fixture"],
            "r_version": "4.6.0",
            "ape_version": "5.0.0",
            "outputs": {{
                "summary_json": str(summary_path),
                "node_depths": str(rows_path),
            }},
        }},
    )
    raise SystemExit(0)

if case_payload["operation"] == "tree-branching-times":
    tree = Phylo.read(case_payload["input_fixture"], "newick")
    depth_lookup = node_depth_lookup(tree)
    tip_depths = [
        depth_lookup[id(terminal)]
        for terminal in tree.get_terminals()
    ]
    rows = branching_time_rows(tree)
    summary = {{
        "internal_node_count": len(rows),
        "rooted": is_rooted_tree(tree),
        "tip_labels": [terminal.name for terminal in tree.get_terminals()],
        "tree_is_ultrametric": (
            abs(max(tip_depths) - min(tip_depths)) <= 1e-12 if tip_depths else True
        ),
        "root_age": max(tip_depths),
        "zero_branch_length_count": sum(
            1
            for clade in tree.find_clades(order="preorder")
            if clade is not tree.root and clade.branch_length == 0.0
        ),
        "minimum_tip_depth": min(tip_depths),
        "maximum_tip_depth": max(tip_depths),
        "max_tip_depth_deviation": max(tip_depths) - min(tip_depths),
        "tolerance": 1e-12,
    }}
    summary.update(SUMMARY_OVERRIDES)
    summary_path = output_root / "summary.json"
    rows_path = output_root / "branching-times.tsv"
    write_json(summary_path, summary)
    write_tsv(rows_path, rows)
    write_json(
        execution_path,
        {{
            "status": "ok",
            "case_id": case_payload["case_id"],
            "function_name": case_payload["function_name"],
            "input_fixture": case_payload["input_fixture"],
            "r_version": "4.6.0",
            "ape_version": "5.0.0",
            "outputs": {{
                "summary_json": str(summary_path),
                "branching_times": str(rows_path),
            }},
        }},
    )
    raise SystemExit(0)

if case_payload["operation"] == "tree-diversification-gamma-statistic":
    tree = Phylo.read(case_payload["input_fixture"], "newick")
    branching_times = sorted(
        row["branching_time"]
        for row in branching_time_rows(tree)
        if row["node_kind"] != "tip"
    )
    intervals = [branching_times[0]]
    intervals.extend(
        branching_times[index] - branching_times[index - 1]
        for index in range(1, len(branching_times))
    )
    waiting_times = list(reversed(intervals))
    tip_count = len(list(tree.get_terminals()))
    total_span = sum(
        multiplier * interval
        for multiplier, interval in zip(range(2, tip_count + 1), waiting_times, strict=True)
    )
    running = 0.0
    cumulative_total = 0.0
    for multiplier, interval in zip(range(2, tip_count), waiting_times[:-1], strict=True):
        running += multiplier * interval
        cumulative_total += running
    gamma_statistic = (
        ((cumulative_total / (tip_count - 2)) - (total_span / 2.0))
        / (total_span * math.sqrt(1.0 / (12.0 * (tip_count - 2))))
    )
    summary = {{
        "tip_count": tip_count,
        "rooted": is_rooted_tree(tree),
        "ultrametric": True,
        "bifurcating": True,
        "root_age": max(branching_times),
        "branching_time_count": len(branching_times),
        "interval_count": len(waiting_times),
        "minimum_branching_time": min(branching_times),
        "maximum_branching_time": max(branching_times),
        "gamma_statistic": gamma_statistic,
    }}
    summary.update(SUMMARY_OVERRIDES)
    rows_path = output_root / "gamma-statistic.tsv"
    summary_path = output_root / "summary.json"
    rows = [summary]
    write_json(summary_path, summary)
    write_tsv(rows_path, rows)
    write_json(
        execution_path,
        {{
            "status": "ok",
            "case_id": case_payload["case_id"],
            "function_name": case_payload["function_name"],
            "input_fixture": case_payload["input_fixture"],
            "r_version": "4.6.0",
            "ape_version": "5.0.0",
            "outputs": {{
                "summary_json": str(summary_path),
                "gamma_statistic": str(rows_path),
            }},
        }},
    )
    raise SystemExit(0)

if case_payload["operation"] == "tree-ultrametricity":
    tree = Phylo.read(case_payload["input_fixture"], "newick")
    rows = ultrametric_rows(tree)
    tip_depths = [row["root_to_tip_depth"] for row in rows]
    minimum_tip_depth = min(tip_depths)
    maximum_tip_depth = max(tip_depths)
    mean_tip_depth = sum(tip_depths) / len(tip_depths)
    max_tip_depth_deviation = maximum_tip_depth - minimum_tip_depth
    option = int(case_payload.get("ultrametric_option", 1))
    if option == 1:
        if math.isclose(maximum_tip_depth, 0.0, abs_tol=1e-15):
            criterion_value = 0.0 if math.isclose(max_tip_depth_deviation, 0.0, abs_tol=1e-15) else math.inf
        else:
            criterion_value = max_tip_depth_deviation / maximum_tip_depth
        criterion_name = "scaled-range"
    else:
        criterion_name = "variance"
        if len(tip_depths) <= 1:
            criterion_value = 0.0
        else:
            criterion_value = sum((depth - mean_tip_depth) ** 2 for depth in tip_depths) / (len(tip_depths) - 1)
    offending_taxa = sorted(
        {{
            row["tip_label"]
            for row in rows
            if math.isclose(row["root_to_tip_depth"], minimum_tip_depth, abs_tol=1e-12)
            or math.isclose(row["root_to_tip_depth"], maximum_tip_depth, abs_tol=1e-12)
        }}
    )
    if math.isclose(max_tip_depth_deviation, 0.0, abs_tol=1e-12):
        offending_taxa = []
    rows = [
        {{
            **row,
            "deviation_from_mean_depth": abs(row["root_to_tip_depth"] - mean_tip_depth),
            "deviation_from_min_depth": row["root_to_tip_depth"] - minimum_tip_depth,
            "deviation_from_max_depth": maximum_tip_depth - row["root_to_tip_depth"],
            "is_offending_taxon": row["tip_label"] in offending_taxa,
        }}
        for row in rows
    ]
    summary = {{
        "tip_count": len(rows),
        "rooted": is_rooted_tree(tree),
        "tip_labels": [terminal.name for terminal in tree.get_terminals()],
        "ultrametric": criterion_value <= case_payload["tolerance"],
        "criterion_name": criterion_name,
        "criterion_value": criterion_value,
        "tolerance": case_payload["tolerance"],
        "option": option,
        "minimum_tip_depth": minimum_tip_depth,
        "maximum_tip_depth": maximum_tip_depth,
        "mean_tip_depth": mean_tip_depth,
        "max_tip_depth_deviation": max_tip_depth_deviation,
        "root_age": maximum_tip_depth,
        "offending_taxa": offending_taxa,
    }}
    summary.update(SUMMARY_OVERRIDES)
    summary_path = output_root / "summary.json"
    rows_path = output_root / "ultrametric-diagnostics.tsv"
    write_json(summary_path, summary)
    write_tsv(rows_path, rows)
    write_json(
        execution_path,
        {{
            "status": "ok",
            "case_id": case_payload["case_id"],
            "function_name": case_payload["function_name"],
            "input_fixture": case_payload["input_fixture"],
            "r_version": "4.6.0",
            "ape_version": "5.0.0",
            "outputs": {{
                "summary_json": str(summary_path),
                "ultrametric_diagnostics": str(rows_path),
            }},
        }},
    )
    raise SystemExit(0)

if case_id in TABULAR_CASES:
    payload = TABULAR_CASES[case_id]
    summary = dict(payload["summary"])
    summary.update(SUMMARY_OVERRIDES)
    summary_path = output_root / "summary.json"
    rows_path = output_root / payload["rows_name"]
    rows = list(payload["rows"])
    if payload["rows_name"] == "distance-matrix.tsv":
        normalized_rows = []
        finite_distance_count = 0
        undefined_distance_count = 0
        infinite_distance_count = 0
        for row in rows:
            distance = row.get("distance", "")
            distance_status = row.get("distance_status")
            if distance_status is None:
                distance_status = "finite"
            if distance_status == "finite":
                finite_distance_count += 1
            elif distance_status == "undefined":
                undefined_distance_count += 1
            elif distance_status == "infinite":
                infinite_distance_count += 1
            normalized_rows.append(
                {{
                    "left_identifier": row["left_identifier"],
                    "right_identifier": row["right_identifier"],
                    "distance": distance,
                    "distance_status": distance_status,
                }}
            )
        rows = normalized_rows
        distance_model = str(case_payload["distance_model"]).lower()
        if distance_model == "raw":
            distance_model = "p-distance"
        elif distance_model == "jc69":
            distance_model = "jukes-cantor"
        elif distance_model == "k80":
            distance_model = "kimura-2-parameter"
        elif distance_model == "f81":
            distance_model = "felsenstein-81"
        elif distance_model == "tn93":
            distance_model = "tamura-nei-93"
        summary.setdefault("distance_model", distance_model)
        summary.setdefault("finite_distance_count", finite_distance_count)
        summary.setdefault("undefined_distance_count", undefined_distance_count)
        summary.setdefault("infinite_distance_count", infinite_distance_count)
    write_json(summary_path, summary)
    write_tsv(rows_path, rows)
    outputs = {{"summary_json": str(summary_path)}}
    if payload["rows_name"] == "base-frequency.tsv":
        outputs["base_frequency"] = str(rows_path)
    elif payload["rows_name"] == "distance-matrix.tsv":
        outputs["distance_matrix"] = str(rows_path)
    else:
        outputs["translation"] = str(rows_path)
    write_json(
        execution_path,
        {{
            "status": "ok",
            "case_id": case_payload["case_id"],
            "function_name": case_payload["function_name"],
            "input_fixture": case_payload["input_fixture"],
            "r_version": "4.6.0",
            "ape_version": "5.0.0",
            "outputs": outputs,
        }},
    )
    raise SystemExit(0)

if case_id in TREE_CASES:
    payload = TREE_CASES[case_id]
    summary = dict(payload["summary"])
    summary.update(SUMMARY_OVERRIDES)
    summary_path = output_root / "summary.json"
    newick_path = output_root / "normalized-tree.nwk"
    write_json(summary_path, summary)
    newick_path.write_text(
        NORMALIZED_TREE_OVERRIDES.get(case_id, payload["newick"]) + "\\n",
        encoding="utf-8",
    )
    write_json(
        execution_path,
        {{
            "status": "ok",
            "case_id": case_payload["case_id"],
            "function_name": case_payload["function_name"],
            "input_fixture": case_payload["input_fixture"],
            "r_version": "4.6.0",
            "ape_version": "5.0.0",
            "outputs": {{
                "summary_json": str(summary_path),
                "normalized_tree": str(newick_path),
            }},
        }},
    )
    raise SystemExit(0)

if case_id in ERROR_CASES:
    payload = ERROR_CASES[case_id]
    write_json(
        execution_path,
        {{
            "status": "failed",
            "mismatch_reason": "reference_execution_failed",
            "error_type": payload["error_type"],
            "message": payload["message"],
            "case_id": case_payload["case_id"],
            "function_name": case_payload["function_name"],
            "input_fixture": case_payload["input_fixture"],
            "r_version": "4.6.0",
            "ape_version": "5.0.0",
        }},
    )
    raise SystemExit(0)

write_json(
    execution_path,
    {{
        "status": "failed",
        "mismatch_reason": "unsupported_operation",
        "message": f"unsupported ape parity operation: {{case_payload['operation']}}",
        "case_id": case_payload["case_id"],
        "function_name": case_payload["function_name"],
        "input_fixture": case_payload["input_fixture"],
        "r_version": "4.6.0",
        "ape_version": "5.0.0",
    }},
)
""",
    )
