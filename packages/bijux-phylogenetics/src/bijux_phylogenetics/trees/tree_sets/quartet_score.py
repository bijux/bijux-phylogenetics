from __future__ import annotations

import csv
from pathlib import Path

from .contracts import CandidateTreeQuartetScoreReport
from .quartet_concordance import compute_gene_tree_quartet_concordance_factors


def compute_candidate_tree_quartet_score(
    candidate_tree_path: Path,
    gene_tree_set_path: Path,
) -> CandidateTreeQuartetScoreReport:
    """Score one candidate species tree by quartet agreement against a gene-tree set."""
    concordance_report = compute_gene_tree_quartet_concordance_factors(
        candidate_tree_path,
        gene_tree_set_path,
    )
    informative_quartet_count = concordance_report.informative_quartet_count
    quartet_score = concordance_report.concordant_quartet_count
    return CandidateTreeQuartetScoreReport(
        candidate_tree_path=candidate_tree_path,
        gene_tree_set_path=gene_tree_set_path,
        tree_count=concordance_report.tree_count,
        processing=concordance_report.processing,
        shared_taxa=concordance_report.shared_taxa,
        branch_count=concordance_report.branch_count,
        total_quartet_count=concordance_report.total_quartet_count,
        concordant_quartet_count=concordance_report.concordant_quartet_count,
        discordant_first_quartet_count=(
            concordance_report.discordant_first_quartet_count
        ),
        discordant_second_quartet_count=(
            concordance_report.discordant_second_quartet_count
        ),
        uninformative_quartet_count=concordance_report.uninformative_quartet_count,
        informative_quartet_count=informative_quartet_count,
        quartet_score=quartet_score,
        normalized_quartet_score=(
            round(quartet_score / informative_quartet_count, 15)
            if informative_quartet_count
            else None
        ),
        rows=concordance_report.rows,
    )


def write_candidate_tree_quartet_score_table(
    path: Path,
    report: CandidateTreeQuartetScoreReport,
) -> Path:
    """Write one candidate-tree quartet score table as TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "branch_id",
                "left_taxa",
                "right_taxa",
                "quartet_count_per_tree",
                "concordant_quartet_count",
                "discordant_first_quartet_count",
                "discordant_second_quartet_count",
                "uninformative_quartet_count",
                "informative_quartet_count",
                "concordance_factor",
                "quartet_score",
                "normalized_quartet_score",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.rows:
            writer.writerow(
                {
                    "branch_id": row.branch_id,
                    "left_taxa": "|".join(row.left_taxa),
                    "right_taxa": "|".join(row.right_taxa),
                    "quartet_count_per_tree": row.quartet_count_per_tree,
                    "concordant_quartet_count": row.concordant_quartet_count,
                    "discordant_first_quartet_count": (
                        row.discordant_first_quartet_count
                    ),
                    "discordant_second_quartet_count": (
                        row.discordant_second_quartet_count
                    ),
                    "uninformative_quartet_count": row.uninformative_quartet_count,
                    "informative_quartet_count": row.informative_quartet_count,
                    "concordance_factor": (
                        ""
                        if row.concordance_factor is None
                        else format(row.concordance_factor, ".15g")
                    ),
                    "quartet_score": report.quartet_score,
                    "normalized_quartet_score": (
                        ""
                        if report.normalized_quartet_score is None
                        else format(report.normalized_quartet_score, ".15g")
                    ),
                }
            )
    return path
